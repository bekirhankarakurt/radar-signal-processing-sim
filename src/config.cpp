#include "radar/config.hpp"

#include "radar/complex_math.hpp"

#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace radar {

double RadarParameters::wavelength() const noexcept { return kSpeedOfLight / carrier_hz; }

std::size_t RadarParameters::pulse_samples() const noexcept {
    return static_cast<std::size_t>(std::llround(pulse_width_s * sample_rate_hz));
}

std::size_t RadarParameters::range_bins() const noexcept {
    return static_cast<std::size_t>(std::ceil(2.0 * max_range_m / kSpeedOfLight * sample_rate_hz)) + 1;
}

std::size_t RadarParameters::fast_time_samples() const noexcept { return range_bins() + pulse_samples(); }

double RadarParameters::range_bin_m() const noexcept { return kSpeedOfLight / (2.0 * sample_rate_hz); }

double RadarParameters::range_resolution_m() const noexcept { return kSpeedOfLight / (2.0 * bandwidth_hz); }

double RadarParameters::unambiguous_range_m() const noexcept { return kSpeedOfLight / (2.0 * prf_hz); }

double RadarParameters::max_unambiguous_velocity_mps() const noexcept { return wavelength() * prf_hz / 4.0; }

double RadarParameters::noise_power_w() const noexcept {
    return kBoltzmann * kStandardTemperature * db_to_linear(noise_figure_db) * sample_rate_hz;
}

namespace {

std::string trim(std::string_view s) {
    std::size_t b = 0;
    std::size_t e = s.size();
    while (b < e && std::isspace(static_cast<unsigned char>(s[b])) != 0) {
        ++b;
    }
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])) != 0) {
        --e;
    }
    return std::string(s.substr(b, e - b));
}

std::string lower(std::string s) {
    for (auto& c : s) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

[[noreturn]] void fail(const IniSection& section, std::string_view key, std::string_view message) {
    std::ostringstream os;
    os << "[" << section.name << "] (line " << section.line << ") " << key << ": " << message;
    throw std::invalid_argument(os.str());
}

/// Typed accessor that records which keys were consumed so that unknown keys
/// (usually typos) can be reported instead of silently ignored.
class SectionReader {
public:
    explicit SectionReader(const IniSection& section) : section_(section) {}

    void read(std::string_view key, double& out) {
        if (const auto* v = find(key)) {
            std::size_t pos = 0;
            try {
                out = std::stod(*v, &pos);
            } catch (const std::exception&) {
                fail(section_, key, "expected a number, got '" + *v + "'");
            }
            if (pos != v->size() || !std::isfinite(out)) {
                fail(section_, key, "expected a number, got '" + *v + "'");
            }
        }
    }

    void read(std::string_view key, std::size_t& out) {
        double d = static_cast<double>(out);
        read(key, d);
        if (d < 0.0 || d != std::floor(d) || d > 1e12) {
            fail(section_, key, "expected a non-negative integer");
        }
        out = static_cast<std::size_t>(d);
    }

    void read_u64(std::string_view key, std::uint64_t& out) {
        if (const auto* v = find(key)) {
            try {
                std::size_t pos = 0;
                out = std::stoull(*v, &pos);
                if (pos != v->size()) {
                    throw std::invalid_argument("trailing characters");
                }
            } catch (const std::exception&) {
                fail(section_, key, "expected an unsigned integer");
            }
        }
    }

    void read(std::string_view key, int& out) {
        double d = static_cast<double>(out);
        read(key, d);
        if (d != std::floor(d) || std::abs(d) > 1e6) {
            fail(section_, key, "expected an integer");
        }
        out = static_cast<int>(d);
    }

    void read(std::string_view key, bool& out) {
        if (const auto* v = find(key)) {
            const auto s = lower(*v);
            if (s == "true" || s == "yes" || s == "on" || s == "1") {
                out = true;
            } else if (s == "false" || s == "no" || s == "off" || s == "0") {
                out = false;
            } else {
                fail(section_, key, "expected true/false");
            }
        }
    }

    void read(std::string_view key, std::string& out) {
        if (const auto* v = find(key)) {
            out = *v;
        }
    }

    void read_window(std::string_view prefix, WindowSpec& spec) {
        const std::string key(prefix);
        if (const auto* v = find(key)) {
            const auto type = parse_window_type(*v);
            if (!type) {
                fail(section_, key, "unknown window '" + *v + "' (none|hann|hamming|taylor)");
            }
            spec.type = *type;
        }
        read(key + "_taylor_nbar", spec.taylor_nbar);
        read(key + "_taylor_sll_db", spec.taylor_sidelobe_db);
    }

    template <typename Enum, std::size_t N>
    void read_enum(std::string_view key, Enum& out, const std::pair<const char*, Enum> (&options)[N]) {
        if (const auto* v = find(key)) {
            const auto s = lower(*v);
            for (const auto& [name, value] : options) {
                if (s == name) {
                    out = value;
                    return;
                }
            }
            std::string allowed;
            for (const auto& option : options) {
                allowed += allowed.empty() ? "" : "|";
                allowed += option.first;
            }
            fail(section_, key, "unknown value '" + *v + "' (" + allowed + ")");
        }
    }

    void finish() const {
        for (const auto& [key, value] : section_.values) {
            if (used_.count(key) == 0) {
                fail(section_, key, "unknown key");
            }
        }
    }

private:
    const std::string* find(std::string_view key) {
        const auto it = section_.values.find(key);
        if (it == section_.values.end()) {
            return nullptr;
        }
        used_.insert(it->first);
        return &it->second;
    }

    const IniSection& section_;
    std::set<std::string, std::less<>> used_;
};

constexpr std::pair<const char*, Swerling> kSwerlingNames[] = {{"0", Swerling::Zero}, {"1", Swerling::One}};
constexpr std::pair<const char*, MtiType> kMtiNames[] = {
    {"none", MtiType::None}, {"two_pulse", MtiType::TwoPulse}, {"three_pulse", MtiType::ThreePulse}};
constexpr std::pair<const char*, CfarType> kCfarNames[] = {{"ca", CfarType::CellAveraging},
                                                           {"os", CfarType::OrderedStatistic}};

void require(bool ok, const std::string& message) {
    if (!ok) {
        throw std::invalid_argument(message);
    }
}

}  // namespace

std::vector<IniSection> parse_ini(std::string_view text) {
    std::vector<IniSection> sections;
    std::size_t line_no = 0;
    std::size_t pos = 0;
    while (pos <= text.size()) {
        const std::size_t eol = text.find('\n', pos);
        std::string_view raw = text.substr(pos, eol == std::string_view::npos ? std::string_view::npos : eol - pos);
        pos = eol == std::string_view::npos ? text.size() + 1 : eol + 1;
        ++line_no;

        const std::size_t comment = raw.find_first_of(";#");
        const std::string line = trim(raw.substr(0, comment));
        if (line.empty()) {
            continue;
        }
        if (line.front() == '[') {
            if (line.back() != ']' || line.size() < 3) {
                throw std::invalid_argument("line " + std::to_string(line_no) + ": malformed section header");
            }
            sections.push_back({lower(trim(std::string_view(line).substr(1, line.size() - 2))), {}, line_no});
            continue;
        }
        const std::size_t eq = line.find('=');
        if (eq == std::string::npos) {
            throw std::invalid_argument("line " + std::to_string(line_no) + ": expected key = value");
        }
        if (sections.empty()) {
            throw std::invalid_argument("line " + std::to_string(line_no) + ": key outside of a section");
        }
        const std::string key = lower(trim(std::string_view(line).substr(0, eq)));
        const std::string value = trim(std::string_view(line).substr(eq + 1));
        if (key.empty()) {
            throw std::invalid_argument("line " + std::to_string(line_no) + ": empty key");
        }
        if (!sections.back().values.emplace(key, value).second) {
            throw std::invalid_argument("line " + std::to_string(line_no) + ": duplicate key '" + key + "'");
        }
    }
    return sections;
}

Scenario scenario_from_ini(std::string_view text) {
    Scenario s;
    for (const auto& section : parse_ini(text)) {
        SectionReader r(section);
        if (section.name == "simulation") {
            r.read("name", s.name);
            r.read("cpis", s.cpis);
            r.read("cpi_interval_s", s.cpi_interval_s);
            r.read_u64("seed", s.seed);
            r.read("noise", s.noise);
        } else if (section.name == "radar") {
            auto& p = s.radar;
            r.read("carrier_hz", p.carrier_hz);
            r.read("bandwidth_hz", p.bandwidth_hz);
            r.read("pulse_width_s", p.pulse_width_s);
            r.read("prf_hz", p.prf_hz);
            r.read("pulses", p.pulses);
            r.read("sample_rate_hz", p.sample_rate_hz);
            r.read("peak_power_w", p.peak_power_w);
            r.read("antenna_gain_db", p.antenna_gain_db);
            r.read("noise_figure_db", p.noise_figure_db);
            r.read("system_loss_db", p.system_loss_db);
            r.read("max_range_m", p.max_range_m);
        } else if (section.name == "processing") {
            auto& p = s.processing;
            r.read_window("range_window", p.range_window);
            r.read_window("doppler_window", p.doppler_window);
            r.read_enum("mti", p.mti, kMtiNames);
            r.read_enum("cfar", p.cfar, kCfarNames);
            r.read("cfar_pfa", p.cfar_pfa);
            r.read("guard_range", p.guard_range);
            r.read("guard_doppler", p.guard_doppler);
            r.read("train_range", p.train_range);
            r.read("train_doppler", p.train_doppler);
            r.read("os_rank", p.os_rank);
            r.read("doppler_fft", p.doppler_fft);
            r.read("clutter_notch_mps", p.clutter_notch_mps);
        } else if (section.name == "clutter") {
            auto& c = s.clutter;
            r.read("enabled", c.enabled);
            r.read("cnr_db", c.cnr_db);
            r.read("min_range_m", c.min_range_m);
            r.read("max_range_m", c.max_range_m);
            r.read("spectral_width_mps", c.spectral_width_mps);
        } else if (section.name == "tracker") {
            auto& t = s.tracker;
            r.read("confirm_hits", t.confirm_hits);
            r.read("confirm_window", t.confirm_window);
            r.read("max_misses", t.max_misses);
            r.read("range_gate_m", t.range_gate_m);
            r.read("velocity_gate_mps", t.velocity_gate_mps);
            r.read("alpha", t.alpha);
            r.read("beta", t.beta);
        } else if (section.name == "monte_carlo") {
            auto& m = s.monte_carlo;
            r.read("snr_min_db", m.snr_min_db);
            r.read("snr_max_db", m.snr_max_db);
            r.read("snr_step_db", m.snr_step_db);
            r.read("trials", m.trials);
            r.read("pfa_trials", m.pfa_trials);
            r.read_enum("swerling", m.swerling, kSwerlingNames);
            r.read("range_m", m.range_m);
            if (section.values.count("velocity_mps") != 0) {
                double v = 0.0;
                r.read("velocity_mps", v);
                m.velocity_mps = v;
            }
        } else if (section.name == "target") {
            Target t;
            r.read("range_m", t.range_m);
            r.read("velocity_mps", t.velocity_mps);
            r.read("rcs_m2", t.rcs_m2);
            r.read_enum("swerling", t.swerling, kSwerlingNames);
            s.targets.push_back(t);
        } else {
            throw std::invalid_argument("line " + std::to_string(section.line) + ": unknown section [" +
                                        section.name + "]");
        }
        r.finish();
    }
    s.validate();
    return s;
}

Scenario load_scenario(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("cannot open scenario file '" + path.string() + "'");
    }
    std::ostringstream buffer;
    buffer << in.rdbuf();
    try {
        return scenario_from_ini(buffer.str());
    } catch (const std::invalid_argument& e) {
        throw std::invalid_argument(path.string() + ": " + e.what());
    }
}

void Scenario::validate() const {
    const auto& r = radar;
    require(r.carrier_hz > 0.0, "radar.carrier_hz must be positive");
    require(r.bandwidth_hz > 0.0, "radar.bandwidth_hz must be positive");
    require(r.sample_rate_hz >= r.bandwidth_hz, "radar.sample_rate_hz must be >= bandwidth_hz (complex sampling)");
    require(r.pulse_width_s > 0.0 && r.pulse_samples() >= 2, "radar.pulse_width_s must span at least 2 samples");
    require(r.prf_hz > 0.0, "radar.prf_hz must be positive");
    require(r.pulse_width_s < r.pri(), "radar.pulse_width_s must be shorter than the PRI");
    require(r.pulses >= 4 && r.pulses <= 65536, "radar.pulses must be in [4, 65536]");
    require(r.peak_power_w > 0.0, "radar.peak_power_w must be positive");
    require(r.max_range_m > 0.0, "radar.max_range_m must be positive");
    require(r.max_range_m + 0.5 * kSpeedOfLight * r.pulse_width_s <= r.unambiguous_range_m(),
            "radar.max_range_m plus the pulse length exceeds the unambiguous range c/(2 PRF)");
    require(r.fast_time_samples() <= (std::size_t{1} << 22U), "fast-time window too long (> 4M samples)");

    const auto& p = processing;
    require(p.cfar_pfa > 0.0 && p.cfar_pfa < 1.0, "processing.cfar_pfa must be in (0, 1)");
    require(p.train_range + p.train_doppler > 0, "processing: CFAR needs training cells");
    require(p.clutter_notch_mps >= 0.0 && p.clutter_notch_mps < r.max_unambiguous_velocity_mps(),
            "processing.clutter_notch_mps must be in [0, unambiguous velocity)");
    require(p.os_rank > 0.0 && p.os_rank <= 1.0, "processing.os_rank must be in (0, 1]");
    require(p.doppler_fft == 0 || is_power_of_two(p.doppler_fft), "processing.doppler_fft must be a power of two");
    require(p.doppler_fft == 0 || p.doppler_fft >= r.pulses, "processing.doppler_fft must be >= radar.pulses");
    const std::size_t doppler_bins = p.doppler_fft == 0 ? next_power_of_two(r.pulses) : p.doppler_fft;
    require(2 * (p.guard_doppler + p.train_doppler) + 1 <= doppler_bins,
            "processing: CFAR Doppler window is wider than the Doppler axis");
    require(2 * (p.guard_range + p.train_range) + 1 <= r.range_bins(),
            "processing: CFAR range window is longer than the range axis");

    require(clutter.max_range_m >= clutter.min_range_m, "clutter.max_range_m must be >= min_range_m");
    require(clutter.spectral_width_mps >= 0.0, "clutter.spectral_width_mps must be non-negative");

    require(tracker.confirm_hits >= 1 && tracker.confirm_hits <= tracker.confirm_window,
            "tracker: need 1 <= confirm_hits <= confirm_window");
    require(tracker.alpha > 0.0 && tracker.alpha <= 1.0, "tracker.alpha must be in (0, 1]");
    require(tracker.beta >= 0.0 && tracker.beta <= 2.0, "tracker.beta must be in [0, 2]");

    require(monte_carlo.snr_step_db > 0.0, "monte_carlo.snr_step_db must be positive");
    require(monte_carlo.snr_max_db >= monte_carlo.snr_min_db, "monte_carlo.snr_max_db must be >= snr_min_db");
    require(monte_carlo.trials >= 1, "monte_carlo.trials must be >= 1");

    require(cpis >= 1 && cpis <= 10000, "simulation.cpis must be in [1, 10000]");
    require(cpi_interval_s == 0.0 || cpi_interval_s >= cpi_duration_s(),
            "simulation.cpi_interval_s must be 0 or at least one CPI (pulses / PRF)");
    for (std::size_t i = 0; i < targets.size(); ++i) {
        const auto& t = targets[i];
        const std::string where = "target " + std::to_string(i + 1) + ": ";
        require(t.range_m >= 0.0 && t.range_m <= r.max_range_m, where + "range_m must be in [0, radar.max_range_m]");
        require(t.rcs_m2 > 0.0, where + "rcs_m2 must be positive");
    }
}

double Scenario::cpi_duration_s() const noexcept { return static_cast<double>(radar.pulses) * radar.pri(); }

double Scenario::cpi_start_s(std::size_t cpi) const noexcept {
    const double interval = cpi_interval_s > 0.0 ? cpi_interval_s : cpi_duration_s();
    return static_cast<double>(cpi) * interval;
}

std::string_view to_string(MtiType type) noexcept {
    switch (type) {
        case MtiType::None:
            return "none";
        case MtiType::TwoPulse:
            return "two_pulse";
        case MtiType::ThreePulse:
            return "three_pulse";
    }
    return "unknown";
}

std::string_view to_string(CfarType type) noexcept {
    return type == CfarType::CellAveraging ? "CA-CFAR" : "OS-CFAR";
}

std::string_view to_string(Swerling model) noexcept { return model == Swerling::Zero ? "Swerling 0" : "Swerling 1"; }

}  // namespace radar
