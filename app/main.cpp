#include "radar/chain.hpp"
#include "radar/complex_math.hpp"
#include "radar/config.hpp"
#include "radar/monte_carlo.hpp"
#include "radar/output.hpp"
#include "radar/signal_model.hpp"

#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr std::string_view kUsage = R"(radar_sim - pulse-Doppler radar signal processing simulator

Usage:
  radar_sim run <scenario.ini> [options]         simulate and process CPIs, write map/detections/tracks
  radar_sim montecarlo <scenario.ini> [options]  Pd vs SNR and empirical Pfa through the full chain
  radar_sim info <scenario.ini>                  print derived radar parameters and link budget
  radar_sim --help | --version

Options:
  -o, --output <dir>   output directory (default: out)
  --cpis <n>           override [simulation] cpis
  --seed <n>           override [simulation] seed
  --trials <n>         override [monte_carlo] trials
  --cfar <ca|os>       override [processing] cfar
  --mti <none|two_pulse|three_pulse>
                       override [processing] mti
  --swerling <0|1>     override [monte_carlo] swerling
  --threads <n>        OpenMP thread count (default: all cores)
  --dynamic-range <dB> PPM image dynamic range (default: 60)
  -q, --quiet          only print the summary and timing
)";

struct Options {
    std::string command;
    std::filesystem::path scenario;
    std::filesystem::path output = "out";
    std::optional<std::size_t> cpis;
    std::optional<std::uint64_t> seed;
    std::optional<std::size_t> trials;
    std::optional<std::string> cfar;
    std::optional<std::string> mti;
    std::optional<std::string> swerling;
    std::optional<int> threads;
    double dynamic_range_db = 60.0;
    bool quiet = false;
};

class UsageError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

unsigned long long parse_unsigned(std::string_view flag, const std::string& value) {
    try {
        std::size_t pos = 0;
        const auto v = std::stoull(value, &pos);
        if (pos == value.size() && value.front() != '-') {
            return v;
        }
    } catch (const std::exception&) {
    }
    throw UsageError(std::string(flag) + " expects a non-negative integer, got '" + value + "'");
}

Options parse_args(const std::vector<std::string>& args) {
    Options o;
    std::size_t i = 0;
    auto next = [&](std::string_view flag) -> const std::string& {
        if (i + 1 >= args.size()) {
            throw UsageError(std::string(flag) + " needs a value");
        }
        return args[++i];
    };
    std::vector<std::string> positional;
    for (; i < args.size(); ++i) {
        const std::string& a = args[i];
        if (a == "-o" || a == "--output") {
            o.output = next(a);
        } else if (a == "--cpis") {
            o.cpis = static_cast<std::size_t>(parse_unsigned(a, next(a)));
        } else if (a == "--seed") {
            o.seed = parse_unsigned(a, next(a));
        } else if (a == "--trials") {
            o.trials = static_cast<std::size_t>(parse_unsigned(a, next(a)));
        } else if (a == "--cfar") {
            o.cfar = next(a);
        } else if (a == "--mti") {
            o.mti = next(a);
        } else if (a == "--swerling") {
            o.swerling = next(a);
        } else if (a == "--threads") {
            o.threads = static_cast<int>(parse_unsigned(a, next(a)));
        } else if (a == "--dynamic-range") {
            o.dynamic_range_db = static_cast<double>(parse_unsigned(a, next(a)));
        } else if (a == "-q" || a == "--quiet") {
            o.quiet = true;
        } else if (!a.empty() && a.front() == '-') {
            throw UsageError("unknown option '" + a + "'");
        } else {
            positional.push_back(a);
        }
    }
    if (positional.size() != 2) {
        throw UsageError("expected a command and a scenario file");
    }
    o.command = positional[0];
    o.scenario = positional[1];
    if (o.command != "run" && o.command != "montecarlo" && o.command != "info") {
        throw UsageError("unknown command '" + o.command + "'");
    }
    return o;
}

void apply_overrides(const Options& o, radar::Scenario& s) {
    if (o.cpis) {
        s.cpis = *o.cpis;
    }
    if (o.seed) {
        s.seed = *o.seed;
    }
    if (o.trials) {
        s.monte_carlo.trials = *o.trials;
    }
    if (o.cfar) {
        if (*o.cfar == "ca") {
            s.processing.cfar = radar::CfarType::CellAveraging;
        } else if (*o.cfar == "os") {
            s.processing.cfar = radar::CfarType::OrderedStatistic;
        } else {
            throw UsageError("--cfar expects ca or os");
        }
    }
    if (o.mti) {
        if (*o.mti == "none") {
            s.processing.mti = radar::MtiType::None;
        } else if (*o.mti == "two_pulse") {
            s.processing.mti = radar::MtiType::TwoPulse;
        } else if (*o.mti == "three_pulse") {
            s.processing.mti = radar::MtiType::ThreePulse;
        } else {
            throw UsageError("--mti expects none, two_pulse or three_pulse");
        }
    }
    if (o.swerling) {
        if (*o.swerling == "0") {
            s.monte_carlo.swerling = radar::Swerling::Zero;
        } else if (*o.swerling == "1") {
            s.monte_carlo.swerling = radar::Swerling::One;
        } else {
            throw UsageError("--swerling expects 0 or 1");
        }
    }
    s.validate();
}

int threads_in_use() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

void print_info(const radar::Scenario& s, std::ostream& os) {
    const auto& r = s.radar;
    os << std::fixed << std::setprecision(2);
    os << "Scenario: " << s.name << "\n"
       << "  wavelength            " << r.wavelength() * 100.0 << " cm\n"
       << "  range resolution      " << r.range_resolution_m() << " m (bin " << r.range_bin_m() << " m)\n"
       << "  unambiguous range     " << r.unambiguous_range_m() / 1e3 << " km (processed to "
       << r.max_range_m / 1e3 << " km)\n"
       << "  unambiguous velocity  +/-" << r.max_unambiguous_velocity_mps() << " m/s\n"
       << "  velocity resolution   " << r.wavelength() * r.prf_hz / (2.0 * static_cast<double>(r.pulses))
       << " m/s\n"
       << "  time-bandwidth        " << r.bandwidth_hz * r.pulse_width_s << "\n"
       << "  CPI duration          " << static_cast<double>(r.pulses) * r.pri() * 1e3 << " ms\n"
       << "  data cube             " << r.pulses << " pulses x " << r.fast_time_samples() << " samples\n";
    os << "\n  " << std::left << std::setw(8) << "target" << std::right << std::setw(10) << "range m" << std::setw(10)
       << "vel m/s" << std::setw(10) << "RCS m2" << std::setw(16) << "SNR/sample dB" << "\n";
    for (std::size_t i = 0; i < s.targets.size(); ++i) {
        const auto& t = s.targets[i];
        os << "  " << std::left << std::setw(8) << i + 1 << std::right << std::setw(10) << std::setprecision(0)
           << t.range_m << std::setw(10) << std::setprecision(1) << t.velocity_mps << std::setw(10)
           << std::setprecision(3) << t.rcs_m2 << std::setw(16) << std::setprecision(1)
           << radar::linear_to_db(radar::single_sample_snr(r, t.range_m, t.rcs_m2)) << "\n";
    }
}

int run(const Options& o) {
    radar::Scenario scenario = radar::load_scenario(o.scenario);
    apply_overrides(o, scenario);
#ifdef _OPENMP
    if (o.threads && *o.threads > 0) {
        omp_set_num_threads(*o.threads);
    }
#endif
    if (o.command == "info") {
        print_info(scenario, std::cout);
        return EXIT_SUCCESS;
    }
    std::filesystem::create_directories(o.output);
    radar::StageTimer timer;

    if (o.command == "montecarlo") {
        std::cout << "Scenario '" << scenario.name << "', " << threads_in_use() << " thread(s)\n";
        const auto result = radar::run_monte_carlo(scenario, timer);
        radar::print_monte_carlo(std::cout, result, scenario);
        radar::write_pd_csv(o.output / "pd_curve.csv", result);
        radar::write_timing_csv(o.output / "timing.csv", timer, 1);
        radar::print_timing(std::cout, timer, 1, 0);
        std::cout << "\nWrote " << (o.output / "pd_curve.csv").string() << "\n";
        return EXIT_SUCCESS;
    }

    std::cout << "Scenario '" << scenario.name << "': " << scenario.targets.size() << " target(s), "
              << scenario.cpis << " CPI(s), " << scenario.radar.pulses << " pulses x "
              << scenario.radar.fast_time_samples() << " samples, " << threads_in_use() << " thread(s)\n";
    const auto result = radar::run_simulation(scenario, timer);
    if (!o.quiet) {
        for (std::size_t c = 0; c < result.cpis.size(); ++c) {
            radar::print_detections(std::cout, result.cpis[c].detections, c);
        }
    }
    radar::print_tracks(std::cout, result.tracks);
    {
        const auto s = timer.scope("output");
        const auto& last = result.cpis.back();
        radar::write_ppm(o.output / "range_doppler.ppm", last.rd_map, last.detections, {o.dynamic_range_db, true});
        radar::write_rd_csv(o.output / "range_doppler.csv", last.rd_map, result.axes);
        radar::write_detections_csv(o.output / "detections.csv", result.detections);
        radar::write_tracks_csv(o.output / "tracks.csv", result.tracks);
        radar::write_html_viewer(o.output / "viewer.html", scenario, result);
    }
    radar::write_timing_csv(o.output / "timing.csv", timer, scenario.cpis);
    radar::print_timing(std::cout, timer, scenario.cpis, result.raw_samples);
    std::cout << "\nOutputs in " << o.output.string()
              << "/: range_doppler.ppm, range_doppler.csv, detections.csv, tracks.csv, timing.csv, viewer.html\n";
    return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char** argv) {
    const std::vector<std::string> args(argv + 1, argv + argc);
    if (args.empty() || args[0] == "-h" || args[0] == "--help") {
        std::cout << kUsage;
        return args.empty() ? EXIT_FAILURE : EXIT_SUCCESS;
    }
    if (args[0] == "--version") {
        std::cout << "radar_sim " << RADAR_VERSION << "\n";
        return EXIT_SUCCESS;
    }
    try {
        return run(parse_args(args));
    } catch (const UsageError& e) {
        std::cerr << "error: " << e.what() << "\n\n" << kUsage;
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }
}
