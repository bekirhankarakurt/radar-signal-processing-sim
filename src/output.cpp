#include "radar/output.hpp"

#include "radar/complex_math.hpp"
#include "radar/viewer_template.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace radar {
namespace {

std::ofstream open_output(const std::filesystem::path& path, std::ios::openmode mode = std::ios::out) {
    std::ofstream out(path, mode);
    if (!out) {
        throw std::runtime_error("cannot write '" + path.string() + "'");
    }
    out << std::setprecision(10);
    return out;
}

/// Google "Turbo" colormap, polynomial approximation. x in [0, 1].
std::array<std::uint8_t, 3> turbo(double x) {
    x = std::clamp(x, 0.0, 1.0);
    const double r = 0.13572138 + x * (4.61539260 + x * (-42.66032258 + x * (132.13108234 + x * (-152.94239396 + x * 59.28637943))));
    const double g = 0.09140261 + x * (2.19418839 + x * (4.84296658 + x * (-14.18503333 + x * (4.27729857 + x * 2.82956604))));
    const double b = 0.10667330 + x * (12.64194608 + x * (-60.58204836 + x * (110.36276771 + x * (-89.90310912 + x * 27.34824973))));
    auto to_byte = [](double v) { return static_cast<std::uint8_t>(std::lround(std::clamp(v, 0.0, 1.0) * 255.0)); };
    return {to_byte(r), to_byte(g), to_byte(b)};
}

double map_peak_db(const RealMatrix& map) {
    const auto it = std::max_element(map.data().begin(), map.data().end());
    return linear_to_db(it == map.data().end() ? 1.0 : *it);
}

std::string json_escape(std::string_view s) {
    std::string out;
    for (const char c : s) {
        switch (c) {
            case '"':
                out += "\\\"";
                break;
            case '\\':
                out += "\\\\";
                break;
            case '<':
                out += "\\u003c";
                break;
            case '\n':
                out += "\\n";
                break;
            default:
                if (static_cast<unsigned char>(c) >= 0x20) {
                    out += c;
                }
        }
    }
    return out;
}

}  // namespace

void write_ppm(const std::filesystem::path& path, const RealMatrix& rd_map, const std::vector<Detection>& detections,
               const ImageOptions& options) {
    if (rd_map.empty()) {
        throw std::invalid_argument("empty range-Doppler map");
    }
    const std::size_t range_bins = rd_map.rows();
    const std::size_t doppler_bins = rd_map.cols();
    // Aim for roughly 512 x 600 pixels: stretch Doppler, max-pool range.
    const std::size_t x_scale = std::max<std::size_t>(1, 512 / doppler_bins);
    const std::size_t pool = std::max<std::size_t>(1, (range_bins + 599) / 600);
    const std::size_t width = doppler_bins * x_scale;
    const std::size_t height = (range_bins + pool - 1) / pool;
    const double top = map_peak_db(rd_map);
    const double bottom = top - options.dynamic_range_db;

    std::vector<std::uint8_t> pixels(width * height * 3);
    for (std::size_t y = 0; y < height; ++y) {
        const std::size_t r0 = (height - 1 - y) * pool;  // near range at the bottom
        const std::size_t r1 = std::min(range_bins, r0 + pool);
        for (std::size_t d = 0; d < doppler_bins; ++d) {
            double p = 0.0;
            for (std::size_t r = r0; r < r1; ++r) {
                p = std::max(p, rd_map(r, d));
            }
            const auto rgb = turbo((linear_to_db(p) - bottom) / options.dynamic_range_db);
            for (std::size_t sx = 0; sx < x_scale; ++sx) {
                std::copy(rgb.begin(), rgb.end(), pixels.begin() + static_cast<std::ptrdiff_t>((y * width + d * x_scale + sx) * 3));
            }
        }
    }
    if (options.mark_detections) {
        for (const auto& det : detections) {
            const long cx = std::lround((det.doppler_bin + 0.5) * static_cast<double>(x_scale));
            const long cy = static_cast<long>(height - 1) - std::lround(det.range_bin / static_cast<double>(pool));
            for (long k = -6; k <= 6; ++k) {
                for (const auto& [x, y] : {std::pair{cx - 6, cy + k}, std::pair{cx + 6, cy + k}, std::pair{cx + k, cy - 6},
                                           std::pair{cx + k, cy + 6}}) {
                    if (x >= 0 && y >= 0 && x < static_cast<long>(width) && y < static_cast<long>(height)) {
                        const auto idx = (static_cast<std::size_t>(y) * width + static_cast<std::size_t>(x)) * 3;
                        pixels[idx] = 255;
                        pixels[idx + 1] = 255;
                        pixels[idx + 2] = 255;
                    }
                }
            }
        }
    }
    auto out = open_output(path, std::ios::out | std::ios::binary);
    out << "P6\n" << width << " " << height << "\n255\n";
    out.write(reinterpret_cast<const char*>(pixels.data()), static_cast<std::streamsize>(pixels.size()));
}

void write_rd_csv(const std::filesystem::path& path, const RealMatrix& rd_map, const RangeDopplerAxes& axes) {
    auto out = open_output(path);
    out << std::setprecision(5) << "range_m";
    for (std::size_t d = 0; d < rd_map.cols(); ++d) {
        out << "," << axes.velocity(static_cast<double>(d));
    }
    out << "\n" << std::fixed << std::setprecision(2);
    for (std::size_t r = 0; r < rd_map.rows(); ++r) {
        out << axes.range(static_cast<double>(r));
        for (std::size_t d = 0; d < rd_map.cols(); ++d) {
            out << "," << linear_to_db(rd_map(r, d), -100.0);
        }
        out << "\n";
    }
}

void write_detections_csv(const std::filesystem::path& path, const std::vector<Detection>& detections) {
    auto out = open_output(path);
    out << "cpi,range_m,velocity_mps,snr_db,peak_power_db,range_bin,doppler_bin,cells\n";
    for (const auto& d : detections) {
        out << d.cpi << "," << d.range_m << "," << d.velocity_mps << "," << d.snr_db << ","
            << linear_to_db(d.peak_power) << "," << d.range_bin << "," << d.doppler_bin << "," << d.cells << "\n";
    }
}

void write_tracks_csv(const std::filesystem::path& path, const std::vector<Track>& tracks) {
    auto out = open_output(path);
    out << "id,status,range_m,velocity_mps,snr_db,hits,first_cpi,last_cpi\n";
    for (const auto& t : tracks) {
        out << t.id << "," << to_string(t.status) << "," << t.range_m << "," << t.velocity_mps << "," << t.snr_db
            << "," << t.hits << "," << t.first_cpi << "," << t.last_cpi << "\n";
    }
}

void write_timing_csv(const std::filesystem::path& path, const StageTimer& timer, std::size_t cpis) {
    auto out = open_output(path);
    out << "stage,total_ms,per_cpi_ms,share_percent\n";
    const double total = timer.total_ms();
    for (const auto& e : timer.entries()) {
        out << e.name << "," << e.total_ms << "," << e.total_ms / static_cast<double>(std::max<std::size_t>(cpis, 1))
            << "," << (total > 0.0 ? 100.0 * e.total_ms / total : 0.0) << "\n";
    }
}

void write_pd_csv(const std::filesystem::path& path, const MonteCarloResult& result) {
    auto out = open_output(path);
    out << "snr_db,trials,detections,pd,pd_ideal,pd_cfar_analytic\n";
    for (const auto& p : result.points) {
        out << p.snr_db << "," << p.trials << "," << p.detections << "," << p.pd << "," << p.pd_ideal << ",";
        if (!std::isnan(p.pd_cfar)) {
            out << p.pd_cfar;
        }
        out << "\n";
    }
}

std::string base64_encode(const std::vector<std::uint8_t>& bytes) {
    static constexpr std::string_view kAlphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    out.reserve((bytes.size() + 2) / 3 * 4);
    std::size_t i = 0;
    for (; i + 2 < bytes.size(); i += 3) {
        const std::uint32_t v = (std::uint32_t{bytes[i]} << 16U) | (std::uint32_t{bytes[i + 1]} << 8U) | bytes[i + 2];
        out += kAlphabet[(v >> 18U) & 63U];
        out += kAlphabet[(v >> 12U) & 63U];
        out += kAlphabet[(v >> 6U) & 63U];
        out += kAlphabet[v & 63U];
    }
    const std::size_t rest = bytes.size() - i;
    if (rest > 0) {
        std::uint32_t v = std::uint32_t{bytes[i]} << 16U;
        if (rest == 2) {
            v |= std::uint32_t{bytes[i + 1]} << 8U;
        }
        out += kAlphabet[(v >> 18U) & 63U];
        out += kAlphabet[(v >> 12U) & 63U];
        out += rest == 2 ? kAlphabet[(v >> 6U) & 63U] : '=';
        out += '=';
    }
    return out;
}

void write_html_viewer(const std::filesystem::path& path, const Scenario& scenario, const SimulationResult& result) {
    if (result.cpis.empty()) {
        throw std::invalid_argument("no CPIs to display");
    }
    const auto& first = result.cpis.front().rd_map;
    double top = -300.0;
    for (const auto& c : result.cpis) {
        top = std::max(top, map_peak_db(c.rd_map));
    }
    const double dynamic_range = 70.0;
    const double bottom = top - dynamic_range;
    // Keep the embedded data around 256k cells per frame: max-pool range bins for large maps.
    constexpr std::size_t kMaxCellsPerFrame = 262144;
    const std::size_t pool = std::max<std::size_t>(1, (first.size() + kMaxCellsPerFrame - 1) / kMaxCellsPerFrame);
    const std::size_t rows = (first.rows() + pool - 1) / pool;
    const std::size_t cols = first.cols();

    std::ostringstream data;
    data << std::setprecision(8);
    const auto& r = scenario.radar;
    data << "{\"name\":\"" << json_escape(scenario.name) << "\",\"rangeBins\":" << rows
         << ",\"dopplerBins\":" << cols << ",\"rangeBinM\":" << result.axes.range_bin_m * static_cast<double>(pool)
         << ",\"velocityBinMps\":" << result.axes.velocity_bin_mps << ",\"dbMin\":" << bottom << ",\"dbMax\":" << top
         << ",\"params\":{\"carrierGHz\":" << r.carrier_hz / 1e9 << ",\"bandwidthMHz\":" << r.bandwidth_hz / 1e6
         << ",\"pulseWidthUs\":" << r.pulse_width_s * 1e6 << ",\"prfKHz\":" << r.prf_hz / 1e3
         << ",\"pulses\":" << r.pulses << ",\"rangeResM\":" << r.range_resolution_m()
         << ",\"vmaxMps\":" << r.max_unambiguous_velocity_mps() << ",\"cfar\":\""
         << to_string(scenario.processing.cfar) << "\",\"pfa\":" << scenario.processing.cfar_pfa << ",\"mti\":\""
         << to_string(scenario.processing.mti) << "\",\"rangeWindow\":\""
         << to_string(scenario.processing.range_window.type) << "\",\"dopplerWindow\":\""
         << to_string(scenario.processing.doppler_window.type) << "\"},\"frames\":[";
    std::vector<std::uint8_t> bytes(rows * cols);
    for (std::size_t c = 0; c < result.cpis.size(); ++c) {
        const auto& map = result.cpis[c].rd_map;
        for (std::size_t row = 0; row < rows; ++row) {
            for (std::size_t d = 0; d < cols; ++d) {
                double p = 0.0;
                for (std::size_t k = row * pool; k < std::min(map.rows(), (row + 1) * pool); ++k) {
                    p = std::max(p, map(k, d));
                }
                const double x = (linear_to_db(p) - bottom) / dynamic_range;
                bytes[row * cols + d] = static_cast<std::uint8_t>(std::lround(std::clamp(x, 0.0, 1.0) * 255.0));
            }
        }
        data << (c == 0 ? "" : ",") << "\"" << base64_encode(bytes) << "\"";
    }
    data << "],\"detections\":[";
    for (std::size_t i = 0; i < result.detections.size(); ++i) {
        const auto& d = result.detections[i];
        data << (i == 0 ? "" : ",") << "[" << d.cpi << "," << d.range_bin / static_cast<double>(pool) << ","
             << d.doppler_bin << "," << d.range_m
             << "," << d.velocity_mps << "," << d.snr_db << "]";
    }
    data << "],\"tracks\":[";
    bool first_track = true;
    for (const auto& t : result.tracks) {
        if (t.status != TrackStatus::Confirmed) {
            continue;
        }
        data << (first_track ? "" : ",") << "[" << t.id << "," << t.range_m << "," << t.velocity_mps << ","
             << t.snr_db << "," << t.hits << "," << t.first_cpi << "," << t.last_cpi << "]";
        first_track = false;
    }
    data << "]}";

    std::string html(kViewerTemplate);
    const std::string marker = "__RADAR_DATA__";
    const auto pos = html.find(marker);
    if (pos == std::string::npos) {
        throw std::logic_error("viewer template has no data marker");
    }
    html.replace(pos, marker.size(), data.str());
    auto out = open_output(path, std::ios::out | std::ios::binary);
    out << html;
}

void print_timing(std::ostream& os, const StageTimer& timer, std::size_t cpis, std::size_t raw_samples) {
    const double total = timer.total_ms();
    const auto n = static_cast<double>(std::max<std::size_t>(cpis, 1));
    os << "\nTiming per stage (" << cpis << " CPI" << (cpis == 1 ? "" : "s") << ")\n";
    os << "  " << std::left << std::setw(20) << "stage" << std::right << std::setw(12) << "total ms" << std::setw(12)
       << "ms / CPI" << std::setw(9) << "share" << "\n";
    os << std::fixed << std::setprecision(2);
    for (const auto& e : timer.entries()) {
        os << "  " << std::left << std::setw(20) << e.name << std::right << std::setw(12) << e.total_ms << std::setw(12)
           << e.total_ms / n << std::setw(8) << (total > 0.0 ? 100.0 * e.total_ms / total : 0.0) << "%\n";
    }
    os << "  " << std::left << std::setw(20) << "total" << std::right << std::setw(12) << total << std::setw(12)
       << total / n << "\n";
    double chain_ms = 0.0;
    for (const auto& e : timer.entries()) {
        if (e.name != "signal generation" && e.name != "output") {
            chain_ms += e.total_ms;
        }
    }
    if (chain_ms > 0.0 && raw_samples > 0) {
        os << "  receive-chain throughput: " << std::setprecision(1) << static_cast<double>(raw_samples) / (chain_ms * 1e3)
           << " Msamples/s (excludes signal generation and file output)\n";
    }
    os.unsetf(std::ios::floatfield);
    os << std::left;
}

void print_detections(std::ostream& os, const std::vector<Detection>& detections, std::size_t cpi) {
    os << "CPI " << cpi << ": " << detections.size() << " detection" << (detections.size() == 1 ? "" : "s") << "\n";
    os << std::fixed;
    for (const auto& d : detections) {
        os << "    range " << std::setw(9) << std::setprecision(1) << d.range_m << " m   velocity " << std::setw(7)
           << std::setprecision(2) << d.velocity_mps << " m/s   SNR " << std::setw(5) << std::setprecision(1)
           << d.snr_db << " dB   cells " << d.cells << "\n";
    }
    os.unsetf(std::ios::floatfield);
}

void print_tracks(std::ostream& os, const std::vector<Track>& tracks) {
    std::size_t confirmed = 0;
    for (const auto& t : tracks) {
        confirmed += t.status == TrackStatus::Confirmed ? 1 : 0;
    }
    os << "\nConfirmed tracks: " << confirmed << "\n";
    os << "  " << std::left << std::setw(5) << "id" << std::right << std::setw(12) << "range m" << std::setw(12)
       << "vel m/s" << std::setw(10) << "SNR dB" << std::setw(7) << "hits" << std::setw(11) << "CPIs" << "\n";
    os << std::fixed;
    for (const auto& t : tracks) {
        if (t.status != TrackStatus::Confirmed) {
            continue;
        }
        std::ostringstream span;
        span << t.first_cpi << "-" << t.last_cpi;
        os << "  " << std::left << std::setw(5) << t.id << std::right << std::setw(12) << std::setprecision(1)
           << t.range_m << std::setw(12) << std::setprecision(2) << t.velocity_mps << std::setw(10)
           << std::setprecision(1) << t.snr_db << std::setw(7) << t.hits << std::setw(11) << span.str() << "\n";
    }
    os.unsetf(std::ios::floatfield);
    os << std::left;
}

void print_monte_carlo(std::ostream& os, const MonteCarloResult& result, const Scenario& scenario) {
    const auto& mc = scenario.monte_carlo;
    os << "Monte Carlo: " << to_string(mc.swerling) << " target at " << std::fixed << std::setprecision(1)
       << result.target_range_m << " m, " << std::setprecision(2) << result.target_velocity_mps << " m/s, "
       << mc.trials << " trials per SNR\n";
    os << to_string(scenario.processing.cfar) << ": N = " << result.training_cells << ", alpha = " << std::setprecision(3)
       << result.cfar_alpha << ", coherent processing gain = " << std::setprecision(1) << result.processing_gain_db
       << " dB\n\n";
    os << "  " << std::right << std::setw(8) << "SNR dB" << std::setw(10) << "Pd sim" << std::setw(12) << "Pd ideal"
       << std::setw(12) << "Pd CFAR" << "\n";
    for (const auto& p : result.points) {
        os << "  " << std::setw(8) << std::setprecision(1) << p.snr_db << std::setw(10) << std::setprecision(3) << p.pd
           << std::setw(12) << p.pd_ideal;
        if (std::isnan(p.pd_cfar)) {
            os << std::setw(12) << "-";
        } else {
            os << std::setw(12) << p.pd_cfar;
        }
        os << "\n";
    }
    os << std::scientific << std::setprecision(2) << "\nPfa design " << result.pfa_design << ", empirical "
       << result.pfa_empirical << " (" << result.false_alarms << " / " << result.cells_tested << " cells)\n";
    os.unsetf(std::ios::floatfield);
    os << std::left;
}

}  // namespace radar
