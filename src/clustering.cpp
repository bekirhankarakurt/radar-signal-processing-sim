#include "radar/clustering.hpp"

#include "radar/complex_math.hpp"

#include <algorithm>
#include <utility>

namespace radar {

std::vector<Detection> cluster_detections(const RealMatrix& power_map, const CfarResult& cfar,
                                          const RangeDopplerAxes& axes, std::size_t cpi) {
    const std::size_t rows = power_map.rows();
    const std::size_t cols = power_map.cols();
    Matrix<std::uint8_t> visited(rows, cols, 0);
    std::vector<std::pair<std::size_t, std::size_t>> stack;
    std::vector<std::pair<std::size_t, std::size_t>> members;
    std::vector<Detection> out;

    for (std::size_t r0 = 0; r0 < rows; ++r0) {
        for (std::size_t d0 = 0; d0 < cols; ++d0) {
            if (cfar.detections(r0, d0) == 0 || visited(r0, d0) != 0) {
                continue;
            }
            // Flood fill the 8-connected component; Doppler wraps around.
            members.clear();
            stack.assign(1, {r0, d0});
            visited(r0, d0) = 1;
            while (!stack.empty()) {
                const auto [r, d] = stack.back();
                stack.pop_back();
                members.emplace_back(r, d);
                for (int dr = -1; dr <= 1; ++dr) {
                    if ((dr < 0 && r == 0) || (dr > 0 && r + 1 >= rows)) {
                        continue;
                    }
                    const std::size_t nr = dr < 0 ? r - 1 : (dr > 0 ? r + 1 : r);
                    for (int dd = -1; dd <= 1; ++dd) {
                        const std::size_t nd = (d + cols + static_cast<std::size_t>(dd + 1) - 1) % cols;
                        if (cfar.detections(nr, nd) != 0 && visited(nr, nd) == 0) {
                            visited(nr, nd) = 1;
                            stack.emplace_back(nr, nd);
                        }
                    }
                }
            }

            // Power-weighted centroid around the peak cell, unwrapping Doppler relative to it.
            const auto peak = *std::max_element(members.begin(), members.end(), [&](const auto& a, const auto& b) {
                return power_map(a.first, a.second) < power_map(b.first, b.second);
            });
            const double peak_power = power_map(peak.first, peak.second);
            double weight = 0.0;
            double range_acc = 0.0;
            double doppler_acc = 0.0;
            const auto half = static_cast<double>(cols / 2);
            for (const auto& [r, d] : members) {
                const double p = power_map(r, d);
                double delta = static_cast<double>(d) - static_cast<double>(peak.second);
                if (delta > half) {
                    delta -= static_cast<double>(cols);
                } else if (delta < -half) {
                    delta += static_cast<double>(cols);
                }
                weight += p;
                range_acc += p * static_cast<double>(r);
                doppler_acc += p * delta;
            }
            Detection det;
            det.cpi = cpi;
            det.range_bin = range_acc / weight;
            double doppler_bin = static_cast<double>(peak.second) + doppler_acc / weight;
            if (doppler_bin < 0.0) {
                doppler_bin += static_cast<double>(cols);
            } else if (doppler_bin >= static_cast<double>(cols)) {
                doppler_bin -= static_cast<double>(cols);
            }
            det.doppler_bin = doppler_bin;
            det.range_m = axes.range(det.range_bin);
            det.velocity_mps = axes.velocity(det.doppler_bin);
            det.peak_power = peak_power;
            const double noise = cfar.noise_estimate(peak.first, peak.second);
            det.snr_db = linear_to_db(noise > 0.0 ? peak_power / noise : 0.0);
            det.cells = members.size();
            out.push_back(det);
        }
    }
    std::sort(out.begin(), out.end(), [](const Detection& a, const Detection& b) { return a.range_m < b.range_m; });
    return out;
}

}  // namespace radar
