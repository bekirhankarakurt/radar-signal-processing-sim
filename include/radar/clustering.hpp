#pragma once

#include "radar/cfar.hpp"
#include "radar/doppler.hpp"
#include "radar/types.hpp"

#include <vector>

namespace radar {

/// A cluster of adjacent CFAR hits reduced to a single power-weighted centroid.
struct Detection {
    std::size_t cpi = 0;
    double range_bin = 0.0;
    double doppler_bin = 0.0;
    double range_m = 0.0;
    double velocity_mps = 0.0;
    double peak_power = 0.0;
    double snr_db = 0.0;  ///< Peak power over the local CFAR noise estimate.
    std::size_t cells = 0;
};

/// Groups 8-connected detections (Doppler wraps around) and centroids each group.
[[nodiscard]] std::vector<Detection> cluster_detections(const RealMatrix& power_map, const CfarResult& cfar,
                                                        const RangeDopplerAxes& axes, std::size_t cpi);

}  // namespace radar
