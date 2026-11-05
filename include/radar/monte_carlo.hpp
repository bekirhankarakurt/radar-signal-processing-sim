#pragma once

#include "radar/config.hpp"
#include "radar/timing.hpp"

#include <vector>

namespace radar {

struct PdPoint {
    double snr_db = 0.0;  ///< SNR in the range-Doppler cell after all processing gain.
    std::size_t trials = 0;
    std::size_t detections = 0;
    double pd = 0.0;
    double pd_ideal = 0.0;  ///< Fixed-threshold detector with perfectly known noise.
    double pd_cfar = 0.0;   ///< Closed form for the configured CFAR (Swerling 1 only, else NaN).
};

struct MonteCarloResult {
    std::vector<PdPoint> points;
    double pfa_design = 0.0;
    double pfa_empirical = 0.0;
    std::size_t false_alarms = 0;
    std::size_t cells_tested = 0;
    double cfar_alpha = 0.0;
    std::size_t training_cells = 0;
    double target_range_m = 0.0;
    double target_velocity_mps = 0.0;
    double processing_gain_db = 0.0;  ///< Coherent gain of compression + Doppler integration.
};

/// Pd versus SNR and empirical Pfa through the complete signal-processing chain.
/// Trials are independent and seeded deterministically, so results do not depend
/// on the thread count.
[[nodiscard]] MonteCarloResult run_monte_carlo(const Scenario& input, StageTimer& timer);

}  // namespace radar
