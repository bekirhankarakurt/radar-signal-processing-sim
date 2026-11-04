#pragma once

#include "radar/config.hpp"

namespace testing {

/// Small, fast scenario used by the chain-level tests: 2 MHz LFM, 4 MHz sampling,
/// 3 km window, 32 pulses, no clutter.
inline radar::Scenario small_scenario() {
    // 5 W peak power keeps a 1 m^2 target at 1-2 km around 35-50 dB in the map, below the
    // level where window sidelobes would cross the CFAR threshold.
    radar::Scenario s;
    s.name = "test";
    s.radar.carrier_hz = 10e9;
    s.radar.bandwidth_hz = 2e6;
    s.radar.pulse_width_s = 10e-6;
    s.radar.prf_hz = 20e3;
    s.radar.pulses = 32;
    s.radar.sample_rate_hz = 4e6;
    s.radar.max_range_m = 3000.0;
    s.radar.peak_power_w = 5.0;
    s.processing.cfar_pfa = 1e-5;
    s.processing.guard_range = 2;
    s.processing.train_range = 8;
    s.processing.guard_doppler = 1;
    s.processing.train_doppler = 3;
    s.seed = 11;
    s.validate();
    return s;
}

}  // namespace testing
