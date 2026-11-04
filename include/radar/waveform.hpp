#pragma once

#include "radar/types.hpp"

namespace radar {

/// Baseband linear-FM chirp s(t) = exp(j pi k t^2), t in [-T/2, T/2), k = B / T.
[[nodiscard]] ComplexVector make_lfm_chirp(double bandwidth_hz, double pulse_width_s, double sample_rate_hz);

/// Value of the same chirp at an arbitrary time t relative to the pulse start (0 outside the pulse).
[[nodiscard]] Complex lfm_sample(double t, double bandwidth_hz, double pulse_width_s) noexcept;

}  // namespace radar
