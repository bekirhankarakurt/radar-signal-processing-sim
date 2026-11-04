#pragma once

#include "radar/config.hpp"
#include "radar/types.hpp"

namespace radar {

/// Binomial MTI canceller coefficients normalised to unit noise gain.
[[nodiscard]] RealVector mti_coefficients(MtiType type);

/// Applies the canceller along slow time to a range-major matrix (rows = range bins,
/// cols = pulses). The output has (pulses - order) columns.
[[nodiscard]] ComplexMatrix apply_mti(const ComplexMatrix& range_major, MtiType type);

/// Power response |H(f)|^2 of the canceller at Doppler frequency f.
[[nodiscard]] double mti_power_response(MtiType type, double doppler_hz, double prf_hz) noexcept;

}  // namespace radar
