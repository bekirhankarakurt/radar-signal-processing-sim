#pragma once

#include <cstddef>

namespace radar {

/// Exponentially scaled modified Bessel function I0(x) * exp(-|x|).
[[nodiscard]] double bessel_i0e(double x) noexcept;

/// Generalised Marcum Q function of order 1, Q1(a, b).
[[nodiscard]] double marcum_q1(double a, double b) noexcept;

/// Pd of a non-fluctuating (Swerling 0) target, square-law detector, known noise level.
[[nodiscard]] double pd_swerling0_ideal(double snr_linear, double pfa) noexcept;

/// Pd of a Swerling 1 target, square-law detector, known noise level: Pfa^(1/(1+SNR)).
[[nodiscard]] double pd_swerling1_ideal(double snr_linear, double pfa) noexcept;

/// Pd of a Swerling 1 target behind a CA-CFAR with N reference cells and scale alpha.
[[nodiscard]] double pd_swerling1_ca_cfar(double snr_linear, std::size_t training_cells, double alpha) noexcept;

/// Pd of a Swerling 1 target behind an OS-CFAR (rank k of N) with scale alpha.
[[nodiscard]] double pd_swerling1_os_cfar(double snr_linear, std::size_t training_cells, std::size_t rank,
                                          double alpha) noexcept;

}  // namespace radar
