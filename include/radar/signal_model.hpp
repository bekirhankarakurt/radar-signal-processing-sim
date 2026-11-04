#pragma once

#include "radar/config.hpp"
#include "radar/fft.hpp"
#include "radar/types.hpp"

#include <cstdint>

namespace radar {

/// Received signal power from the radar range equation,
/// P_r = P_t G^2 lambda^2 sigma / ((4 pi)^3 R^4 L).
[[nodiscard]] double received_power_w(const RadarParameters& radar, double range_m, double rcs_m2) noexcept;

/// Single-pulse, single-sample SNR (noise bandwidth = sampling rate) of a target.
[[nodiscard]] double single_sample_snr(const RadarParameters& radar, double range_m, double rcs_m2) noexcept;

/// A target as seen during one CPI: its range at the first pulse and its complex
/// amplitude in units of the per-sample noise standard deviation.
struct TargetEcho {
    double range_m = 0.0;
    double velocity_mps = 0.0;
    Complex amplitude{0.0, 0.0};
};

/// Generates the raw fast-time x slow-time data cube for one coherent processing interval.
///
/// Rows are pulses, columns are fast-time samples. Noise is complex white Gaussian
/// with unit variance, so every amplitude in the model is a voltage SNR.
/// Each pulse draws from its own random stream, which keeps the output identical
/// regardless of the number of threads.
class SignalModel {
public:
    explicit SignalModel(const Scenario& scenario);

    /// Target echoes for CPI `cpi` (targets move and Swerling-1 targets decorrelate between CPIs).
    [[nodiscard]] std::vector<TargetEcho> echoes_for_cpi(std::size_t cpi) const;

    [[nodiscard]] ComplexMatrix generate(std::size_t cpi) const;
    [[nodiscard]] ComplexMatrix generate(const std::vector<TargetEcho>& echoes, std::uint64_t stream,
                                         bool add_noise) const;

private:
    void add_clutter(ComplexMatrix& cube, std::uint64_t stream) const;

    const Scenario& scenario_;
    ComplexVector chirp_;
    FftPlan clutter_plan_;            // linear convolution of clutter reflectivity with the chirp
    ComplexVector chirp_spectrum_;
    RealMatrix clutter_coloring_;     // Cholesky factor of the slow-time clutter correlation matrix
};

/// Lower-triangular Cholesky factor L of a symmetric positive-definite matrix A = L L^T.
[[nodiscard]] RealMatrix cholesky(const RealMatrix& a);

/// Slow-time correlation of clutter with a Gaussian Doppler spectrum of standard deviation
/// sigma_v: rho(m) = exp(-2 (pi sigma_f m T)^2), sigma_f = 2 sigma_v / lambda.
[[nodiscard]] double clutter_correlation(double spectral_width_mps, double wavelength_m, double pri_s,
                                         std::size_t lag) noexcept;

/// Mixes a 64-bit seed with stream identifiers (splitmix64), giving independent RNG seeds.
[[nodiscard]] std::uint64_t mix_seed(std::uint64_t seed, std::uint64_t a, std::uint64_t b = 0) noexcept;

}  // namespace radar
