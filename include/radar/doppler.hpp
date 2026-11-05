#pragma once

#include "radar/fft.hpp"
#include "radar/types.hpp"
#include "radar/window.hpp"

#include <cmath>

namespace radar {

/// Slow-time FFT processor. Produces a range-Doppler power map with rows = range
/// bins and cols = Doppler bins, zero Doppler in the centre column.
class DopplerProcessor {
public:
    DopplerProcessor(std::size_t pulses, std::size_t fft_size, const WindowSpec& window);

    /// `range_major` holds one row of slow-time samples per range bin.
    [[nodiscard]] RealMatrix process(const ComplexMatrix& range_major) const;

    [[nodiscard]] std::size_t fft_size() const noexcept { return plan_.size(); }
    [[nodiscard]] double snr_gain() const noexcept { return snr_gain_; }

private:
    std::size_t pulses_;
    FftPlan plan_;
    RealVector window_;  // normalised to unit noise gain
    double snr_gain_;
};

/// Physical axes of a range-Doppler map.
struct RangeDopplerAxes {
    double range_bin_m = 0.0;
    double velocity_bin_mps = 0.0;
    std::size_t doppler_bins = 0;

    [[nodiscard]] double range(double bin) const noexcept { return bin * range_bin_m; }
    [[nodiscard]] double velocity(double bin) const noexcept {
        return (bin - std::floor(static_cast<double>(doppler_bins) / 2.0)) * velocity_bin_mps;
    }
};

}  // namespace radar
