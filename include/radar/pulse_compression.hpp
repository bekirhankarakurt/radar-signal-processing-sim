#pragma once

#include "radar/fft.hpp"
#include "radar/types.hpp"
#include "radar/window.hpp"

namespace radar {

/// Frequency-domain matched filter for one transmitted waveform.
///
/// The reference is the windowed chirp normalised to unit energy, so white
/// noise keeps unit variance after compression and the peak SNR gain equals the
/// windowed time-bandwidth product.
class PulseCompressor {
public:
    PulseCompressor(const ComplexVector& waveform, const WindowSpec& window, std::size_t input_samples,
                    std::size_t output_samples);

    /// Compresses every row (pulse) of `cube`. Output row n is the range bin with delay n / fs.
    [[nodiscard]] ComplexMatrix compress(const ComplexMatrix& cube) const;

    [[nodiscard]] std::size_t fft_size() const noexcept { return plan_.size(); }
    [[nodiscard]] double snr_gain() const noexcept { return snr_gain_; }

private:
    std::size_t input_samples_;
    std::size_t output_samples_;
    FftPlan plan_;
    ComplexVector reference_spectrum_;  // conj(FFT(windowed reference))
    double snr_gain_;
};

}  // namespace radar
