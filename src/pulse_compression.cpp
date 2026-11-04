#include "radar/pulse_compression.hpp"

#include "radar/complex_math.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace radar {

PulseCompressor::PulseCompressor(const ComplexVector& waveform, const WindowSpec& window, std::size_t input_samples,
                                 std::size_t output_samples)
    : input_samples_(input_samples),
      output_samples_(output_samples),
      plan_(next_power_of_two(input_samples + waveform.size() - 1)),
      reference_spectrum_(plan_.size()) {
    if (waveform.empty() || output_samples > input_samples) {
        throw std::invalid_argument("invalid pulse compressor dimensions");
    }
    const RealVector w = make_window(window, waveform.size());
    ComplexVector reference(plan_.size(), Complex{0.0, 0.0});
    for (std::size_t n = 0; n < waveform.size(); ++n) {
        reference[n] = waveform[n] * w[n];
    }
    const double norm = std::sqrt(total_energy(reference));
    for (auto& z : reference) {
        z /= norm;
    }
    // Peak gain for a unit-amplitude echo: |sum s* s w|^2 / sum |s w|^2.
    Complex peak{0.0, 0.0};
    for (std::size_t n = 0; n < waveform.size(); ++n) {
        peak += waveform[n] * std::conj(reference[n]);
    }
    snr_gain_ = power(peak);

    plan_.execute(reference, FftDirection::Forward);
    for (std::size_t k = 0; k < reference.size(); ++k) {
        reference_spectrum_[k] = std::conj(reference[k]);
    }
}

ComplexMatrix PulseCompressor::compress(const ComplexMatrix& cube) const {
    if (cube.cols() != input_samples_) {
        throw std::invalid_argument("pulse compressor input length mismatch");
    }
    ComplexMatrix out(cube.rows(), output_samples_);
    const std::size_t nfft = plan_.size();
    const auto pulses = static_cast<std::ptrdiff_t>(cube.rows());
#pragma omp parallel
    {
        ComplexVector buffer(nfft);
#pragma omp for schedule(static)
        for (std::ptrdiff_t pi = 0; pi < pulses; ++pi) {
            const auto m = static_cast<std::size_t>(pi);
            const Complex* in = cube.row(m);
            std::copy(in, in + input_samples_, buffer.begin());
            std::fill(buffer.begin() + static_cast<std::ptrdiff_t>(input_samples_), buffer.end(), Complex{0.0, 0.0});
            plan_.execute(buffer.data(), FftDirection::Forward);
            for (std::size_t k = 0; k < nfft; ++k) {
                const Complex& a = buffer[k];
                const Complex& b = reference_spectrum_[k];
                buffer[k] = {a.real() * b.real() - a.imag() * b.imag(), a.real() * b.imag() + a.imag() * b.real()};
            }
            plan_.execute(buffer.data(), FftDirection::Inverse);
            std::copy(buffer.begin(), buffer.begin() + static_cast<std::ptrdiff_t>(output_samples_), out.row(m));
        }
    }
    return out;
}

}  // namespace radar
