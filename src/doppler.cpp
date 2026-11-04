#include "radar/doppler.hpp"

#include "radar/complex_math.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace radar {

DopplerProcessor::DopplerProcessor(std::size_t pulses, std::size_t fft_size, const WindowSpec& window)
    : pulses_(pulses), plan_(fft_size), window_(make_window(window, pulses)), snr_gain_(window_snr_gain(window_)) {
    if (pulses > fft_size) {
        throw std::invalid_argument("Doppler FFT shorter than the number of pulses");
    }
    double energy = 0.0;
    for (const double w : window_) {
        energy += w * w;
    }
    const double norm = 1.0 / std::sqrt(energy);
    for (auto& w : window_) {
        w *= norm;
    }
}

RealMatrix DopplerProcessor::process(const ComplexMatrix& range_major) const {
    if (range_major.cols() != pulses_) {
        throw std::invalid_argument("Doppler processor pulse count mismatch");
    }
    const std::size_t nfft = plan_.size();
    RealMatrix map(range_major.rows(), nfft);
    const auto rows = static_cast<std::ptrdiff_t>(range_major.rows());
#pragma omp parallel
    {
        ComplexVector buffer(nfft);
#pragma omp for schedule(static)
        for (std::ptrdiff_t ri = 0; ri < rows; ++ri) {
            const auto r = static_cast<std::size_t>(ri);
            const Complex* in = range_major.row(r);
            for (std::size_t m = 0; m < pulses_; ++m) {
                buffer[m] = in[m] * window_[m];
            }
            std::fill(buffer.begin() + static_cast<std::ptrdiff_t>(pulses_), buffer.end(), Complex{0.0, 0.0});
            plan_.execute(buffer.data(), FftDirection::Forward);
            fft_shift(buffer.data(), nfft);
            double* out = map.row(r);
            for (std::size_t k = 0; k < nfft; ++k) {
                out[k] = power(buffer[k]);
            }
        }
    }
    return map;
}

}  // namespace radar
