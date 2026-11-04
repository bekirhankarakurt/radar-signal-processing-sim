#include "radar/fft.hpp"

#include "radar/complex_math.hpp"

#include <stdexcept>
#include <utility>

namespace radar {

FftPlan::FftPlan(std::size_t size)
    : size_(size), bit_reverse_(size), twiddles_(size > 0 ? size - 1 : 0), inverse_twiddles_(twiddles_.size()) {
    if (!is_power_of_two(size)) {
        throw std::invalid_argument("FFT size must be a power of two, got " + std::to_string(size));
    }
    std::size_t bits = 0;
    while ((std::size_t{1} << bits) < size) {
        ++bits;
    }
    for (std::size_t i = 0; i < size; ++i) {
        std::size_t r = 0;
        for (std::size_t b = 0; b < bits; ++b) {
            if ((i >> b) & 1U) {
                r |= std::size_t{1} << (bits - 1 - b);
            }
        }
        bit_reverse_[i] = r;
    }
    // Twiddles stored stage by stage (offset half - 1 for butterfly span 2 * half) so the
    // inner loop reads them contiguously; the inverse set is kept separately to avoid a branch.
    for (std::size_t half = 1; half < size; half <<= 1U) {
        for (std::size_t k = 0; k < half; ++k) {
            const Complex w = cis(-kPi * static_cast<double>(k) / static_cast<double>(half));
            twiddles_[half - 1 + k] = w;
            inverse_twiddles_[half - 1 + k] = std::conj(w);
        }
    }
}

void FftPlan::execute(Complex* data, FftDirection direction) const noexcept {
    const std::size_t n = size_;
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j = bit_reverse_[i];
        if (i < j) {
            std::swap(data[i], data[j]);
        }
    }
    const bool inverse = direction == FftDirection::Inverse;
    const Complex* stage_twiddles = inverse ? inverse_twiddles_.data() : twiddles_.data();
    for (std::size_t half = 1; half < n; half <<= 1U) {
        const Complex* tw = stage_twiddles + (half - 1);
        for (std::size_t start = 0; start < n; start += 2 * half) {
            Complex* lo = data + start;
            Complex* hi = lo + half;
            for (std::size_t k = 0; k < half; ++k) {
                const Complex& w = tw[k];
                // Explicit complex multiply avoids the NaN/Inf handling of operator*.
                const double tr = hi[k].real() * w.real() - hi[k].imag() * w.imag();
                const double ti = hi[k].real() * w.imag() + hi[k].imag() * w.real();
                const Complex t{tr, ti};
                hi[k] = lo[k] - t;
                lo[k] += t;
            }
        }
    }
    if (inverse) {
        const double scale = 1.0 / static_cast<double>(n);
        for (std::size_t i = 0; i < n; ++i) {
            data[i] *= scale;
        }
    }
}

void FftPlan::execute(ComplexVector& data, FftDirection direction) const {
    if (data.size() != size_) {
        throw std::invalid_argument("FFT input length does not match the plan");
    }
    execute(data.data(), direction);
}

ComplexVector naive_dft(const ComplexVector& input, FftDirection direction) {
    const std::size_t n = input.size();
    const double sign = direction == FftDirection::Forward ? -1.0 : 1.0;
    ComplexVector out(n);
    for (std::size_t k = 0; k < n; ++k) {
        Complex acc{0.0, 0.0};
        for (std::size_t t = 0; t < n; ++t) {
            const auto kt = static_cast<double>((k * t) % n);
            acc += input[t] * cis(sign * 2.0 * kPi * kt / static_cast<double>(n));
        }
        out[k] = direction == FftDirection::Inverse ? acc / static_cast<double>(n) : acc;
    }
    return out;
}

void fft_shift(Complex* data, std::size_t size) noexcept {
    const std::size_t half = size / 2;
    for (std::size_t i = 0; i < half; ++i) {
        std::swap(data[i], data[i + half]);
    }
}

}  // namespace radar
