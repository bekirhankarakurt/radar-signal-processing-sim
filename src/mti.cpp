#include "radar/mti.hpp"

#include "radar/complex_math.hpp"

#include <cmath>
#include <stdexcept>

namespace radar {

RealVector mti_coefficients(MtiType type) {
    switch (type) {
        case MtiType::None:
            return {1.0};
        case MtiType::TwoPulse:
            return {1.0 / std::sqrt(2.0), -1.0 / std::sqrt(2.0)};
        case MtiType::ThreePulse:
            return {1.0 / std::sqrt(6.0), -2.0 / std::sqrt(6.0), 1.0 / std::sqrt(6.0)};
    }
    throw std::invalid_argument("unknown MTI type");
}

ComplexMatrix apply_mti(const ComplexMatrix& range_major, MtiType type) {
    const RealVector h = mti_coefficients(type);
    const std::size_t order = h.size() - 1;
    if (range_major.cols() <= order) {
        throw std::invalid_argument("not enough pulses for the MTI canceller");
    }
    const std::size_t out_cols = range_major.cols() - order;
    ComplexMatrix out(range_major.rows(), out_cols);
    const auto rows = static_cast<std::ptrdiff_t>(range_major.rows());
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t ri = 0; ri < rows; ++ri) {
        const auto r = static_cast<std::size_t>(ri);
        const Complex* in = range_major.row(r);
        Complex* o = out.row(r);
        for (std::size_t m = 0; m < out_cols; ++m) {
            Complex acc{0.0, 0.0};
            for (std::size_t k = 0; k <= order; ++k) {
                acc += h[k] * in[m + order - k];
            }
            o[m] = acc;
        }
    }
    return out;
}

double mti_power_response(MtiType type, double doppler_hz, double prf_hz) noexcept {
    const double s = std::sin(kPi * doppler_hz / prf_hz);
    switch (type) {
        case MtiType::None:
            return 1.0;
        case MtiType::TwoPulse:
            return 2.0 * s * s;  // |1 - z^-1|^2 / 2
        case MtiType::ThreePulse:
            return 16.0 * s * s * s * s / 6.0;  // |1 - z^-1|^4 / 6
    }
    return 1.0;
}

}  // namespace radar
