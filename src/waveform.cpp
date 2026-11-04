#include "radar/waveform.hpp"

#include "radar/complex_math.hpp"

#include <cmath>
#include <stdexcept>

namespace radar {

Complex lfm_sample(double t, double bandwidth_hz, double pulse_width_s) noexcept {
    if (t < 0.0 || t >= pulse_width_s) {
        return {0.0, 0.0};
    }
    const double k = bandwidth_hz / pulse_width_s;
    const double tc = t - 0.5 * pulse_width_s;
    return cis(kPi * k * tc * tc);
}

ComplexVector make_lfm_chirp(double bandwidth_hz, double pulse_width_s, double sample_rate_hz) {
    if (bandwidth_hz <= 0.0 || pulse_width_s <= 0.0 || sample_rate_hz <= 0.0) {
        throw std::invalid_argument("chirp parameters must be positive");
    }
    const auto samples = static_cast<std::size_t>(std::llround(pulse_width_s * sample_rate_hz));
    if (samples == 0) {
        throw std::invalid_argument("pulse shorter than one sample");
    }
    ComplexVector chirp(samples);
    for (std::size_t n = 0; n < samples; ++n) {
        chirp[n] = lfm_sample(static_cast<double>(n) / sample_rate_hz, bandwidth_hz, pulse_width_s);
    }
    return chirp;
}

}  // namespace radar
