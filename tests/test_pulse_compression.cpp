#include "radar/complex_math.hpp"
#include "radar/pulse_compression.hpp"
#include "radar/waveform.hpp"
#include "test_framework.hpp"

#include <random>

using namespace radar;

namespace {

constexpr double kBandwidth = 5e6;
constexpr double kPulse = 20e-6;
constexpr double kFs = 40e6;  // 8x oversampled so the main lobe is well resolved

ComplexMatrix single_echo(std::size_t samples, double delay_samples, Complex amplitude) {
    ComplexMatrix cube(1, samples);
    for (std::size_t n = 0; n < samples; ++n) {
        const double t = (static_cast<double>(n) - delay_samples) / kFs;
        cube(0, n) = amplitude * lfm_sample(t, kBandwidth, kPulse);
    }
    return cube;
}

RealVector power_profile(const ComplexMatrix& m) {
    RealVector p(m.cols());
    for (std::size_t n = 0; n < m.cols(); ++n) {
        p[n] = power(m(0, n));
    }
    return p;
}

/// -3 dB width of the main lobe around `peak`, in samples, with linear interpolation.
double half_power_width(const RealVector& p, std::size_t peak) {
    const double half = 0.5 * p[peak];
    std::size_t lo = peak;
    while (lo > 0 && p[lo - 1] > half) {
        --lo;
    }
    std::size_t hi = peak;
    while (hi + 1 < p.size() && p[hi + 1] > half) {
        ++hi;
    }
    const double left = static_cast<double>(lo) - (p[lo] - half) / (p[lo] - p[lo - 1]);
    const double right = static_cast<double>(hi) + (p[hi] - half) / (p[hi] - p[hi + 1]);
    return right - left;
}

double peak_sidelobe_db(const RealVector& p, std::size_t peak) {
    std::size_t lo = peak;
    while (lo > 0 && p[lo - 1] < p[lo]) {
        --lo;
    }
    std::size_t hi = peak;
    while (hi + 1 < p.size() && p[hi + 1] < p[hi]) {
        ++hi;
    }
    double side = 0.0;
    for (std::size_t n = 0; n < p.size(); ++n) {
        if (n < lo || n > hi) {
            side = std::max(side, p[n]);
        }
    }
    return linear_to_db(side / p[peak]);
}

}  // namespace

TEST(peak_position_and_gain_match_theory) {
    const ComplexVector chirp = make_lfm_chirp(kBandwidth, kPulse, kFs);
    const std::size_t samples = 2400;
    const PulseCompressor pc(chirp, {}, samples, samples - chirp.size());
    const ComplexMatrix out = pc.compress(single_echo(samples, 500.0, {1.0, 0.0}));
    const RealVector p = power_profile(out);
    const auto peak = static_cast<std::size_t>(std::max_element(p.begin(), p.end()) - p.begin());
    CHECK(peak == 500);
    // Unit-energy matched filter: peak SNR gain = number of pulse samples (time-bandwidth x oversampling).
    CHECK_NEAR(p[peak], static_cast<double>(chirp.size()), 1e-6 * static_cast<double>(chirp.size()));
    CHECK_NEAR(pc.snr_gain(), static_cast<double>(chirp.size()), 1e-6);
}

TEST(range_resolution_matches_c_over_2b) {
    const ComplexVector chirp = make_lfm_chirp(kBandwidth, kPulse, kFs);
    const std::size_t samples = 2400;
    const PulseCompressor rect(chirp, {}, samples, samples - chirp.size());
    const RealVector p = power_profile(rect.compress(single_echo(samples, 700.0, {1.0, 0.0})));
    // A sinc main lobe has a -3 dB width of 0.886 / B seconds, i.e. 0.886 c / (2B) metres.
    const double width_m = half_power_width(p, 700) * kSpeedOfLight / (2.0 * kFs);
    CHECK_NEAR(width_m, 0.886 * kSpeedOfLight / (2.0 * kBandwidth), 0.05 * 0.886 * kSpeedOfLight / (2.0 * kBandwidth));

    const PulseCompressor taylor(chirp, {WindowType::Taylor, 4, -35.0}, samples, samples - chirp.size());
    const RealVector pt = power_profile(taylor.compress(single_echo(samples, 700.0, {1.0, 0.0})));
    CHECK(half_power_width(pt, 700) > half_power_width(p, 700) * 1.15);  // tapering broadens the main lobe
}

TEST(sidelobes_follow_the_window) {
    const ComplexVector chirp = make_lfm_chirp(kBandwidth, kPulse, kFs);  // TB = 100
    const std::size_t samples = 2400;
    const PulseCompressor rect(chirp, {}, samples, samples - chirp.size());
    const PulseCompressor taylor(chirp, {WindowType::Taylor, 4, -35.0}, samples, samples - chirp.size());
    const auto echo = single_echo(samples, 900.0, {1.0, 0.0});
    CHECK_NEAR(peak_sidelobe_db(power_profile(rect.compress(echo)), 900), -13.3, 0.8);
    CHECK(peak_sidelobe_db(power_profile(taylor.compress(echo)), 900) < -30.0);
}

TEST(white_noise_keeps_unit_variance) {
    const ComplexVector chirp = make_lfm_chirp(kBandwidth, kPulse, kFs);
    const std::size_t samples = 4096;
    const PulseCompressor pc(chirp, {WindowType::Hamming, 4, -35.0}, samples, samples - chirp.size());
    std::mt19937_64 rng(5);
    std::normal_distribution<double> g(0.0, std::sqrt(0.5));
    ComplexMatrix noise(8, samples);
    for (auto& z : noise.data()) {
        z = {g(rng), g(rng)};
    }
    const ComplexMatrix out = pc.compress(noise);
    double mean_power = 0.0;
    for (const auto& z : out.data()) {
        mean_power += power(z);
    }
    mean_power /= static_cast<double>(out.size());
    CHECK_NEAR(mean_power, 1.0, 0.05);
}

TEST(two_targets_two_resolution_cells_apart_are_resolved) {
    const ComplexVector chirp = make_lfm_chirp(kBandwidth, kPulse, kFs);
    const std::size_t samples = 2400;
    const PulseCompressor pc(chirp, {WindowType::Hamming, 4, -35.0}, samples, samples - chirp.size());
    ComplexMatrix cube = single_echo(samples, 600.0, {1.0, 0.0});
    const ComplexMatrix second = single_echo(samples, 616.0, {0.0, 1.0});  // 2 x c/(2B) = 16 samples
    for (std::size_t n = 0; n < samples; ++n) {
        cube(0, n) += second(0, n);
    }
    const RealVector p = power_profile(pc.compress(cube));
    CHECK(p[608] < 0.5 * std::min(p[600], p[616]));
}

TEST(rejects_mismatched_input) {
    const ComplexVector chirp = make_lfm_chirp(kBandwidth, kPulse, kFs);
    const PulseCompressor pc(chirp, {}, 1000, 900);
    CHECK_THROWS(pc.compress(ComplexMatrix(1, 999)));
    CHECK_THROWS(PulseCompressor(chirp, {}, 100, 200));
    CHECK_THROWS(make_lfm_chirp(0.0, kPulse, kFs));
}

TEST_MAIN()
