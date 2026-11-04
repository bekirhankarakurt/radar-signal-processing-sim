#include "radar/complex_math.hpp"
#include "radar/doppler.hpp"
#include "radar/mti.hpp"
#include "test_framework.hpp"

#include <random>

using namespace radar;

TEST(tone_lands_in_expected_velocity_bin) {
    const std::size_t pulses = 64;
    const double prf = 10e3;
    const double lambda = 0.03;
    const DopplerProcessor dp(pulses, 64, {WindowType::Hamming, 4, -35.0});
    const RangeDopplerAxes axes{15.0, prf * lambda / (2.0 * 64.0), 64};
    for (const int k : {-20, -1, 0, 5, 31}) {
        const double fd = static_cast<double>(k) * prf / 64.0;
        ComplexMatrix slow(1, pulses);
        for (std::size_t m = 0; m < pulses; ++m) {
            slow(0, m) = cis(2.0 * kPi * fd * static_cast<double>(m) / prf);
        }
        const RealMatrix map = dp.process(slow);
        const auto peak = static_cast<std::size_t>(std::max_element(map.row(0), map.row(0) + 64) - map.row(0));
        CHECK(peak == static_cast<std::size_t>(32 + k));
        CHECK_NEAR(axes.velocity(static_cast<double>(peak)), fd * lambda / 2.0, 1e-9);
        // Coherent gain of the windowed FFT for an on-bin tone.
        CHECK_NEAR(map(0, peak), dp.snr_gain(), 1e-6 * dp.snr_gain());
    }
}

TEST(doppler_processing_preserves_noise_power) {
    const DopplerProcessor dp(48, 64, {WindowType::Taylor, 4, -40.0});
    std::mt19937_64 rng(9);
    std::normal_distribution<double> g(0.0, std::sqrt(0.5));
    ComplexMatrix noise(2000, 48);
    for (auto& z : noise.data()) {
        z = {g(rng), g(rng)};
    }
    const RealMatrix map = dp.process(noise);
    double mean = 0.0;
    for (const double v : map.data()) {
        mean += v;
    }
    CHECK_NEAR(mean / static_cast<double>(map.size()), 1.0, 0.03);
    CHECK_THROWS(dp.process(ComplexMatrix(1, 47)));
    CHECK_THROWS(DopplerProcessor(65, 64, {}));
}

TEST(mti_cancels_stationary_clutter_and_keeps_unit_noise_gain) {
    for (const auto type : {MtiType::TwoPulse, MtiType::ThreePulse}) {
        double energy = 0.0;
        for (const double c : mti_coefficients(type)) {
            energy += c * c;
        }
        CHECK_NEAR(energy, 1.0, 1e-12);
        ComplexMatrix dc(3, 16, Complex{5.0, -2.0});
        const ComplexMatrix out = apply_mti(dc, type);
        CHECK(out.cols() == 16 - (type == MtiType::TwoPulse ? 1U : 2U));
        for (const auto& z : out.data()) {
            CHECK(std::abs(z) < 1e-12);
        }
        CHECK_NEAR(mti_power_response(type, 0.0, 1e3), 0.0, 1e-12);
    }
    // Maximum response at PRF / 2: |1 - (-1)|^2 / 2 = 2 and |1 + 2 + 1|^2 / 6 = 8/3.
    CHECK_NEAR(mti_power_response(MtiType::TwoPulse, 500.0, 1e3), 2.0, 1e-12);
    CHECK_NEAR(mti_power_response(MtiType::ThreePulse, 500.0, 1e3), 8.0 / 3.0, 1e-12);
    CHECK_NEAR(mti_power_response(MtiType::None, 123.0, 1e3), 1.0, 0.0);
    CHECK(apply_mti(ComplexMatrix(2, 8), MtiType::None).cols() == 8);
    CHECK_THROWS(apply_mti(ComplexMatrix(2, 2), MtiType::ThreePulse));
}

TEST(mti_response_matches_measured_gain) {
    const double prf = 1e3;
    const double fd = 137.0;
    ComplexMatrix tone(1, 64);
    for (std::size_t m = 0; m < 64; ++m) {
        tone(0, m) = cis(2.0 * kPi * fd * static_cast<double>(m) / prf);
    }
    const ComplexMatrix out = apply_mti(tone, MtiType::ThreePulse);
    CHECK_NEAR(power(out(0, 10)), mti_power_response(MtiType::ThreePulse, fd, prf), 1e-9);
}

TEST_MAIN()
