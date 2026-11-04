#include "radar/complex_math.hpp"
#include "radar/fft.hpp"
#include "test_framework.hpp"

#include <random>

using namespace radar;

namespace {

ComplexVector random_vector(std::size_t n, unsigned seed) {
    std::mt19937 rng(seed);
    std::normal_distribution<double> g(0.0, 1.0);
    ComplexVector v(n);
    for (auto& z : v) {
        z = {g(rng), g(rng)};
    }
    return v;
}

double max_error(const ComplexVector& a, const ComplexVector& b) {
    double e = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        e = std::max(e, std::abs(a[i] - b[i]));
    }
    return e;
}

}  // namespace

TEST(fft_matches_naive_dft) {
    for (const std::size_t n : {1U, 2U, 4U, 8U, 64U, 256U}) {
        const ComplexVector x = random_vector(n, static_cast<unsigned>(n));
        ComplexVector y = x;
        FftPlan(n).execute(y, FftDirection::Forward);
        CHECK(max_error(y, naive_dft(x, FftDirection::Forward)) < 1e-9 * static_cast<double>(n));
    }
}

TEST(inverse_fft_restores_input) {
    const ComplexVector x = random_vector(1024, 7);
    ComplexVector y = x;
    const FftPlan plan(1024);
    plan.execute(y, FftDirection::Forward);
    plan.execute(y, FftDirection::Inverse);
    CHECK(max_error(x, y) < 1e-12);
}

TEST(parseval_energy_is_preserved) {
    const ComplexVector x = random_vector(512, 3);
    ComplexVector y = x;
    FftPlan(512).execute(y, FftDirection::Forward);
    CHECK_NEAR(total_energy(y) / 512.0, total_energy(x), 1e-9 * total_energy(x));
}

TEST(complex_tone_lands_in_its_bin) {
    const std::size_t n = 128;
    const std::size_t k = 17;
    ComplexVector x(n);
    for (std::size_t t = 0; t < n; ++t) {
        x[t] = cis(2.0 * kPi * static_cast<double>(k * t) / static_cast<double>(n));
    }
    FftPlan(n).execute(x, FftDirection::Forward);
    CHECK_NEAR(std::abs(x[k]), static_cast<double>(n), 1e-9);
    for (std::size_t i = 0; i < n; ++i) {
        if (i != k) {
            CHECK(std::abs(x[i]) < 1e-9);
        }
    }
}

TEST(non_power_of_two_is_rejected) {
    CHECK_THROWS(FftPlan(0));
    CHECK_THROWS(FftPlan(12));
    ComplexVector wrong(8);
    CHECK_THROWS(FftPlan(16).execute(wrong, FftDirection::Forward));
}

TEST(fft_shift_centres_zero_frequency) {
    ComplexVector v{{0, 0}, {1, 0}, {2, 0}, {3, 0}};
    fft_shift(v.data(), v.size());
    CHECK_NEAR(v[0].real(), 2.0, 0.0);
    CHECK_NEAR(v[2].real(), 0.0, 0.0);
}

TEST(power_of_two_helpers) {
    CHECK(is_power_of_two(1024));
    CHECK(!is_power_of_two(1000));
    CHECK(next_power_of_two(1000) == 1024);
    CHECK(next_power_of_two(1) == 1);
}

TEST_MAIN()
