#include "radar/cfar.hpp"
#include "radar/complex_math.hpp"
#include "radar/statistics.hpp"
#include "test_framework.hpp"

using namespace radar;

namespace {

/// Independent reference: Q1(a, b) = exp(-(a^2 + b^2) / 2) sum_k (a / b)^k I_k(a b), valid for a < b.
double marcum_series(double a, double b) {
    const double x = a * b;
    double sum = 0.0;
    for (int k = 0; k < 200; ++k) {
        // I_k(x) by its power series.
        double term = std::pow(x / 2.0, k) / std::tgamma(k + 1.0);
        double ik = term;
        for (int m = 1; m < 200; ++m) {
            term *= (x * x / 4.0) / (m * (m + k));
            ik += term;
        }
        sum += std::pow(a / b, k) * ik;
    }
    return std::exp(-(a * a + b * b) / 2.0) * sum;
}

}  // namespace

TEST(bessel_i0e_reference_values) {
    CHECK_NEAR(bessel_i0e(0.0), 1.0, 1e-15);
    CHECK_NEAR(bessel_i0e(1.0), 1.2660658777520082 * std::exp(-1.0), 1e-12);
    CHECK_NEAR(bessel_i0e(30.0), 7.816722978239774e11 * std::exp(-30.0), 1e-9);
    CHECK_NEAR(bessel_i0e(-2.0), bessel_i0e(2.0), 0.0);
}

TEST(marcum_q_special_cases_and_series) {
    for (const double b : {0.5, 1.0, 3.0}) {
        CHECK_NEAR(marcum_q1(0.0, b), std::exp(-b * b / 2.0), 1e-7);
    }
    CHECK_NEAR(marcum_q1(2.0, 0.0), 1.0, 0.0);
    CHECK_NEAR(marcum_q1(1.0, 2.0), marcum_series(1.0, 2.0), 1e-7);
    CHECK_NEAR(marcum_q1(2.0, 4.0), marcum_series(2.0, 4.0), 1e-7);
    CHECK_NEAR(marcum_q1(3.0, 5.25), marcum_series(3.0, 5.25), 1e-7);
}

TEST(classic_detection_numbers) {
    // Single-pulse, Pfa = 1e-6, Pd = 0.9: about 13.2 dB for Swerling 0 and 21.1 dB for Swerling 1.
    CHECK_NEAR(pd_swerling0_ideal(db_to_linear(13.2), 1e-6), 0.9, 0.02);
    CHECK_NEAR(pd_swerling1_ideal(db_to_linear(21.14), 1e-6), 0.9, 0.005);
    CHECK_NEAR(pd_swerling0_ideal(0.0, 1e-4), 1e-4, 1e-6);
}

TEST(cfar_pd_reduces_to_pfa_and_approaches_ideal) {
    const std::size_t n = 32;
    const double alpha = ca_cfar_alpha(n, 1e-6);
    CHECK_NEAR(pd_swerling1_ca_cfar(0.0, n, alpha), 1e-6, 1e-10);
    const double big_alpha = ca_cfar_alpha(100000, 1e-6);
    CHECK_NEAR(pd_swerling1_ca_cfar(100.0, 100000, big_alpha), pd_swerling1_ideal(100.0, 1e-6), 1e-4);
    // CFAR loss: the finite-N detector needs more SNR for the same Pd.
    CHECK(pd_swerling1_ca_cfar(100.0, n, alpha) < pd_swerling1_ideal(100.0, 1e-6));

    const double os_alpha = os_cfar_alpha(48, 36, 1e-4);
    CHECK_NEAR(pd_swerling1_os_cfar(0.0, 48, 36, os_alpha), 1e-4, 1e-9);
    CHECK(pd_swerling1_os_cfar(100.0, 48, 36, os_alpha) > 0.85);
}

TEST_MAIN()
