#include "radar/statistics.hpp"

#include "radar/types.hpp"

#include <algorithm>
#include <cmath>

namespace radar {

double bessel_i0e(double x) noexcept {
    const double ax = std::abs(x);
    if (ax < 20.0) {
        // Power series sum (x^2/4)^k / (k!)^2, scaled by exp(-|x|).
        const double q = 0.25 * ax * ax;
        double term = 1.0;
        double sum = 1.0;
        for (int k = 1; k < 200; ++k) {
            term *= q / (static_cast<double>(k) * static_cast<double>(k));
            sum += term;
            if (term < sum * 1e-17) {
                break;
            }
        }
        return sum * std::exp(-ax);
    }
    // Asymptotic expansion for large arguments.
    double term = 1.0;
    double sum = 1.0;
    for (int k = 1; k < 12; ++k) {
        const double odd = 2.0 * static_cast<double>(k) - 1.0;
        term *= odd * odd / (8.0 * static_cast<double>(k) * ax);
        sum += term;
    }
    return sum / std::sqrt(2.0 * kPi * ax);
}

double marcum_q1(double a, double b) noexcept {
    if (b <= 0.0) {
        return 1.0;
    }
    // Q1(a, b) = integral_b^inf x exp(-(x^2 + a^2) / 2) I0(a x) dx, written with the scaled
    // Bessel function so the integrand is exp(-(x - a)^2 / 2) x I0e(a x). Simpson's rule.
    const double upper = std::max(a, b) + 40.0;
    const int steps = 20000;
    const double h = (upper - b) / steps;
    auto f = [a](double x) { return x * std::exp(-0.5 * (x - a) * (x - a)) * bessel_i0e(a * x); };
    double sum = f(b) + f(upper);
    for (int i = 1; i < steps; ++i) {
        sum += f(b + h * i) * ((i % 2 == 1) ? 4.0 : 2.0);
    }
    const double q = sum * h / 3.0;
    return std::min(1.0, std::max(0.0, q));
}

double pd_swerling0_ideal(double snr_linear, double pfa) noexcept {
    return marcum_q1(std::sqrt(2.0 * snr_linear), std::sqrt(-2.0 * std::log(pfa)));
}

double pd_swerling1_ideal(double snr_linear, double pfa) noexcept { return std::pow(pfa, 1.0 / (1.0 + snr_linear)); }

double pd_swerling1_ca_cfar(double snr_linear, std::size_t training_cells, double alpha) noexcept {
    const auto n = static_cast<double>(training_cells);
    return std::pow(1.0 + alpha / (n * (1.0 + snr_linear)), -n);
}

double pd_swerling1_os_cfar(double snr_linear, std::size_t training_cells, std::size_t rank, double alpha) noexcept {
    // A Swerling 1 cell is exponential with mean (1 + SNR), so Pd is the OS-CFAR Pfa
    // expression evaluated with the scale alpha / (1 + SNR).
    const double scaled = alpha / (1.0 + snr_linear);
    double log_pd = 0.0;
    for (std::size_t i = 0; i < rank; ++i) {
        const auto ni = static_cast<double>(training_cells - i);
        log_pd += std::log(ni / (ni + scaled));
    }
    return std::exp(log_pd);
}

}  // namespace radar
