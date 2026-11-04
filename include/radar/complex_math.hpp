#pragma once

#include "radar/types.hpp"

#include <cmath>

namespace radar {

/// Unit phasor exp(j * phase).
[[nodiscard]] inline Complex cis(double phase) noexcept { return {std::cos(phase), std::sin(phase)}; }

/// |z|^2 without the square root of std::abs.
[[nodiscard]] inline double power(const Complex& z) noexcept { return z.real() * z.real() + z.imag() * z.imag(); }

[[nodiscard]] inline double db_to_linear(double db) noexcept { return std::pow(10.0, db / 10.0); }

/// Power ratio to decibels with a floor so that zero maps to a finite value.
[[nodiscard]] inline double linear_to_db(double linear, double floor_db = -300.0) noexcept {
    return linear > 0.0 ? 10.0 * std::log10(linear) : floor_db;
}

[[nodiscard]] inline double total_energy(const ComplexVector& v) noexcept {
    double e = 0.0;
    for (const auto& z : v) {
        e += power(z);
    }
    return e;
}

[[nodiscard]] constexpr bool is_power_of_two(std::size_t n) noexcept { return n != 0 && (n & (n - 1)) == 0; }

[[nodiscard]] constexpr std::size_t next_power_of_two(std::size_t n) noexcept {
    std::size_t p = 1;
    while (p < n) {
        p <<= 1U;
    }
    return p;
}

}  // namespace radar
