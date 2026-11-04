#pragma once

#include "radar/types.hpp"

#include <array>
#include <cmath>
#include <cstdint>

namespace radar {

/// xoshiro256** pseudo-random generator (Blackman & Vigna), seeded through splitmix64.
/// Much faster than std::mt19937_64 and statistically strong enough for Monte Carlo work.
class Xoshiro256 {
public:
    explicit Xoshiro256(std::uint64_t seed) noexcept {
        for (auto& s : state_) {
            seed += 0x9E37'79B9'7F4A'7C15ULL;
            std::uint64_t z = seed;
            z = (z ^ (z >> 30U)) * 0xBF58'476D'1CE4'E5B9ULL;
            z = (z ^ (z >> 27U)) * 0x94D0'49BB'1331'11EBULL;
            s = z ^ (z >> 31U);
        }
    }

    [[nodiscard]] std::uint64_t next() noexcept {
        const std::uint64_t result = rotl(state_[1] * 5U, 7) * 9U;
        const std::uint64_t t = state_[1] << 17U;
        state_[2] ^= state_[0];
        state_[3] ^= state_[1];
        state_[1] ^= state_[2];
        state_[0] ^= state_[3];
        state_[2] ^= t;
        state_[3] = rotl(state_[3], 45);
        return result;
    }

    /// Uniform double in [0, 1).
    [[nodiscard]] double uniform() noexcept { return static_cast<double>(next() >> 11U) * 0x1.0p-53; }

    /// Circular complex Gaussian sample with E|z|^2 = variance (Marsaglia polar method).
    [[nodiscard]] Complex complex_gaussian(double variance = 1.0) noexcept {
        double u = 0.0;
        double v = 0.0;
        double s = 0.0;
        do {
            u = 2.0 * uniform() - 1.0;
            v = 2.0 * uniform() - 1.0;
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);
        // Each component has variance / 2.
        const double scale = std::sqrt(-variance * std::log(s) / s);
        return {u * scale, v * scale};
    }

    /// Unit-mean exponential sample.
    [[nodiscard]] double exponential() noexcept { return -std::log1p(-uniform()); }

private:
    [[nodiscard]] static constexpr std::uint64_t rotl(std::uint64_t x, int k) noexcept {
        return (x << static_cast<unsigned>(k)) | (x >> static_cast<unsigned>(64 - k));
    }

    std::array<std::uint64_t, 4> state_{};
};

}  // namespace radar
