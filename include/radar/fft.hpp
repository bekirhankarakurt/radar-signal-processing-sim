#pragma once

#include "radar/types.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace radar {

enum class FftDirection : std::uint8_t { Forward, Inverse };

/// Precomputed iterative radix-2 decimation-in-time FFT.
///
/// The plan owns the twiddle factors and bit-reversal permutation for one
/// power-of-two length, so repeated transforms (one per pulse or per range
/// bin) cost only the butterflies. A plan is immutable after construction and
/// may be shared between threads.
class FftPlan {
public:
    explicit FftPlan(std::size_t size);

    [[nodiscard]] std::size_t size() const noexcept { return size_; }

    /// In-place transform of `size()` contiguous samples. The inverse is scaled by 1/N.
    void execute(Complex* data, FftDirection direction) const noexcept;
    void execute(ComplexVector& data, FftDirection direction) const;

private:
    std::size_t size_;
    std::vector<std::size_t> bit_reverse_;
    ComplexVector twiddles_;          // per stage: exp(-j pi k / half), k = 0 .. half - 1
    ComplexVector inverse_twiddles_;  // complex conjugates of twiddles_
};

/// Reference O(N^2) DFT, used to validate the FFT.
[[nodiscard]] ComplexVector naive_dft(const ComplexVector& input, FftDirection direction);

/// Swaps the two halves so that the zero-frequency bin is centred.
void fft_shift(Complex* data, std::size_t size) noexcept;

}  // namespace radar
