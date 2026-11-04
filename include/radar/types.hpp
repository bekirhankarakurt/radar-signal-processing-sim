#pragma once

#include <algorithm>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace radar {

using Complex = std::complex<double>;
using ComplexVector = std::vector<Complex>;
using RealVector = std::vector<double>;

inline constexpr double kSpeedOfLight = 299'792'458.0;  // m/s
inline constexpr double kBoltzmann = 1.380'649e-23;     // J/K
inline constexpr double kStandardTemperature = 290.0;   // K
inline constexpr double kPi = 3.141'592'653'589'793'238'46;

/// Dense row-major 2-D array. Rows are contiguous so that a row can be handed
/// to an FFT or a filter without copying.
template <typename T>
class Matrix {
public:
    Matrix() = default;
    Matrix(std::size_t rows, std::size_t cols, T value = T{})
        : rows_(rows), cols_(cols), data_(rows * cols, value) {}

    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t cols() const noexcept { return cols_; }
    [[nodiscard]] std::size_t size() const noexcept { return data_.size(); }
    [[nodiscard]] bool empty() const noexcept { return data_.empty(); }

    [[nodiscard]] T& operator()(std::size_t r, std::size_t c) noexcept { return data_[r * cols_ + c]; }
    [[nodiscard]] const T& operator()(std::size_t r, std::size_t c) const noexcept {
        return data_[r * cols_ + c];
    }

    [[nodiscard]] T* row(std::size_t r) noexcept { return data_.data() + r * cols_; }
    [[nodiscard]] const T* row(std::size_t r) const noexcept { return data_.data() + r * cols_; }

    [[nodiscard]] std::vector<T>& data() noexcept { return data_; }
    [[nodiscard]] const std::vector<T>& data() const noexcept { return data_; }

private:
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    std::vector<T> data_;
};

using ComplexMatrix = Matrix<Complex>;
using RealMatrix = Matrix<double>;

/// Transposes `in` into `out` using cache-sized tiles (the radar "corner turn").
template <typename T>
void transpose(const Matrix<T>& in, Matrix<T>& out) {
    constexpr std::size_t kTile = 32;
    out = Matrix<T>(in.cols(), in.rows());
    const auto row_tiles = static_cast<std::ptrdiff_t>((in.rows() + kTile - 1) / kTile);
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t tile = 0; tile < row_tiles; ++tile) {
        const std::size_t r0 = static_cast<std::size_t>(tile) * kTile;
        for (std::size_t c0 = 0; c0 < in.cols(); c0 += kTile) {
            const std::size_t r1 = std::min(r0 + kTile, in.rows());
            const std::size_t c1 = std::min(c0 + kTile, in.cols());
            for (std::size_t r = r0; r < r1; ++r) {
                for (std::size_t c = c0; c < c1; ++c) {
                    out(c, r) = in(r, c);
                }
            }
        }
    }
}

}  // namespace radar
