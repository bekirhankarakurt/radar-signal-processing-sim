#pragma once

#include "radar/config.hpp"
#include "radar/types.hpp"

#include <cstdint>

namespace radar {

struct CfarWindow {
    std::size_t guard_range = 2;
    std::size_t guard_doppler = 1;
    std::size_t train_range = 8;
    std::size_t train_doppler = 4;

    /// Number of training cells in the 2-D annulus around the cell under test.
    [[nodiscard]] std::size_t training_cells() const noexcept;
};

struct CfarResult {
    Matrix<std::uint8_t> detections;  ///< 1 where the cell exceeded its threshold.
    RealMatrix noise_estimate;        ///< Local interference estimate per cell (0 where untested).
    std::size_t tested_cells = 0;
    std::size_t detection_count = 0;
};

/// CA-CFAR scale factor for square-law detected exponential noise:
/// alpha = N (Pfa^(-1/N) - 1), threshold = alpha * mean.
[[nodiscard]] double ca_cfar_alpha(std::size_t training_cells, double pfa);

/// OS-CFAR scale factor solving Pfa = prod_{i=0}^{k-1} (N - i) / (N - i + alpha), threshold = alpha * x_(k).
[[nodiscard]] double os_cfar_alpha(std::size_t training_cells, std::size_t rank, double pfa);

/// Analytic false-alarm probability of an OS-CFAR for a given scale factor.
[[nodiscard]] double os_cfar_pfa(std::size_t training_cells, std::size_t rank, double alpha) noexcept;

/// 2-D CFAR over a range-Doppler power map. The Doppler axis is treated as circular;
/// range cells whose window would leave the map are not tested.
class CfarDetector {
public:
    CfarDetector(CfarType type, const CfarWindow& window, double pfa, double os_rank_fraction = 0.75);

    [[nodiscard]] CfarResult detect(const RealMatrix& power_map) const;

    /// Decision for a single cell (used by the Monte Carlo engine where only the target cell matters).
    [[nodiscard]] bool test_cell(const RealMatrix& power_map, std::size_t r, std::size_t d) const;

    [[nodiscard]] double alpha() const noexcept { return alpha_; }
    [[nodiscard]] std::size_t os_rank() const noexcept { return rank_; }
    [[nodiscard]] CfarType type() const noexcept { return type_; }
    [[nodiscard]] const CfarWindow& window() const noexcept { return window_; }

private:
    [[nodiscard]] CfarResult detect_ca(const RealMatrix& map) const;
    [[nodiscard]] CfarResult detect_os(const RealMatrix& map) const;

    CfarType type_;
    CfarWindow window_;
    std::size_t rank_ = 0;
    double alpha_ = 0.0;
};

}  // namespace radar
