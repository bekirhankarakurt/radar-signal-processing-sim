#include "radar/cfar.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace radar {
namespace {

/// Copies the training cells around (r, d) into `cells` (size = training_cells()).
void gather_training(const RealMatrix& map, const CfarWindow& w, std::size_t r, std::size_t d, RealVector& cells) {
    const std::size_t cols = map.cols();
    const std::size_t hr = w.guard_range + w.train_range;
    const std::size_t hd = w.guard_doppler + w.train_doppler;
    std::size_t k = 0;
    for (std::size_t rr = r - hr; rr <= r + hr; ++rr) {
        const double* row = map.row(rr);
        const bool guard_row = rr + w.guard_range >= r && rr <= r + w.guard_range;
        for (std::size_t off = 0; off <= 2 * hd; ++off) {
            if (guard_row && off + w.guard_doppler >= hd && off <= hd + w.guard_doppler) {
                continue;
            }
            cells[k++] = row[(d + cols + off - hd) % cols];
        }
    }
}

}  // namespace

std::size_t CfarWindow::training_cells() const noexcept {
    const std::size_t outer = (2 * (guard_range + train_range) + 1) * (2 * (guard_doppler + train_doppler) + 1);
    const std::size_t inner = (2 * guard_range + 1) * (2 * guard_doppler + 1);
    return outer - inner;
}

double ca_cfar_alpha(std::size_t training_cells, double pfa) {
    if (training_cells == 0 || pfa <= 0.0 || pfa >= 1.0) {
        throw std::invalid_argument("invalid CA-CFAR parameters");
    }
    const auto n = static_cast<double>(training_cells);
    return n * (std::pow(pfa, -1.0 / n) - 1.0);
}

double os_cfar_pfa(std::size_t training_cells, std::size_t rank, double alpha) noexcept {
    double log_pfa = 0.0;
    for (std::size_t i = 0; i < rank; ++i) {
        const auto ni = static_cast<double>(training_cells - i);
        log_pfa += std::log(ni / (ni + alpha));
    }
    return std::exp(log_pfa);
}

double os_cfar_alpha(std::size_t training_cells, std::size_t rank, double pfa) {
    if (rank == 0 || rank > training_cells || pfa <= 0.0 || pfa >= 1.0) {
        throw std::invalid_argument("invalid OS-CFAR parameters");
    }
    double lo = 0.0;
    double hi = 1.0;
    while (os_cfar_pfa(training_cells, rank, hi) > pfa) {
        hi *= 2.0;
        if (hi > 1e12) {
            throw std::invalid_argument("OS-CFAR Pfa not reachable with this window");
        }
    }
    for (int it = 0; it < 200; ++it) {
        const double mid = 0.5 * (lo + hi);
        (os_cfar_pfa(training_cells, rank, mid) > pfa ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

CfarDetector::CfarDetector(CfarType type, const CfarWindow& window, double pfa, double os_rank_fraction)
    : type_(type), window_(window) {
    const std::size_t n = window.training_cells();
    if (n == 0) {
        throw std::invalid_argument("CFAR window has no training cells");
    }
    if (type == CfarType::CellAveraging) {
        alpha_ = ca_cfar_alpha(n, pfa);
    } else {
        const auto k = static_cast<std::size_t>(std::ceil(os_rank_fraction * static_cast<double>(n)));
        rank_ = std::clamp<std::size_t>(k, 1, n);
        alpha_ = os_cfar_alpha(n, rank_, pfa);
    }
}

bool CfarDetector::test_cell(const RealMatrix& power_map, std::size_t r, std::size_t d) const {
    const std::size_t hr = window_.guard_range + window_.train_range;
    if (r < hr || r + hr >= power_map.rows() || d >= power_map.cols()) {
        throw std::out_of_range("CFAR test cell too close to the map edge");
    }
    RealVector cells(window_.training_cells());
    gather_training(power_map, window_, r, d, cells);
    double statistic = 0.0;
    if (type_ == CfarType::CellAveraging) {
        for (const double v : cells) {
            statistic += v;
        }
        statistic /= static_cast<double>(cells.size());
    } else {
        const auto kth = cells.begin() + static_cast<std::ptrdiff_t>(rank_ - 1);
        std::nth_element(cells.begin(), kth, cells.end());
        statistic = *kth;
    }
    return power_map(r, d) > alpha_ * statistic;
}

CfarResult CfarDetector::detect(const RealMatrix& power_map) const {
    const std::size_t hr = window_.guard_range + window_.train_range;
    const std::size_t hd = window_.guard_doppler + window_.train_doppler;
    if (power_map.rows() < 2 * hr + 1 || power_map.cols() < 2 * hd + 1) {
        throw std::invalid_argument("CFAR window larger than the map");
    }
    return type_ == CfarType::CellAveraging ? detect_ca(power_map) : detect_os(power_map);
}

CfarResult CfarDetector::detect_ca(const RealMatrix& map) const {
    const std::size_t rows = map.rows();
    const std::size_t cols = map.cols();
    const std::size_t gr = window_.guard_range;
    const std::size_t gd = window_.guard_doppler;
    const std::size_t hr = gr + window_.train_range;
    const std::size_t hd = gd + window_.train_doppler;
    const std::size_t padded_cols = cols + 2 * hd;

    // Summed-area table over the map padded circularly along Doppler: every
    // rectangle sum is then four look-ups, so the cost per cell is O(1).
    RealMatrix sat(rows + 1, padded_cols + 1, 0.0);
    for (std::size_t r = 0; r < rows; ++r) {
        double row_sum = 0.0;
        const double* in = map.row(r);
        const double* above = sat.row(r);
        double* out = sat.row(r + 1);
        for (std::size_t c = 0; c < padded_cols; ++c) {
            row_sum += in[(c + cols - hd) % cols];
            out[c + 1] = above[c + 1] + row_sum;
        }
    }
    auto rect = [&sat](std::size_t r0, std::size_t r1, std::size_t c0, std::size_t c1) {
        return sat(r1, c1) - sat(r0, c1) - sat(r1, c0) + sat(r0, c0);  // [r0, r1) x [c0, c1)
    };

    CfarResult result;
    result.detections = Matrix<std::uint8_t>(rows, cols, 0);
    result.noise_estimate = RealMatrix(rows, cols, 0.0);
    const auto n = static_cast<double>(window_.training_cells());
    std::size_t hits = 0;
    const auto last_row = static_cast<std::ptrdiff_t>(rows - hr);
#pragma omp parallel for schedule(static) reduction(+ : hits)
    for (auto ri = static_cast<std::ptrdiff_t>(hr); ri < last_row; ++ri) {
        const auto r = static_cast<std::size_t>(ri);
        const double* cut = map.row(r);
        for (std::size_t d = 0; d < cols; ++d) {
            const std::size_t pc = d + hd;  // CUT column in padded coordinates
            const double outer = rect(r - hr, r + hr + 1, pc - hd, pc + hd + 1);
            const double inner = rect(r - gr, r + gr + 1, pc - gd, pc + gd + 1);
            const double noise = std::max(outer - inner, 0.0) / n;
            result.noise_estimate(r, d) = noise;
            if (cut[d] > alpha_ * noise) {
                result.detections(r, d) = 1;
                ++hits;
            }
        }
    }
    result.tested_cells = (rows - 2 * hr) * cols;
    result.detection_count = hits;
    return result;
}

CfarResult CfarDetector::detect_os(const RealMatrix& map) const {
    const std::size_t rows = map.rows();
    const std::size_t cols = map.cols();
    const std::size_t hr = window_.guard_range + window_.train_range;
    const std::size_t n = window_.training_cells();

    // E[x_(k)] for unit-mean exponential noise, used to turn the order statistic into a power estimate.
    double expected_order = 0.0;
    for (std::size_t i = 0; i < rank_; ++i) {
        expected_order += 1.0 / static_cast<double>(n - i);
    }

    CfarResult result;
    result.detections = Matrix<std::uint8_t>(rows, cols, 0);
    result.noise_estimate = RealMatrix(rows, cols, 0.0);
    std::size_t hits = 0;
    const auto last_row = static_cast<std::ptrdiff_t>(rows - hr);
#pragma omp parallel reduction(+ : hits)
    {
        RealVector cells(n);
#pragma omp for schedule(static)
        for (auto ri = static_cast<std::ptrdiff_t>(hr); ri < last_row; ++ri) {
            const auto r = static_cast<std::size_t>(ri);
            for (std::size_t d = 0; d < cols; ++d) {
                gather_training(map, window_, r, d, cells);
                const auto kth = cells.begin() + static_cast<std::ptrdiff_t>(rank_ - 1);
                std::nth_element(cells.begin(), kth, cells.end());
                result.noise_estimate(r, d) = *kth / expected_order;
                if (map(r, d) > alpha_ * *kth) {
                    result.detections(r, d) = 1;
                    ++hits;
                }
            }
        }
    }
    result.tested_cells = (rows - 2 * hr) * cols;
    result.detection_count = hits;
    return result;
}

}  // namespace radar
