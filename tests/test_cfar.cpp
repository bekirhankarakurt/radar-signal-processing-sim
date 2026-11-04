#include "radar/cfar.hpp"
#include "test_framework.hpp"

#include <random>

using namespace radar;

namespace {

RealMatrix exponential_noise(std::size_t rows, std::size_t cols, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::exponential_distribution<double> e(1.0);
    RealMatrix m(rows, cols);
    for (auto& v : m.data()) {
        v = e(rng);
    }
    return m;
}

double empirical_pfa(const CfarDetector& det, std::size_t maps) {
    std::size_t hits = 0;
    std::size_t tested = 0;
    for (std::size_t i = 0; i < maps; ++i) {
        const CfarResult r = det.detect(exponential_noise(400, 128, 100 + i));
        hits += r.detection_count;
        tested += r.tested_cells;
    }
    return static_cast<double>(hits) / static_cast<double>(tested);
}

}  // namespace

TEST(ca_alpha_closed_form) {
    CHECK_NEAR(ca_cfar_alpha(1, 1e-3), 999.0, 1e-9);
    // Large N tends to the fixed-threshold value -ln(Pfa).
    CHECK_NEAR(ca_cfar_alpha(100000, 1e-6), -std::log(1e-6), 0.01);
    CHECK_THROWS(ca_cfar_alpha(0, 1e-3));
    CHECK_THROWS(ca_cfar_alpha(10, 1.0));
}

TEST(os_alpha_inverts_the_pfa_expression) {
    for (const double pfa : {1e-3, 1e-6}) {
        const double a = os_cfar_alpha(48, 36, pfa);
        CHECK_NEAR(std::log10(os_cfar_pfa(48, 36, a)), std::log10(pfa), 1e-6);
    }
    // k = N = 1 reduces to the single-cell CA case.
    CHECK_NEAR(os_cfar_alpha(1, 1, 1e-3), ca_cfar_alpha(1, 1e-3), 1e-6);
    CHECK_THROWS(os_cfar_alpha(10, 11, 1e-3));
    CHECK_THROWS(os_cfar_alpha(10, 0, 1e-3));
}

TEST(window_cell_count) {
    const CfarWindow w{2, 1, 8, 4};
    CHECK(w.training_cells() == 21 * 11 - 5 * 3);
}

TEST(ca_cfar_false_alarm_rate_matches_design) {
    const CfarDetector det(CfarType::CellAveraging, {2, 1, 8, 3}, 1e-3);
    CHECK_NEAR(empirical_pfa(det, 20) / 1e-3, 1.0, 0.15);
}

TEST(os_cfar_false_alarm_rate_matches_design) {
    const CfarDetector det(CfarType::OrderedStatistic, {2, 1, 8, 3}, 1e-3, 0.75);
    CHECK_NEAR(empirical_pfa(det, 20) / 1e-3, 1.0, 0.15);
}

TEST(detects_strong_target_and_estimates_noise) {
    for (const auto type : {CfarType::CellAveraging, CfarType::OrderedStatistic}) {
        RealMatrix map = exponential_noise(200, 64, 1);
        map(100, 20) = 1000.0;  // 30 dB
        const CfarDetector det(type, {2, 1, 8, 3}, 1e-6);
        const CfarResult r = det.detect(map);
        CHECK(r.detections(100, 20) == 1);
        CHECK_NEAR(r.noise_estimate(100, 20), 1.0, 0.35);
        CHECK(r.tested_cells == (200 - 20) * 64);
        CHECK(det.test_cell(map, 100, 20));
    }
}

TEST(single_cell_test_agrees_with_full_map) {
    const RealMatrix map = exponential_noise(120, 32, 77);
    for (const auto type : {CfarType::CellAveraging, CfarType::OrderedStatistic}) {
        const CfarDetector det(type, {1, 1, 6, 2}, 0.05);
        const CfarResult r = det.detect(map);
        for (std::size_t row = 7; row < 113; row += 3) {
            for (std::size_t d = 0; d < 32; ++d) {
                CHECK(det.test_cell(map, row, d) == (r.detections(row, d) == 1));
            }
        }
        CHECK_THROWS(static_cast<void>(det.test_cell(map, 2, 0)));
    }
}

TEST(doppler_axis_wraps_around) {
    RealMatrix map(60, 16, 1.0);
    map(30, 0) = 100.0;
    map(30, 15) = 50.0;  // neighbour across the wrap, inside the guard band
    const CfarDetector det(CfarType::CellAveraging, {1, 1, 4, 2}, 1e-4);
    const CfarResult r = det.detect(map);
    CHECK(r.detections(30, 0) == 1);
    CHECK(r.detections(30, 15) == 1);
    CHECK(r.detection_count == 2);
}

TEST(os_cfar_resists_masking_by_an_interfering_target) {
    // Flat unit noise, a 30 dB target inside the training window of a 22 dB target.
    RealMatrix map(100, 32, 1.0);
    map(50, 10) = 1000.0;
    map(54, 10) = 158.0;
    const CfarWindow w{1, 0, 8, 1};
    const CfarResult ca = CfarDetector(CfarType::CellAveraging, w, 1e-6).detect(map);
    const CfarResult os = CfarDetector(CfarType::OrderedStatistic, w, 1e-6, 0.75).detect(map);
    CHECK(ca.detections(50, 10) == 1);
    CHECK(ca.detections(54, 10) == 0);  // masked
    CHECK(os.detections(50, 10) == 1);
    CHECK(os.detections(54, 10) == 1);
}

TEST(rejects_window_larger_than_map) {
    const CfarDetector det(CfarType::CellAveraging, {2, 1, 8, 3}, 1e-3);
    CHECK_THROWS(det.detect(RealMatrix(10, 64, 1.0)));
    CHECK_THROWS(CfarDetector(CfarType::CellAveraging, {1, 1, 0, 0}, 1e-3));
}

TEST_MAIN()
