#include "radar/clustering.hpp"
#include "test_framework.hpp"

#include <cmath>

using namespace radar;

namespace {

CfarResult hits_above(const RealMatrix& map, double threshold) {
    CfarResult r;
    r.detections = Matrix<std::uint8_t>(map.rows(), map.cols(), 0);
    r.noise_estimate = RealMatrix(map.rows(), map.cols(), 1.0);
    for (std::size_t i = 0; i < map.rows(); ++i) {
        for (std::size_t j = 0; j < map.cols(); ++j) {
            if (map(i, j) > threshold) {
                r.detections(i, j) = 1;
                ++r.detection_count;
            }
        }
    }
    return r;
}

}  // namespace

TEST(centroid_recovers_sub_bin_position) {
    RealMatrix map(64, 32, 1.0);
    const double r0 = 20.3;
    const double d0 = 10.6;
    for (std::size_t r = 0; r < 64; ++r) {
        for (std::size_t d = 0; d < 32; ++d) {
            const double dr = static_cast<double>(r) - r0;
            const double dd = static_cast<double>(d) - d0;
            map(r, d) += 1000.0 * std::exp(-(dr * dr + dd * dd) / 1.5);
        }
    }
    const RangeDopplerAxes axes{15.0, 0.5, 32};
    const auto dets = cluster_detections(map, hits_above(map, 20.0), axes, 4);
    CHECK(dets.size() == 1);
    CHECK_NEAR(dets[0].range_bin, r0, 0.1);
    CHECK_NEAR(dets[0].doppler_bin, d0, 0.1);
    CHECK_NEAR(dets[0].range_m, r0 * 15.0, 1.5);
    CHECK_NEAR(dets[0].velocity_mps, (d0 - 16.0) * 0.5, 0.05);
    CHECK(dets[0].cpi == 4);
    CHECK(dets[0].cells > 1);
    CHECK(dets[0].snr_db > 25.0);
}

TEST(separate_blobs_and_doppler_wrap) {
    RealMatrix map(40, 16, 1.0);
    map(10, 0) = 100.0;
    map(10, 15) = 100.0;  // same cluster across the Doppler wrap
    map(30, 8) = 50.0;    // isolated second cluster
    const RangeDopplerAxes axes{10.0, 1.0, 16};
    const auto dets = cluster_detections(map, hits_above(map, 10.0), axes, 0);
    CHECK(dets.size() == 2);
    CHECK_NEAR(dets[0].range_bin, 10.0, 1e-12);
    const double wrapped = dets[0].doppler_bin < 8.0 ? dets[0].doppler_bin + 16.0 : dets[0].doppler_bin;
    CHECK_NEAR(wrapped, 15.5, 1e-9);
    CHECK(dets[0].cells == 2);
    CHECK_NEAR(dets[1].range_m, 300.0, 1e-9);
}

TEST(no_hits_no_detections) {
    const RealMatrix map(10, 8, 1.0);
    CHECK(cluster_detections(map, hits_above(map, 5.0), {1.0, 1.0, 8}, 0).empty());
}

TEST_MAIN()
