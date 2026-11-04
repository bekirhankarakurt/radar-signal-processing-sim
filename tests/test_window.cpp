#include "radar/complex_math.hpp"
#include "radar/fft.hpp"
#include "radar/window.hpp"
#include "test_framework.hpp"

using namespace radar;

namespace {

/// Peak sidelobe level (dB below the main lobe) of a real taper via a zero-padded FFT.
double peak_sidelobe_db(const RealVector& w) {
    const std::size_t n = 16384;
    ComplexVector x(n);
    for (std::size_t i = 0; i < w.size(); ++i) {
        x[i] = w[i];
    }
    FftPlan(n).execute(x, FftDirection::Forward);
    std::size_t k = 1;
    while (k < n / 2 && std::abs(x[k + 1]) < std::abs(x[k])) {
        ++k;  // walk down the main lobe to the first null
    }
    double sidelobe = 0.0;
    for (; k < n / 2; ++k) {
        sidelobe = std::max(sidelobe, std::abs(x[k]));
    }
    return 20.0 * std::log10(sidelobe / std::abs(x[0]));
}

}  // namespace

TEST(hamming_has_textbook_shape) {
    const RealVector w = make_window({WindowType::Hamming, 4, -35.0}, 65);
    CHECK_NEAR(w.front(), 0.08, 1e-12);
    CHECK_NEAR(w[32], 1.0, 1e-12);
    for (std::size_t i = 0; i < w.size(); ++i) {
        CHECK_NEAR(w[i], w[w.size() - 1 - i], 1e-12);
    }
    CHECK_NEAR(peak_sidelobe_db(w), -42.7, 1.0);
}

TEST(rectangular_sidelobes_are_minus_13_db) {
    CHECK_NEAR(peak_sidelobe_db(make_window({}, 64)), -13.26, 0.3);
}

TEST(taylor_meets_design_sidelobe_level) {
    for (const double sll : {-30.0, -35.0, -40.0}) {
        const RealVector w = make_window({WindowType::Taylor, 5, sll}, 128);
        const double psl = peak_sidelobe_db(w);
        CHECK(psl <= sll + 0.5);
        CHECK(psl >= sll - 3.0);
        CHECK_NEAR(*std::max_element(w.begin(), w.end()), 1.0, 1e-12);
    }
    CHECK_THROWS(make_window({WindowType::Taylor, 0, -35.0}, 16));
    CHECK_THROWS(make_window({WindowType::Taylor, 4, 10.0}, 16));
}

TEST(snr_gain_reflects_taper_loss) {
    CHECK_NEAR(window_snr_gain(make_window({}, 100)), 100.0, 1e-9);
    // Hamming processing loss is 1.34 dB.
    const double loss_db = 10.0 * std::log10(100.0 / window_snr_gain(make_window({WindowType::Hamming, 4, -35.0}, 100)));
    CHECK_NEAR(loss_db, 1.34, 0.05);
}

TEST(window_names_round_trip) {
    for (const auto t : {WindowType::Rectangular, WindowType::Hann, WindowType::Hamming, WindowType::Taylor}) {
        CHECK(parse_window_type(to_string(t)) == t);
    }
    CHECK(parse_window_type("NONE") == WindowType::Rectangular);
    CHECK(!parse_window_type("kaiser").has_value());
    CHECK_THROWS(make_window({}, 0));
    CHECK_NEAR(make_window({WindowType::Hann, 4, -35.0}, 1)[0], 1.0, 0.0);
}

TEST_MAIN()
