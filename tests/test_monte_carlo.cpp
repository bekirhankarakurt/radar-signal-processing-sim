#include "radar/complex_math.hpp"
#include "radar/monte_carlo.hpp"
#include "test_helpers.hpp"
#include "test_framework.hpp"

using namespace radar;

TEST(swerling1_pd_follows_ca_cfar_closed_form) {
    Scenario s = testing::small_scenario();
    s.processing.cfar_pfa = 1e-3;
    s.monte_carlo.swerling = Swerling::One;
    s.monte_carlo.snr_min_db = 5.0;
    s.monte_carlo.snr_max_db = 25.0;
    s.monte_carlo.snr_step_db = 10.0;
    s.monte_carlo.trials = 400;
    s.monte_carlo.pfa_trials = 30;
    s.monte_carlo.velocity_mps = 20.0;
    StageTimer timer;
    const MonteCarloResult r = run_monte_carlo(s, timer);
    CHECK(r.points.size() == 3);
    for (const auto& p : r.points) {
        CHECK_NEAR(p.pd, p.pd_cfar, 0.08);
        CHECK(p.pd_cfar <= p.pd_ideal + 1e-12);
    }
    CHECK(r.points.front().pd < r.points.back().pd);
    // Empirical Pfa over ~59k tested cells. Window-correlated reference cells push it above the
    // design value (see README), so the bound is a factor of three.
    CHECK(r.pfa_empirical > 0.5e-3);
    CHECK(r.pfa_empirical < 3.0e-3);
    CHECK(r.cells_tested > 50000);
}

TEST(swerling0_pd_is_monotonic_and_reaches_one) {
    Scenario s = testing::small_scenario();
    s.processing.cfar = CfarType::OrderedStatistic;
    s.processing.cfar_pfa = 1e-4;
    s.monte_carlo.swerling = Swerling::Zero;
    s.monte_carlo.snr_min_db = 0.0;
    s.monte_carlo.snr_max_db = 24.0;
    s.monte_carlo.snr_step_db = 8.0;
    s.monte_carlo.trials = 100;
    s.monte_carlo.pfa_trials = 2;
    StageTimer timer;
    const MonteCarloResult r = run_monte_carlo(s, timer);
    CHECK(r.points.size() == 4);
    CHECK(r.points.front().pd < 0.1);
    CHECK_NEAR(r.points.back().pd, 1.0, 1e-12);
    CHECK(std::isnan(r.points.front().pd_cfar));
    // Coherent gain = compression gain x Doppler gain (Taylor + Hamming tapers).
    CHECK(r.processing_gain_db > 25.0);
}

TEST(rejects_test_target_outside_the_map) {
    Scenario s = testing::small_scenario();
    s.monte_carlo.range_m = 10.0;
    StageTimer timer;
    CHECK_THROWS(run_monte_carlo(s, timer));
    s.monte_carlo.range_m = 0.0;
    s.monte_carlo.velocity_mps = 500.0;
    CHECK_THROWS(run_monte_carlo(s, timer));
}

TEST_MAIN()
