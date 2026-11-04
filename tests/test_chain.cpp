#include "radar/chain.hpp"
#include "radar/complex_math.hpp"
#include "radar/mti.hpp"
#include "radar/signal_model.hpp"
#include "test_helpers.hpp"
#include "test_framework.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace radar;

namespace {

Target target(double range, double velocity, double rcs) {
    Target t;
    t.range_m = range;
    t.velocity_mps = velocity;
    t.rcs_m2 = rcs;
    return t;
}

}  // namespace

TEST(radar_equation_scaling) {
    const RadarParameters r;
    const double p1 = received_power_w(r, 1000.0, 1.0);
    CHECK_NEAR(linear_to_db(received_power_w(r, 2000.0, 1.0) / p1), -12.041, 1e-3);
    CHECK_NEAR(received_power_w(r, 1000.0, 10.0) / p1, 10.0, 1e-9);
    // Hand calculation for the default radar at 1 km, 1 m^2.
    const double g = db_to_linear(32.0);
    const double lambda = kSpeedOfLight / 9.5e9;
    const double expected = 2e3 * g * g * lambda * lambda / (std::pow(4.0 * kPi, 3.0) * 1e12 * db_to_linear(4.0));
    CHECK_NEAR(p1 / expected, 1.0, 1e-12);
    CHECK_NEAR(r.noise_power_w(), kBoltzmann * 290.0 * db_to_linear(3.0) * 10e6, 1e-30);
}

TEST(targets_appear_at_their_range_and_velocity) {
    Scenario s = testing::small_scenario();
    s.targets = {target(900.0, 20.0, 0.1), target(2100.0, -35.0, 1.0)};
    const SignalModel model(s);
    const ProcessingChain chain(s);
    const CpiResult r = chain.process(model.generate(0), 0, nullptr);
    CHECK(r.detections.size() == 2);
    const double range_tol = s.radar.range_resolution_m();
    const double vel_tol = chain.axes().velocity_bin_mps;
    for (std::size_t i = 0; i < std::min<std::size_t>(r.detections.size(), 2); ++i) {
        CHECK_NEAR(r.detections[i].range_m, s.targets[i].range_m, range_tol);
        CHECK_NEAR(r.detections[i].velocity_mps, s.targets[i].velocity_mps, vel_tol);
    }
}

TEST(on_bin_peak_snr_equals_radar_equation_plus_processing_gain) {
    Scenario s = testing::small_scenario();
    s.noise = false;
    const ProcessingChain chain(s);
    const auto& axes = chain.axes();
    const double range = 50.0 * axes.range_bin_m;
    const double velocity = 5.0 * axes.velocity_bin_mps;
    s.targets = {target(range, velocity, 1.0)};
    const SignalModel model(s);
    const RealMatrix map = chain.range_doppler_map(model.generate(model.echoes_for_cpi(0), 0, false), nullptr);
    const double expected_db = linear_to_db(single_sample_snr(s.radar, range, 1.0)) +
                               linear_to_db(chain.compressor().snr_gain()) + linear_to_db(chain.doppler().snr_gain());
    CHECK_NEAR(linear_to_db(map(50, 16 + 5)), expected_db, 0.2);
}

TEST(closing_targets_have_positive_doppler) {
    Scenario s = testing::small_scenario();
    s.noise = false;
    s.targets = {target(1500.0, 30.0, 10.0)};
    const SignalModel model(s);
    const ProcessingChain chain(s);
    const RealMatrix map = chain.range_doppler_map(model.generate(0), nullptr);
    const auto peak = static_cast<std::size_t>(std::max_element(map.data().begin(), map.data().end()) -
                                               map.data().begin());
    CHECK(chain.axes().velocity(static_cast<double>(peak % map.cols())) > 25.0);
}

TEST(mti_suppresses_ground_clutter) {
    Scenario s = testing::small_scenario();
    s.clutter.enabled = true;
    s.clutter.cnr_db = 40.0;
    s.clutter.min_range_m = 300.0;
    s.clutter.max_range_m = 1500.0;
    s.clutter.spectral_width_mps = 0.1;
    auto zero_doppler_power = [](const Scenario& sc) {
        const SignalModel model(sc);
        const ProcessingChain chain(sc);
        const RealMatrix map = chain.range_doppler_map(model.generate(0), nullptr);
        double p = 0.0;
        for (std::size_t r = 20; r < 36; ++r) {  // 0.75 - 1.35 km, inside the clutter patch
            p += map(r, map.cols() / 2);
        }
        return p / 16.0;
    };
    const double without = zero_doppler_power(s);
    s.processing.mti = MtiType::ThreePulse;
    const double with = zero_doppler_power(s);
    CHECK(linear_to_db(without) > 50.0);
    CHECK(linear_to_db(without / with) > 30.0);
}

TEST(cholesky_reconstructs_the_clutter_correlation) {
    const std::size_t n = 32;
    RealMatrix a(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            a(i, j) = clutter_correlation(0.5, 0.03, 1e-4, i > j ? i - j : j - i) + (i == j ? 1e-8 : 0.0);
        }
    }
    const RealMatrix l = cholesky(a);
    double max_err = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            double v = 0.0;
            for (std::size_t k = 0; k < n; ++k) {
                v += l(i, k) * l(j, k);
            }
            max_err = std::max(max_err, std::abs(v - a(i, j)));
        }
    }
    CHECK(max_err < 1e-10);
    CHECK_NEAR(clutter_correlation(0.5, 0.03, 1e-4, 0), 1.0, 0.0);
    // sigma_f = 2 * 0.5 / 0.03 = 33.3 Hz; rho(1) = exp(-2 (pi * 33.3 * 1e-4)^2).
    CHECK_NEAR(clutter_correlation(0.5, 0.03, 1e-4, 1), std::exp(-2.0 * std::pow(kPi * 33.333333 * 1e-4, 2.0)), 1e-9);
    RealMatrix bad(2, 2, 1.0);
    bad(1, 1) = 0.5;
    CHECK_THROWS(cholesky(bad));
    CHECK_THROWS(cholesky(RealMatrix(2, 3)));
}

TEST(clutter_notch_discards_zero_velocity_hits) {
    Scenario s = testing::small_scenario();
    s.targets = {target(1000.0, 0.0, 0.05), target(2000.0, 40.0, 0.5)};
    const SignalModel model(s);
    const ProcessingChain open_chain(s);
    CHECK(open_chain.process(model.generate(0), 0, nullptr).detections.size() == 2);
    s.processing.clutter_notch_mps = 15.0;  // covers the Hamming main lobe (bins -1..1)
    const ProcessingChain notched(s);
    const CpiResult r = notched.process(model.generate(0), 0, nullptr);
    CHECK(r.detections.size() == 1);
    CHECK(!r.detections.empty() && std::abs(r.detections[0].velocity_mps - 40.0) < 5.0);
}

TEST(swerling1_amplitude_fluctuates_between_cpis) {
    Scenario s = testing::small_scenario();
    Target t = target(1000.0, 0.0, 1.0);
    t.swerling = Swerling::One;
    s.targets = {t, target(1200.0, 0.0, 1.0)};
    const SignalModel model(s);
    double min_ratio = 1e9;
    double max_ratio = 0.0;
    for (std::size_t c = 0; c < 20; ++c) {
        const auto e = model.echoes_for_cpi(c);
        const double ratio = power(e[0].amplitude) / single_sample_snr(s.radar, e[0].range_m, 1.0);
        min_ratio = std::min(min_ratio, ratio);
        max_ratio = std::max(max_ratio, ratio);
        CHECK_NEAR(power(e[1].amplitude), single_sample_snr(s.radar, e[1].range_m, 1.0), 1e-9);
    }
    CHECK(max_ratio / min_ratio > 3.0);
}

TEST(output_is_deterministic_across_runs_and_threads) {
    Scenario s = testing::small_scenario();
    s.targets = {target(1000.0, 10.0, 1.0)};
    s.clutter.enabled = true;
    const SignalModel model(s);
    const ComplexMatrix a = model.generate(3);
#ifdef _OPENMP
    const int saved = omp_get_max_threads();
    omp_set_num_threads(1);
#endif
    const ComplexMatrix b = model.generate(3);
#ifdef _OPENMP
    omp_set_num_threads(saved);
#endif
    CHECK(a.data() == b.data());
    CHECK(!(model.generate(4).data() == a.data()));
}

TEST(full_simulation_confirms_tracks) {
    Scenario s = testing::small_scenario();
    s.targets = {target(800.0, 25.0, 2.0), target(2200.0, -15.0, 2.0)};
    s.cpis = 6;
    StageTimer timer;
    const SimulationResult r = run_simulation(s, timer);
    CHECK(r.cpis.size() == 6);
    std::size_t confirmed = 0;
    for (const auto& t : r.tracks) {
        if (t.status == TrackStatus::Confirmed) {
            ++confirmed;
            const bool near_first = std::abs(t.velocity_mps - 25.0) < 2.0;
            const bool near_second = std::abs(t.velocity_mps + 15.0) < 2.0;
            CHECK(near_first || near_second);
        }
    }
    CHECK(confirmed == 2);
    CHECK(r.raw_samples == 6 * s.radar.pulses * s.radar.fast_time_samples());
    CHECK(timer.entries().size() >= 6);
    CHECK(timer.total_ms() > 0.0);
}

TEST_MAIN()
