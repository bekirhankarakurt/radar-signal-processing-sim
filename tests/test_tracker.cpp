#include "radar/tracker.hpp"
#include "test_framework.hpp"

using namespace radar;

namespace {

Detection det(double range, double velocity) {
    Detection d;
    d.range_m = range;
    d.velocity_mps = velocity;
    d.snr_db = 20.0;
    return d;
}

TrackerParameters params() {
    TrackerParameters p;
    p.confirm_hits = 3;
    p.confirm_window = 4;
    p.max_misses = 2;
    p.range_gate_m = 50.0;
    p.velocity_gate_mps = 3.0;
    return p;
}

}  // namespace

TEST(constant_velocity_target_is_confirmed_and_followed) {
    const double dt = 0.5;
    Tracker tracker(params(), dt);
    for (std::size_t k = 0; k < 6; ++k) {
        const double range = 5000.0 - 40.0 * dt * static_cast<double>(k);  // closing at 40 m/s
        tracker.update({det(range, 40.0)}, k);
        if (k == 1) {
            CHECK(tracker.confirmed().empty());
        }
    }
    const auto confirmed = tracker.confirmed();
    CHECK(confirmed.size() == 1);
    CHECK(confirmed[0].hits == 6);
    CHECK_NEAR(confirmed[0].range_m, 5000.0 - 40.0 * dt * 5.0, 1.0);
    CHECK_NEAR(confirmed[0].velocity_mps, 40.0, 1e-9);
    CHECK(confirmed[0].first_cpi == 0);
    CHECK(confirmed[0].last_cpi == 5);
}

TEST(m_of_n_tolerates_one_miss) {
    Tracker tracker(params(), 1.0);
    tracker.update({det(1000.0, 0.0)}, 0);
    tracker.update({}, 1);
    tracker.update({det(1000.0, 0.0)}, 2);
    tracker.update({det(1000.0, 0.0)}, 3);
    CHECK(tracker.confirmed().size() == 1);
}

TEST(isolated_false_alarm_is_dropped) {
    Tracker tracker(params(), 1.0);
    tracker.update({det(1000.0, 0.0)}, 0);
    tracker.update({}, 1);
    tracker.update({}, 2);
    CHECK(tracker.tracks().size() == 1);
    CHECK(tracker.tracks()[0].status == TrackStatus::Deleted);
    CHECK(tracker.confirmed().empty());
}

TEST(confirmed_track_deleted_after_consecutive_misses) {
    Tracker tracker(params(), 1.0);
    for (std::size_t k = 0; k < 3; ++k) {
        tracker.update({det(2000.0, 0.0)}, k);
    }
    CHECK(tracker.confirmed().size() == 1);
    tracker.update({}, 3);
    CHECK(tracker.confirmed().size() == 1);
    tracker.update({}, 4);
    CHECK(tracker.confirmed().empty());
}

TEST(gates_separate_crossing_targets) {
    Tracker tracker(params(), 1.0);
    for (std::size_t k = 0; k < 4; ++k) {
        const double t = static_cast<double>(k);
        tracker.update({det(3000.0 - 10.0 * t, 10.0), det(3010.0 + 20.0 * t, -20.0)}, k);
    }
    CHECK(tracker.confirmed().size() == 2);
    CHECK(to_string(TrackStatus::Confirmed) == "confirmed");
}

TEST_MAIN()
