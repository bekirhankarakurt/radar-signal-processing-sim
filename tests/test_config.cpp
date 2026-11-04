#include "radar/config.hpp"
#include "test_framework.hpp"

#include <filesystem>

using namespace radar;

namespace {

const char* const kScenario = R"(
# comment line
[simulation]
name = unit ; trailing comment
cpis = 3
seed = 99

[radar]
carrier_hz = 10e9
bandwidth_hz = 2e6
pulse_width_s = 20e-6
prf_hz = 5e3
pulses = 32
sample_rate_hz = 4e6
max_range_m = 10000

[processing]
range_window = hamming
doppler_window = taylor
doppler_window_taylor_sll_db = -40
mti = three_pulse
cfar = os
os_rank = 0.5

[target]
range_m = 1000
velocity_mps = 10
rcs_m2 = 2

[target]
range_m = 2000
velocity_mps = -5
rcs_m2 = 0.5
swerling = 1
)";

}  // namespace

TEST(parses_all_sections_and_repeated_targets) {
    const Scenario s = scenario_from_ini(kScenario);
    CHECK(s.name == "unit");
    CHECK(s.cpis == 3);
    CHECK(s.seed == 99);
    CHECK_NEAR(s.radar.carrier_hz, 10e9, 0.0);
    CHECK(s.radar.pulses == 32);
    CHECK(s.processing.range_window.type == WindowType::Hamming);
    CHECK(s.processing.doppler_window.type == WindowType::Taylor);
    CHECK_NEAR(s.processing.doppler_window.taylor_sidelobe_db, -40.0, 0.0);
    CHECK(s.processing.mti == MtiType::ThreePulse);
    CHECK(s.processing.cfar == CfarType::OrderedStatistic);
    CHECK(s.targets.size() == 2);
    CHECK(s.targets[1].swerling == Swerling::One);
    CHECK_NEAR(s.targets[1].velocity_mps, -5.0, 0.0);
}

TEST(derived_radar_quantities) {
    const Scenario s = scenario_from_ini(kScenario);
    CHECK_NEAR(s.radar.wavelength(), 0.0299792458, 1e-12);
    CHECK_NEAR(s.radar.range_resolution_m(), 74.948, 1e-3);
    CHECK_NEAR(s.radar.unambiguous_range_m(), 29979.2458, 1e-3);
    CHECK_NEAR(s.radar.max_unambiguous_velocity_mps(), 37.474, 1e-3);
    CHECK(s.radar.pulse_samples() == 80);
    CHECK(s.radar.range_bins() == 268);
}

TEST(rejects_bad_input_with_context) {
    CHECK_THROWS(scenario_from_ini("[radar]\nbogus = 1\n"));
    CHECK_THROWS(scenario_from_ini("[nonsense]\n"));
    CHECK_THROWS(scenario_from_ini("[radar]\nprf_hz = fast\n"));
    CHECK_THROWS(scenario_from_ini("[radar]\nprf_hz = 1\nprf_hz = 2\n"));
    CHECK_THROWS(scenario_from_ini("prf_hz = 1\n"));
    CHECK_THROWS(scenario_from_ini("[radar\n"));
    CHECK_THROWS(scenario_from_ini("[radar]\nno equals sign\n"));
    CHECK_THROWS(scenario_from_ini("[radar]\npulses = 2.5\n"));
    CHECK_THROWS(scenario_from_ini("[processing]\nmti = four_pulse\n"));
    CHECK_THROWS(scenario_from_ini("[simulation]\nnoise = maybe\n"));
    try {
        static_cast<void>(scenario_from_ini("\n\n[radar]\nbandwith_hz = 1e6\n"));
        CHECK(false);
    } catch (const std::invalid_argument& e) {
        const std::string msg = e.what();
        CHECK(msg.find("line 3") != std::string::npos);
        CHECK(msg.find("bandwith_hz") != std::string::npos);
    }
}

TEST(validation_catches_inconsistent_parameters) {
    CHECK_THROWS(scenario_from_ini("[radar]\nsample_rate_hz = 1e6\nbandwidth_hz = 5e6\n"));
    CHECK_THROWS(scenario_from_ini("[radar]\nprf_hz = 20e3\nmax_range_m = 9000\n"));
    CHECK_THROWS(scenario_from_ini("[target]\nrange_m = 50000\n"));
    CHECK_THROWS(scenario_from_ini("[processing]\ndoppler_fft = 48\n"));
    CHECK_THROWS(scenario_from_ini("[processing]\ncfar_pfa = 1.5\n"));
    CHECK_THROWS(scenario_from_ini("[tracker]\nconfirm_hits = 5\nconfirm_window = 4\n"));
}

TEST(shipped_scenarios_load) {
    const std::filesystem::path dir = std::filesystem::path(RADAR_SOURCE_DIR) / "scenarios";
    std::size_t count = 0;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        if (entry.path().extension() == ".ini") {
            const Scenario s = load_scenario(entry.path());
            CHECK(s.radar.pulses >= 4);
            ++count;
        }
    }
    CHECK(count >= 1);
    CHECK_THROWS(load_scenario(dir / "does_not_exist.ini"));
}

TEST(enum_names) {
    CHECK(to_string(MtiType::TwoPulse) == "two_pulse");
    CHECK(to_string(CfarType::OrderedStatistic) == "OS-CFAR");
    CHECK(to_string(Swerling::One) == "Swerling 1");
}

TEST_MAIN()
