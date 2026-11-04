#include "radar/monte_carlo.hpp"

#include "radar/chain.hpp"
#include "radar/complex_math.hpp"
#include "radar/random.hpp"
#include "radar/signal_model.hpp"
#include "radar/statistics.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace radar {

MonteCarloResult run_monte_carlo(const Scenario& input, StageTimer& timer) {
    Scenario scenario = input;
    scenario.targets.clear();
    const auto& mc = scenario.monte_carlo;
    const SignalModel model(scenario);
    const ProcessingChain chain(scenario);
    const auto& axes = chain.axes();
    const auto& cfar = chain.cfar();

    // Put the test target exactly on a range-bin and Doppler-bin centre so that the
    // configured SNR is the SNR of the cell under test (no straddle loss).
    const std::size_t range_bins = scenario.radar.range_bins();
    const double requested_range = mc.range_m > 0.0 ? mc.range_m : 0.5 * scenario.radar.max_range_m;
    const auto range_bin = static_cast<std::size_t>(std::llround(requested_range / axes.range_bin_m));
    const std::size_t guard = cfar.window().guard_range + cfar.window().train_range;
    if (range_bin < guard || range_bin + guard >= range_bins) {
        throw std::invalid_argument("monte_carlo.range_m too close to the edge of the range window");
    }
    const double requested_velocity = mc.velocity_mps.value_or(0.5 * scenario.radar.max_unambiguous_velocity_mps());
    const auto doppler_offset = std::llround(requested_velocity / axes.velocity_bin_mps);
    const auto half = static_cast<long long>(axes.doppler_bins / 2);
    if (doppler_offset <= -half || doppler_offset >= half) {
        throw std::invalid_argument("monte_carlo.velocity_mps outside the unambiguous velocity interval");
    }
    const auto doppler_bin = static_cast<std::size_t>(half + doppler_offset);

    MonteCarloResult result;
    result.target_range_m = static_cast<double>(range_bin) * axes.range_bin_m;
    result.target_velocity_mps = static_cast<double>(doppler_offset) * axes.velocity_bin_mps;
    result.pfa_design = scenario.processing.cfar_pfa;
    result.cfar_alpha = cfar.alpha();
    result.training_cells = cfar.window().training_cells();

    // Noise-free response of the whole chain to a unit-amplitude echo: this is the exact
    // processing gain into the test cell, including window, MTI and range-migration losses.
    double unit_gain = 0.0;
    {
        const auto s = timer.scope("calibration");
        const std::vector<TargetEcho> unit{{result.target_range_m, result.target_velocity_mps, Complex{1.0, 0.0}}};
        Scenario quiet = scenario;
        quiet.clutter.enabled = false;
        const SignalModel quiet_model(quiet);
        const RealMatrix map = chain.range_doppler_map(quiet_model.generate(unit, 0, false), nullptr);
        unit_gain = map(range_bin, doppler_bin);
    }
    result.processing_gain_db = linear_to_db(unit_gain);

    const std::size_t steps =
        static_cast<std::size_t>(std::floor((mc.snr_max_db - mc.snr_min_db) / mc.snr_step_db + 1e-9)) + 1;
    {
        const auto s = timer.scope("Pd trials");
        for (std::size_t si = 0; si < steps; ++si) {
            const double snr_db = mc.snr_min_db + static_cast<double>(si) * mc.snr_step_db;
            const double snr = db_to_linear(snr_db);
            const double amplitude = std::sqrt(snr / unit_gain);
            std::size_t hits = 0;
            const auto trials = static_cast<std::ptrdiff_t>(mc.trials);
#pragma omp parallel for schedule(dynamic) reduction(+ : hits)
            for (std::ptrdiff_t ti = 0; ti < trials; ++ti) {
                const auto trial = static_cast<std::uint64_t>(ti);
                const std::uint64_t stream = mix_seed(scenario.seed, 0x5D00ULL + si, trial);
                Xoshiro256 rng(stream);
                const double a = mc.swerling == Swerling::One ? amplitude * std::sqrt(rng.exponential()) : amplitude;
                const std::vector<TargetEcho> echo{
                    {result.target_range_m, result.target_velocity_mps, a * cis(2.0 * kPi * rng.uniform())}};
                const RealMatrix map = chain.range_doppler_map(model.generate(echo, stream, true), nullptr);
                if (cfar.test_cell(map, range_bin, doppler_bin)) {
                    ++hits;
                }
            }
            PdPoint p;
            p.snr_db = snr_db;
            p.trials = mc.trials;
            p.detections = hits;
            p.pd = static_cast<double>(hits) / static_cast<double>(mc.trials);
            const double pfa = scenario.processing.cfar_pfa;
            if (mc.swerling == Swerling::Zero) {
                p.pd_ideal = pd_swerling0_ideal(snr, pfa);
                p.pd_cfar = std::numeric_limits<double>::quiet_NaN();
            } else {
                p.pd_ideal = pd_swerling1_ideal(snr, pfa);
                p.pd_cfar = cfar.type() == CfarType::CellAveraging
                                ? pd_swerling1_ca_cfar(snr, result.training_cells, cfar.alpha())
                                : pd_swerling1_os_cfar(snr, result.training_cells, cfar.os_rank(), cfar.alpha());
            }
            result.points.push_back(p);
        }
    }

    {
        const auto s = timer.scope("Pfa trials");
        std::size_t false_alarms = 0;
        std::size_t tested = 0;
        const auto trials = static_cast<std::ptrdiff_t>(mc.pfa_trials);
#pragma omp parallel for schedule(dynamic) reduction(+ : false_alarms, tested)
        for (std::ptrdiff_t ti = 0; ti < trials; ++ti) {
            const std::uint64_t stream = mix_seed(scenario.seed, 0xFA00ULL, static_cast<std::uint64_t>(ti));
            const RealMatrix map = chain.range_doppler_map(model.generate({}, stream, true), nullptr);
            const CfarResult r = cfar.detect(map);
            false_alarms += r.detection_count;
            tested += r.tested_cells;
        }
        result.false_alarms = false_alarms;
        result.cells_tested = tested;
        result.pfa_empirical = tested > 0 ? static_cast<double>(false_alarms) / static_cast<double>(tested) : 0.0;
    }
    return result;
}

}  // namespace radar
