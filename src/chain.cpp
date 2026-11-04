#include "radar/chain.hpp"

#include "radar/complex_math.hpp"
#include "radar/mti.hpp"
#include "radar/waveform.hpp"

#include <cmath>
#include <optional>

namespace radar {
namespace {

std::size_t mti_order(MtiType type) { return mti_coefficients(type).size() - 1; }

std::size_t doppler_fft_size(const Scenario& s) {
    const std::size_t slow_time = s.radar.pulses - mti_order(s.processing.mti);
    return s.processing.doppler_fft == 0 ? next_power_of_two(slow_time) : s.processing.doppler_fft;
}

CfarWindow cfar_window(const ProcessingParameters& p) {
    return {p.guard_range, p.guard_doppler, p.train_range, p.train_doppler};
}

/// Optional stage timing: a no-op when no timer is supplied.
class MaybeScope {
public:
    MaybeScope(StageTimer* timer, const char* name) {
        if (timer != nullptr) {
            scope_.emplace(*timer, name);
        }
    }

private:
    std::optional<StageTimer::Scope> scope_;
};

}  // namespace

ProcessingChain::ProcessingChain(const Scenario& scenario)
    : scenario_(scenario),
      compressor_(make_lfm_chirp(scenario.radar.bandwidth_hz, scenario.radar.pulse_width_s,
                                 scenario.radar.sample_rate_hz),
                  scenario.processing.range_window, scenario.radar.fast_time_samples(), scenario.radar.range_bins()),
      doppler_(scenario.radar.pulses - mti_order(scenario.processing.mti), doppler_fft_size(scenario),
               scenario.processing.doppler_window),
      cfar_(scenario.processing.cfar, cfar_window(scenario.processing), scenario.processing.cfar_pfa,
            scenario.processing.os_rank) {
    const auto& r = scenario.radar;
    axes_.range_bin_m = r.range_bin_m();
    axes_.doppler_bins = doppler_.fft_size();
    axes_.velocity_bin_mps = r.prf_hz * r.wavelength() / (2.0 * static_cast<double>(doppler_.fft_size()));
}

RealMatrix ProcessingChain::range_doppler_map(const ComplexMatrix& cube, StageTimer* timer) const {
    ComplexMatrix compressed;
    {
        const MaybeScope s(timer, "pulse compression");
        compressed = compressor_.compress(cube);
    }
    ComplexMatrix range_major;
    {
        const MaybeScope s(timer, "corner turn");
        transpose(compressed, range_major);
    }
    if (scenario_.processing.mti != MtiType::None) {
        const MaybeScope s(timer, "MTI canceller");
        range_major = apply_mti(range_major, scenario_.processing.mti);
    }
    const MaybeScope s(timer, "Doppler FFT");
    return doppler_.process(range_major);
}

CpiResult ProcessingChain::process(const ComplexMatrix& cube, std::size_t cpi, StageTimer* timer) const {
    CpiResult result;
    result.rd_map = range_doppler_map(cube, timer);
    {
        const MaybeScope s(timer, scenario_.processing.cfar == CfarType::CellAveraging ? "CA-CFAR" : "OS-CFAR");
        result.cfar = cfar_.detect(result.rd_map);
        const double notch = scenario_.processing.clutter_notch_mps;
        if (notch > 0.0) {
            for (std::size_t d = 0; d < result.rd_map.cols(); ++d) {
                if (std::abs(axes_.velocity(static_cast<double>(d))) >= notch) {
                    continue;
                }
                for (std::size_t r = 0; r < result.rd_map.rows(); ++r) {
                    result.cfar.detection_count -= result.cfar.detections(r, d);
                    result.cfar.detections(r, d) = 0;
                }
            }
        }
    }
    {
        const MaybeScope s(timer, "clustering");
        result.detections = cluster_detections(result.rd_map, result.cfar, axes_, cpi);
    }
    return result;
}

SimulationResult run_simulation(const Scenario& scenario, StageTimer& timer) {
    const SignalModel model(scenario);
    const ProcessingChain chain(scenario);
    Tracker tracker(scenario.tracker, scenario.cpi_start_s(1));

    SimulationResult result;
    result.axes = chain.axes();
    for (std::size_t cpi = 0; cpi < scenario.cpis; ++cpi) {
        ComplexMatrix cube;
        {
            const auto s = timer.scope("signal generation");
            cube = model.generate(cpi);
        }
        result.raw_samples += cube.size();
        CpiResult cpi_result = chain.process(cube, cpi, &timer);
        {
            const auto s = timer.scope("tracking");
            tracker.update(cpi_result.detections, cpi);
        }
        result.detections.insert(result.detections.end(), cpi_result.detections.begin(),
                                 cpi_result.detections.end());
        result.cpis.push_back(std::move(cpi_result));
    }
    result.tracks = tracker.tracks();
    return result;
}

}  // namespace radar
