#pragma once

#include "radar/cfar.hpp"
#include "radar/clustering.hpp"
#include "radar/config.hpp"
#include "radar/doppler.hpp"
#include "radar/pulse_compression.hpp"
#include "radar/signal_model.hpp"
#include "radar/timing.hpp"
#include "radar/tracker.hpp"

#include <vector>

namespace radar {

struct CpiResult {
    RealMatrix rd_map;  ///< Range-Doppler power, noise normalised to 1.
    CfarResult cfar;
    std::vector<Detection> detections;
};

/// The full receive chain for one scenario:
/// pulse compression -> corner turn -> MTI -> Doppler FFT -> CFAR -> clustering.
/// All members are immutable after construction, so one chain can process
/// independent CPIs concurrently.
class ProcessingChain {
public:
    explicit ProcessingChain(const Scenario& scenario);

    [[nodiscard]] CpiResult process(const ComplexMatrix& cube, std::size_t cpi, StageTimer* timer) const;

    /// Range-Doppler map only (no detection), used by the Monte Carlo engine.
    [[nodiscard]] RealMatrix range_doppler_map(const ComplexMatrix& cube, StageTimer* timer) const;

    [[nodiscard]] const RangeDopplerAxes& axes() const noexcept { return axes_; }
    [[nodiscard]] const CfarDetector& cfar() const noexcept { return cfar_; }
    [[nodiscard]] const PulseCompressor& compressor() const noexcept { return compressor_; }
    [[nodiscard]] const DopplerProcessor& doppler() const noexcept { return doppler_; }

private:
    const Scenario& scenario_;
    PulseCompressor compressor_;
    DopplerProcessor doppler_;
    CfarDetector cfar_;
    RangeDopplerAxes axes_;
};

struct SimulationResult {
    std::vector<CpiResult> cpis;
    std::vector<Detection> detections;  ///< All CPIs, in order.
    std::vector<Track> tracks;          ///< Final track table (all states).
    RangeDopplerAxes axes;
    std::size_t raw_samples = 0;        ///< Complex samples processed across all CPIs.
};

[[nodiscard]] SimulationResult run_simulation(const Scenario& scenario, StageTimer& timer);

}  // namespace radar
