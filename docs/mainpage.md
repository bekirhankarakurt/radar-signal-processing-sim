# radar-signal-processing-sim API {#mainpage}

C++17 pulse-Doppler radar simulator: scenario model, LFM pulse compression, MTI,
Doppler processing, 2-D CA/OS-CFAR, clustering, tracking and Monte Carlo
detection analysis. No third-party dependencies.

## Processing chain

| Stage | Entry point |
|---|---|
| Scenario file | radar::load_scenario(), radar::Scenario |
| Raw data cube | radar::SignalModel (targets, thermal noise, Gaussian-spectrum clutter) |
| Pulse compression | radar::PulseCompressor (FFT matched filter, Hamming/Taylor taper) |
| Corner turn | radar::transpose() |
| MTI | radar::apply_mti() |
| Doppler | radar::DopplerProcessor |
| Detection | radar::CfarDetector (CA via summed-area table, OS via `nth_element`) |
| Clustering | radar::cluster_detections() |
| Tracking | radar::Tracker (M-of-N confirmation, alpha-beta smoothing) |
| Whole chain | radar::ProcessingChain, radar::run_simulation() |
| Pd / Pfa study | radar::run_monte_carlo(), radar/statistics.hpp |
| Outputs | radar/output.hpp (PPM, CSV, HTML viewer) |

Building blocks: radar::FftPlan (radix-2), radar::make_window(),
radar::Matrix (row-major, contiguous rows), radar::Xoshiro256.

Source and README: https://github.com/bekirhankarakurt/radar-signal-processing-sim
