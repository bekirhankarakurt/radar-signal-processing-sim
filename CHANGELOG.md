# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2026-11-04

### Added

- INI scenario format with radar, processing, clutter, tracker, Monte Carlo and repeated target sections; unknown keys and inconsistent parameters are rejected with line numbers.
- Signal model: baseband LFM chirp, stop-and-hop target motion, radar-range-equation amplitudes, Swerling 0/1 targets, complex white thermal noise and ground clutter with a Gaussian Doppler spectrum.
- Processing chain: FFT matched filter with Hamming/Hann/Taylor tapers, tiled corner turn, two- and three-pulse MTI cancellers, windowed Doppler FFT and a zero-velocity notch.
- 2-D CA-CFAR (summed-area table, O(1) per cell) and OS-CFAR with circular Doppler handling.
- Detection clustering with power-weighted centroids and an M-of-N track-before-display tracker.
- Monte Carlo mode for Pd vs SNR through the full chain and empirical Pfa, with closed-form Swerling 0/1 and CFAR references.
- Outputs: PPM range-Doppler image, CSV map/detections/tracks/timing/Pd curve and a self-contained HTML canvas viewer.
- Per-stage timing report and optional OpenMP parallelism with thread-count-independent results.
- Dependency-free unit tests (CTest), CI on GCC and Clang with ASan/UBSan, cppcheck, clang-tidy, CodeQL and coverage.
- Static Linux release archives for x86_64 and aarch64.

[Unreleased]: https://github.com/bekirhankarakurt/radar-signal-processing-sim/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/bekirhankarakurt/radar-signal-processing-sim/releases/tag/v1.0.0
