#pragma once

#include "radar/window.hpp"

#include <cstdint>
#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace radar {

enum class Swerling : std::uint8_t { Zero, One };
enum class MtiType : std::uint8_t { None, TwoPulse, ThreePulse };
enum class CfarType : std::uint8_t { CellAveraging, OrderedStatistic };

struct RadarParameters {
    double carrier_hz = 9.5e9;
    double bandwidth_hz = 5e6;
    double pulse_width_s = 10e-6;
    double prf_hz = 10e3;
    std::size_t pulses = 64;
    double sample_rate_hz = 10e6;
    double peak_power_w = 2e3;
    double antenna_gain_db = 32.0;
    double noise_figure_db = 3.0;
    double system_loss_db = 4.0;
    double max_range_m = 12e3;

    [[nodiscard]] double wavelength() const noexcept;
    [[nodiscard]] double pri() const noexcept { return 1.0 / prf_hz; }
    [[nodiscard]] std::size_t pulse_samples() const noexcept;
    [[nodiscard]] std::size_t range_bins() const noexcept;
    [[nodiscard]] std::size_t fast_time_samples() const noexcept;
    [[nodiscard]] double range_bin_m() const noexcept;
    [[nodiscard]] double range_resolution_m() const noexcept;
    [[nodiscard]] double unambiguous_range_m() const noexcept;
    [[nodiscard]] double max_unambiguous_velocity_mps() const noexcept;
    [[nodiscard]] double noise_power_w() const noexcept;
};

struct ProcessingParameters {
    WindowSpec range_window{WindowType::Taylor, 4, -35.0};
    WindowSpec doppler_window{WindowType::Hamming, 4, -35.0};
    MtiType mti = MtiType::None;
    CfarType cfar = CfarType::CellAveraging;
    double cfar_pfa = 1e-6;
    std::size_t guard_range = 2;
    std::size_t guard_doppler = 1;
    std::size_t train_range = 8;
    std::size_t train_doppler = 4;
    double os_rank = 0.75;          ///< OS-CFAR order statistic as a fraction of training cells.
    std::size_t doppler_fft = 0;    ///< 0 selects the next power of two above the pulse count.
    double clutter_notch_mps = 0.0;  ///< Discard CFAR hits with |velocity| below this (zero-velocity notch).
};

struct ClutterParameters {
    bool enabled = false;
    double cnr_db = 30.0;              ///< Clutter-to-noise ratio per raw sample.
    double min_range_m = 0.0;
    double max_range_m = 5e3;
    double spectral_width_mps = 0.25;  ///< Standard deviation of the clutter velocity spectrum.
};

struct TrackerParameters {
    std::size_t confirm_hits = 3;       ///< M in the M-of-N confirmation rule.
    std::size_t confirm_window = 4;     ///< N in the M-of-N confirmation rule.
    std::size_t max_misses = 3;         ///< Consecutive misses before a track is deleted.
    double range_gate_m = 60.0;
    double velocity_gate_mps = 4.0;
    double alpha = 0.6;
    double beta = 0.2;
};

struct MonteCarloParameters {
    double snr_min_db = 0.0;
    double snr_max_db = 20.0;
    double snr_step_db = 1.0;
    std::size_t trials = 200;
    std::size_t pfa_trials = 20;
    Swerling swerling = Swerling::One;
    double range_m = 0.0;      ///< 0 places the test target in the middle of the range window.
    /// Snapped to the nearest Doppler bin centre. Unset selects half the unambiguous velocity
    /// (Doppler = PRF / 4), clear of the zero-Doppler MTI notch.
    std::optional<double> velocity_mps;
};

struct Target {
    double range_m = 0.0;
    double velocity_mps = 0.0;  ///< Radial velocity, positive when closing.
    double rcs_m2 = 1.0;
    Swerling swerling = Swerling::Zero;
};

struct Scenario {
    std::string name = "scenario";
    RadarParameters radar;
    ProcessingParameters processing;
    ClutterParameters clutter;
    TrackerParameters tracker;
    MonteCarloParameters monte_carlo;
    std::vector<Target> targets;
    std::size_t cpis = 1;
    double cpi_interval_s = 0.0;  ///< Revisit time between CPI starts; 0 means back-to-back CPIs.
    std::uint64_t seed = 1;
    bool noise = true;

    /// Throws std::invalid_argument describing the first inconsistent parameter.
    void validate() const;

    [[nodiscard]] double cpi_duration_s() const noexcept;
    [[nodiscard]] double cpi_start_s(std::size_t cpi) const noexcept;
};

/// A parsed INI document. Repeated sections (e.g. several `[target]`) are kept in order.
struct IniSection {
    std::string name;
    std::map<std::string, std::string, std::less<>> values;
    std::size_t line = 0;
};

[[nodiscard]] std::vector<IniSection> parse_ini(std::string_view text);
[[nodiscard]] Scenario scenario_from_ini(std::string_view text);
[[nodiscard]] Scenario load_scenario(const std::filesystem::path& path);

[[nodiscard]] std::string_view to_string(MtiType type) noexcept;
[[nodiscard]] std::string_view to_string(CfarType type) noexcept;
[[nodiscard]] std::string_view to_string(Swerling model) noexcept;

}  // namespace radar
