#pragma once

#include "radar/chain.hpp"
#include "radar/monte_carlo.hpp"
#include "radar/timing.hpp"

#include <filesystem>
#include <iosfwd>

namespace radar {

struct ImageOptions {
    double dynamic_range_db = 60.0;
    bool mark_detections = true;
};

/// Range-Doppler map as a binary PPM (P6): Doppler on x, range on y (near range at the bottom).
void write_ppm(const std::filesystem::path& path, const RealMatrix& rd_map, const std::vector<Detection>& detections,
               const ImageOptions& options);

void write_rd_csv(const std::filesystem::path& path, const RealMatrix& rd_map, const RangeDopplerAxes& axes);
void write_detections_csv(const std::filesystem::path& path, const std::vector<Detection>& detections);
void write_tracks_csv(const std::filesystem::path& path, const std::vector<Track>& tracks);
void write_timing_csv(const std::filesystem::path& path, const StageTimer& timer, std::size_t cpis);
void write_pd_csv(const std::filesystem::path& path, const MonteCarloResult& result);

/// Self-contained HTML page (no external resources) with a canvas range-Doppler viewer,
/// a CPI slider, hover read-out, detection overlay and the confirmed track table.
void write_html_viewer(const std::filesystem::path& path, const Scenario& scenario, const SimulationResult& result);

void print_timing(std::ostream& os, const StageTimer& timer, std::size_t cpis, std::size_t raw_samples);
void print_detections(std::ostream& os, const std::vector<Detection>& detections, std::size_t cpi);
void print_tracks(std::ostream& os, const std::vector<Track>& tracks);
void print_monte_carlo(std::ostream& os, const MonteCarloResult& result, const Scenario& scenario);

/// Base64 encoding used to embed binary map data in the HTML viewer.
[[nodiscard]] std::string base64_encode(const std::vector<std::uint8_t>& bytes);

}  // namespace radar
