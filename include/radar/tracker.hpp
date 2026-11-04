#pragma once

#include "radar/clustering.hpp"
#include "radar/config.hpp"

#include <vector>

namespace radar {

enum class TrackStatus { Tentative, Confirmed, Deleted };

struct Track {
    std::size_t id = 0;
    TrackStatus status = TrackStatus::Tentative;
    double range_m = 0.0;
    double velocity_mps = 0.0;  ///< Positive when closing.
    double snr_db = 0.0;
    std::size_t hits = 0;
    std::size_t misses = 0;
    std::size_t first_cpi = 0;
    std::size_t last_cpi = 0;
    std::vector<bool> history;  ///< Hit/miss record for the M-of-N rule, newest last.
};

/// Track-before-display list: nearest-neighbour association with range and velocity
/// gates, alpha-beta smoothing of range, M-of-N confirmation and miss-count deletion.
/// Only confirmed tracks are meant for display.
class Tracker {
public:
    Tracker(const TrackerParameters& params, double cpi_duration_s);

    void update(const std::vector<Detection>& detections, std::size_t cpi);

    [[nodiscard]] const std::vector<Track>& tracks() const noexcept { return tracks_; }
    [[nodiscard]] std::vector<Track> confirmed() const;

private:
    TrackerParameters params_;
    double dt_;
    std::vector<Track> tracks_;
    std::size_t next_id_ = 1;
};

[[nodiscard]] std::string_view to_string(TrackStatus status) noexcept;

}  // namespace radar
