#include "radar/tracker.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <tuple>

namespace radar {

Tracker::Tracker(const TrackerParameters& params, double cpi_duration_s) : params_(params), dt_(cpi_duration_s) {}

void Tracker::update(const std::vector<Detection>& detections, std::size_t cpi) {
    // Predict: closing velocity is positive, so range decreases.
    std::vector<double> predicted(tracks_.size());
    for (std::size_t i = 0; i < tracks_.size(); ++i) {
        predicted[i] = tracks_[i].range_m - tracks_[i].velocity_mps * dt_;
    }

    // Global nearest neighbour: all gated pairs sorted by normalised distance, assigned greedily.
    std::vector<std::tuple<double, std::size_t, std::size_t>> pairs;
    for (std::size_t i = 0; i < tracks_.size(); ++i) {
        if (tracks_[i].status == TrackStatus::Deleted) {
            continue;
        }
        for (std::size_t j = 0; j < detections.size(); ++j) {
            const double dr = (detections[j].range_m - predicted[i]) / params_.range_gate_m;
            const double dv = (detections[j].velocity_mps - tracks_[i].velocity_mps) / params_.velocity_gate_mps;
            if (std::abs(dr) <= 1.0 && std::abs(dv) <= 1.0) {
                pairs.emplace_back(dr * dr + dv * dv, i, j);
            }
        }
    }
    std::sort(pairs.begin(), pairs.end());
    std::vector<bool> track_used(tracks_.size(), false);
    std::vector<bool> det_used(detections.size(), false);
    for (const auto& [cost, i, j] : pairs) {
        if (track_used[i] || det_used[j]) {
            continue;
        }
        track_used[i] = true;
        det_used[j] = true;
        auto& t = tracks_[i];
        const auto& d = detections[j];
        t.range_m = predicted[i] + params_.alpha * (d.range_m - predicted[i]);
        t.velocity_mps += params_.beta * (d.velocity_mps - t.velocity_mps);
        t.snr_db = d.snr_db;
        ++t.hits;
        t.misses = 0;
        t.last_cpi = cpi;
        t.history.push_back(true);
    }

    for (std::size_t i = 0; i < tracks_.size(); ++i) {
        auto& t = tracks_[i];
        if (t.status == TrackStatus::Deleted || track_used[i]) {
            continue;
        }
        t.range_m = predicted[i];
        ++t.misses;
        t.history.push_back(false);
    }

    for (auto& t : tracks_) {
        if (t.status == TrackStatus::Deleted) {
            continue;
        }
        if (t.history.size() > params_.confirm_window) {
            t.history.erase(t.history.begin());
        }
        const auto recent_hits = static_cast<std::size_t>(std::count(t.history.begin(), t.history.end(), true));
        if (t.status == TrackStatus::Tentative) {
            if (recent_hits >= params_.confirm_hits) {
                t.status = TrackStatus::Confirmed;
            } else if (t.history.size() >= params_.confirm_window ||
                       recent_hits + (params_.confirm_window - t.history.size()) < params_.confirm_hits) {
                t.status = TrackStatus::Deleted;
            }
        } else if (t.misses >= params_.max_misses) {
            t.status = TrackStatus::Deleted;
        }
    }

    for (std::size_t j = 0; j < detections.size(); ++j) {
        if (det_used[j]) {
            continue;
        }
        Track t;
        t.id = next_id_++;
        t.range_m = detections[j].range_m;
        t.velocity_mps = detections[j].velocity_mps;
        t.snr_db = detections[j].snr_db;
        t.hits = 1;
        t.first_cpi = cpi;
        t.last_cpi = cpi;
        t.history.push_back(true);
        if (params_.confirm_hits <= 1) {
            t.status = TrackStatus::Confirmed;
        }
        tracks_.push_back(std::move(t));
    }
}

std::vector<Track> Tracker::confirmed() const {
    std::vector<Track> out;
    std::copy_if(tracks_.begin(), tracks_.end(), std::back_inserter(out),
                 [](const Track& t) { return t.status == TrackStatus::Confirmed; });
    return out;
}

std::string_view to_string(TrackStatus status) noexcept {
    switch (status) {
        case TrackStatus::Tentative:
            return "tentative";
        case TrackStatus::Confirmed:
            return "confirmed";
        case TrackStatus::Deleted:
            return "deleted";
    }
    return "unknown";
}

}  // namespace radar
