#pragma once

#include "radar/types.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string_view>

namespace radar {

enum class WindowType : std::uint8_t { Rectangular, Hann, Hamming, Taylor };

struct WindowSpec {
    WindowType type = WindowType::Rectangular;
    int taylor_nbar = 4;               ///< Number of nearly constant-level sidelobes.
    double taylor_sidelobe_db = -35.0;  ///< Design peak sidelobe level (negative dB).
};

[[nodiscard]] RealVector make_window(const WindowSpec& spec, std::size_t length);

/// Coherent gain (sum w)^2 / sum w^2: the SNR gain of a windowed coherent integration.
[[nodiscard]] double window_snr_gain(const RealVector& w) noexcept;

[[nodiscard]] std::optional<WindowType> parse_window_type(std::string_view name) noexcept;
[[nodiscard]] std::string_view to_string(WindowType type) noexcept;

}  // namespace radar
