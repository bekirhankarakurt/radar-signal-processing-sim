#include "radar/window.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>
#include <string>

namespace radar {
namespace {

/// Taylor window (Carrara/Goodman/Majewski formulation), normalised to a peak of 1.
RealVector taylor_window(std::size_t length, int nbar, double sidelobe_db) {
    if (nbar < 1) {
        throw std::invalid_argument("Taylor nbar must be >= 1");
    }
    if (sidelobe_db >= 0.0) {
        throw std::invalid_argument("Taylor sidelobe level must be negative dB");
    }
    const double eta = std::pow(10.0, -sidelobe_db / 20.0);
    const double a = std::acosh(eta) / kPi;
    const auto nb = static_cast<double>(nbar);
    const double sigma2 = nb * nb / (a * a + (nb - 0.5) * (nb - 0.5));

    std::vector<double> fm(static_cast<std::size_t>(nbar - 1));
    for (int m = 1; m < nbar; ++m) {
        const auto md = static_cast<double>(m);
        double num = 1.0;
        double den = 1.0;
        for (int n = 1; n < nbar; ++n) {
            const auto nd = static_cast<double>(n);
            num *= 1.0 - (md * md / sigma2) / (a * a + (nd - 0.5) * (nd - 0.5));
            if (n != m) {
                den *= 1.0 - (md * md) / (nd * nd);
            }
        }
        const double sign = (m % 2 == 1) ? 1.0 : -1.0;
        fm[static_cast<std::size_t>(m - 1)] = sign * num / (2.0 * den);
    }

    RealVector w(length);
    const auto len = static_cast<double>(length);
    double peak = 0.0;
    for (std::size_t i = 0; i < length; ++i) {
        const double x = (static_cast<double>(i) + 0.5) / len - 0.5;
        double v = 1.0;
        for (std::size_t m = 0; m < fm.size(); ++m) {
            v += 2.0 * fm[m] * std::cos(2.0 * kPi * static_cast<double>(m + 1) * x);
        }
        w[i] = v;
        peak = std::max(peak, v);
    }
    for (auto& v : w) {
        v /= peak;
    }
    return w;
}

RealVector cosine_window(std::size_t length, double a0) {
    RealVector w(length);
    if (length == 1) {
        w[0] = 1.0;
        return w;
    }
    const auto denom = static_cast<double>(length - 1);
    for (std::size_t i = 0; i < length; ++i) {
        w[i] = a0 - (1.0 - a0) * std::cos(2.0 * kPi * static_cast<double>(i) / denom);
    }
    return w;
}

}  // namespace

RealVector make_window(const WindowSpec& spec, std::size_t length) {
    if (length == 0) {
        throw std::invalid_argument("window length must be positive");
    }
    switch (spec.type) {
        case WindowType::Rectangular:
            return RealVector(length, 1.0);  // NOLINT(modernize-return-braced-init-list): count, value
        case WindowType::Hann:
            return cosine_window(length, 0.5);
        case WindowType::Hamming:
            return cosine_window(length, 0.54);
        case WindowType::Taylor:
            return taylor_window(length, spec.taylor_nbar, spec.taylor_sidelobe_db);
    }
    throw std::invalid_argument("unknown window type");
}

double window_snr_gain(const RealVector& w) noexcept {
    double sum = 0.0;
    double sum_sq = 0.0;
    for (const double v : w) {
        sum += v;
        sum_sq += v * v;
    }
    return sum_sq > 0.0 ? sum * sum / sum_sq : 0.0;
}

std::optional<WindowType> parse_window_type(std::string_view name) noexcept {
    std::string lower;
    lower.reserve(name.size());
    for (const char c : name) {
        lower.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    if (lower == "none" || lower == "rect" || lower == "rectangular") {
        return WindowType::Rectangular;
    }
    if (lower == "hann" || lower == "hanning") {
        return WindowType::Hann;
    }
    if (lower == "hamming") {
        return WindowType::Hamming;
    }
    if (lower == "taylor") {
        return WindowType::Taylor;
    }
    return std::nullopt;
}

std::string_view to_string(WindowType type) noexcept {
    switch (type) {
        case WindowType::Rectangular:
            return "rectangular";
        case WindowType::Hann:
            return "hann";
        case WindowType::Hamming:
            return "hamming";
        case WindowType::Taylor:
            return "taylor";
    }
    return "unknown";
}

}  // namespace radar
