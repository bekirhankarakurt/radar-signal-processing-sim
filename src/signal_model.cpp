#include "radar/signal_model.hpp"

#include "radar/complex_math.hpp"
#include "radar/random.hpp"
#include "radar/waveform.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace radar {

std::uint64_t mix_seed(std::uint64_t seed, std::uint64_t a, std::uint64_t b) noexcept {
    auto splitmix = [](std::uint64_t x) {
        x += 0x9E37'79B9'7F4A'7C15ULL;
        x = (x ^ (x >> 30U)) * 0xBF58'476D'1CE4'E5B9ULL;
        x = (x ^ (x >> 27U)) * 0x94D0'49BB'1331'11EBULL;
        return x ^ (x >> 31U);
    };
    return splitmix(splitmix(splitmix(seed) ^ a) ^ b);
}

RealMatrix cholesky(const RealMatrix& a) {
    const std::size_t n = a.rows();
    if (a.cols() != n) {
        throw std::invalid_argument("cholesky: matrix must be square");
    }
    RealMatrix l(n, n, 0.0);
    for (std::size_t j = 0; j < n; ++j) {
        double diag = a(j, j);
        for (std::size_t k = 0; k < j; ++k) {
            diag -= l(j, k) * l(j, k);
        }
        if (diag <= 0.0) {
            throw std::invalid_argument("cholesky: matrix is not positive definite");
        }
        l(j, j) = std::sqrt(diag);
        for (std::size_t i = j + 1; i < n; ++i) {
            double v = a(i, j);
            for (std::size_t k = 0; k < j; ++k) {
                v -= l(i, k) * l(j, k);
            }
            l(i, j) = v / l(j, j);
        }
    }
    return l;
}

double clutter_correlation(double spectral_width_mps, double wavelength_m, double pri_s, std::size_t lag) noexcept {
    const double sigma_f = 2.0 * spectral_width_mps / wavelength_m;
    const double x = kPi * sigma_f * pri_s * static_cast<double>(lag);
    return std::exp(-2.0 * x * x);
}

double received_power_w(const RadarParameters& radar, double range_m, double rcs_m2) noexcept {
    const double gain = db_to_linear(radar.antenna_gain_db);
    const double lambda = radar.wavelength();
    const double loss = db_to_linear(radar.system_loss_db);
    const double four_pi_cubed = std::pow(4.0 * kPi, 3.0);
    const double r = std::max(range_m, 1.0);
    return radar.peak_power_w * gain * gain * lambda * lambda * rcs_m2 / (four_pi_cubed * std::pow(r, 4.0) * loss);
}

double single_sample_snr(const RadarParameters& radar, double range_m, double rcs_m2) noexcept {
    return received_power_w(radar, range_m, rcs_m2) / radar.noise_power_w();
}

SignalModel::SignalModel(const Scenario& scenario)
    : scenario_(scenario),
      chirp_(make_lfm_chirp(scenario.radar.bandwidth_hz, scenario.radar.pulse_width_s, scenario.radar.sample_rate_hz)),
      clutter_plan_(next_power_of_two(scenario.radar.fast_time_samples() + chirp_.size() - 1)),
      chirp_spectrum_(clutter_plan_.size(), Complex{0.0, 0.0}) {
    std::copy(chirp_.begin(), chirp_.end(), chirp_spectrum_.begin());
    clutter_plan_.execute(chirp_spectrum_, FftDirection::Forward);
    if (scenario.clutter.enabled) {
        // Exact Gaussian-spectrum clutter: colour white noise with the Cholesky factor of the
        // slow-time correlation matrix. A -80 dB diagonal load keeps the factorisation stable
        // and sets the spectral floor of the clutter model.
        const std::size_t pulses = scenario.radar.pulses;
        RealMatrix correlation(pulses, pulses);
        for (std::size_t i = 0; i < pulses; ++i) {
            for (std::size_t j = 0; j < pulses; ++j) {
                const std::size_t lag = i > j ? i - j : j - i;
                correlation(i, j) = clutter_correlation(scenario.clutter.spectral_width_mps,
                                                        scenario.radar.wavelength(), scenario.radar.pri(), lag) +
                                    (i == j ? 1e-8 : 0.0);
            }
        }
        clutter_coloring_ = cholesky(correlation);
    }
}

std::vector<TargetEcho> SignalModel::echoes_for_cpi(std::size_t cpi) const {
    const auto& radar = scenario_.radar;
    const double cpi_start = scenario_.cpi_start_s(cpi);
    std::vector<TargetEcho> echoes;
    echoes.reserve(scenario_.targets.size());
    for (std::size_t i = 0; i < scenario_.targets.size(); ++i) {
        const auto& t = scenario_.targets[i];
        Xoshiro256 rng(mix_seed(scenario_.seed, 0xEC40ULL + i, cpi));
        // Swerling 1: exponentially distributed RCS, constant within a CPI, independent between CPIs.
        const double rcs = t.swerling == Swerling::One ? t.rcs_m2 * rng.exponential() : t.rcs_m2;
        const double range = t.range_m - t.velocity_mps * cpi_start;
        const double amplitude = std::sqrt(single_sample_snr(radar, range, rcs));
        echoes.push_back({range, t.velocity_mps, amplitude * cis(2.0 * kPi * rng.uniform())});
    }
    return echoes;
}

ComplexMatrix SignalModel::generate(std::size_t cpi) const {
    return generate(echoes_for_cpi(cpi), cpi, scenario_.noise);
}

ComplexMatrix SignalModel::generate(const std::vector<TargetEcho>& echoes, std::uint64_t stream,
                                    bool add_noise) const {
    const auto& radar = scenario_.radar;
    const std::size_t pulses = radar.pulses;
    const std::size_t samples = radar.fast_time_samples();
    const double fs = radar.sample_rate_hz;
    const double lambda = radar.wavelength();
    ComplexMatrix cube(pulses, samples);

    const auto pulse_count = static_cast<std::ptrdiff_t>(pulses);
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t pi = 0; pi < pulse_count; ++pi) {
        const auto m = static_cast<std::size_t>(pi);
        Complex* row = cube.row(m);
        if (add_noise) {
            Xoshiro256 rng(mix_seed(scenario_.seed, stream, 0x0015E000ULL + m));
            for (std::size_t n = 0; n < samples; ++n) {
                row[n] = rng.complex_gaussian();
            }
        }
        const double t_pulse = static_cast<double>(m) * radar.pri();
        for (const auto& echo : echoes) {
            // Stop-and-hop model: the target is frozen during each pulse and moves between pulses.
            const double range = echo.range_m - echo.velocity_mps * t_pulse;
            const double delay = 2.0 * range / kSpeedOfLight;
            const Complex carrier = echo.amplitude * cis(-4.0 * kPi * range / lambda);
            const auto first = static_cast<std::size_t>(std::max(0.0, std::ceil(delay * fs)));
            const auto last = static_cast<std::size_t>(std::ceil((delay + radar.pulse_width_s) * fs));
            for (std::size_t n = first; n < std::min(last, samples); ++n) {
                const double t = static_cast<double>(n) / fs - delay;
                row[n] += carrier * lfm_sample(t, radar.bandwidth_hz, radar.pulse_width_s);
            }
        }
    }
    if (scenario_.clutter.enabled) {
        add_clutter(cube, stream);
    }
    return cube;
}

void SignalModel::add_clutter(ComplexMatrix& cube, std::uint64_t stream) const {
    const auto& radar = scenario_.radar;
    const auto& clutter = scenario_.clutter;
    const std::size_t pulses = cube.rows();
    const std::size_t samples = cube.cols();
    const double bin = radar.range_bin_m();
    const auto first = static_cast<std::size_t>(std::ceil(clutter.min_range_m / bin));
    const std::size_t last = std::min(samples, static_cast<std::size_t>(std::floor(clutter.max_range_m / bin)) + 1);
    if (first >= last) {
        return;
    }

    // Reflectivity per range cell with the configured Doppler spectrum (see constructor).
    // Each echo sample sums `chirp_.size()` reflectivity cells, so this gives the configured CNR per sample.
    const double variance = db_to_linear(clutter.cnr_db) / static_cast<double>(chirp_.size());
    const std::size_t cells = last - first;
    ComplexMatrix reflectivity(cells, pulses);  // cell-major so each slow-time sequence is contiguous

    const auto cell_count = static_cast<std::ptrdiff_t>(cells);
#pragma omp parallel
    {
        ComplexVector white(pulses);
#pragma omp for schedule(static)
        for (std::ptrdiff_t ci = 0; ci < cell_count; ++ci) {
            const auto c = static_cast<std::size_t>(ci);
            Xoshiro256 rng(mix_seed(scenario_.seed, stream, 0xC1A7'0000ULL + c));
            for (auto& w : white) {
                w = rng.complex_gaussian(variance);
            }
            Complex* seq = reflectivity.row(c);
            for (std::size_t m = 0; m < pulses; ++m) {
                const double* l = clutter_coloring_.row(m);
                Complex acc{0.0, 0.0};
                for (std::size_t j = 0; j <= m; ++j) {
                    acc += l[j] * white[j];
                }
                seq[m] = acc;
            }
        }
    }

    // Echo = reflectivity convolved with the transmitted chirp, done per pulse in the frequency domain.
    const std::size_t nfft = clutter_plan_.size();
    const auto pulse_count = static_cast<std::ptrdiff_t>(pulses);
#pragma omp parallel
    {
        ComplexVector buffer(nfft);
#pragma omp for schedule(static)
        for (std::ptrdiff_t pi = 0; pi < pulse_count; ++pi) {
            const auto m = static_cast<std::size_t>(pi);
            std::fill(buffer.begin(), buffer.end(), Complex{0.0, 0.0});
            for (std::size_t c = 0; c < cells; ++c) {
                buffer[first + c] = reflectivity(c, m);
            }
            clutter_plan_.execute(buffer.data(), FftDirection::Forward);
            for (std::size_t k = 0; k < nfft; ++k) {
                buffer[k] *= chirp_spectrum_[k];
            }
            clutter_plan_.execute(buffer.data(), FftDirection::Inverse);
            Complex* row = cube.row(m);
            for (std::size_t n = first; n < samples; ++n) {
                row[n] += buffer[n];
            }
        }
    }
}

}  // namespace radar
