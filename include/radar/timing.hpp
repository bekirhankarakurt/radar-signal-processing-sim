#pragma once

#include <chrono>
#include <string>
#include <utility>
#include <vector>

namespace radar {

/// Accumulates wall-clock time per named processing stage, preserving first-seen order.
class StageTimer {
public:
    struct Entry {
        std::string name;
        double total_ms = 0.0;
        std::size_t calls = 0;
    };

    /// RAII scope that adds its lifetime to a stage when destroyed.
    class Scope {
    public:
        Scope(StageTimer& timer, std::string name)
            : timer_(timer), name_(std::move(name)), start_(std::chrono::steady_clock::now()) {}
        Scope(const Scope&) = delete;
        Scope& operator=(const Scope&) = delete;
        Scope(Scope&&) = delete;
        Scope& operator=(Scope&&) = delete;
        ~Scope() {
            const auto elapsed = std::chrono::steady_clock::now() - start_;
            timer_.add(name_, std::chrono::duration<double, std::milli>(elapsed).count());
        }

    private:
        StageTimer& timer_;
        std::string name_;
        std::chrono::steady_clock::time_point start_;
    };

    [[nodiscard]] Scope scope(std::string name) { return Scope(*this, std::move(name)); }

    void add(const std::string& name, double ms) {
        for (auto& e : entries_) {
            if (e.name == name) {
                e.total_ms += ms;
                ++e.calls;
                return;
            }
        }
        entries_.push_back({name, ms, 1});
    }

    [[nodiscard]] const std::vector<Entry>& entries() const noexcept { return entries_; }

    [[nodiscard]] double total_ms() const noexcept {
        double t = 0.0;
        for (const auto& e : entries_) {
            t += e.total_ms;
        }
        return t;
    }

private:
    std::vector<Entry> entries_;
};

}  // namespace radar
