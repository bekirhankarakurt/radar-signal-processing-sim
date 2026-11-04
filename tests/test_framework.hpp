#pragma once

// Minimal dependency-free test harness: TEST registers a function, CHECK* record
// failures without aborting the test, and main() runs everything and returns
// non-zero if any check failed.

#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace testing {

struct TestCase {
    std::string name;
    std::function<void()> body;
};

inline std::vector<TestCase>& registry() {
    static std::vector<TestCase> tests;
    return tests;
}

inline int& failures() {
    static int count = 0;
    return count;
}

struct Registrar {
    Registrar(const char* name, std::function<void()> body) { registry().push_back({name, std::move(body)}); }
};

inline void report(const char* file, int line, const std::string& message) {
    ++failures();
    std::cerr << "  " << file << ":" << line << ": " << message << "\n";
}

inline int run_all() {
    int failed_tests = 0;
    for (const auto& test : registry()) {
        const int before = failures();
        try {
            test.body();
        } catch (const std::exception& e) {
            report(__FILE__, __LINE__, std::string("unexpected exception: ") + e.what());
        }
        const bool ok = failures() == before;
        failed_tests += ok ? 0 : 1;
        std::cout << (ok ? "[ PASS ] " : "[ FAIL ] ") << test.name << "\n";
    }
    std::cout << registry().size() - static_cast<std::size_t>(failed_tests) << "/" << registry().size()
              << " tests passed\n";
    return failed_tests == 0 ? 0 : 1;
}

}  // namespace testing

#define TEST_CONCAT_INNER(a, b) a##b
#define TEST_CONCAT(a, b) TEST_CONCAT_INNER(a, b)

#define TEST(name)                                                                                   \
    static void name();                                                                              \
    static const testing::Registrar TEST_CONCAT(registrar_, name){#name, name};                      \
    static void name()

#define CHECK(cond)                                                                                  \
    do {                                                                                             \
        if (!(cond)) {                                                                               \
            testing::report(__FILE__, __LINE__, "CHECK(" #cond ") failed");                          \
        }                                                                                            \
    } while (false)

#define CHECK_NEAR(actual, expected, tolerance)                                                      \
    do {                                                                                             \
        const double check_a = static_cast<double>(actual);                                          \
        const double check_e = static_cast<double>(expected);                                        \
        if (!(std::abs(check_a - check_e) <= static_cast<double>(tolerance))) {                      \
            testing::report(__FILE__, __LINE__,                                                      \
                            std::string(#actual " = ") + std::to_string(check_a) + ", expected " +   \
                                std::to_string(check_e) + " +/- " + std::to_string(tolerance));      \
        }                                                                                            \
    } while (false)

#define CHECK_THROWS(expr)                                                                           \
    do {                                                                                             \
        bool check_threw = false;                                                                    \
        try {                                                                                        \
            static_cast<void>(expr);                                                                 \
        } catch (const std::exception&) {                                                            \
            check_threw = true;                                                                      \
        }                                                                                            \
        if (!check_threw) {                                                                          \
            testing::report(__FILE__, __LINE__, "expected exception from " #expr);                   \
        }                                                                                            \
    } while (false)

#define TEST_MAIN()                                                                                  \
    int main() { return testing::run_all(); }
