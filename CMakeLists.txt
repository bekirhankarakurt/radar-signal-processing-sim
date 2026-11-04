cmake_minimum_required(VERSION 3.16)

project(radar_signal_processing_sim
    VERSION 1.0.0
    DESCRIPTION "Pulse-Doppler radar signal processing simulator"
    LANGUAGES CXX)

option(RADAR_ENABLE_OPENMP "Parallelise the processing chain with OpenMP" ON)
option(RADAR_BUILD_TESTS "Build the unit tests" ON)
option(RADAR_ENABLE_SANITIZERS "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(RADAR_ENABLE_COVERAGE "Instrument the build for gcov coverage" OFF)
option(RADAR_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(RADAR_NATIVE "Optimise for the build machine (-march=native)" OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(radar_options INTERFACE)
target_compile_options(radar_options INTERFACE
    -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
    $<$<BOOL:${RADAR_WARNINGS_AS_ERRORS}>:-Werror>
    $<$<CONFIG:Release>:-O3>
    $<$<BOOL:${RADAR_NATIVE}>:-march=native>)

if(RADAR_ENABLE_SANITIZERS)
    target_compile_options(radar_options INTERFACE
        -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all)
    target_link_options(radar_options INTERFACE -fsanitize=address,undefined)
endif()

if(RADAR_ENABLE_COVERAGE)
    target_compile_options(radar_options INTERFACE --coverage -O0 -g)
    target_link_options(radar_options INTERFACE --coverage)
endif()

add_library(radar_core
    src/cfar.cpp
    src/chain.cpp
    src/clustering.cpp
    src/config.cpp
    src/doppler.cpp
    src/fft.cpp
    src/monte_carlo.cpp
    src/mti.cpp
    src/output.cpp
    src/pulse_compression.cpp
    src/signal_model.cpp
    src/statistics.cpp
    src/tracker.cpp
    src/waveform.cpp
    src/window.cpp)
target_include_directories(radar_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(radar_core PUBLIC radar_options)
target_compile_definitions(radar_core PUBLIC RADAR_VERSION="${PROJECT_VERSION}")

if(RADAR_ENABLE_OPENMP)
    find_package(OpenMP COMPONENTS CXX)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(radar_core PUBLIC OpenMP::OpenMP_CXX)
        message(STATUS "OpenMP enabled")
    else()
        message(STATUS "OpenMP not found; building single-threaded")
        target_compile_options(radar_options INTERFACE -Wno-unknown-pragmas)
    endif()
else()
    target_compile_options(radar_options INTERFACE -Wno-unknown-pragmas)
endif()

add_executable(radar_sim app/main.cpp)
target_link_libraries(radar_sim PRIVATE radar_core)

install(TARGETS radar_sim RUNTIME DESTINATION bin)

if(RADAR_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
