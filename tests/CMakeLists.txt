set(RADAR_TESTS
    test_fft
    test_window
    test_config
    test_pulse_compression
    test_doppler
    test_cfar
    test_clustering
    test_tracker
    test_statistics
    test_chain
    test_monte_carlo
    test_output)

foreach(name IN LISTS RADAR_TESTS)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE radar_core)
    target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_definitions(${name} PRIVATE RADAR_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
    add_test(NAME ${name} COMMAND ${name})
    set_tests_properties(${name} PROPERTIES TIMEOUT 300)
endforeach()
