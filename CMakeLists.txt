cmake_minimum_required(VERSION 3.14)
project(SchwartzTwoFactorADE VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
    add_compile_options(/O2 /W4)
else()
    add_compile_options(-O3 -Wall -Wextra -Wpedantic)
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)

set(CORE_SOURCES
    src/ADESolver.cpp
    src/MonteCarloLSMC.cpp
)

# Main benchmark and paper replication executable
add_executable(schwartz_ade_suite
    src/main.cpp
    ${CORE_SOURCES}
)

# Automated unit test executable
add_executable(schwartz_unit_tests
    tests/unit_tests.cpp
    ${CORE_SOURCES}
)

enable_testing()
add_test(NAME UnitTests COMMAND schwartz_unit_tests)
