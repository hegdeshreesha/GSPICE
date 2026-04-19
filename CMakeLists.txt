cmake_minimum_required(VERSION 3.10)
project(GSPICE VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(include)

# Find and enable OpenMP
find_package(OpenMP REQUIRED)

# Core source files
set(GSPICE_CORE_SOURCES
    src/core/main.cpp
    src/core/parser.cpp
)

add_executable(gspice ${GSPICE_CORE_SOURCES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(gspice PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(gspice PRIVATE include)
