cmake_minimum_required(VERSION 3.15)
project(sudoku_gen)

include(FetchContent)

FetchContent_Declare(
    minisat
    URL http://minisat.se/downloads/MiniSat_v1.14.2006-Aug-29.src.zip
    URL_HASH MD5=69787b2942a8d2b2ab4eefe22d3dd262
)

FetchContent_MakeAvailable(minisat)

add_library(minisat STATIC)
target_sources(minisat PRIVATE
    "${minisat_SOURCE_DIR}/Solver.C"
)
target_compile_features(minisat PRIVATE cxx_std_14)
target_include_directories(minisat PUBLIC "${minisat_SOURCE_DIR}")

add_executable(sudoku_gen)
target_sources(sudoku_gen PRIVATE
    sudoku.cpp
    sudokuGen.cpp
)
target_compile_features(sudoku_gen PRIVATE cxx_std_20)
set_target_properties(sudoku_gen PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
target_link_libraries(sudoku_gen PRIVATE minisat)

if (NOT MSVC)
    #target_compile_options(sudoku_gen PRIVATE -Wall -Wextra -Wpedantic -Wconversion)
endif()