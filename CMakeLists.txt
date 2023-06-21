cmake_minimum_required(VERSION 3.15)
project(sudoku_gen)

add_executable(sudoku_gen)
target_sources(sudoku_gen PRIVATE
    File.cpp
    Proof.cpp
    Solver.cpp
    sudoku.cpp
    sudokuGen.cpp
)
target_compile_features(sudoku_gen PRIVATE cxx_std_20)
set_target_properties(sudoku_gen PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")