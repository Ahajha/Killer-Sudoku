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