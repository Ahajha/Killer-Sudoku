cmake_minimum_required(VERSION 3.15)
project(sudoku_gen)

find_package(minisat REQUIRED)

add_executable(sudoku_gen)
target_sources(sudoku_gen PRIVATE
    sudoku.cpp
    sudokuGen.cpp
)
target_compile_features(sudoku_gen PRIVATE cxx_std_20)
set_target_properties(sudoku_gen PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
target_link_libraries(sudoku_gen PRIVATE minisat::minisat)

if (NOT MSVC)
    #target_compile_options(sudoku_gen PRIVATE -Wall -Wextra -Wpedantic -Wconversion)
endif()