#include <iostream>
#include <string>
#include <filesystem>

#include "sudoku.h"

// START: The main function
int main([[maybe_unused]] int argc, char const *argv[])
{
  // Initialising seed for random number generation
  srand(time(NULL));

  // Creating an instance of Sudoku
  std::cout << "======================================" << std::endl;
  std::cout << "random generating puzzle by SAT...." << std::endl;
  std::cout << "======================================" << std::endl;
  Sudoku *puzzle = new Sudoku();

  // Creating a seed for puzzle generation
  puzzle->createSeed();

  // Generating the puzzle cage
  puzzle->genPuzzle();

  // testing by printing the grid
  // puzzle->printGrid();

  // Printing the grid into SVG file
  const std::string path = std::filesystem::path(argv[0]).remove_filename();
  puzzle->printSVG(path);
  puzzle->printSVG(path, "images/puzzles_sol.svg", true);

  std::cout << "======================================" << std::endl;
  std::cout << "trying to solve puzzle by SAT...." << std::endl;
  std::cout << "======================================" << std::endl;

  puzzle->solveBySAT();
  puzzle->printSVG(path, "images/puzzles_solbySAT.svg", true);
  // cout<<"The above sudoku puzzle has been stored in puzzles.svg in current folder\n";
  // freeing the memory
  // puzzle->printGrid();
  delete puzzle;

  return 0;
}
// END: The main function
