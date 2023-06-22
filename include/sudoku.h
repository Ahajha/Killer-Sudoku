#ifndef SUDOKU_H
#define SUDOKU_H

#include <array>
#include <string>
#include <vector>

#include "sat.h"

using Gate = Var;

struct Position {
  size_t x, y;
  int val;
};

bool sortGrid(const Position &a, const Position &b);
const size_t boxSize = 3;
const size_t gridSize = boxSize * boxSize;
const size_t maxCageSizeMinus1 = 5; // gridSize - 4;

class Cage {
public:
  Cage(std::vector<Position> e) : _eles{std::move(e)} {}
  Cage(int i, int s, std::vector<Position> e)
      : _cageId{i}, _sum{s}, _eles{std::move(e)} {}

  int getSum() const { return _sum; }
  int getID() const { return _cageId; }
  size_t getCageSize() const { return _eles.size(); }
  size_t getPox(size_t i) const { return _eles[i].x; }
  size_t getPoy(size_t i) const { return _eles[i].y; }
  Gate &getGate() { return _g; }
  void setGate(SatSolver &s) {
    Var v = s.newVar();
    _g = v;
  }
  bool addEle(Position e, int i) {
    if (getCageSize() > maxCageSizeMinus1) {
      return false;
    }
    for (Position x : _eles) {
      if (x.val == i) {
        return false;
      }
    }
    _eles.push_back(e);
    _sum += i;
    return true;
  }

private:
  int _cageId{-1};
  int _sum{0};
  std::vector<Position> _eles;
  Gate _g;
};

class Sudoku {
private:
  std::array<std::array<int, gridSize>, gridSize> _grid;
  // std::array<std::array<int, gridSize>, gridSize> _solnGrid;
  std::array<std::array<int, gridSize>, gridSize> _cageId;
  std::vector<Cage> _cages;

  using Gates =
      std::array<std::array<std::array<Gate, gridSize>, gridSize>, gridSize>;

  void gateInitial(SatSolver &solver, Gates &gates);
  void genProofModel(SatSolver &solver, Gates &gates);

public:
  Sudoku();
  void createSeed();
  void printGrid();
  std::string getGrid();
  void genPuzzle();
  void printSVG(std::string path = "",
                std::string svgName = "images/puzzle.svg",
                bool printSol = false);
  void solveBySAT();
};

#endif