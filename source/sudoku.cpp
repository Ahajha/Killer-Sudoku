#include "sudoku.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string_view>

namespace {

constexpr std::array<std::string_view, 12> colors{
    "LightCoral", "DarkMagenta", "DarkOrange",  "Yellow",
    "DarkKhaki",  "Lavender",    "LightPink",   "LawnGreen",
    "SteelBlue",  "Teal",        "SaddleBrown", "Cyan",
};

float opacity(int n) { return n % 2 ? 0.3f : 1.0f; }

void subsetSum(const int minValue, const std::size_t target_size,
               const int target_sum, std::vector<int> &partial,
               std::vector<std::vector<int>> &answer) {
  int sum = std::accumulate(partial.begin(), partial.end(), 0);

  if (sum == target_sum && partial.size() == target_size) {
    answer.push_back(partial);
  }

  if (sum >= target_sum || partial.size() >= target_size) {
    return;
  }

  for (auto value = minValue; value <= static_cast<int>(gridSize); ++value) {
    partial.push_back(value);
    subsetSum(value + 1, target_size, target_sum, partial, answer);
    partial.pop_back();
  }
}

} // namespace

using Minisat::mkLit;

// START: Get grid as string in row major order
std::string Sudoku::getGrid() {
  std::string s = "";
  for (std::size_t row_num = 0; row_num < gridSize; ++row_num) {
    for (std::size_t col_num = 0; col_num < gridSize; ++col_num) {
      s = s + std::to_string(_grid[row_num][col_num]);
    }
  }

  return s;
}
// END: Get grid as string in row major order

// START: Create seed grid
void Sudoku::createSeed() {
  this->solveBySAT();

  // Saving the solution grid
  //_solnGrid = _grid;
}
// END: Create seed grid

// START: Intialising
Sudoku::Sudoku() {
  // Initialising the grid
  for (std::size_t i = 0; i < gridSize; i++) {
    for (std::size_t j = 0; j < gridSize; j++) {
      _grid[i][j] = 0;
      _cageId[i][j] = -1;
    }
  }
}
// END: Initialising

// START: Printing the grid
void Sudoku::printGrid() {
  for (std::size_t i = 0; i < gridSize; i++) {
    for (std::size_t j = 0; j < gridSize; j++) {
      if (_grid[i][j] == 0)
        std::cout << ".";
      else
        std::cout << _grid[i][j];
      std::cout << "|";
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;
}
// END: Printing the grid

// START: Gneerate puzzle
void Sudoku::genPuzzle() {
  size_t dir;
  std::vector<int> cageAppeared;
  for (std::size_t j = 0; j < gridSize; ++j) {
    for (std::size_t i = 0; i < gridSize; ++i) {
      if (_cageId[i][j] > -1) {
        continue;
      }

      int deadLock = 0;
      cageAppeared.clear();
      std::size_t sizeOfCage =
          static_cast<std::size_t>(rand()) % maxCageSizeMinus1 + 2;
      int currentID = static_cast<int>(_cages.size());

      _cageId[i][j] = currentID;

      int sum = _grid[i][j];

      std::vector<Position> temp;

      temp.push_back({i, j, static_cast<int>(sum)});
      cageAppeared.push_back(sum);

      while (temp.size() < sizeOfCage && deadLock < 1) {
        std::size_t extPo = static_cast<std::size_t>(rand()) % temp.size();
        std::size_t pox = temp[extPo].x;
        std::size_t poy = temp[extPo].y;

        int dir = rand() % 2;
        switch (dir) {
        case 0:
          if (pox < 8 && _cageId[pox + 1][poy] < 0) {
            int gridValue = _grid[pox + 1][poy];
            if (std::find(cageAppeared.begin(), cageAppeared.end(),
                          gridValue) == cageAppeared.end()) {
              temp.push_back({pox + 1, poy, gridValue});
              _cageId[pox + 1][poy] = currentID;
              cageAppeared.push_back(gridValue);
              sum += gridValue;
              break;
            }
          }
        case 1:
          if (poy < 8 && _cageId[pox][poy + 1] < 0) {
            int gridValue = _grid[pox][poy + 1];
            if (std::find(cageAppeared.begin(), cageAppeared.end(),
                          gridValue) == cageAppeared.end()) {
              temp.push_back({pox, poy + 1, gridValue});
              _cageId[pox][poy + 1] = currentID;
              cageAppeared.push_back(gridValue);
              sum += gridValue;
              break;
            }
          }
        case 2:
          if (pox > 1 && _cageId[pox - 1][poy] < 0) {
            int gridValue = _grid[pox - 1][poy];
            if (std::find(cageAppeared.begin(), cageAppeared.end(),
                          gridValue) == cageAppeared.end()) {
              temp.push_back({pox - 1, poy, gridValue});
              _cageId[pox - 1][poy] = currentID;
              cageAppeared.push_back(gridValue);
              sum += gridValue;
              break;
            }
          }
        case 3:
          if (poy > 1 && _cageId[pox][poy - 1] < 0) {
            int gridValue = _grid[pox][poy - 1];
            if (std::find(cageAppeared.begin(), cageAppeared.end(),
                          gridValue) == cageAppeared.end()) {
              temp.push_back({pox, poy - 1, gridValue});
              _cageId[pox][poy - 1] = currentID;
              cageAppeared.push_back(gridValue);
              sum += gridValue;
              break;
            }
          }
        default:
          ++deadLock;
          break;
        }
      }

      if (temp.size() == 1) {
        switch (dir) {
        case 0:
          if (i > 1 &&
              _cages[static_cast<std::size_t>(_cageId[i - 1][j])].addEle(
                  temp[0], sum)) {
            _cageId[i][j] = _cageId[i - 1][j];
            continue;
          }
        default:
          if (j > 1 &&
              _cages[static_cast<std::size_t>(_cageId[i][j - 1])].addEle(
                  temp[0], sum)) {
            _cageId[i][j] = _cageId[i][j - 1];
            continue;
          }
        }
      }

      std::sort(temp.begin(), temp.end(), sortGrid);
      _cages.push_back(Cage(currentID, sum, temp));
    }
  }
}
// END: Generate puzzle

// START: Printing into SVG file
void Sudoku::printSVG(std::string path, std::string svgName, bool printSol) {
  std::string fileName = path + "svgHead.txt";
  std::ifstream file1(fileName.c_str());
  std::stringstream svgHead;
  svgHead << file1.rdbuf();

  std::ofstream outFile(svgName);
  outFile << svgHead.rdbuf();

  file1.close();

  int width = 50 * (gridSize + 2);

  std::stringstream head;
  head << width << "\" height=\"" << width << "\">\n";
  outFile << head.rdbuf();

  for (std::size_t j = 0; j < gridSize; j++) {
    for (std::size_t i = 0; i < gridSize; i++) {
      if (this->_grid[i][j] != 0) {
        std::size_t x = 50 * j;
        std::size_t y = 50 * i;

        std::stringstream text;
        text << "<rect x=\"" << x << "\" y=\"" << y
             << "\" width=\"50\" height=\"50\" style=\"fill:"
             << colors[static_cast<std::size_t>(_cageId[i][j]) % colors.size()]
             << ";opacity:" << opacity(_cageId[i][j]) << "\"/>\n";
        if (printSol) {
          text << "<text x=\"" << x + 16 << "\" y=\"" << y + 35
               << "\" style=\"font-weight:bold\" font-size=\"30px\">"
               << _grid[i][j] << "</text>\n";
        }
        outFile << text.rdbuf();
      }
    }
  }

  for (auto it = _cages.begin(); it != _cages.end(); ++it) {
    std::size_t x = 50 * it->getPoy(0) + 8;
    std::size_t y = 50 * it->getPox(0) + 18;

    std::stringstream text;
    text << "<text x=\"" << x << "\" y=\"" << y
         << "\" style=\"font-weight:bold\" fill=\"red\" font-size=\"15px\">"
         << it->getSum() << "</text>\n";
    outFile << text.rdbuf();
  }

  for (std::size_t i = 0; i <= gridSize; ++i) {
    std::stringstream text, text2;
    text << "<polyline points=\"" << 50 * i << ",0 " << 50 * i << ","
         << width - 100 << "\" style=\"fill:none; stroke:black ; stroke-width:";
    text2 << "<polyline points=\"0," << 50 * i << " " << width - 100 << ","
          << 50 * i << "\" style=\"fill:none; stroke:black ; stroke-width:";
    if (i % boxSize == 0) {
      text << 5;
      text2 << 5;
    } else {
      text << 1;
      text2 << 1;
    }
    text << "\" />\n";
    text2 << "\" />\n";

    outFile << text.rdbuf();
    outFile << text2.rdbuf();
  }

  outFile << "</svg>";

  outFile.close();
}
// END: Printing into SVG file

bool sortGrid(const Position &a, const Position &b) {
  return (a.y != b.y) ? (a.y < b.y) : (a.x < b.x);
}

void Sudoku::gateInitial(SatSolver &solver, Gates &gates) {
  for (std::size_t i = 0; i < gridSize; ++i) {
    for (std::size_t j = 0; j < gridSize; ++j) {
      for (std::size_t k = 0; k < gridSize; ++k) {
        Var v = solver.newVar();
        gates[i][j][k] = v;
      }
    }
  }
}

void Sudoku::genProofModel(SatSolver &solver, Gates &gates) {
  // Ahajha:
  // From my basic initial understanding, my guess is that each
  // solver variable represents "Can this cell be this value",
  // which gives an initial set of gridSize * gridSize * gridSize variables.
  // These are stored in `gates`.

  // Throughout, the variable names 'row', 'col', and 'num' refer to row
  // indexes, column indexes, and possible value for a cell (minus one).
  vec<Lit> lits;

  // entry condition

  // Together, these two conditions ensure that each cell has exactly one
  // value. For a grid size of 9, and representing "can this cell be N" as
  // "LN", for a single cell we have:
  //
  // (L1 | L2 | L3 | L4 | L5 | L6 | L7 | L8 | L9) &
  // (~L1 | ~L2) & (~L1 | ~L3) & (~L1 | ~L4) & (~L1 | ~L5) & (~L1 | ~L6) &
  // (~L1 | ~L7) & (~L1 | ~L8) & (~L1 | ~L9) &
  // (~L2 | ~L3) & (~L2 | ~L4) & (~L2 | ~L5) & (~L2 | ~L6) & (~L2 | ~L7) &
  // (~L2 | ~L8) & (~L2 | ~L9) &
  // (~L3 | ~L4) & (~L3 | ~L5) & (~L3 | ~L6) & (~L3 | ~L7) & (~L3 | ~L8) &
  // (~L3 | ~L9) &
  // (~L4 | ~L5) & (~L4 | ~L6) & (~L4 | ~L7) & (~L4 | ~L8) & (~L4 | ~L9) &
  // (~L5 | ~L6) & (~L5 | ~L7) & (~L5 | ~L8) & (~L5 | ~L9) &
  // (~L6 | ~L7) & (~L6 | ~L8) & (~L6 | ~L9) &
  // (~L7 | ~L8) & (~L7 | ~L59) &
  // (~L8 | ~L9)
  // (Could this be simplified? Probably some interesting research around
  // minimal XOR of N variables)
  for (std::size_t row = 0; row < gridSize; ++row) {
    for (std::size_t col = 0; col < gridSize; ++col) {
      // Each cell has gridSize gates, at least one of which must be true
      for (std::size_t num = 0; num < gridSize; ++num) {
        lits.push(mkLit(gates[col][row][num]));
      }
      solver.addClause_(lits);
      lits.clear();
    }
  }
  for (std::size_t row = 0; row < gridSize; ++row) {
    for (std::size_t col = 0; col < gridSize; ++col) {
      // Each cell has gridSize gates, no two of which can be true
      for (std::size_t num1 = 0; num1 < gridSize - 1; ++num1) {
        for (std::size_t num2 = num1 + 1; num2 < gridSize; ++num2) {
          lits.push(~mkLit(gates[col][row][num1]));
          lits.push(~mkLit(gates[col][row][num2]));
          solver.addCNF(lits);
          lits.clear();
        }
      }
    }
  }

  // row
  for (std::size_t row = 0; row < gridSize; ++row) {
    // Each pair of values in a row must be different
    for (std::size_t num = 0; num < gridSize; ++num) {
      for (std::size_t col1 = 0; col1 < gridSize - 1; ++col1) {
        for (std::size_t col2 = col1 + 1; col2 < gridSize; ++col2) {
          lits.push(~mkLit(gates[col1][row][num]));
          lits.push(~mkLit(gates[col2][row][num]));
          solver.addCNF(lits);
          lits.clear();
        }
      }
    }
  }

  // column
  for (std::size_t col = 0; col < gridSize; ++col) {
    // Each pair of values in a column must be different
    for (std::size_t num = 0; num < gridSize; ++num) {
      for (std::size_t row1 = 0; row1 < gridSize - 1; ++row1) {
        for (std::size_t row2 = row1 + 1; row2 < gridSize; ++row2) {
          lits.push(~mkLit(gates[col][row1][num]));
          lits.push(~mkLit(gates[col][row2][num]));
          solver.addCNF(lits);
          lits.clear();
        }
      }
    }
  }

  // box
  for (std::size_t num = 0; num < gridSize; ++num) {
    for (std::size_t col_of_box = 0; col_of_box < boxSize; ++col_of_box) {
      for (std::size_t row_of_box = 0; row_of_box < boxSize; ++row_of_box) {
        for (std::size_t col1_in_box = 0; col1_in_box < boxSize;
             ++col1_in_box) {
          for (std::size_t row1_in_box = 0; row1_in_box < boxSize;
               ++row1_in_box) {
            const std::size_t col1 = boxSize * col_of_box + col1_in_box;
            const std::size_t row1 = boxSize * row_of_box + row1_in_box;
            for (std::size_t col2_in_box = 0; col2_in_box < boxSize;
                 ++col2_in_box) {
              for (std::size_t row2_in_box = 0; row2_in_box < boxSize;
                   ++row2_in_box) {
                const std::size_t col2 = boxSize * col_of_box + col2_in_box;
                const std::size_t row2 = boxSize * row_of_box + row2_in_box;
                if (row1 != row2 || col1 != col2) {
                  lits.push(~mkLit(gates[col1][row1][num]));
                  lits.push(~mkLit(gates[col2][row2][num]));
                  solver.addCNF(lits);
                  lits.clear();
                }
              }
            }
          }
        }
      }
    }
  }

  std::vector<int> partial;

  // sum
  for (auto it = _cages.begin(); it != _cages.end(); ++it) {
    /** this is concerned already bellow
     // each num appears at most once in the cage
     for(size_t s=it->getCageSize()-1, i=0; i<s; ++i){
        for(int k=0; k<9; ++k){
            lits.push(~mkLit(gates[it->getPox(i)][it->getPoy(i)][k]));
            lits.push(~mkLit(gates[it->getPox(i+1)][it->getPoy(i+1)][k]));
            solver.addCNF(lits); lits.clear();
        }
     }
    **/

    size_t cs = it->getCageSize();
    int sum = it->getSum();
    std::vector<std::vector<int>> answers;
    subsetSum(1, cs, sum, partial, answers);

    vec<Lit> validSols;
    for (auto &ans : answers) {
      do {
        for (std::size_t i = 0; i < cs; ++i) {
          lits.push(mkLit(gates[it->getPox(i)][it->getPoy(i)]
                               [static_cast<std::size_t>(ans[i] - 1)]));
        }
        Var v = solver.newVar();
        solver.addAND(v, lits);
        lits.clear();

        validSols.push(mkLit(v));
      } while (std::next_permutation(ans.begin(), ans.end()));
    }
    solver.addOR(it->getGate(), validSols);
  }
}

void Sudoku::solveBySAT() {
  Gates gates;

  SatSolver solver;

  for (size_t s = _cages.size(), i = 0; i < s; ++i) {
    _cages[i].setGate(solver);
  }

  gateInitial(solver, gates);

  clock_t start, end;
  start = clock();
  genProofModel(solver, gates);
  end = clock();

  double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
  std::cout << "Time taken by genProofModel() is : " << std::fixed << time_taken
            << std::setprecision(5);
  std::cout << " sec " << std::endl;

  // k = Solve(Gate(5) ^ !Gate(8))
  // Var newV = solver.newVar();
  // solver.addXorCNF(newV, gates[5]->getVar(), false, gates[8]->getVar(),
  // true);

  vec<Lit> assumptions;
  // letting all sum condition be true
  for (std::size_t i = 0, s = _cages.size(); i < s; ++i) {
    assumptions.push(mkLit(_cages[i].getGate()));
  }
  // solver.assumeProperty(newV, true);  // k = 1
  start = clock();
  bool result = solver.solve(assumptions);
  end = clock();
  time_taken = double(end - start) / double(CLOCKS_PER_SEC);
  std::cout << "Time taken by solver.assumpSolve() is : " << std::fixed
            << time_taken << std::setprecision(5);
  std::cout << " sec " << std::endl;
  solver.printStats();
  std::cout << (result ? "SAT" : "UNSAT") << std::endl;
  if (result) {
    for (std::size_t i = 0, n = gridSize; i < n; ++i) {
      for (std::size_t j = 0; j < gridSize; ++j) {
        for (std::size_t k = 0; k < gridSize; ++k) {
          if (solver.getValue(gates[i][j][k])) {
            _grid[i][j] = k + 1;
            break;
          }
        }
      }
    }
  } else {
    for (std::size_t j = 0; j < gridSize; ++j) {
      for (std::size_t i = 0; i < gridSize; ++i) {
        _grid[i][j] = 0;
      }
    }
  }
}
