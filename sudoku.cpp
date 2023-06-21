#include "sudoku.h"
#include "color.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <ios>
#include <iostream>
#include <numeric>
#include <iomanip>

// START: Get grid as string in row major order
std::string Sudoku::getGrid()
{
  std::string s = "";
  for(std::size_t row_num=0; row_num<gridSize; ++row_num)
  {
    for (std::size_t col_num = 0; col_num < gridSize; ++col_num) {
      s = s + std::to_string(_grid[row_num][col_num]);
    }
  }

  return s;
}
// END: Get grid as string in row major order

// START: Create seed grid
void Sudoku::createSeed()
{ 
  this->solveBySAT();
  
  // Saving the solution grid
  //_solnGrid = _grid;
}
// END: Create seed grid


// START: Intialising
Sudoku::Sudoku()
{
  // Initialising the grid
  for (std::size_t i = 0; i < gridSize; i++) {
    for (std::size_t j = 0; j < gridSize; j++) {
      _grid[i][j]=0;
      _cageId[i][j] = -1;
    }
  }
}
// END: Initialising


// START: Printing the grid
void Sudoku::printGrid()
{
  for (std::size_t i = 0; i < gridSize; i++) {
    for (std::size_t j = 0; j < gridSize; j++) {
      if(_grid[i][j] == 0)
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
void Sudoku::genPuzzle()
{
  size_t dir;
  std::vector<int> cageAppeared;
  for (std::size_t j = 0; j < gridSize; ++j) {
    for (std::size_t i = 0; i < gridSize; ++i) {
      if (_cageId[i][j] > -1) {
        continue;
      }

      int deadLock = 0;
      cageAppeared.clear();
      std::size_t sizeOfCage = rand() % maxCageSizeMinus1 + 2;
      int currentID = static_cast<int>(_cages.size());

      _cageId[i][j] = currentID;

      int sum = _grid[i][j];

      std::vector<Position> temp;

      temp.push_back({i, j, static_cast<int>(sum)});
      cageAppeared.push_back(sum);

      while (temp.size() < sizeOfCage && deadLock < 1) {
        std::size_t extPo = rand() % temp.size();
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
          if (i > 1 && _cages[_cageId[i - 1][j]].addEle(temp[0], sum)) {
            _cageId[i][j] = _cageId[i - 1][j];
            continue;
          }
        default:
          if (j > 1 && _cages[_cageId[i][j - 1]].addEle(temp[0], sum)) {
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

  int width = 50 * (gridSize+2);

  std::stringstream head;
  head << width << "\" height=\"" << width << "\">\n" ;
  outFile << head.rdbuf();

  for (std::size_t j = 0; j < gridSize; j++) {
    for (std::size_t i = 0; i < gridSize; i++) {
      if(this->_grid[i][j]!=0)
      {
        std::size_t x = 50*j;
        std::size_t y = 50*i;

        std::stringstream text;
        text << "<rect x=\""<<x<<"\" y=\""<<y<<"\" width=\"50\" height=\"50\" style=\"fill:" << colors[_cageId[i][j] % colors.size()] <<";opacity:"<< opacity(_cageId[i][j]) << "\"/>\n";
        if(printSol){
            text<<"<text x=\""<<x+16<<"\" y=\""<<y+35<<"\" style=\"font-weight:bold\" font-size=\"30px\">"<<_grid[i][j]<<"</text>\n";
        }
        outFile << text.rdbuf();
      }
    }
  }

  for(auto it=_cages.begin(); it!=_cages.end(); ++it){
      std::size_t x = 50* it->getPoy(0) + 8;
      std::size_t y = 50* it->getPox(0) + 18;

      std::stringstream text;
      text<<"<text x=\""<<x<<"\" y=\""<<y<<"\" style=\"font-weight:bold\" fill=\"red\" font-size=\"15px\">"<<it->getSum()<<"</text>\n";
      outFile << text.rdbuf();
  }

  for (std::size_t i = 0; i <= gridSize; ++i) {
      std::stringstream text, text2;
      text<<  "<polyline points=\"" << 50*i << ",0 " << 50*i << "," << width-100 << "\" style=\"fill:none; stroke:black ; stroke-width:";
      text2<<  "<polyline points=\"0," << 50*i << " " << width-100 << "," << 50*i<< "\" style=\"fill:none; stroke:black ; stroke-width:";
      if(i%boxSize == 0){
          text << 5;
          text2 << 5;
      }
      else{
          text << 1;
          text2 << 1;
      }
      text <<"\" />\n";
      text2 <<"\" />\n";
      
      outFile << text.rdbuf();
      outFile << text2.rdbuf();
  }

  outFile << "</svg>";

  outFile.close();
}
// END: Printing into SVG file

bool sortGrid(const Position& a, const Position& b){
    return (a.y != b.y)?(a.y < b.y):(a.x < b.x);
}

void Sudoku::gateInitial(SatSolver& solver, Gates& gates){
    for (std::size_t i = 0; i < gridSize; ++i) {
      for (std::size_t j = 0; j < gridSize; ++j) {
          for (std::size_t k = 0; k < gridSize; ++k) {
            Var v = solver.newVar();
            gates[i][j][k].setVar(v);
          }
      }
    }
}

void Sudoku::genProofModel(SatSolver& solver, Gates& gates){
    vec<Lit> lits;
    // entry condition
    for (std::size_t j = 0; j < gridSize; ++j) {
        for (std::size_t i = 0; i < gridSize; ++i) {
          for (std::size_t k = 0; k < gridSize; ++k) {
        lits.push(Lit(gates[i][j][k].getVar()));
          }
          solver.addCNF(lits);
          lits.clear();
        }
    }
    for (std::size_t j = 0; j < gridSize; ++j) {
        for (std::size_t i = 0; i < gridSize; ++i) {
          for (std::size_t k = 0; k < gridSize - 1; ++k) {
            for (std::size_t z = k + 1; z < gridSize; ++z) {
                lits.push(~Lit(gates[i][j][k].getVar()));
                lits.push(~Lit(gates[i][j][z].getVar()));
                solver.addCNF(lits);
                lits.clear();
            }
          }
        }
    }

    // row
    for (std::size_t j = 0; j < gridSize; ++j) {
        for (std::size_t k = 0; k < gridSize; ++k) {
          for (std::size_t i = 0; i < gridSize - 1; ++i) {
            for (std::size_t q = i + 1; q < gridSize; ++q) {
                lits.push(~Lit(gates[i][j][k].getVar()));
                lits.push(~Lit(gates[q][j][k].getVar()));
                solver.addCNF(lits);
                lits.clear();
            }
          }
        }
    }

    // column
    for (std::size_t i = 0; i < gridSize; ++i) {
        for (std::size_t k = 0; k < gridSize; ++k) {
          for (std::size_t j = 0; j < gridSize - 1; ++j) {
            for (std::size_t q = j + 1; q < gridSize; ++q) {
                lits.push(~Lit(gates[i][j][k].getVar()));
                lits.push(~Lit(gates[i][q][k].getVar()));
                solver.addCNF(lits);
                lits.clear();
            }
          }
        }
    }

    // box
    for (std::size_t k = 0; k < gridSize; ++k) {
        for (std::size_t q = 0; q < boxSize; ++q) {
          for (std::size_t r = 0; r < boxSize; ++r) {
            for (std::size_t i = 0; i < boxSize; ++i) {
                for (std::size_t j = 0; j < boxSize; ++j) {
                    for (std::size_t s = j + 1; s < boxSize; ++s) {
                      lits.push(
                          ~Lit(gates[boxSize * q + i][boxSize * r + j][k].getVar()));
                      lits.push(
                          ~Lit(gates[boxSize * q + i][boxSize * r + s][k].getVar()));
                      solver.addCNF(lits);
                      lits.clear();
                    }
                }
            }
          }
        }
    }

    for (std::size_t k = 0; k < gridSize; ++k) {
        for (std::size_t q = 0; q < boxSize; ++q) {
          for (std::size_t r = 0; r < boxSize; ++r) {
            for (std::size_t i = 0; i < boxSize; ++i) {
                for (std::size_t j = 0; j < boxSize; ++j) {
                    for (std::size_t s = i + 1; s < boxSize; ++s) {
                      for (std::size_t t = 0; t < boxSize; ++t) {
                        lits.push(
                            ~Lit(gates[boxSize * q + i][boxSize * r + j][k]
                                     .getVar()));
                        lits.push(
                            ~Lit(gates[boxSize * q + s][boxSize * r + t][k]
                                     .getVar()));
                        solver.addCNF(lits);
                        lits.clear();
                      }
                    }
                }
            }
          }
        }
    }

    std::vector<int> numbers;
    std::vector<int> partial;

    for (std::size_t i = 1; i <= gridSize; ++i) {
        numbers.push_back(static_cast<int>(i));
    }

    // sum
     for(auto it=_cages.begin(); it!=_cages.end(); ++it){
        /** this is concerned already bellow
         // each num appears at most once in the cage
         for(size_t s=it->getCageSize()-1, i=0; i<s; ++i){
            for(int k=0; k<9; ++k){
                lits.push(~Lit(gates[it->getPox(i)][it->getPoy(i)][k].getVar()));
                lits.push(~Lit(gates[it->getPox(i+1)][it->getPoy(i+1)][k].getVar()));
                solver.addCNF(lits); lits.clear();
            }
         }
        **/ 

         size_t cs = it->getCageSize();
         int sum = it->getSum();
         std::vector<std::vector<int>> answers;
         subsetSum(numbers, cs, sum, partial, answers);
         vec<Lit> validSols;
         for(auto ans=answers.begin(); ans!=answers.end(); ++ans){
             if(ans->size() == cs){
                 do{
                    for(size_t i=0; i<cs; ++i){
                        lits.push(Lit(gates[it->getPox(i)][it->getPoy(i)][(*ans)[i]-1].getVar()));
                    }
                    Var v = solver.newVar();
                    solver.addAND(v, lits); lits.clear();

                    validSols.push(Lit(v));
                 }while(next_permutation(ans->begin(), ans->end()));
             }
         }
         solver.addOR(it->getGate().getVar(), validSols);
     }

}

void Sudoku::solveBySAT(){
    Gates gates;

    SatSolver solver;
    solver.initialize();

    for (size_t s = _cages.size(), i = 0; i < s; ++i) {
         _cages[i].setGate(solver);
    }

    gateInitial(solver, gates);

    clock_t start, end;
    start = clock();
    genProofModel(solver, gates);
    end = clock();

    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    std::cout << "Time taken by genProofModel() is : " << std::fixed
              << time_taken << std::setprecision(5);
    std::cout << " sec " << std::endl;

    bool result;
    // k = Solve(Gate(5) ^ !Gate(8))
    // Var newV = solver.newVar();
    // solver.addXorCNF(newV, gates[5]->getVar(), false, gates[8]->getVar(), true);
    solver.assumeRelease();  // Clear assumptions

    // letting all sum condition be true
    for (std::size_t i = 0, s = _cages.size(); i < s; ++i) {
         solver.assumeProperty(_cages[i].getGate().getVar(), true);
    }
    // solver.assumeProperty(newV, true);  // k = 1
    start = clock();
    result = solver.assumpSolve();
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
                    if(solver.getValue(gates[i][j][k].getVar())){
                        _grid[i][j] = k+1;
                        break;
                    }
                 }
             }
         }
    }
    else{
         for (std::size_t j = 0; j < gridSize; ++j) {
             for (std::size_t i = 0; i < gridSize; ++i) {
                 _grid[i][j] = 0;
             }
         }
    }
}

void Sudoku::subsetSum(std::vector<int> numbers, const int &s,
                       const int &target, std::vector<int> partial,
                       std::vector<std::vector<int>> &answer) {
    int sum = std::accumulate(partial.begin(), partial.end(), 0);

    if(sum == target){
        answer.push_back(partial);
    }

    if(sum >= target || partial.size() >= s){
        return;
    }

    for(auto it=numbers.begin(); it!=numbers.end(); ++it){
        partial.push_back(*it);
        std::vector<int> remainig(it + 1, numbers.end());
        subsetSum(remainig, s, target, partial, answer);
        partial.pop_back();
    }
}