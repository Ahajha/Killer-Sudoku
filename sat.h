/****************************************************************************
  FileName     [ sat.h ]
  PackageName  [ sat ]
  Synopsis     [ Define miniSat solver interface functions ]
  Author       [ Chung-Yang (Ric) Huang, Cheng-Yin Wu ]
  Copyright    [ Copyleft(c) 2010-present LaDs(III), GIEE, NTU, Taiwan ]
****************************************************************************/

#ifndef SAT_H
#define SAT_H

#include "core/Solver.h"

using Minisat::lbool;
using Minisat::Lit;
using Minisat::mkLit;
using Minisat::Solver;
using Minisat::Var;
using Minisat::vec;

/********** MiniSAT_Solver **********/
class SatSolver : public Solver {
public:
  // Constructing proof model

  void addCNF(const vec<Lit> &ps);
  // fa/fb = true if it is inverted
  void addAigCNF(Var vf, Var va, bool fa, Var vb, bool fb);
  // f <-> abc
  void addAND(Var f, const vec<Lit> &fanin);
  // f <-> a or b or c
  void addOR(Var f, const vec<Lit> &fanin);
  // fa/fb = true if it is inverted
  void addXorCNF(Var vf, Var va, bool fa, Var vb, bool fb);

  // Functions about Reporting
  // Return 1/0/-1; -1 means unknown value
  int getValue(Var v) const;

  void printStats() const;
};

#endif // SAT_H
