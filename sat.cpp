#include "sat.h"

void SatSolver::addAigCNF(Var vf, Var va, bool fa, Var vb, bool fb) {
  vec<Lit> lits;
  Lit lf = mkLit(vf);
  Lit la = fa ? ~mkLit(va) : mkLit(va);
  Lit lb = fb ? ~mkLit(vb) : mkLit(vb);
  lits.push(la);
  lits.push(~lf);
  addClause(lits);
  lits.clear();
  lits.push(lb);
  lits.push(~lf);
  addClause(lits);
  lits.clear();
  lits.push(~la);
  lits.push(~lb);
  lits.push(lf);
  addClause(lits);
  lits.clear();
}

void SatSolver::addCNF(const vec<Lit> &ps) { addClause(ps); }

void SatSolver::addAND(Var f, const vec<Lit> &fanin) {
  vec<Lit> rEq, lEq;
  lEq.push(mkLit(f));
  for (size_t i = 0, s = fanin.size(); i < s; ++i) {
    rEq.push(~mkLit(f));
    rEq.push(fanin[i]);
    addClause(rEq);
    rEq.clear();
    lEq.push(~fanin[i]);
  }
  addClause(lEq);
}

void SatSolver::addOR(Var f, const vec<Lit> &fanin) {
  vec<Lit> rEq, lEq;
  lEq.push(~mkLit(f));
  for (size_t i = 0, s = fanin.size(); i < s; ++i) {
    rEq.push(mkLit(f));
    rEq.push(~fanin[i]);
    addClause(rEq);
    rEq.clear();
    lEq.push(fanin[i]);
  }
  addClause(lEq);
}

void SatSolver::addXorCNF(Var vf, Var va, bool fa, Var vb, bool fb) {
  vec<Lit> lits;
  Lit lf = mkLit(vf);
  Lit la = fa ? ~mkLit(va) : mkLit(va);
  Lit lb = fb ? ~mkLit(vb) : mkLit(vb);
  lits.push(~la);
  lits.push(lb);
  lits.push(lf);
  addClause(lits);
  lits.clear();
  lits.push(la);
  lits.push(~lb);
  lits.push(lf);
  addClause(lits);
  lits.clear();
  lits.push(la);
  lits.push(lb);
  lits.push(~lf);
  addClause(lits);
  lits.clear();
  lits.push(~la);
  lits.push(~lb);
  lits.push(~lf);
  addClause(lits);
  lits.clear();
}

int SatSolver::getValue(Var v) const {
  return (model[v] == l_True ? 1 : (model[v] == l_False ? 0 : -1));
}

void SatSolver::printStats() const {
  printf("==============================[MINISAT]");
  printf("===============================\n");
  printf("| Conflicts |     ORIGINAL     |          ");
  printf("LEARNT          | Progress |\n");
  printf("|           | Clauses Literals | Clauses ");
  printf("Literals  Lit/Cl |          |\n");
  printf("=======================================");
  printf("===============================\n");
  printf("| %9d | %7d %8d | %7d %8d %7.1f | %6.3f %% |\n", (int)conflicts,
         nClauses(), (int)clauses_literals, nLearnts(), (int)learnts_literals,
         (double)learnts_literals / nLearnts(), progress_estimate * 100);
  printf("=======================================");
  printf("===============================\n");
}