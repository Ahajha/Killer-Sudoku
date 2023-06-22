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

using Minisat::vec;
using Minisat::Lit;
using Minisat::mkLit;
using Minisat::Var;
using Minisat::Solver;
using Minisat::lbool;

/********** MiniSAT_Solver **********/
class SatSolver : public Solver {
   public: 
      // Constructing proof model

      // fa/fb = true if it is inverted
      void addAigCNF(Var vf, Var va, bool fa, Var vb, bool fb) {
         vec<Lit> lits;
         Lit lf = mkLit(vf);
         Lit la = fa ? ~mkLit(va) : mkLit(va);
         Lit lb = fb ? ~mkLit(vb) : mkLit(vb);
         lits.push(la); lits.push(~lf);
         addClause(lits); lits.clear();
         lits.push(lb); lits.push(~lf);
         addClause(lits); lits.clear();
         lits.push(~la); lits.push(~lb); lits.push(lf);
         addClause(lits); lits.clear();
      }
      void addCNF(const vec<Lit>& ps){
         addClause(ps);
      }
      // f <-> abc 
      void addAND(Var f, const vec<Lit>& fanin){
         vec<Lit> rEq, lEq;
         lEq.push(mkLit(f));
         for(size_t i=0, s=fanin.size(); i<s; ++i){
            rEq.push(~mkLit(f)); rEq.push(fanin[i]);
            addClause(rEq); rEq.clear();
            lEq.push(~fanin[i]);
         }
         addClause(lEq);
      }
      // f <-> a or b or c
      void addOR(Var f, const vec<Lit>& fanin){
         vec<Lit> rEq, lEq;
         lEq.push(~mkLit(f));
         for(size_t i=0, s=fanin.size(); i<s; ++i){
            rEq.push(mkLit(f)); rEq.push(~fanin[i]);
            addClause(rEq); rEq.clear();
            lEq.push(fanin[i]);
         }
         addClause(lEq);
      }
      // fa/fb = true if it is inverted
      void addXorCNF(Var vf, Var va, bool fa, Var vb, bool fb) {
         vec<Lit> lits;
         Lit lf = mkLit(vf);
         Lit la = fa ? ~mkLit(va) : mkLit(va);
         Lit lb = fb ? ~mkLit(vb) : mkLit(vb);
         lits.push(~la); lits.push( lb); lits.push( lf);
         addClause(lits); lits.clear();
         lits.push( la); lits.push(~lb); lits.push( lf);
         addClause(lits); lits.clear();
         lits.push( la); lits.push( lb); lits.push(~lf);
         addClause(lits); lits.clear();
         lits.push(~la); lits.push(~lb); lits.push(~lf);
         addClause(lits); lits.clear();
      }

      // Functions about Reporting
      // Return 1/0/-1; -1 means unknown value
      int getValue(Var v) const {
         return (model[v]==l_True?1:
                (model[v]==l_False?0:-1));
      }

      void printStats() const {
         printf("==============================[MINISAT]");
         printf("===============================\n");
         printf("| Conflicts |     ORIGINAL     |          ");
         printf("LEARNT          | Progress |\n");
         printf("|           | Clauses Literals | Clauses ");
         printf("Literals  Lit/Cl |          |\n");
         printf("=======================================");
         printf("===============================\n");
         printf("| %9d | %7d %8d | %7d %8d %7.1f | %6.3f %% |\n",
                 (int)conflicts, nClauses(), (int)clauses_literals,
                 nLearnts(), (int)learnts_literals,
                 (double)learnts_literals / nLearnts(),
                 progress_estimate * 100);
         printf("=======================================");
         printf("===============================\n");
      }
};

#endif  // SAT_H

