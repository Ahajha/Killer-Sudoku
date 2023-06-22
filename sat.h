/****************************************************************************
  FileName     [ sat.h ]
  PackageName  [ sat ]
  Synopsis     [ Define miniSat solver interface functions ]
  Author       [ Chung-Yang (Ric) Huang, Cheng-Yin Wu ]
  Copyright    [ Copyleft(c) 2010-present LaDs(III), GIEE, NTU, Taiwan ]
****************************************************************************/

#ifndef SAT_H
#define SAT_H

#include "Solver.h"

/********** MiniSAT_Solver **********/
class SatSolver : public Solver {
   public: 
      // Constructing proof model

      // fa/fb = true if it is inverted
      void addAigCNF(Var vf, Var va, bool fa, Var vb, bool fb) {
         vec<Lit> lits;
         Lit lf = Lit{vf};
         Lit la = fa ? ~Lit{va} : Lit{va};
         Lit lb = fb ? ~Lit{vb} : Lit{vb};
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
         lEq.push(Lit{f});
         for(size_t i=0, s=fanin.size(); i<s; ++i){
            rEq.push(~Lit{f}); rEq.push(fanin[i]);
            addClause(rEq); rEq.clear();
            lEq.push(~fanin[i]);
         }
         addClause(lEq);
      }
      // f <-> a or b or c
      void addOR(Var f, const vec<Lit>& fanin){
         vec<Lit> rEq, lEq;
         lEq.push(~Lit{f});
         for(size_t i=0, s=fanin.size(); i<s; ++i){
            rEq.push(Lit{f}); rEq.push(~fanin[i]);
            addClause(rEq); rEq.clear();
            lEq.push(fanin[i]);
         }
         addClause(lEq);
      }
      // fa/fb = true if it is inverted
      void addXorCNF(Var vf, Var va, bool fa, Var vb, bool fb) {
         vec<Lit> lits;
         Lit lf = Lit{vf};
         Lit la = fa ? ~Lit{va} : Lit{va};
         Lit lb = fb ? ~Lit{vb} : Lit{vb};
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

      void printStats() /* const */ {
         reportf("==============================[MINISAT]");
         reportf("===============================\n");
         reportf("| Conflicts |     ORIGINAL     |          ");
         reportf("LEARNT          | Progress |\n");
         reportf("|           | Clauses Literals | Clauses ");
         reportf("Literals  Lit/Cl |          |\n");
         reportf("=======================================");
         reportf("===============================\n");
         reportf("| %9d | %7d %8d | %7d %8d %7.1f | %6.3f %% |\n",
                 (int)stats.conflicts, nClauses(), (int)stats.clauses_literals,
                 nLearnts(), (int)stats.learnts_literals,
                 (double)stats.learnts_literals / nLearnts(),
                 progress_estimate * 100);
         reportf("=======================================");
         reportf("===============================\n");
      }
};

#endif  // SAT_H

