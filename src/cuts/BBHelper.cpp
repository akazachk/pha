//============================================================================
// Name        : BBHelper.cpp
// Author      : akazachk
// Version     : 2017.03.16
// Description : Helper functions for branch-and-bound
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#include <cstdio> // for tmpnam

// For SHOULD
#include "BBHelper.hpp"
#include "GlobalConstants.hpp"
//#include "Utility.hpp"
//#include "Output.hpp"
//#include "CutHelper.hpp"
//#include "CoinTime.hpp"
//#include "CglSIC.hpp"

// Gurobi
#ifdef SHOULD_USE_GUROBI
#include "gurobi_c++.h"

void setStrategyForBBTestGurobi(GRBModel& model, 
    //const int seed = GlobalVariables::random_seed,
    const double best_bound = GlobalVariables::bestObjValue) {
  // Parameters that should always be set
//  model.set(GRB_DoubleParam_TimeLimit, GlobalVariables::param.getBB_TIMELIMIT()); // time limit
//  model.set(GRB_IntParam_Threads, 1); // single-threaded
//  model.set(GRB_IntParam_Seed, seed); // random seed
//  model.set(GRB_DoubleParam_MIPGap, param.getEPS()); // I guess the default 1e-4 is okay, though it messes up for large objective values
  if (!isInfinity(std::abs(best_bound))) {
//    model.set(GRB_DoubleParam_BestObjStop, best_bound + 1e-3); // give the solver the best IP objective value (it is a minimization problem) with a tolerance
    model.set(GRB_DoubleParam_BestBdStop, best_bound + 1e-3); // give the solver the best IP objective value (it is a minimization problem) with a tolerance
  }

#ifndef TRACE
  model.set(GRB_IntParam_OutputFlag, 0); // turn off output
#endif
} /* setStrategyForBBTestGurobi */

void doBranchAndBoundWithGurobi(GRBModel& model, BBInfo& info,
    std::vector<double>* const solution = NULL) {
//#ifdef TRACE
  printf("\n## Running B&B with Gurobi. ##\n"); 
  //Strategy: %d. Random seed: %d. ##\n", param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND), GlobalVariables::random_seed);
//#endif
  try {
    setStrategyForBBTestGurobi(model);

    model.optimize();

    int optimstatus = model.get(GRB_IntAttr_Status);

    if (optimstatus == GRB_INF_OR_UNBD) {
      const int presolve_flag = model.get(GRB_IntParam_Presolve);
      if (presolve_flag) {
        model.set(GRB_IntParam_Presolve, 0);
        model.optimize();
        optimstatus = model.get(GRB_IntAttr_Status);
      }
    }

    if (optimstatus == GRB_INFEASIBLE || optimstatus == GRB_UNBOUNDED || optimstatus == GRB_INF_OR_UNBD) {
      error_msg(errorstring, "Gurobi: Failed to optimize MIP.\n");
      writeErrorToII(errorstring, GlobalVariables::log_file);
      exit(1);
    }

    switch (optimstatus) {
      case GRB_CUTOFF:
      case GRB_ITERATION_LIMIT:
      case GRB_NODE_LIMIT:
      case GRB_TIME_LIMIT: {
        info.obj = model.get(GRB_DoubleAttr_ObjVal);
        info.bound = model.get(GRB_DoubleAttr_ObjBound);
        break;
      }
      case GRB_USER_OBJ_LIMIT: {
        info.obj = model.get(GRB_DoubleAttr_ObjVal);
        info.bound = model.get(GRB_DoubleAttr_ObjBound);
        break;
      }
      case GRB_OPTIMAL: {
        info.obj = model.get(GRB_DoubleAttr_ObjVal);
        info.bound = model.get(GRB_DoubleAttr_ObjBound);
        break;
      }
      default: {
        error_msg(errorstring, "Gurobi: Other status after solve: %d.\n", optimstatus);
        writeErrorToII(errorstring, GlobalVariables::log_file);
        exit(1);
      }
    } /* switch optimistatus */
    info.iters = (long) model.get(GRB_DoubleAttr_IterCount);
    info.nodes = (long) model.get(GRB_DoubleAttr_NodeCount);
    info.time = model.get(GRB_DoubleAttr_Runtime);

#ifdef TRACE
    printf("Gurobi: Solution value: %1.6f.\n", info.obj);
    printf("Gurobi: Best bound: %1.6f.\n", info.bound);
    printf("Gurobi: Number iterations: %ld.\n", info.iters);
    printf("Gurobi: Number nodes: %ld.\n", info.nodes);
    printf("Gurobi: Time: %f.\n", info.time);
#endif

   // Save the solution if needed
   if (solution) {
     // Get variables
     GRBVar* vars = model.getVars();
     const int num_vars = model.get(GRB_IntAttr_NumVars);
     (*solution).resize(num_vars);
     for (int i = 0; i < num_vars; i++) {
       (*solution)[i] = vars[i].get(GRB_DoubleAttr_X);
     }
     if (vars) {
       delete[] vars;
     }
   }
  } catch(GRBException& e) {
    error_msg(errorstring, "Gurobi: Exception caught: %s\n", e.getMessage().c_str());
    writeErrorToII(errorstring, GlobalVariables::log_file);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "Gurobi: Unknown exception caught.\n");
    writeErrorToII(errorstring, GlobalVariables::log_file);
    exit(1);
  }
} /* doBranchAndBoundWithGurobi (GRBModel) */

void doBranchAndBoundWithGurobi(const char* f_name, BBInfo& info,
    std::vector<double>* const solution) {
#ifdef TRACE
  printf("\n## Reading from file into Gurobi. ##\n");
#endif
  try {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env, f_name);
    doBranchAndBoundWithGurobi(model, info, solution);
  } catch (GRBException& e) {
    error_msg(errorstring, "Gurobi: Exception caught: %s\n", e.getMessage().c_str());
    writeErrorToII(errorstring, GlobalVariables::log_file);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "Gurobi: Unknown exception caught.\n");
    writeErrorToII(errorstring, GlobalVariables::log_file);
    exit(1);
  }
} /* doBranchAndBoundWithGurobi (filename) */

#endif /* SHOULD_USE_GUROBI */
