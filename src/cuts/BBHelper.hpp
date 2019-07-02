//============================================================================
// Name        : BBHelper.hpp
// Author      : akazachk
// Version     : 2018-09-18
// Description : Helper functions for branch-and-bound
//============================================================================

#pragma once

#include "typedefs.hpp"
#include "AdvCut.hpp"

struct BBInfo {
  double obj;
  double bound;
  long iters;
  long nodes;
  long root_passes;
  double first_cut_pass;
  double last_cut_pass;
  double root_time;
  double last_sol_time;
  double time;
}; /* BBInfo */

enum BBInfoEnum {
  OBJ_BB_INFO_IND,
  BOUND_BB_INFO_IND,
  ITERS_BB_INFO_IND,
  NODES_BB_INFO_IND,
  ROOT_PASSES_BB_INFO_IND,
  FIRST_CUT_PASS_BB_INFO_IND,
  LAST_CUT_PASS_BB_INFO_IND,
  ROOT_TIME_BB_INFO_IND,
  LAST_SOL_TIME_BB_INFO_IND,
  TIME_BB_INFO_IND,
  NUM_BB_INFO
};

const std::vector<std::string> BB_INFO_CONTENTS = {
    "OBJ", "BOUND", "ITERS", "NODES", "ROOT_PASSES", "FIRST_CUT_PASS", "LAST_CUT_PASS", "ROOT_TIME", "LAST_SOL_TIME", "TIME"
};

inline void initializeBBInfo(BBInfo& info, double opt = 0.) {
  info.obj = opt;
  info.bound = opt;
  info.iters = 0;
  info.nodes = 0;
  info.root_passes = 0;
  info.first_cut_pass = 0.;
  info.last_cut_pass = 0.;
  info.root_time = 0.;
  info.last_sol_time = 0.;
  info.time = 0.;
}

void updateBestBBInfo(BBInfo& min_info, const BBInfo& curr_info, const bool first);
void averageBBInfo(BBInfo& avg_info, const std::vector<BBInfo>& info);
void printBBInfo(const BBInfo& info, FILE* myfile, const bool print_blanks =
    false, const char SEP = ',');
void printBBInfo(const BBInfo& info_mycuts, const BBInfo& info_allcuts,
    FILE* myfile, const bool print_blanks = false, const char SEP = ',');
void createStringFromBBInfoVec(const std::vector<BBInfo>& vec_info,
    std::vector<std::string>& vec_str);

#ifdef SHOULD_USE_GUROBI
// Gurobi stuff
void doBranchAndBoundWithGurobi(const char* f_name, BBInfo& info, std::vector<double>* const solution = NULL);
void doBranchAndBoundWithGurobi(const OsiSolverInterface* const solver,
    BBInfo& info, std::vector<double>* const solution = NULL);
#endif /* SHOULD_USE_GUROBI */
