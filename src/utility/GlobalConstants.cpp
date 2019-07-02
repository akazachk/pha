//============================================================================
// Name        : GlobalConstants.cpp
// Author      : akazachk
// Version     : 0.2014.01.10
// Copyright   : Your copyright notice
// Description : Miscellaneous useful definitions and functions used globally.
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#include "GlobalConstants.hpp"

namespace GlobalVariables {
  std::chrono::time_point<std::chrono::system_clock> start_time;
  std::FILE* log_file = stdout;
  std::string LOG_FILE;
  std::string prob_name, out_f_name_stub;

  Stats timeStats;
  std::vector<std::string> time_name;
  CglGICParam param;

  std::vector<int> numCutSolverFails(GlobalConstants::CutSolverFails::NUM_CUTSOLVER_FAILS, 0);
  std::vector<int> numCutsFromHeur(GlobalConstants::CutHeuristics::NUM_CUT_HEUR, 0);
  double bestObjValue = std::numeric_limits<double>::max();
  double bestSplitObjValue = std::numeric_limits<double>::max();
  double worstSplitObjValue = std::numeric_limits<double>::max();

#ifdef SHOULD_CALC_ACTIVE_FINAL_POINTS
  int totalNumTightFinalPointsGICs = 0;
  int totalNumTightFinalPointsSICs = 0;
  std::vector<int> numTightFinalPointsGICs;
  std::vector<int> numTightFinalPointsSICs;

  int totalNumTightFinalPointsActiveGICs = 0;
  int totalNumTightFinalPointsActiveSICs = 0;
  std::vector<int> numTightFinalPointsActiveGICs;
  std::vector<int> numTightFinalPointsActiveSICs;

  int totalNumTightFinalRaysGICs = 0;
  int totalNumTightFinalRaysSICs = 0;
  std::vector<int> numTightFinalRaysGICs;
  std::vector<int> numTightFinalRaysSICs;

  int totalNumTightFinalRaysActiveGICs = 0;
  int totalNumTightFinalRaysActiveSICs = 0;
  std::vector<int> numTightFinalRaysActiveGICs;
  std::vector<int> numTightFinalRaysActiveSICs;

  std::vector<std::vector<double> > cutCoeffsNB;
  std::vector<double> cutRHSNB;
#endif

  // Instantiate timers
  void initialize() {
    if (LOG_FILE.empty()) {
      LOG_FILE = GlobalConstants::INSTANCES_INFO_DEFAULT;
    }

    time_name.push_back(GlobalConstants::INIT_SOLVE_TIME);
    time_name.push_back(GlobalConstants::GEN_MSICS_TIME);
    time_name.push_back(GlobalConstants::APPLY_MSICS_TIME);
    time_name.push_back(GlobalConstants::CHOOSE_HPLANES_TIME);
    time_name.push_back(GlobalConstants::GEN_INTERSECTION_POINTS_TIME);
    time_name.push_back(GlobalConstants::GEN_VERTICES_TIME);
    time_name.push_back(GlobalConstants::GEN_PHA_TIME);
    time_name.push_back(GlobalConstants::APPLY_PHA_TIME);
    time_name.push_back(GlobalConstants::GEN_MPHA_TIME);
    time_name.push_back(GlobalConstants::APPLY_MPHA_TIME);
    time_name.push_back(GlobalConstants::COLLECT_STATS_TIME);
    for (int heur_ind = 0;
        heur_ind < GlobalConstants::CutHeuristics::NUM_CUT_HEUR; heur_ind++) {
      time_name.push_back(
          GlobalConstants::CutHeuristicsName[heur_ind] + "_TIME");
    }
    time_name.push_back(GlobalConstants::TOTAL_TIME);
    for (int t = 0; t < (int) time_name.size(); t++) {
      timeStats.register_name(time_name[t]);
    }
  } /* initialize */
} /* namespace GlobalVariables */
