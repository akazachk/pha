//============================================================================
// Name        : GlobalConstants.hpp
// Author      : akazachk
// Version     : 0.2014.01.10
// Copyright   : Your copyright notice
// Description : Header for GlobalConstants.hpp, which has miscellaneous useful definitions and functions used globally.
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#pragma once
#include <iostream>
#include <fstream>
#include <Stats.hpp>
#include <cmath>
#include <chrono>
#include "CglGICParam.hpp"

#ifdef TRACE
  #define SHOULD_WRITE_BRIEF_SOLN // Should brief solution be printed
  #define SHOULD_WRITE_SOLN // Should full solution be printed
  #define SHOULD_WRITE_INTPOINTS // Print intersection points and rays
  #define SHOULD_WRITE_HPLANES // Should hplanes activated be printed
  #define SHOULD_WRITE_SICS// Should SICs be printed
  #define SHOULD_WRITE_GICS // Should GICs be printed
  #define WRITE_LP_FILE // All possible LPs generated should be saved
  #define SHOULD_CALC_POINT_OBJ_DEPTH // Should we calculate point's objective violation?
  //#define SHOULD_CALC_ACTIVE_FINAL_POINTS // Should we calculate how many total final points are tight across all cuts?
  #define SHOULD_WRITE_PTSUMMARY // Print how good the points are (depth, etc.) into pointSummary
  #define SHOULD_WRITE_CUTSOLVER_ERRORS // Print errors encountered during creating of cuts
  //#define SHOULD_WRITE_SPLIT_FEAS_INFO // Print feasibility information about splits
  #define SHOULD_WRITE_TIME_INFO // Print times into a separate file
  #define SHOULD_WRITE_PARAMS // Print parameters used to a separate file
#else
  //#define SHOULD_CALC_ACTIVE_FINAL_POINTS // Should we calculate how many total final points are tight across all cuts?
#endif

#define error_msg(str, fmt, ...) \
  char str[500]; \
  snprintf(str, sizeof(str)/sizeof(char), "*** ERROR: %s:%d (%.1fs): " fmt, __FILE__, __LINE__, std::chrono::duration<double>(std::chrono::system_clock::now() - GlobalVariables::start_time).count(), ##__VA_ARGS__); \
  std::cerr << str

#define warning_msg(str, fmt, ...) \
  char str[500]; \
  snprintf(str, sizeof(str)/sizeof(char), "*** WARNING: %s:%d (%.1fs): " fmt, __FILE__, __LINE__, std::chrono::duration<double>(std::chrono::system_clock::now() - GlobalVariables::start_time).count(), ##__VA_ARGS__); \
  std::cout << str

namespace GlobalVariables {
  extern std::chrono::time_point<std::chrono::system_clock> start_time;
  extern std::FILE* log_file;
  extern std::string LOG_FILE;
  extern std::string prob_name;
  extern std::string out_f_name_stub;

  extern Stats timeStats;
  extern std::vector<std::string> time_name; // TODO Make const string vec 
  extern CglGICParam param;

  extern std::vector<int> numCutSolverFails;
  extern std::vector<int> numCutsFromHeur;
  extern double bestObjValue;
  extern double bestSplitObjValue;
  extern double worstSplitObjValue;

#ifdef SHOULD_CALC_ACTIVE_FINAL_POINTS
  extern int totalNumTightFinalPointsGICs;
  extern int totalNumTightFinalPointsSICs;
  extern std::vector<int> numTightFinalPointsGICs;
  extern std::vector<int> numTightFinalPointsSICs;

  extern int totalNumTightFinalPointsActiveGICs;
  extern int totalNumTightFinalPointsActiveSICs;
  extern std::vector<int> numTightFinalPointsActiveGICs;
  extern std::vector<int> numTightFinalPointsActiveSICs;

  extern int totalNumTightFinalRaysGICs;
  extern int totalNumTightFinalRaysSICs;
  extern std::vector<int> numTightFinalRaysGICs;
  extern std::vector<int> numTightFinalRaysSICs;

  extern int totalNumTightFinalRaysActiveGICs;
  extern int totalNumTightFinalRaysActiveSICs;
  extern std::vector<int> numTightFinalRaysActiveGICs;
  extern std::vector<int> numTightFinalRaysActiveSICs;

  // Below are needed to store NB cut coeffs and rhs, as we don't keep those otherwise
  extern std::vector<std::vector<double> > cutCoeffsNB;
  extern std::vector<double> cutRHSNB;
#endif

  void initialize();
} /* namespace GlobalVariables */

using namespace GlobalVariables;

namespace GlobalConstants {
  const std::string II_OUT_EXT = ".csv";
  const std::string INSTANCES_INFO_DEFAULT = "pha-info";

  // cut heuristic numbers
  enum CutHeuristics {
    SIC_CUT_GEN,
    CUT_VERTICES_CUT_HEUR,
    DUMMY_OBJ_CUT_HEUR,
    ITER_BILINEAR_CUT_HEUR,
    UNIT_VECTORS_CUT_HEUR,
    TIGHT_POINTS_CUT_HEUR,
    SPLIT_SHARE_CUT_HEUR,
    NUM_CUT_HEUR
  };

  const std::vector<std::string> CutHeuristicsName = { "SIC_CUT_GEN",
      "CUT_VERTICES_CUT_HEUR", "DUMMY_OBJ_CUT_HEUR", "ITER_BILINEAR_CUT_HEUR",
      "UNIT_VECTORS_CUT_HEUR", "TIGHT_POINTS_CUT_HEUR", "SPLIT_SHARE_CUT_HEUR"
      };

  // cut solver fails
  enum CutSolverFails {
    DUAL_CUTSOLVER_FAIL_IND,
    DUPLICATE_SIC_CUTSOLVER_FAIL_IND,
    DUPLICATE_GIC_CUTSOLVER_FAIL_IND,
    ORTHOGONALITY_CUTSOLVER_FAIL_IND,
    TIMELIMIT_CUTSOLVER_FAIL_IND,
    ITERATION_CUTSOLVER_FAIL_IND,
    ABANDONED_CUTSOLVER_FAIL_IND,
    NON_CUTTING_CUTSOLVER_FAIL_IND,
    SCALING_CUTSOLVER_FAIL_IND,
    UNSPECIFIED_CUTSOLVER_FAIL_IND,
    CUT_LIMIT_FAIL_IND,
    PRIMAL_CUTSOLVER_FAIL_IND,
    PRIMAL_CUTSOLVER_NO_OBJ_FAIL_IND,
    NUMERICAL_ISSUES_WARNING_NO_OBJ_FAIL_IND,
    NUMERICAL_ISSUES_NO_OBJ_FAIL_IND,
    NUM_CUTSOLVER_FAILS
  };

  const std::vector<std::string> CutSolverFailName = { "DUAL_CUT_SOLVER_FAILS",
      "DUPLICATE_SIC_FAILS", "DUPLICATE_GIC_FAILS", "ORTHOGONALITY_FAILS",
      "TIMELIMIT_FAILS", "ITERATION_FAILS", "ABANDONED_FAILS",
      "NON_CUTTING_FAILS", "SCALING_FAILS", "UNSPEC_FAILS", "CUT_LIMIT_FAILS",
      "PRIMAL_CUT_SOLVER_FAILS", "PRIMAL_CUT_SOLVER_NO_OBJ_FAILS",
      "NUMERICAL_ISSUES_WARNING_NO_OBJ", "NUMERICAL_ISSUES_NO_OBJ" };

  // For timing
  const std::string TOTAL_TIME = "TOTAL_TIME";
  const std::string INIT_SOLVE_TIME = "INIT_SOLVE_TIME";
  const std::string GEN_SUB_TIME = "GEN_SUB_TIME";
  const std::string GEN_MSICS_TIME = "GEN_MSICS_TIME";
  const std::string APPLY_MSICS_TIME = "APPLY_MSICS_TIME";
  const std::string LIFT_SICS_TIME = "LIFT_SICS_TIME";
  const std::string GEN_PHA_TIME = "GEN_PHA_TIME";
  const std::string LIFT_PHA_TIME = "LIFT_PHA_TIME";
  const std::string APPLY_PHA_TIME = "APPLY_PHA_TIME";
  const std::string CHOOSE_HPLANES_TIME = "CHOOSE_HPLANES_TIME";
  const std::string GEN_INTERSECTION_POINTS_TIME = "GEN_INTERSECTION_POINTS_TIME";
  const std::string GEN_VERTICES_TIME = "GEN_VERTICES_TIME";
  const std::string GEN_MPHA_TIME = "GEN_MPHA_TIME";
  const std::string APPLY_MPHA_TIME = "APPLY_MPHA_TIME";
  const std::string COLLECT_STATS_TIME = "COLLECT_STATS_TIME";

  const int SPACE = 15;
  const char SEP = ',';

  const double MIN_RAYS_RATIO = .1;

  // Should we use lopsided splits or not?
  const bool YESLOP = true;
  const bool USE_EMPTY = true; // Are we okay when both sides are empty? This overrides YESLOP.

  // This is for what we multiply by in the big-M constraint for the split LPs, to make them work on lopsided splits
  const int MMULT = 0;

  const bool ACTIVATE_BASIC_BOUNDS = true;
  const bool ACTIVATE_NB_BOUNDS = true;
  const bool SECOND_ACTIVATE_BASIC_BOUNDS = true;
  const bool SECOND_ACTIVATE_NB_BOUNDS = true;
  const bool CUT_PARALLEL_RAYS = false;
  const bool CHECK_FOR_DUP_RAYS = true;

  const bool NO_GEN_CUTS = false;

  const bool RANK2_ACTIVATION = false;

  /***********************************************************************/
  /**
   * @brief Writes error and closes file myfile.
   */
  inline void writeErrorToII(std::string text, FILE *myfile) {
    if (myfile == NULL)
      return;
    fprintf(myfile, "||%c%s", ',', text.c_str());
    fclose(myfile);
  }

  // Nearest integer
  inline double nearestInt(const double x) {
    return std::floor((x) + 0.5);
  }

  /***********************************************************************/
  /**
   * @brief Checks for equivalence to arbitrary value
   * @param val1 :: value to check
   * @param val2 :: value to check against
   * @param eps :: zero tolerance
   */
  inline bool isVal(const double val1, const double val2,
      const double eps = param.getEPS()) {
    return (std::abs((val1) - (val2)) <= eps);
  }

  /***********************************************************************/
  /**
   * @brief Checks for equivalence to zero
   */
  inline bool isZero(const double val, const double eps = param.getEPS()) {
    return isVal(val, 0.0, eps);
  }

  /***********************************************************************/
  /**
   * @brief Checks for equivalence to one
   * @param val :: value to check
   */
  inline bool isOne(const double val) {
    return isVal(val, 1.0);
  }

  /***********************************************************************/
  /**
   * @brief Checks for less than arbitrary value
   * @param val1 :: value to check
   * @param val2 :: value to check against
   * @param eps :: zero tolerance
   */
  inline bool lessThanVal(const double val1, const double val2, const double eps =
      param.getEPS()) {
    return (val1 < (val2 - eps));
  }

  /***********************************************************************/
  /**
   * @brief Checks for greater than arbitrary value
   * @param val1 :: value to check
   * @param val2 :: value to check against
   * @param eps :: zero tolerance
   */
  inline bool greaterThanVal(const double val1, const double val2,
      const double eps = param.getEPS()) {
    return (val1 > (val2 + eps));
  }

  /***********************************************************************/
  /**
   * @brief Checks for greater than infinity - eps
   * @param val :: value to check
   * @param infinity :: value of infinity to use (default same as OsiClp)
   * @param eps :: zero tolerance
   */
  inline bool isInfinity(const double val, const double infinity =
      __DBL_MAX__, const double eps = param.getEPS()) {
    return (val >= (infinity - eps));
  }

  /***********************************************************************/
  /**
   * @brief Checks for less than -infinity + eps
   * @param val :: value to check
   * @param infinity :: value of infinity to use (default same as OsiClp)
   * @param eps :: zero tolerance
   */
  inline bool isNegInfinity(const double val, const double infinity =
      __DBL_MAX__, const double eps = param.getEPS()) {
    return (val <= (-infinity + eps));
  }

  // Function for checking equality with user tolerance
  /*
  inline bool areEqual(const double x, const double y, const double epsAbs =
      EQUALITY_EPS_ABS, const double epsRel = EQUALITY_EPS_REL) {
    return isVal(x, y, epsAbs) || isVal(x, y, epsRel * std::abs(x))
      || isVal(x, y, epsRel * std::abs(y));
  }
  // Function for checking if a number is integer
  inline bool isIntegerValue(const double x) {
    return isVal(x, nearestInt(x), INTEGRALITY_EPS_ABS)
      || isVal(x, nearestInt(x), INTEGRALITY_EPS_REL * std::abs(x));
  }
  */

  inline const std::string stringValue(const int value,
      const char* format = "%d") {
    if (!isInfinity(std::abs(value))) {
      char temp[500];
      snprintf(temp, sizeof(temp) / sizeof(char), format, value);
      std::string tmp(temp);
      return tmp;
    } else {
      if (greaterThanVal(value, 0.0)) {
        const std::string infty = "\'inf\'";
        return infty;
      } else {
        const std::string neg_infty = "\'-inf\'";
        return neg_infty;
      }
    }
  }

  inline const std::string stringValue(const double value, const char* format =
      "%f") {
    if (!isInfinity(std::abs(value))) {
      char temp[500];
      snprintf(temp, sizeof(temp) / sizeof(char), format, value);
      std::string tmp(temp);
      return tmp;
    } else {
      if (greaterThanVal(value, 0.0)) {
        const std::string infty = "\'+inf\'";
        return infty;
      } else {
        const std::string neg_infty = "\'-inf\'";
        return neg_infty;
      }
    }
  }
} /* namespace GlobalConstants */

using namespace GlobalConstants;
