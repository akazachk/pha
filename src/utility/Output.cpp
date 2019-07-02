//============================================================================
// Name        : Output.cpp
// Author      : akazachk
// Version     : 0.2013.08.15
// Copyright   : Your copyright notice
// Description : 
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

//#include <iostream>
#include <fstream>
#include <cmath> // isnan
#include <algorithm> // For max_element, min_element
#include "Utility.hpp"
#include "Output.hpp"

/*************/
/* VARIABLES */
/*************/
const int countParamInfoEntries = ParamIndices::NUM_PARAMS + 4; // +1 for ortho, +1 for eps, +1 for rayeps, +1 for total time
const int countLPRelaxEntries = 23;
const int countGICPointInfoEntries = 38;
const int countRootBoundEntries = 13;
const int countGICCutHeurInfoEntries = 2 + 1 + 2 * CutHeuristics::NUM_CUT_HEUR
    + CutSolverFails::NUM_CUTSOLVER_FAILS; // The +1 for "total num cuts" and "total obj tried"
const int countRoundInfoEntries = 18;

/***********************************************************************/
/**
 * @brief Adds the header lines to the file
 */
void writeHeaderToII(FILE *myfile) {
  if (myfile == NULL)
    return;
  std::string tmpstring = "";
  fprintf(myfile, "%c", SEP);
  fprintf(myfile, "%s", "PARAM INFO");
  tmpstring.assign(countParamInfoEntries, SEP);
  fprintf(myfile, "%s", tmpstring.c_str());
  fprintf(myfile, "%s", "ORIG PROB");
  tmpstring.assign(countLPRelaxEntries, SEP);
  fprintf(myfile, "%s", tmpstring.c_str());
  if (param.paramVal[ParamIndices::PHA_PARAM_IND]) {
    for (int i = 0; i < param.getParamVal(ParamIndices::ROUNDS_PARAM_IND); i++) {
      fprintf(myfile, "%s ROUND %d", "GIC POINT INFO", i+1);
      tmpstring.assign(countGICPointInfoEntries, SEP);
      fprintf(myfile, "%s", tmpstring.c_str());
    }
  }
  fprintf(myfile, "%s", "GIC CUT HEUR INFO");
  tmpstring.assign(countGICCutHeurInfoEntries, SEP);
  fprintf(myfile, "%s", tmpstring.c_str());
  fprintf(myfile, "%s", "ROOT BOUND");
  tmpstring.assign(countRootBoundEntries, SEP);
  fprintf(myfile, "%s", tmpstring.c_str());
  fprintf(myfile, "%s", "ROUND INFO");
  tmpstring.assign(countRoundInfoEntries, SEP);
  fprintf(myfile, "%s", tmpstring.c_str());
#ifdef SHOULD_USE_GUROBI
  if (param.getTEMP() == 1) {
    fprintf(myfile, "%s,,", "ERRORS");
  }
#endif
  fprintf(myfile, "%s", "TIME INFO");
  fprintf(myfile, "%c", SEP);
  fprintf(myfile, "\n");

  fprintf(myfile, "%s%c", "INSTANCE", SEP);
  ////////////////////// PARAM INFO
  for (int param_ind = 0; param_ind < ParamIndices::NUM_PARAMS; param_ind++) {
    fprintf(myfile, "%s%c", ParamName[param_ind].c_str(), SEP);
  }
  fprintf(myfile, "%s%c", "MIN_ORTHOGONALITY", SEP);
  fprintf(myfile, "%s%c", "EPS", SEP);
  fprintf(myfile, "%s%c", "RAYEPS", SEP);
  fprintf(myfile, "%s%c", "TIMELIMIT", SEP);
  ////////////////////// ORIG PROB
  fprintf(myfile, "%s%c", "NUM ROWS", SEP); // 1
  fprintf(myfile, "%s%c", "NUM COLS", SEP); // 2
  fprintf(myfile, "%s%c", "NUM EQ ROWS", SEP); // 3
  fprintf(myfile, "%s%c", "NUM INEQ ROWS", SEP); // 4
  fprintf(myfile, "%s%c", "NUM BOUND ROWS", SEP); // 5
  fprintf(myfile, "%s%c", "NUM ASSIGN ROWS", SEP); // 6
  fprintf(myfile, "%s%c", "INTEGERS NON-BIN", SEP); // 7
  fprintf(myfile, "%s%c", "BINARIES", SEP); // 8
  fprintf(myfile, "%s%c", "CONTINUOUS", SEP); // 9
  fprintf(myfile, "%s%c", "NUM BASIC", SEP); // 10
  fprintf(myfile, "%s%c", "NUM NON-BASIC", SEP); // 11
//  fprintf(myfile, "%s%c", "NUM ORIG BASIC", SEP);
  fprintf(myfile, "%s%c", "NUM ORIG NON-BASIC", SEP); // 13
  fprintf(myfile, "%s%c", "NUM PRIMAL DEGEN", SEP); // 20
  fprintf(myfile, "%s%c", "NUM DUAL DEGEN", SEP); // 21 // 20 (removed 17) *
  fprintf(myfile, "%s%c", "FRAC CORE", SEP); // 14
  fprintf(myfile, "%s%c", "A NONZERO", SEP); // 15
  fprintf(myfile, "%s%c", "A-DENSITY", SEP); // 16
  fprintf(myfile, "%s%c", "MIN RAYVAL NO EPS", SEP); // 12
  fprintf(myfile, "%s%c", "MIN RAYVAL", SEP); // 19
  fprintf(myfile, "%s%c", "LP OBJ", SEP); // 18 // 21 above
  fprintf(myfile, "%s%c", "WORST SPLIT OBJ", SEP); // 18 // 21 above
  fprintf(myfile, "%s%c", "BEST SPLIT OBJ", SEP); // 18 // 21 above
  fprintf(myfile, "%s%c", "IP OBJ", SEP); // 17
  if (param.paramVal[ParamIndices::PHA_PARAM_IND] > 0) {
    for (int i = 0; i < param.getParamVal(ParamIndices::ROUNDS_PARAM_IND); i++) {
    ////////////////////// GIC POINT INFO
    fprintf(myfile, "%s%c", "HPLANE HEUR", SEP); // 1
    fprintf(myfile, "%s%c", "NUM HPLANES PER RAY", SEP); // 2
    fprintf(myfile, "%s%c", "NUM RAYS C1", SEP); // 19
    fprintf(myfile, "%s%c", "AVG NUM PARALLEL RAYS C1", SEP); // 26
    fprintf(myfile, "%s%c", "AVG NUM SIC FINAL POINTS", SEP); // 18 // 19 above
    fprintf(myfile, "%s%c", "AVG NUM SIC FINAL RAYS", SEP); // 18 // 19 above
    fprintf(myfile, "%s%c", "NUM RAYS TO CUT", SEP); // 3
    fprintf(myfile, "%s%c", "NUM CUT GEN SETS", SEP); // 4 // 5 missing
    fprintf(myfile, "%s%c", "NUM POINTS (TOTAL)", SEP); // 6
    fprintf(myfile, "%s%c", "NUM POINTS (AVG)", SEP); // 7
    fprintf(myfile, "%s%c", "NUM POINTS (MIN)", SEP); // 8
    fprintf(myfile, "%s%c", "NUM POINTS (MAX)", SEP); // 9
    fprintf(myfile, "%s%c", "NUM FINAL POINTS (TOTAL)", SEP); // 10
    fprintf(myfile, "%s%c", "NUM FINAL POINTS (AVG)", SEP); // 11
    fprintf(myfile, "%s%c", "NUM FINAL POINTS (MIN)", SEP); // 12
    fprintf(myfile, "%s%c", "NUM FINAL POINTS (MAX)", SEP); // 13
    fprintf(myfile, "%s%c", "SIC AVG TIGHT FINAL POINTS", SEP); // 27
    fprintf(myfile, "%s%c", "GIC AVG TIGHT FINAL POINTS", SEP); // 28
      fprintf(myfile, "%s%c", "ACTIVE SIC AVG TIGHT FINAL POINTS", SEP); // 31
      fprintf(myfile, "%s%c", "ACTIVE GIC AVG TIGHT FINAL POINTS", SEP); // 32
    fprintf(myfile, "%s%c", "NUM RAYS (TOTAL)", SEP); // 14
    fprintf(myfile, "%s%c", "NUM RAYS (AVG)", SEP); // 15
    fprintf(myfile, "%s%c", "NUM RAYS (MIN)", SEP); // 16
    fprintf(myfile, "%s%c", "NUM RAYS (MAX)", SEP); // 17 // 18 above
    fprintf(myfile, "%s%c", "NUM FINAL RAYS (TOTAL)", SEP); // 10
    fprintf(myfile, "%s%c", "NUM FINAL RAYS (AVG)", SEP); // 11
    fprintf(myfile, "%s%c", "NUM FINAL RAYS (MIN)", SEP); // 12
    fprintf(myfile, "%s%c", "NUM FINAL RAYS (MAX)", SEP); // 13
    fprintf(myfile, "%s%c", "SIC AVG TIGHT FINAL RAYS", SEP); // 29
    fprintf(myfile, "%s%c", "GIC AVG TIGHT FINAL RAYS", SEP); // 30
      fprintf(myfile, "%s%c", "ACTIVE SIC AVG TIGHT FINAL RAYS", SEP); // 33
      fprintf(myfile, "%s%c", "ACTIVE GIC AVG TIGHT FINAL RAYS", SEP); // 34 *
    fprintf(myfile, "%s%c", "SIC DEPTH POINTS (AVG)", SEP); // 20
    fprintf(myfile, "%s%c", "SIC DEPTH POINTS (MIN)", SEP); // 21
    fprintf(myfile, "%s%c", "SIC DEPTH POINTS (MAX)", SEP); // 22
    fprintf(myfile, "%s%c", "OBJ DEPTH POINTS (AVG)", SEP); // 23
    fprintf(myfile, "%s%c", "OBJ DEPTH POINTS (MIN)", SEP); // 24
    fprintf(myfile, "%s%c", "OBJ DEPTH POINTS (MAX)", SEP); // 25 // 26 27 above
    } /* print this as many times as there are rounds */
  } /* only print GIC point info when pha is used */
  //////////////////// CUT HEUR INFO
  fprintf(myfile, "%s%c", "NUM CUTS", SEP); // 1
  fprintf(myfile, "%s%c", "NUM CGS", SEP); // 2
  // Cut heuristic stats
  for (int heur_ind = 0; heur_ind < CutHeuristics::NUM_CUT_HEUR; heur_ind++) { // 3 * # cutHeur
    fprintf(myfile, "%s%c", CutHeuristicsName[heur_ind].c_str(), SEP);
    fprintf(myfile, "ACTIVE %s%c", CutHeuristicsName[heur_ind].c_str(),
        SEP);
  }
  // Cut solver fails
  fprintf(myfile, "%s%c", "OBJ", SEP); // 5
  for (int fail_ind = 0; fail_ind < CutSolverFails::NUM_CUTSOLVER_FAILS;
      fail_ind++) {
    fprintf(myfile, "%s%c", CutSolverFailName[fail_ind].c_str(), SEP);
  }
  ////////////////////// ROOT NODE BOUNDS
  fprintf(myfile, "%s%c", "INIT BOUND", SEP); // 1
  fprintf(myfile, "%s%c", "NUM SIC", SEP); // 2
  fprintf(myfile, "%s%c", "SIC BOUND", SEP); // 3
  fprintf(myfile, "%s%c", "NUM PHA", SEP); // 4
  fprintf(myfile, "%s%c", "PHA BOUND", SEP); // 5
  fprintf(myfile, "%s%c", "ALL CUTS (WITH SIC)", SEP); // 6
  fprintf(myfile, "%s%c", "GIC OVER SIC", SEP); // 7
  fprintf(myfile, "%s%c", "ALL OVER SIC", SEP); // 8
  fprintf(myfile, "%s%c", "ACTIVE SIC", SEP); // 9
  fprintf(myfile, "%s%c", "ACTIVE PHA", SEP); // 10
  fprintf(myfile, "%s%c", "SIC % GAP CLOSED", SEP); // 11
  fprintf(myfile, "%s%c", "GIC % GAP CLOSED", SEP); // 12
  fprintf(myfile, "%s%c", "SIC+GIC % GAP CLOSED", SEP); // 13

  //////////////////// ROUND INFO
  fprintf(myfile, "%s%c", "TOTAL NUM SICS", SEP); // 1
  fprintf(myfile, "%s%c", "NUM SIC ROUNDS", SEP); // 2
  fprintf(myfile, "%s%c", "FINAL SIC BOUND", SEP); // 3
  fprintf(myfile, "%s%c", "PHA+ BOUND RD 1", SEP); // 4
  fprintf(myfile, "%s%c", "PHA+ BOUND RD 2", SEP); // 5
  fprintf(myfile, "%s%c", "PHA+ BOUND RD 5", SEP); // 6
//  fprintf(myfile, "%s%c", "LP BASIS COND", SEP); // 10
//  fprintf(myfile, "%s%c", "SIC BASIS COND INIT", SEP); // 11
//  fprintf(myfile, "%s%c", "SIC BASIS COND FINAL", SEP); // 12
//  fprintf(myfile, "%s%c", "PHA BASIS COND RD 1", SEP); // 21
//  fprintf(myfile, "%s%c", "PHA BASIS COND RD 2", SEP); // 22
//  fprintf(myfile, "%s%c", "PHA BASIS COND RD 5", SEP); // 23
//  fprintf(myfile, "%s%c", "PHA BASIS COND FINAL", SEP); // 24
//  fprintf(myfile, "%s%c", "PHA+ BASIS COND RD 1", SEP); // 25
//  fprintf(myfile, "%s%c", "PHA+ BASIS COND RD 2", SEP); // 26
//  fprintf(myfile, "%s%c", "PHA+ BASIS COND RD 5", SEP); // 27
//  fprintf(myfile, "%s%c", "PHA+ BASIS COND FINAL", SEP); // 28
  fprintf(myfile, "%s%c", "MIN SIC ORTHOGONALITY INIT", SEP); // 29
  fprintf(myfile, "%s%c", "MAX SIC ORTHOGONALITY INIT", SEP); // 30
  fprintf(myfile, "%s%c", "AVG SIC ORTHOGONALITY INIT", SEP); // 31
  fprintf(myfile, "%s%c", "MIN PHA ORTHOGONALITY INIT", SEP); // 35
  fprintf(myfile, "%s%c", "MAX PHA ORTHOGONALITY INIT", SEP); // 36
  fprintf(myfile, "%s%c", "AVG PHA ORTHOGONALITY INIT", SEP); // 37
  fprintf(myfile, "%s%c", "MIN SIC ORTHOGONALITY FINAL", SEP); // 38
  fprintf(myfile, "%s%c", "MAX SIC ORTHOGONALITY FINAL", SEP); // 39
  fprintf(myfile, "%s%c", "AVG SIC ORTHOGONALITY FINAL", SEP); // 40
  fprintf(myfile, "%s%c", "MIN PHA ORTHOGONALITY FINAL", SEP); // 44
  fprintf(myfile, "%s%c", "MAX PHA ORTHOGONALITY FINAL", SEP); // 45
  fprintf(myfile, "%s%c", "AVG PHA ORTHOGONALITY FINAL", SEP); // 46 - 17 = 29 - 11 = 18*

  ////////////////////// ERRORS
#ifdef SHOULD_USE_GUROBI
  if (param.getTEMP() == 1) {
    fprintf(myfile, "%s%c", "NUM ERRORS", SEP);
    fprintf(myfile, "%s%c", "MAX ERROR", SEP);
  }
#endif

  ////////////////////// TIME
  for (int t = 0; t < (int) GlobalVariables::time_name.size(); t++) {
    fprintf(myfile, "%s%c", GlobalVariables::time_name[t].c_str(), SEP);
  }
  ////////////////////// newline
  fprintf(myfile, "\n");
  fflush(myfile);
} /* writeHeaderToII */

/***********************************************************************/
/**
 * @brief Writes problem name to start a new instance info line
 */
void writeProbNameToII(std::string &prob_name, FILE *myfile) {
  if (myfile == NULL)
    return;
  // Instance name
  fprintf(myfile, "%s%c", prob_name.c_str(), SEP);
  fflush(myfile);
} /* writeProbNameToInstanceInfo */

/***********************************************************************/
/**
 * Writes which settings we are using to instance info
 */
void writeParamInfoToII(const CglGICParam& param, FILE* myfile) {
  if (myfile == NULL)
    return;
  for (int param_ind = 0; param_ind < ParamIndices::NUM_PARAMS; param_ind++) {
    fprintf(myfile, "%d%c", param.paramVal[param_ind], SEP);
  }
  fprintf(myfile, "%e%c", param.getMINORTHOGONALITY(), SEP);
  fprintf(myfile, "%e%c", param.getEPS(), SEP);
  fprintf(myfile, "%e%c", param.getRAYEPS(), SEP);
  fprintf(myfile, "%e%c", param.getTIMELIMIT(), SEP);
  fflush(myfile);
} /* writeParamInfoToII */

/***********************************************************************/
/**
 * @brief Appends information about LP relaxation solution to prob
 */
void writeSolnInfoToII(const SolutionInfo& probSolnInfo,
    FILE *myfile) {
  if (myfile == NULL)
    return;
  // Print original problem info:
  //  NUM ROWS, NUM COLS, NUM EQ ROWS, NUM INEQ ROWS, NUM BOUND ROWS, NUM ASSIGN ROWS,
  //  INTEGERS NON-BIN, BINARIES, CONTINUOUS, NUM BASIC, NUM NON-BASIC,
  //Deleted //NUM ORIG BASIC,
  //  NUM ORIG NON-BASIC, FRAC CORE, A NONZERO, A-DENSITY, INTEGER OBJ, LP OBJ
//  if (probSolnInfo != NULL) {
  fprintf(myfile, "%d%c", probSolnInfo.numRows, SEP);
  fprintf(myfile, "%d%c", probSolnInfo.numCols, SEP);
  fprintf(myfile, "%d%c", probSolnInfo.numEqRows, SEP);
  fprintf(myfile, "%d%c", probSolnInfo.numIneqRows, SEP);
  fprintf(myfile, "%d%c", probSolnInfo.numBoundRows, SEP);
  fprintf(myfile, "%d%c", probSolnInfo.numAssignRows, SEP);
  fprintf(myfile, "%d%c", probSolnInfo.numIntegerNonBinary, SEP);
  fprintf(myfile, "%d%c", probSolnInfo.numBinary, SEP);
  fprintf(myfile, "%d%c", probSolnInfo.numContinuous, SEP);
  fprintf(myfile, "%d%c", probSolnInfo.numBasicCols, SEP);
  fprintf(myfile, "%d%c", probSolnInfo.numNB, SEP);
//    fprintf(myfile, "%d%c", (int) probSolnInfo.basicOrigVarIndex.size(),
//        SEP);
  fprintf(myfile, "%d%c", (int) probSolnInfo.nonBasicOrigVarIndex.size(),
      SEP);
  fprintf(myfile, "%d%c", probSolnInfo.numPrimalDegeneratePivots, SEP);
  fprintf(myfile, "%d%c", probSolnInfo.numDualDegeneratePivots, SEP);
  fprintf(myfile, "%d%c", (int) probSolnInfo.fractionalCore.size(), SEP);
  fprintf(myfile, "%d%c", probSolnInfo.numANonzero, SEP);
  fprintf(myfile, "%2.3f%c",
      (double) probSolnInfo.numANonzero
          / (probSolnInfo.numCols * probSolnInfo.numRows), SEP);
  fprintf(myfile, "%1.10f%c", probSolnInfo.minAbsRayValNoEPS, SEP);
  fprintf(myfile, "%1.10f%c", probSolnInfo.minAbsRayVal, SEP);
  fprintf(myfile, "%.20f%c", probSolnInfo.LPopt, SEP);
  fprintf(myfile, "%s%c", stringValue(GlobalVariables::worstSplitObjValue, "%.20f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(GlobalVariables::bestSplitObjValue, "%.20f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(GlobalVariables::bestObjValue, "%.20f").c_str(), SEP);
  fflush(myfile);
} /* writeSolnInfoToInstanceInfo */

/***********************************************************************/
/**
 * Writes which settings we are using to instance info,
 * as well as the number of points and rays generated
 */
void writeGICPointInfoToII(const int hplaneHeur,
    const int numHplanesPerRay, const int numRaysOfC1,
    const double avgNumParallelRays, 
    const double avg_num_SIC_final_points,
    const double avg_num_SIC_final_rays,
    const int numRaysToBeCut, const int numCutGenSets,
    std::vector<int> &numPoints, std::vector<int> &numFinalPoints,
    std::vector<int> &numRays, std::vector<int>& numFinalRays, 
    const int totalNumPoints, const int totalNumFinalPoints, 
    const int totalNumRays, const int totalNumFinalRays,
    const double avgNumTightFinalPointsSICs,
    const double avgNumTightFinalPointsGICs,
    const double avgNumTightFinalRaysSICs,
    const double avgNumTightFinalRaysGICs,
    const double avgNumTightFinalPointsActiveSICs,
    const double avgNumTightFinalPointsActiveGICs,
    const double avgNumTightFinalRaysActiveSICs,
    const double avgNumTightFinalRaysActiveGICs,
    const std::vector<double>& avgSICDepth,
    const std::vector<double>& minSICDepth,
    const std::vector<double>& maxSICDepth,
    const std::vector<double>& avgObjDepth,
    const std::vector<double>& minObjDepth,
    const std::vector<double>& maxObjDepth, FILE* myfile) {
  if (myfile == NULL)
    return;
  fprintf(myfile, "%d%c", hplaneHeur, SEP); // 1
  fprintf(myfile, "%d%c", numHplanesPerRay, SEP); // 2
  fprintf(myfile, "%d%c", numRaysOfC1, SEP); // 3
  fprintf(myfile, "%2.3f%c", avgNumParallelRays, SEP); // 4
  fprintf(myfile, "%2.3f%c", avg_num_SIC_final_points, SEP); // 5
  fprintf(myfile, "%2.3f%c", avg_num_SIC_final_rays, SEP); // 6
  fprintf(myfile, "%d%c", numRaysToBeCut, SEP); // 7
  fprintf(myfile, "%d%c", numCutGenSets, SEP); // 8
  fprintf(myfile, "%d%c", totalNumPoints, SEP); // 9
  fprintf(myfile, "%2.3f%c", computeAverage(numPoints), SEP); // 10
  fprintf(myfile, "%d%c", *computeMin(numPoints), SEP); // 11
  fprintf(myfile, "%d%c", *computeMax(numPoints), SEP); // 12
  fprintf(myfile, "%d%c", totalNumFinalPoints, SEP); // 13
  fprintf(myfile, "%2.3f%c", computeAverage(numFinalPoints), SEP); // 14
  fprintf(myfile, "%d%c", *computeMin(numFinalPoints), SEP); // 15
  fprintf(myfile, "%d%c", *computeMax(numFinalPoints), SEP); // 16
  fprintf(myfile, "%2.3f%c", avgNumTightFinalPointsSICs, SEP); // 17
  fprintf(myfile, "%2.3f%c", avgNumTightFinalPointsGICs, SEP); // 18
  fprintf(myfile, "%2.3f%c", avgNumTightFinalPointsActiveSICs, SEP); // 17
  fprintf(myfile, "%2.3f%c", avgNumTightFinalPointsActiveGICs, SEP); // 18
  fprintf(myfile, "%d%c", totalNumRays, SEP); // 19
  fprintf(myfile, "%2.3f%c", computeAverage(numRays), SEP); // 20
  fprintf(myfile, "%d%c", *computeMin(numRays), SEP); // 21
  fprintf(myfile, "%d%c", *computeMax(numRays), SEP); // 22
  fprintf(myfile, "%d%c", totalNumFinalRays, SEP); // 23
  fprintf(myfile, "%2.3f%c", computeAverage(numFinalRays), SEP); // 24
  fprintf(myfile, "%d%c", *computeMin(numFinalRays), SEP); // 25
  fprintf(myfile, "%d%c", *computeMax(numFinalRays), SEP); // 26
  fprintf(myfile, "%2.3f%c", avgNumTightFinalRaysSICs, SEP); // 27
  fprintf(myfile, "%2.3f%c", avgNumTightFinalRaysGICs, SEP); // 28
  fprintf(myfile, "%2.3f%c", avgNumTightFinalRaysActiveSICs, SEP); // 27
  fprintf(myfile, "%2.3f%c", avgNumTightFinalRaysActiveGICs, SEP); // 28
  fprintf(myfile, "%2.3f%c", computeAverage(avgSICDepth), SEP); // 29
  fprintf(myfile, "%2.3f%c", computeAverage(minSICDepth), SEP); // 30
  fprintf(myfile, "%2.3f%c", computeAverage(maxSICDepth), SEP); // 31
  fprintf(myfile, "%2.3f%c", computeAverage(avgObjDepth), SEP); // 32
  fprintf(myfile, "%2.3f%c", computeAverage(minObjDepth), SEP); // 33
  fprintf(myfile, "%2.3f%c", computeAverage(maxObjDepth), SEP); // 34
  fflush(myfile);
} /* writeGICPointInfoToII */

void writeRootBoundToII(const double lp_opt, const int numSICs,
    const double sic_opt, 
    const int numPHA, const double pha_opt, 
    const double sic_pha_opt, 
    const int num_active_sics_sicsolver,
    const int num_active_pha_sicsolver,
    FILE *myfile) {
  if (myfile == NULL)
    return;

  const bool GICsOverSICs = greaterThanVal(pha_opt, sic_opt);
  const bool allOverSICs = greaterThanVal(sic_pha_opt, sic_opt);

  //////////
  // Root node bound
  //////////
  fprintf(myfile, "%.20f%c", lp_opt, SEP);  // 1
  fprintf(myfile, "%d%c", numSICs, SEP);  // 2
  fprintf(myfile, "%s%c", stringValue(sic_opt, "%2.20f").c_str(), SEP);  // 3
  fprintf(myfile, "%d%c", numPHA, SEP);  // 4
  fprintf(myfile, "%s%c", stringValue(pha_opt, "%2.20f").c_str(), SEP);  // 5
  fprintf(myfile, "%s%c", stringValue(sic_pha_opt, "%2.20f").c_str(), SEP);  // 6
  fprintf(myfile, "%s%c", stringValue(
      GICsOverSICs * (pha_opt - sic_opt), "%2.6f").c_str(), SEP);  // 7
  fprintf(myfile, "%s%c", stringValue(
      allOverSICs * (sic_pha_opt - sic_opt), "%2.6f").c_str(), SEP);  // 8
  fprintf(myfile, "%d%c", num_active_sics_sicsolver, SEP);  // 9
  fprintf(myfile, "%d%c", num_active_pha_sicsolver, SEP);  // 10
  if (!isInfinity(std::abs(GlobalVariables::bestObjValue))) {
    const double ip_opt = GlobalVariables::bestObjValue;
    fprintf(myfile, "%2.6f%c", 100. * (sic_opt - lp_opt) / (ip_opt - lp_opt), SEP);  // 11
    fprintf(myfile, "%2.6f%c", 100. * (pha_opt - lp_opt) / (ip_opt - lp_opt), SEP);  // 12
    fprintf(myfile, "%2.6f%c", 100. * (sic_pha_opt - lp_opt) / (ip_opt - lp_opt), SEP);  // 13
  } else {
    fprintf(myfile, "%c",SEP);
    fprintf(myfile, "%c",SEP);
    fprintf(myfile, "%c",SEP);
  }
  fflush(myfile);
} /* writeRootBoundToII */

void writeRoundInfoToII(const double sic_pha_opt,
    const int total_num_sics, const int numRoundsSIC, const double final_sic_bound,
    const std::vector<double>& boundByRoundPHAPlus,
//    const double LPBasisCond,
//    const std::vector<double>& basisCondByRoundSIC,
//    const std::vector<double>& basisCondByRoundPHA,
//    const std::vector<double>& basisCondByRoundPHAPlus,
    const double minSICOrthogonalityInit, const double maxSICOrthogonalityInit,
    const double avgSICOrthogonalityInit, 
    const double minPHAOrthogonalityInit, const double maxPHAOrthogonalityInit,
    const double avgPHAOrthogonalityInit, const double minSICOrthogonalityFinal,
    const double maxSICOrthogonalityFinal,
    const double avgSICOrthogonalityFinal,
    const double minPHAOrthogonalityFinal,
    const double maxPHAOrthogonalityFinal,
    const double avgPHAOrthogonalityFinal, FILE *myfile) {
  if (myfile == NULL)
    return;
//  const int numRoundsSICCond = basisCondByRoundSIC.size();
  const int numRoundsPHA = boundByRoundPHAPlus.size();
  // SICs
//  const double SICBasisCondInit =
//      (numRoundsSICCond > 0) ? basisCondByRoundSIC[0] : 0.;
//  const double SICBasisCondFinal =
//      (numRoundsSICCond > 0) ? basisCondByRoundSIC[numRoundsSICCond - 1] : 0.;
  // PHA
  const double PHAPlusBoundRd1 =
      (numRoundsPHA > 0) ? boundByRoundPHAPlus[0] : sic_pha_opt;
  const double PHAPlusBoundRd2 =
      (numRoundsPHA > 1) ? boundByRoundPHAPlus[1] : sic_pha_opt;
  const double PHAPlusBoundRd5 =
      (numRoundsPHA > 4) ? boundByRoundPHAPlus[4] : sic_pha_opt;
//  const double PHABasisCondRd1 =
//      (numRoundsPHA > 0) ? basisCondByRoundPHA[0] : 0.;
//  const double PHABasisCondRd2 =
//      (numRoundsPHA > 1) ? basisCondByRoundPHA[1] :
//      (numRoundsPHA > 0) ? basisCondByRoundPHA[numRoundsPHA - 1] : 0.;
//  const double PHABasisCondRd5 =
//      (numRoundsPHA > 4) ? basisCondByRoundPHA[4] :
//      (numRoundsPHA > 0) ? basisCondByRoundPHA[numRoundsPHA - 1] : 0.;
//  const double PHABasisCondFinal =
//      (numRoundsPHA > 0) ? basisCondByRoundPHA[numRoundsPHA - 1] : 0.;
//  const double PHAPlusBasisCondRd1 =
//      (numRoundsPHA > 0) ? basisCondByRoundPHAPlus[0] : 0.;
//  const double PHAPlusBasisCondRd2 =
//      (numRoundsPHA > 1) ? basisCondByRoundPHAPlus[1] :
//      (numRoundsPHA > 0) ? basisCondByRoundPHAPlus[numRoundsPHA - 1] : 0.;
//  const double PHAPlusBasisCondRd5 =
//      (numRoundsPHA > 4) ? basisCondByRoundPHAPlus[4] :
//      (numRoundsPHA > 0) ? basisCondByRoundPHAPlus[numRoundsPHA - 1] : 0.;
//  const double PHAPlusBasisCondFinal =
//      (numRoundsPHA > 0) ? basisCondByRoundPHAPlus[numRoundsPHA - 1] : 0.;
  fprintf(myfile, "%d%c", total_num_sics, SEP);
  fprintf(myfile, "%d%c", numRoundsSIC, SEP);
  fprintf(myfile, "%s%c", stringValue(final_sic_bound, "%2.6f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(PHAPlusBoundRd1, "%2.6f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(PHAPlusBoundRd2, "%2.6f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(PHAPlusBoundRd5, "%2.6f").c_str(), SEP);
//  fprintf(myfile, "%s%c", stringValue(LPBasisCond, "%2.3f").c_str(), SEP);
//  fprintf(myfile, "%s%c", stringValue(SICBasisCondInit, "%2.3f").c_str(), SEP);
//  fprintf(myfile, "%s%c", stringValue(SICBasisCondFinal, "%2.3f").c_str(), SEP);
//  fprintf(myfile, "%s%c", stringValue(PHABasisCondRd1, "%2.3f").c_str(), SEP);
//  fprintf(myfile, "%s%c", stringValue(PHABasisCondRd2, "%2.3f").c_str(), SEP);
//  fprintf(myfile, "%s%c", stringValue(PHABasisCondRd5, "%2.3f").c_str(), SEP);
//  fprintf(myfile, "%s%c", stringValue(PHABasisCondFinal, "%2.3f").c_str(), SEP);
//  fprintf(myfile, "%s%c", stringValue(PHAPlusBasisCondRd1, "%2.3f").c_str(), SEP);
//  fprintf(myfile, "%s%c", stringValue(PHAPlusBasisCondRd2, "%2.3f").c_str(), SEP);
//  fprintf(myfile, "%s%c", stringValue(PHAPlusBasisCondRd5, "%2.3f").c_str(), SEP);
//  fprintf(myfile, "%s%c", stringValue(PHAPlusBasisCondFinal, "%2.3f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(minSICOrthogonalityInit, "%2.3f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(maxSICOrthogonalityInit, "%2.3f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(avgSICOrthogonalityInit, "%2.3f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(minPHAOrthogonalityInit, "%2.3f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(maxPHAOrthogonalityInit, "%2.3f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(avgPHAOrthogonalityInit, "%2.3f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(minSICOrthogonalityFinal, "%2.3f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(maxSICOrthogonalityFinal, "%2.3f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(avgSICOrthogonalityFinal, "%2.3f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(minPHAOrthogonalityFinal, "%2.3f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(maxPHAOrthogonalityFinal, "%2.3f").c_str(), SEP);
  fprintf(myfile, "%s%c", stringValue(avgPHAOrthogonalityFinal, "%2.3f").c_str(), SEP);
  fflush(myfile);
} /* writeRoundInfoToII */

/***********************************************************************/
/**
 * @brief Appends information about GICs (cut heuristics that worked and fails)
 */
void writeCutHeurInfoToII(const int totalNumCuts, const int totalNumObj,
    const std::vector<int>& active_cut_heur_sicsolver,
    FILE *myfile) {
  if (myfile == NULL)
    return;
  // How many cuts total
  fprintf(myfile, "%d%c", totalNumCuts, SEP);
  // Number cgs
  fprintf(myfile, "%d%c",
      GlobalVariables::param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND), SEP);
  // Cut heuristic stats
  for (int heur_ind = 0; heur_ind < CutHeuristics::NUM_CUT_HEUR; heur_ind++) {
    fprintf(myfile, "%d%c", GlobalVariables::numCutsFromHeur[heur_ind], SEP);
    fprintf(myfile, "%d%c", active_cut_heur_sicsolver[heur_ind], SEP);
  }
  // How many objectives tried
  fprintf(myfile, "%d%c", totalNumObj, SEP);
  // Cut solver fails
  for (int cut_fails_ind = 0;
      cut_fails_ind < CutSolverFails::NUM_CUTSOLVER_FAILS; cut_fails_ind++) {
    fprintf(myfile, "%d%c", GlobalVariables::numCutSolverFails[cut_fails_ind],
        SEP);
  }
  fflush(myfile);
} /* writeCutHeurInfoToII */

/***********************************************************************/
/**
 * @brief Appends information about time to solve various parts of the prob
 */
void writeTimeInfoToII(Stats &timeStats, std::vector<std::string> &times,
    const char* end_time_string, FILE *myfile) {
  if (myfile == NULL)
    return;
  for (int t = 0; t < (int) times.size(); t++)
    fprintf(myfile, "%2.3f%c", timeStats.get_time(times[t]), SEP);
  fprintf(myfile, "%s%c", end_time_string, SEP);
  fflush(myfile);
} /* writeTimeInfoToII */

/***********************************************************************/
/**
 * @brief Writes entry to II
 */
void writeEntryToII(const int x, FILE *myfile) {
  if (myfile == NULL)
    return;
  fprintf(myfile, "%d%c", x, SEP);
}

/***********************************************************************/
/**
 * @brief Writes entry to II
 */
void writeEntryToII(const double x, FILE *myfile) {
  if (myfile == NULL)
    return;
  fprintf(myfile, "%f%c", x, SEP);
}

/***********************************************************************/
/**
 * @brief Writes entry to II
 */
void writeEntryToII(const std::string x, FILE *myfile, const char SEPCHAR) {
  if (myfile == NULL)
    return;
  fprintf(myfile, "%s%c", x.c_str(), SEPCHAR);
}

/***********************************************************************/
/**
 * @brief Writes time info to separate file from instances info.
 */
void writeTimeInfo(Stats &timeStats, std::vector<std::string> &times,
    char* filename, const char* end_time_string) {
  FILE* TimeOut = fopen(filename, "w");
  if (TimeOut == NULL) {
    error_msg(errorstring,
        "Unable to open file to write time info for filename %s.\n",
        filename);
    writeErrorToII(errorstring, GlobalVariables::log_file);
    exit(1);
  }
  fprintf(TimeOut, "\n## Writing time info. ##\n");
  for (int t = 0; t < (int) times.size(); t++)
    fprintf(TimeOut, "%s%c", times[t].c_str(), SEP);
  fprintf(TimeOut, "\n");
  for (int t = 0; t < (int) times.size(); t++)
    fprintf(TimeOut, "%2.3f%c", timeStats.get_time(times[t]), SEP);
  fprintf(TimeOut, "%s", end_time_string);
  fclose(TimeOut);
}
