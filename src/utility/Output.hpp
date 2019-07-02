//============================================================================
// Name        : Output.hpp
// Author      : akazachk
// Version     : 0.2013.mm.dd
// Copyright   : Your copyright notice
// Description : 
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#pragma once

#include <string>

#include "AdvCut.hpp"
#include "GlobalConstants.hpp"

#include "CoinPackedMatrix.hpp"
#include "typedefs.hpp"
#include "SolutionInfo.hpp"
#include "Stats.hpp"


/***********************************************************************/
/**
 * @brief Adds the header lines to the file
 */
void writeHeaderToII(FILE *myfile);

/***********************************************************************/
/**
 * @brief Writes problem name to start a new instance info line
 */
void writeProbNameToII(std::string &prob_name, FILE *myfile);

/***********************************************************************/
/**
 * Writes which settings we are using to instance info
 */
void writeParamInfoToII(const CglGICParam& param, FILE* myfile);

/***********************************************************************/
/**
 * @brief Appends information about solution to prob of given type
 * @param type  :: May be either lpRelax, subprob, or afterCut
 */
void writeSolnInfoToII(const SolutionInfo& probSolnInfo,
    FILE *myfile);

/***********************************************************************/
/**
 * Writes which settings we are using to instance info,
 * as well as the number of points and rays generated
 */
void writeGICPointInfoToII(const int hplaneHeur, const int numHplanesPerRay,
    const int numRaysOfC1, const double avgNumParallelRays,
    const double avg_num_SIC_final_points, 
    const double avg_num_SIC_final_rays, 
    const int numRaysToBeCut,
    const int numCutGenSets, std::vector<int> &numPoints,
    std::vector<int> &numFinalPoints, std::vector<int> &numRays,
    std::vector<int>& numFinalRays,
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
    const std::vector<double>& maxObjDepth, FILE* myfile);

/***********************************************************************/
/**
 * @brief Appends information about GICs (cut heuristics that worked and fails)
 */
void writeCutHeurInfoToII(const int totalNumCuts, const int totalNumObj,
    const std::vector<int>& active_cut_heur_sicsolver,
    FILE *myfile);

/***********************************************************************/
/**
 * @brief Appends information about time to solve various parts of the prob
 */
void writeTimeInfoToII(Stats &timeStats, std::vector<std::string> &times,
    const char* end_time_string, FILE *myfile);

/***********************************************************************/
/**
 * @brief Writes entry to II
 */
void writeEntryToII(const int x, FILE *myfile);

/***********************************************************************/
/**
 * @brief Writes entry to II
 */
void writeEntryToII(const double x, FILE *myfile);

/***********************************************************************/
/**
 * @brief Writes entry to II
 */
void writeEntryToII(const std::string x, FILE *myfile, const char SEPCHAR = SEP);

/***********************************************************************/
/**
 * @brief Writes info about root node bound after all cuts added
 */
void writeRootBoundToII(const double LPBound, const int numSICs,
    const double sic_opt, 
    const int numPHA, const double pha_opt, 
    const double sic_pha_opt,
    const int num_active_sics_sicsolver, 
    const int num_active_pha_sicsolver,
    FILE *myfile);

/***********************************************************************/
/**
 * @brief Write bound by round and other round info
 */
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
    const double avgPHAOrthogonalityFinal, FILE *myfile);

/***********************************************************************/
/**
 * @brief Writes time info to separate file from instances info.
 */
void writeTimeInfo(Stats &timeStats, std::vector<std::string> &times,
    char* filename, const char* end_time_string);
