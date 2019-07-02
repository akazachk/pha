#pragma once

#include <string>

#include "AdvCut.hpp"
#include "GlobalConstants.hpp"

#include "CoinPackedMatrix.hpp"
#include "typedefs.hpp"
#include "SolutionInfo.hpp"
#include "Stats.hpp"

void writeCondNumToII(double LGICBasisCondNumInit, double LGICBasisCondNumAll,
    double LSICBasisCondNumInit, double LSICBasisCondNumAll,
    double OSICBasisCondNum, double LGICAndLSICBasisCondNumInit,
    double LGICAndLSICBasisCondNumAll, double LGICAndOSICBasisCondNumInit,
    double LGICAndOSICBasisCondNumAll, FILE* myfile);

/***********************************************************************/
void printSplitInfo(const OsiSolverInterface* const solverMain,
    SolutionInfo& solnInfoMain);

/***********************************************************************/
/**
 * @brief Writes euclidean distance by which initial vertex is cut,
 * and the objective function after adding the single cut
 */
void writeCutSummary(std::string headerText, AdvCuts &structCut, char* filename,
    FILE* inst_info_out);

/***********************************************************************/
/**
 * @brief Writes euclidean distance by which initial vertex is cut,
 * and the objective function after adding the single cut
 */
void writeCutSummary(std::string headerText, std::vector<AdvCuts> &structCut, char* filename,
    FILE* inst_info_out);

/***********************************************************************/
/**
 * @brief Writes sets of cuts that are stored in structural space.
 */
void writeCutsInStructSpace(std::string headerText, AdvCuts &structCut,
    SolutionInfo &solnInfo, const PHASolverInterface* const solver, char* filename,
    FILE* inst_info_out);

/***********************************************************************/
/**
 * @brief Writes sets of cuts that are stored in structural space.
 */
void writeCutsInStructSpace(std::string headerText,
    std::vector<AdvCuts> &structCut, SolutionInfo &solnInfo,
    const PHASolverInterface* const solver, char* filename, FILE* inst_info_out);

/***********************************************************************/
/**
 * @brief Writes set of cuts that are stored in NB space.
 */
void writeCutsInNBSpace(std::string headerText, AdvCuts &NBcutinfo,
    SolutionInfo &solnInfo, const PHASolverInterface* const solver, char* filename,
    FILE* inst_info_out);

/***********************************************************************/
/**
 * @brief Writes set of cuts that are stored in NB space.
 */
void writeCutsInNBSpace(std::string headerText, std::vector<AdvCuts> &NBcutinfo,
    SolutionInfo &solnInfo, const PHASolverInterface* const solver, char* filename,
    FILE* inst_info_out);

/***********************************************************************/
/**
 * @brief Writes how far along each NB ray we go for each cut.
 */
void writeRayDistInfo(std::string headerText, const PHASolverInterface* const solver,
    SolutionInfo &solnInfo, AdvCuts &structCut,
    std::vector<std::vector<double> > &distAlongRay, char* filename,
    FILE* inst_info_out);

/***********************************************************************/
/**
 * @brief Writes how far along each NB ray we go for each cut.
 */
void writeRayDistInfo(std::string headerText, const PHASolverInterface* const solver,
    SolutionInfo &solnInfo, std::vector<AdvCuts> &structCut,
    std::vector<std::vector<std::vector<double> > > &distAlongRay,
    char* filename, FILE* inst_info_out);

void writeInterPtsAndRaysInJSpace(FILE* fptr, const SolutionInfo &solnInfo,
    const std::vector<IntersectionInfo>& interPtsAndRays,
    const std::vector<Hplane>& hplaneStore);

void writeInterPtsAndRaysInJSpace(FILE* fptr, const SolutionInfo &solnInfo,
    const std::vector<std::vector<IntersectionInfo> >& interPtsAndRays,
    const std::vector<std::string>& cgsName,
    const std::vector<Hplane>& hplaneStore);

void writeInterPtsAndRaysInJSpace(FILE* fptr, const SolutionInfo &solnInfo,
    const std::vector<IntersectionInfo>& interPtsAndRays,
    const std::vector<Hplane>& hplaneStore,
    const std::vector<int>& num_SIC_points,
    const std::vector<int>& num_SIC_final_points,
    const std::vector<int>& num_SIC_rays);

/***********************************************************************/
/**
 * @brief Writes cut solver errors to separate file from instances info.
 */
void writeCutSolverErrors(int num_obj_tried, char* filename);

/***********************************************************************/
/**
 * Prints our the basis corresponding the COIN instance solver
 */
void printSimplexTableauWithNames(FILE* fptr, PHASolverInterface* const solver);

/**
 * Prints our the basis corresponding the COIN instance solver; only non-basic non-artificial columns printed.
 */
void printNBSimplexTableauWithNames(FILE* fptr, PHASolverInterface* const solver);

/***********************************************************************/
/**
 * @brief Writes solution (both primal and dual) to file.
 * @param solver          ::  The solver whose solution we are writing.
 * @param soln_name       ::  The name we are outputting to.
 * @param out_stream      ::  Out stream.
 * @param width_varname   ::  Width of variable name field.
 * @param width_val       ::  Width of variable value field.
 * @param width_stat      ::  Width of cstat/rstat field.
 * @param width_rc        ::  Width of reduced cost field.
 * @param width_activity  ::  Width of row activity field.
 * @param width_rhs       ::  Width of rhs field.
 * @param after_dec       ::  Number of digits after the decimal.
 * @param EXT             ::  Extension added to end of file output.
 */
void writeSoln(const PHASolverInterface* const solver, std::string soln_name,
    int width_varname = 16, int width_val = 16, int width_stat = 16,
    int width_rc = 16, int width_activity = 16, int width_rhs = 16,
    int after_dec = 6, const std::string EXT = "-optSoln.alex");

double rayPostPivot(const int component, const int nb_var_in1,
    const int struct_var_out1, const int nb_var_in2,
   // const PHASolverInterface* const solver, 
    const SolutionInfo& solnInfo);

double vertexPostSinglePivot(const int component, const double rayDirnToHplane1,
    const double init_val_out1, const double final_val_out1);

double pointPostDoublePivot(const int component, const double rayDirnToHplane1,
    const double init_val_out1, const double final_val_out1,
    const double rayDirnToSplit1, const double rayDirnToSplit2,
    const double rayDirnToHplane2, const double init_val_out2,
    const double final_val_out2);

double pointPostPivot(const int component, const int nb_var_in1,
    const int struct_var_out1, const int nb_var_in2,
    const int struct_var_out2, const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const bool computeVertexFlag = false);

void packedShallowVectorSum(CoinPackedVector& sum, const double mult1,
    const CoinShallowPackedVector& vec1, const double mult2,
    const CoinShallowPackedVector& vec2);

double packedSum(const CoinPackedVector& vec1, const CoinPackedVector& vec2);

double packedDotProduct(const CoinPackedVector& vec1,
    const CoinPackedVector& vec2);

/********************************************************************************/
/**
 * Check if two solvers have the same everything:
 * 1. colSolution
 * 2. rowPrice
 * 3. slackVals
 * 4. cstat and rstat
 * 5. reducedCosts
 */
bool solversEqual(const PHASolverInterface* const clp1, const PHASolverInterface* const clp2);
