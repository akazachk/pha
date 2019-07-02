//============================================================================
// Name        : CutHelper.hpp
// Author      : akazachk
// Version     : 0.2014.03.20
// Copyright   : Your copyright notice
// Description : Helper functions of cuts
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#pragma once

#include "Point.hpp"
#include "Vertex.hpp"
#include "Ray.hpp"
#include "AdvCut.hpp"
#include "SolutionInfo.hpp"

#include "typedefs.hpp"

#include "IntersectionInfo.hpp"

/***********************************************************************/
/***********************************************************************/
/** Methods related to generating cuts from non-basic space */
/***********************************************************************/
/***********************************************************************/

/**
 * Generates cuts based on points and rays, in the *non-basic* space
 * @return Number of cuts generated
 */
int genCutsFromPointsRaysNB(PHASolverInterface* const cutSolver,
    const double beta, AdvCuts & cuts, int& num_cuts_this_cgs,
    int& num_cuts_total, int& num_obj_tried,
    const PHASolverInterface* const solver, const SolutionInfo & solnInfo,
//    const IntersectionInfo& interPtsAndRays,
    const int cgsIndex, const std::string& cgsName, const AdvCuts& structSICs);

/**
 * Set the cut solver objective based on option, in the *non-basic* space
 * Options are:
 *  1. max_eff
 */
void setCutSolverObjectiveFromVertexNB(const int vert_ind,
    PHASolverInterface* const cutSolver,
    const std::vector<Vertex>& vertexStore);

/**
 * Set up solver from points and rays provided
 */
bool setupCutSolver(PHASolverInterface* const cutSolver,
    const IntersectionInfo& interPtsAndRays,
    const int cgsIndex, const std::string& cgsName);

/***********************************************************************/
/***********************************************************************/
/** Cut heuristics **/
/***********************************************************************/
/***********************************************************************/

/**
 * @brief Try being tight on each of the intersection points/rays
 */
int tightOnPointsRaysCutGeneration(PHASolverInterface* const cutSolver,
    const double beta, AdvCuts& cuts, int& num_cuts_this_cgs,
    int& num_cuts_total, int& num_obj_tried,
    const PHASolverInterface* const solver, const SolutionInfo & probData,
    //const IntersectionInfo& interPtsAndRays,
    const int cgsIndex, const std::string& cgsName, const AdvCuts& structSICs,
    const bool inNBSpace = true);

/**
 * Generate a cut, if possible, from the solver
 * Rhs is given as beta
 */
int genCutsFromCutSolver(AdvCuts & pcuts,
    PHASolverInterface* const cutSolver,
    const PHASolverInterface* const origSolver,
    const SolutionInfo & origProbData, const double beta, const int cgsIndex,
    const std::string& cgsName, const AdvCuts& structSICs,
    int& num_cuts_this_cgs, int& num_cuts_total, const bool inNBSpace,
    const CutHeuristics cutHeur, const bool tryExtraHard = false);
int genCutsHelper(AdvCuts & pcuts, PHASolverInterface* const cutSolver,
    const PHASolverInterface* const origSolver,
    const SolutionInfo & origProbData, const double beta, const int cgsIndex,
    const std::string& cgsName, const AdvCuts& structSICs,
    int& num_cuts_this_cgs, int& num_cuts_total, const bool inNBSpace,
    const CutHeuristics cutHeur);
int exitGenCutsFromCutSolver(PHASolverInterface* const cutSolver,
    const int return_status);

/**
 * @brief Try cutting off each of the vertices generated in the intermediate stages
 */
int cutVerticesCutGeneration(PHASolverInterface* const cutSolver,
    const double beta, AdvCuts & cuts, int& num_cuts_this_cgs,
    int& num_cuts_total, int& num_obj_tried,
    const PHASolverInterface* const solver, const SolutionInfo & probData,
    const std::vector<Vertex>& vertexStore, const int cgsIndex,
    const std::string& cgsName, const AdvCuts& structSICs);

/**
 * For each ray, look at set of intersection points (across all splits) generated by this ray
 * Call this set P(r), corresponding to ray r
 * Let p(r,s) be the point generated by ray r for split s
 * Let d(r,s) = distance along ray r to split s
 * Then for any split s, we can try to cut points p(r,s') for all s' : d(r,s') < d(r,s)
 */
void splitSharingCutGeneration(std::vector<PHASolverInterface*>& cutSolver,
    AdvCuts & cuts, std::vector<int>& num_cuts_per_split, int& num_cuts_total,
    int& num_obj_tried, const std::vector<IntersectionInfo>& interPtsAndRays,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
    const std::vector<Vertex>& vertexStore,
    const std::vector<Hplane>& hplaneStore, const AdvCuts& structSICs,
    const std::vector<bool>& isCutSolverPrimalFeasible);

/**
 * @brief For the ray that generated the given point, find the deepest split the ray intersects
 */
int findDeepestSplit(const IntersectionInfo& interPtsAndRays, const int pt_ind,
    const int split_ind, const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo,
//    const std::vector<Ray>& rayStore,
    const std::vector<Vertex>& vertexStore,
    const std::vector<Hplane>& hplaneStore);

/**
 * We want to solve the bilinear program
 * min <\alpha, x>
 * s.t.
 * <\alpha, p> \ge 1 \forall p \in \pointset
 * <\alpha, r> \ge 0 \forall r \in \rayset
 * x \in Pbar (P + all SICs)
 */
int iterateDeepestCutPostMSIC(AdvCuts & cuts, int& num_cuts_this_cgs,
    int& num_cuts_total, int& num_obj_tried,
    PHASolverInterface* const cutSolver,
    const PHASolverInterface* const origSolver,
    const SolutionInfo & origSolnInfo, const double beta, const int cgsIndex,
    const std::string& cgsName, const AdvCuts& structSICs,
    const bool inNBSpace);

int exitFromIterateDeepestCutPostMSIC(const int num_cuts_gen,
    PHASolverInterface* const cutSolver,
    PHASolverInterface* const MSICSolver);

/**
 * Scales ray generated.
 * If max is not too big, then we do not bother rescaling
 * If min / norm too small, then we do not do that scaling
 *
 * @return status: -1 is bad ray, don't add it; 0 is no scaling; 1 is okay scaling
 */
int scalePrimalRay(int numCols, double** rayToScale, const double SENSIBLE_MAX =
    100.0);

/**
 * Compute score as in SCIP from initialized orthogonalities
 */
void computeScore(AdvCuts& cs, const int pos);

/**
 * As in SCIP
 */
void updateOrthogonalities(AdvCuts& cs, const OsiRowCut& cut,
    const double min_cut_orthogononality,
    const std::vector<int>& roundWhenCutApplied);

/**
 * As in SCIP, returns the position of the best non-forced cut in the cuts array
 * Currently not forcing any cuts
 * */
int getBestCutIndex(const AdvCuts& cs,
    const std::vector<int>& roundWhenCutApplied);

/**
 * Generic way to apply a set of cuts to a solver
 * @return Solver value after adding the cuts, if they were added successfully
 */
double applyCuts(PHASolverInterface* solver, const OsiCuts& cs,
    const int startCutIndex, const int numCutsOverload);
double applyCuts(PHASolverInterface* solver, const OsiCuts& cs,
    const int numCutsOverload = -1);
double applyCutsInRounds(PHASolverInterface* solver, AdvCuts& cs,
    std::vector<int>& numCutsByRound, std::vector<double>& boundByRound,
    std::vector<double>& basisCond, const int numCutsToAddPerRound,
    std::vector<int>* const oldRoundWhenCutApplied =
    NULL, bool populateRoundWhenCutApplied = false, const int maxRounds = 5, 
    double compare_obj = std::numeric_limits<double>::lowest());

/**
 * Get bound from applying
 * 1. SICs
 * 2. PHAs
 * 3. SIC+PHA
 */
void compareCuts(const int num_obj_tried,
    const OsiSolverInterface* const solver, AdvCuts& structSICs,
    AdvCuts& structPHAs, const int numCutsToAddPerRound, std::string& output);