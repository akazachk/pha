// Name:     CglPHA.cpp
// Author:   A.M. Kazachkov
//           Tepper School of Business
//-----------------------------------------------------------------------------

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cstdlib>
#include <cstdio>
#include <cmath> // log2
#include <cfloat>
#include <cassert>
#include <iostream>
#include <fenv.h>
#include <climits>

#include "OsiSolverInterface.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CoinFactorization.hpp"
#include "CglPHA.hpp"
#include "CoinFinite.hpp"

#include "GlobalConstants.hpp"

#include "Utility.hpp"
#include "CutHelper.hpp"

// Useful classes and such for GICs
#include "Output.hpp"
#include "CglSIC.hpp"

#ifdef TRACE
#include "Debug.hpp"
#endif

/***************************************************************************/
CglPHA::CglPHA() :
    CglCutGenerator(), solver(NULL), castsolver(NULL), xlp(NULL), byRow(NULL), byCol(NULL), nb_advcs(true), struct_advcs(
        false), structSICs(false), NBSICs(true), num_obj_tried(0) {
  copyOurStuff();
  param.setParams();
} /* empty constructor */

/***************************************************************************/
CglPHA::CglPHA(const CglGICParam& parameters) :
    CglCutGenerator(), param(parameters), solver(NULL), castsolver(NULL), xlp(
    NULL), byRow(NULL), byCol(NULL), nb_advcs(true), struct_advcs(false), structSICs(
        false), NBSICs(true), num_obj_tried(0) {
  copyOurStuff();
} /* constructor copying params */

/***************************************************************************/
CglPHA::CglPHA(const CglPHA& rhs) :
    CglCutGenerator(rhs) {
  copyOurStuff(&rhs);
} /* copy constructor */

/***************************************************************************/
CglPHA & CglPHA::operator=(const CglPHA& rhs) {
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    copyOurStuff(&rhs);
  }
  return *this;
} /* assignment operator */

/***************************************************************************/
CglPHA::~CglPHA() {} /* destructor */

/***************************************************************************/
void CglPHA::copyOurStuff(const CglPHA* const rhs) {
  if (rhs) {
    param = rhs->param;
    solver = rhs->solver;
    castsolver = rhs->castsolver;
    interPtsAndRays = rhs->interPtsAndRays;
    xlp = rhs->xlp;
    byRow = rhs->byRow;
    byCol = rhs->byCol;
    nb_advcs = rhs->nb_advcs;
    struct_advcs = rhs->struct_advcs;
    structSICs = rhs->structSICs;
    NBSICs = rhs->NBSICs;
    min_nb_obj_val = rhs->min_nb_obj_val;
    num_obj_tried = rhs->num_obj_tried;
    num_total_gics = rhs->num_total_gics;
  } else {
    this->interPtsAndRays.resize(0);
    this->min_nb_obj_val = std::numeric_limits<double>::max(); //std::numeric_limits<double>::lowest();
    this->num_obj_tried = 0;
    this->num_total_gics = 0;
  }
} /* copyOurStuff */

/*********************************************************************/
CglCutGenerator *
CglPHA::clone() const {
  return new CglPHA(*this);
} /* clone */

/************************************************************************/
void CglPHA::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
    const CglTreeInfo info) const {
  // kludge to be able to modify the CglPHA object if it is const
  CglPHA temp(*this);
  temp.generateCuts(si, cs, info);
} /* generateCuts (required by CglCutGenerator) */

/************************************************************************/
void CglPHA::generateCuts(const OsiSolverInterface &si, OsiCuts & cs,
    const CglTreeInfo info) {
  // Create SolutionInfo
  solver = const_cast<OsiSolverInterface *>(&si); // Removes const-ness
  SolutionInfo tmpSolnInfo(solver,
      param.paramVal[ParamIndices::MAX_FRAC_VAR_PARAM_IND], false);
  if (tmpSolnInfo.numSplits == 0) {
    return;
  }
  if (!param.finalizeParams(tmpSolnInfo.numNB, tmpSolnInfo.numSplits)) {
    error_msg(errstr,
        "Issue finalizing parameters (no rays to cut or issue with parameter combinations).\n");
    writeErrorToII(errstr, GlobalVariables::log_file);
    exit(1);
  }
  if (param.paramVal[MAX_FRAC_VAR_PARAM_IND] > tmpSolnInfo.numFeasSplits) {
    warning_msg(warnstring,
        "Number of requested cut generating sets is %d, but the maximum number we can choose is %d. Reducing the requested number of splits.\n",
        param.paramVal[MAX_FRAC_VAR_PARAM_IND], tmpSolnInfo.numFeasSplits);
    param.paramVal[MAX_FRAC_VAR_PARAM_IND] = tmpSolnInfo.numFeasSplits;
  }

  if (param.paramVal[ParamIndices::SICS_PARAM_IND] && (NBSICs.size() == 0)) {
    generateSICs(NBSICs, structSICs, tmpSolnInfo);
  }
//  GlobalVariables::param = this->param;
  generateCuts(si, cs, tmpSolnInfo, structSICs, NBSICs);
} /* generateCuts */

/************************************************************************/
void CglPHA::generateCuts(const OsiSolverInterface &si, OsiCuts & cs,
    const SolutionInfo& probData,
    const AdvCuts& structMSICs, const AdvCuts& NBMSICs) {
  solver = const_cast<OsiSolverInterface *>(&si); // Removes const-ness
  if (solver == NULL) {
    printf("*** WARNING: CglPHA::generateCuts(): no solver available.\n");
    return;
  }

  if (!solver->optimalBasisIsAvailable()) {
    printf(
        "*** WARNING: CglPHA::generateCuts(): no optimal basis available.\n");
    return;
  }

  xlp = solver->getColSolution();
  try {
    startPHAGeneration(dynamic_cast<AdvCuts&>(cs), probData, structMSICs, NBMSICs);
    this->struct_advcs = dynamic_cast<AdvCuts&>(cs);
  } catch (std::exception& e) {
    warning_msg(warnstr,
        "Not able to cast OsiCuts as AdvCuts, used to keep extra information.\n");

    startPHAGeneration(this->struct_advcs, probData, structMSICs, NBMSICs);
    cs = struct_advcs;
  }
} /* generateCuts */

/************************************************************************/
/**
 * This is the main PHA routine. It sets things up and sends calls to Algorithms 2 and 3.
 */
void CglPHA::startPHAGeneration(AdvCuts &structPHA, const SolutionInfo& solnInfo,
    const AdvCuts &structSICs, const AdvCuts &NBSICs) {
  castsolver = dynamic_cast<OsiClpSolverInterface*>(this->solver->clone());
  castsolver->enableFactorization(); // Need factorization for chooseHplanes and splitShare

  /***********************************************************************************
   * Save objective
   ***********************************************************************************/
  // In order to calculate how good the points we find are
  // with respect to the objective function (which should be in NB space when appropriate)
  // Clearly a valid inequality is c^T x \ge c^\T v,
  // where v is the LP optimum, since we are minimizing
  AdvCut NBObj(true, 1), structObj(false, 1);
  structObj.setOsiRowCut(solver->getNumCols(), solver->getObjCoefficients(),
      solver->getObjValue(), CoinMin(solnInfo.EPS, param.getEPS()));

  // Either we have the objective, or we have the reduced costs
  // of the nonbasic variables (which is exactly the objective
  // in the nonbasic space).
  // Recall that the reduced cost of slack variables is exactly
  // the negative of the row price (dual variable) for that row.
  NBObj.NBSpace = true;
  NBObj.num_coeff = solnInfo.numNB;
  NBObj.setOsiRowCut(solnInfo.nonBasicReducedCost, 0.0,
      CoinMin(solnInfo.EPS, param.getEPS()));
  const double objNorm = NBObj.getTwoNorm();

  /***********************************************************************************
   * Set up for generating GICs
   ***********************************************************************************/
  printf(
      "\n## Instance info:\n== Splits chosen on the following variables:\n");
  for (int s = 0; s < solnInfo.numFeasSplits; s++) {
    const int curr_var = solnInfo.feasSplitVar[s];
    printf(
        ":::: Var name: %s. Var value: %2.3f. Var index: %d. Row index for var: %d.\n",
        solver->getColName(curr_var).c_str(),
        solver->getColSolution()[curr_var],
        solnInfo.feasSplitVar[s],
        solnInfo.rowOfVar[curr_var]);
  }
  printf("====\n");

  printf("== Number rays to be cut set to %d. ==\n",
      param.paramVal[NUM_RAYS_CUT_PARAM_IND]);
  printf("== Number effective cut-generating sets is %d. ==\n",
      solnInfo.numFeasSplits);

  // For tracking time limit
  Stats phaTimeStats;
  phaTimeStats.register_name("PHA_TOTAL");
  phaTimeStats.start_timer("PHA_TOTAL");

  /***********************************************************************************
   * Generate initial points and rays (Alg 2 & 3 steps 2-3)
   ***********************************************************************************/
  std::vector<Hplane> hplaneStore; // activated hyperplanes
  std::vector<Vertex> vertexStore; // vertices generated along the way
  std::vector<Ray> rayStore; // rays generated along the way

  // Set up Vertex and Ray storage
  vertexStore.resize(1); // Initially just the origin
  vertexStore[0].rayIndices.resize(solnInfo.numNB);
  vertexStore[0].rayActiveFlag.resize(solnInfo.numNB, true);
  rayStore.resize(solnInfo.numNB); // Initially the rays of C1
  for (int i = 0; i < solnInfo.numNB; i++) {
    const int rayIndices[1] = { i };
    const double rayValues[1] = { 1.0 };
    rayStore[i].setVector(1, rayIndices, rayValues);
    rayStore[i].cutFlag.resize(solnInfo.numFeasSplits, false);

    // This ray emanates from the original vertex
    vertexStore[0].rayIndices[i] = i;
  }

  // The initial constraint matrices for each split and initial points/rays
  std::vector<int> num_SIC_points(solnInfo.numFeasSplits, 0);
  std::vector<int> num_SIC_final_points(solnInfo.numFeasSplits, 0);
  std::vector<int> num_SIC_rays(solnInfo.numFeasSplits, 0);
  std::vector<int> num_SIC_final_rays(solnInfo.numFeasSplits, 0);
  std::vector<std::vector<double> > distToSplitAlongRay(solnInfo.numFeasSplits);
//  std::vector<std::vector<bool> > SICPointIsFinal(solnInfo.numFeasSplits);
  this->interPtsAndRays.resize(solnInfo.numFeasSplits);
  for (int split_ind = 0; split_ind < solnInfo.numFeasSplits; split_ind++) {
    // Reverse ordering so that the packed matrix is row ordered
    // and set number of columns
    interPtsAndRays[split_ind].reverseOrdering();
    interPtsAndRays[split_ind].setDimensions(0, solnInfo.numNB);

//    SICPointIsFinal[split_ind].resize(solnInfo.numNB);
    const bool recalc_dist = NBSICs.sizeCuts() != solnInfo.numFeasSplits;
    if (recalc_dist) {
      distToSplitAlongRay[split_ind].resize(solnInfo.numNB);
    } else {
      distToSplitAlongRay[split_ind] = NBSICs.cuts[split_ind].distAlongRay;
    }
    // This adds all the original intersection points
    for (int ray_ind = 0; ray_ind < solnInfo.numNB; ray_ind++) {
      if (!addInterPtFromUncutRayOfC1(interPtsAndRays[split_ind],
          distToSplitAlongRay[split_ind][ray_ind], solnInfo, split_ind,
          rayStore, ray_ind, vertexStore, recalc_dist)) {
        continue;
      }

      // Gather information about number of final initial points
      GlobalVariables::timeStats.start_timer(COLLECT_STATS_TIME);

      // If it's a point get it, o/w count the ray
      const int row_ind = interPtsAndRays[split_ind].getNumRows() - 1;
      if (isZero(interPtsAndRays[split_ind].RHS[row_ind])) {
        num_SIC_rays[split_ind]++;
        if (interPtsAndRays[split_ind].finalFlag[row_ind]) {
          num_SIC_final_rays[split_ind]++;
        }
      } else {
        num_SIC_points[split_ind]++;
        // TODO Check if the ray value is equal to the distance to split; should be the same
        if (interPtsAndRays[split_ind].finalFlag[row_ind]) {
//          SICPointIsFinal[split_ind][row_ind] = true;
          num_SIC_final_points[split_ind]++;
        }
      }
      GlobalVariables::timeStats.end_timer(COLLECT_STATS_TIME);
    }
  } /* iterate over splits, collecting initial points and rays */

  /***********************************************************************************
   * Alg 3, step 4: sort rays by decreasing reduced cost
   ***********************************************************************************/
  //Choose which of the initial rays to be cut
  //(based on weakest average cost coefficients among Gomory cuts over the chosen splits)
  const int num_rays_to_be_cut = param.getParamVal(ParamIndices::NUM_RAYS_CUT_PARAM_IND);
  std::vector<double> avgCoeffAcrossCuts(solnInfo.numNB, 0.0);
  std::vector<int> sortedRaysIndex(solnInfo.numNB);
  for (int j = 0; j < solnInfo.numNB; j++) {
    sortedRaysIndex[j] = j;
    for (int s = 0; s < (int) NBSICs.size(); s++) {
      avgCoeffAcrossCuts[j] += NBSICs[s][j];
    }
    avgCoeffAcrossCuts[j] /= (double) solnInfo.numFeasSplits;
  }

  //Sort average cut coefficients in descending order
  //Larger cut coefficients in the nonbasic space correspond to weaker cuts
  sort(sortedRaysIndex.begin(), sortedRaysIndex.end(),
      index_cmp<double*>(&avgCoeffAcrossCuts[0]));

  // Resize arrays and put default values
//  std::vector<int> allRayIndices(solnInfo.numNB);
//  for (int j = 0; j < solnInfo.numNB; j++) {
//    allRayIndices[j] = j;
//  }

  // Number of first activations (if any) is ceil(log(num_rays_to_be_cut + 1)),
  // as in each round we activate double the number of ray as in the prior round
  const int num_alg3_rounds =
      param.getParamVal(ParamIndices::NUM_ALG2_ROUNDS_PARAM_IND) < 0 ?
          0 : std::ceil(std::log2(num_rays_to_be_cut + 1));
  const int num_alg2_rounds = std::abs(param.getParamVal(ParamIndices::NUM_ALG2_ROUNDS_PARAM_IND));
  const int num_total_act_rounds = num_alg3_rounds + num_alg2_rounds;

  for (int split_ind = 0; split_ind < solnInfo.numFeasSplits; split_ind++) {
    interPtsAndRays[split_ind].resizeInfoByActivationRound(num_total_act_rounds);
  }

  // Across all splits, all rays that have been cut
  // [act][ray]
  std::vector<int> allRaysCutAcrossSplits;

  // For each split, which rays are cut by at least one hplane
  // [act][split][rays]
  std::vector<std::vector<int> > raysCutForSplit(solnInfo.numFeasSplits);

  std::vector<int> num_gics_per_split(solnInfo.numFeasSplits, 0);

  if (num_alg3_rounds > 0) {
    /***********************************************************************************
     * Alg 3
     ***********************************************************************************/
    algorithm3(structPHA, solnInfo, structSICs, NBSICs,
        phaTimeStats, num_gics_per_split, NBObj, objNorm, distToSplitAlongRay,
        hplaneStore, vertexStore, rayStore, sortedRaysIndex,
//        allRayIndices,
        allRaysCutAcrossSplits, raysCutForSplit);
  }

  /***********************************************************************************
   * Alg 2
   ***********************************************************************************/
  algorithm2(structPHA, solnInfo, structSICs, NBSICs,
      phaTimeStats, num_gics_per_split, NBObj, objNorm, distToSplitAlongRay,
      hplaneStore, vertexStore, rayStore, sortedRaysIndex,
//      allRayIndices,
      allRaysCutAcrossSplits, raysCutForSplit);

  /***********************************************************************************
   * Finished cut generation, now do processing
   ***********************************************************************************/
  structPHA.setOsiCuts();

#ifdef SHOULD_WRITE_GICS
  structPHA.printCuts("PHA", solver, solnInfo);
#endif

  // Count up number of cut limit failures
  for (int s = 0; s < solnInfo.numFeasSplits; s++) {
    if (reachedCutLimit(num_gics_per_split[s], num_total_gics)) {
      GlobalVariables::numCutSolverFails[CutSolverFails::CUT_LIMIT_FAIL_IND]++;
    }
  }

  // Analyze resulting point and ray collection
  runAnalyses(&solnInfo, hplaneStore, num_SIC_points, num_SIC_final_points,
      num_SIC_final_rays, num_SIC_rays, structSICs, NBSICs, structPHA);

  // Clear things
  if (castsolver) {
    delete castsolver;
  }
  phaTimeStats.end_timer("PHA_TOTAL");

  printf("\n## Finished round of PHA generation with %d cuts. ##\n", num_total_gics);
} /* startPHAGeneration */

/**
 * ALGORITHM 1
 * (SINGLE) HYPERPLANE ACTIVATION
 * Description: Modifies the given point-ray collection by activating a single hyperplane
 *
 * Input:
 * polyhedron P defined by hyperplane set H;
 * P_I-free convex set S;
 * hyperplane H_h valid for P;
 * set of rays R_A \subseteq \bar{R};
 * proper point-ray collection (\pointset,\rayset) contained in \bar{C}.
 */
void CglPHA::algorithm1(int& totalNumPointsRemoved, double& avgRemovedDepth,
    int& totalNumPointsAdded, double& avgAddedDepth,
    IntersectionInfo& interPtsAndRays, const int split_ind,
    const int hplane_ind, const int hplane_ind_in_split,
    const std::vector<Hplane>& hplaneStore, const Hplane& hplane,
    std::vector<Ray>& rayStore, const std::vector<Vertex>& vertexStore,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
    const AdvCut* NBSIC, const double NBSICnorm, const AdvCut& objCut,
    const double objNorm, const int act_ind) {
  const int split_var = solnInfo.feasSplitVar[split_ind];

  // We we want to track the average depth of the removed intersection points
  double sumRemovedDepth = 0.0, sumAddedDepth = 0.0;
  int numNewPointsRemoved = 0, numNewPointsAdded = 0;

  std::vector<int> newRayColIndices;

  // Mark rays that will be cut by a hyperplane
//  // In PHA 1.1., all the rays for which the first hplane is this one
  std::vector<bool> raysCutFlag(solnInfo.numNB, false);
  const int start_ind = getStartIndex(hplane, split_ind, act_ind);
  const int end_ind = getEndIndex(hplane, split_ind, act_ind);
  for (int j = start_ind; j < end_ind; j++) {
    const int ray_ind = hplane.rayToBeCutByHplane[split_ind][j];
    raysCutFlag[ray_ind] = true;
    rayStore[ray_ind].cutFlag[split_ind] = true;
//    // Removed depth is zero for these points, since they are from C1
//    // However, don't count any rays that have been previously cut
//    if (!rayStore[ray_ind].cutFlag[split_ind]) {
//      rayStore[ray_ind].cutFlag[split_ind] = true;
//      // Only count it if it has not been added;
//      // otherwise the point should be removed in removePtsAndRaysCutByHplane
//      if (!rayStore[ray_ind].addedPointFlag[split_ind])
//        numNewPointsRemoved++;
//    }
  }

  /***********************************************************************************
   * Alg 1, step 3
   ***********************************************************************************/
  //Remove all intersection points cut off by this hyperplane
  numNewPointsRemoved += removePtsAndRaysCutByHplane(solver, solnInfo,
      raysCutFlag, split_var, hplane, hplane_ind_in_split, vertexStore,
      rayStore, interPtsAndRays, NBSIC, NBSICnorm, sumRemovedDepth);

  // Calculate average depth of the removed points
  if (numNewPointsRemoved > 0) {
    avgRemovedDepth += sumRemovedDepth / numNewPointsRemoved;
    totalNumPointsRemoved += numNewPointsRemoved;
  }

  /***********************************************************************************
   * Alg 1, step 4 (find indices of new edges emanating from activated hyperplane)
   ***********************************************************************************/
  //For each hyperplane, create the set J \setminus K
  storeNewRayNBIndices(solnInfo.numNB, newRayColIndices,
      hplane.allRaysIntByHplane[split_ind], raysCutFlag);
  const int numNewRays = newRayColIndices.size();

  //For ray intersected by hyperplane H generate new rays corresponding to each element in J \setminus K
  //This is of course available in closed form from the tableau
  for (int j = start_ind; j < end_ind; j++) {
    const int ray_to_cut = hplane.rayToBeCutByHplane[split_ind][j];
    const int vertex_ind = hplane.vertCreatedByRayToBeCutByHplane[split_ind][j];

    for (int r = 0; r < numNewRays; r++) {
      Ray newRay;
      Point tmpPoint;
      calcNewRayCoordinatesRank1(newRay, solnInfo, ray_to_cut,
          newRayColIndices[r], hplane, CoinMin(param.getRAYEPS(), solnInfo.EPS));
      bool ptFlagToAdd = calcIntersectionPointWithSplit(tmpPoint, solver,
          solnInfo, vertexStore[vertex_ind], newRay, split_var);

      if (ptFlagToAdd) { // if the new ray intersects the split bd, it creates a new point
        if (!duplicatePoint(interPtsAndRays, interPtsAndRays.RHS, tmpPoint)) {
          // Check that the point lies on bd S, as it should
          const int split_row_ind = solnInfo.rowOfVar[split_var];
          const double split_val = solnInfo.a0[split_row_ind];
          double newValOnSplit = dotProductWithOptTableau(tmpPoint,
              split_row_ind, solnInfo);
          double newValOnHplane = 0.0;

          bool notOnSplit = (!isVal(newValOnSplit, std::floor(split_val),
              param.getEPS()))
              && (!isVal(newValOnSplit, std::ceil(split_val), param.getEPS()));

          // Check that the point is on the hyperplane it is supposed to lie on
          // If hplane.row < 0, meaning the hplane activate was the "other" bound
          // on a nb variable, then there is no point checking anything;
          // the only nb var that intersects this hplane is the one it bounds,
          // and we force the variable to have the appropriate value (ub - lb)
          bool notOnHplane = false;
          double colLB = 0.0, colUB = solver->getInfinity();
          if (hplane.row >= 0) {
            newValOnHplane = dotProductWithOptTableau(tmpPoint, hplane.row,
                solnInfo);

            // The hplane is either structural or slack
            if (hplane.var < solver->getNumCols()) {
              colLB = getVarLB(solver, hplane.var);
              colUB = getVarUB(solver, hplane.var);
            } else {
              colLB = 0.0;
              colUB = solver->getInfinity();
            }
            notOnHplane = !isVal(newValOnHplane, colLB, param.getEPS())
                && !isVal(newValOnHplane, colUB, param.getEPS());
          }

          if (notOnSplit || notOnHplane) {
            // Sometimes worth redoing things with a higher eps tolerance
            const double prevRAYEPS = param.getRAYEPS();
            param.setRAYEPS(prevRAYEPS * prevRAYEPS);
            calcNewRayCoordinatesRank1(newRay, solnInfo, ray_to_cut,
                newRayColIndices[r], hplane, param.getRAYEPS());
            calcIntersectionPointWithSplit(tmpPoint, solver, solnInfo,
                vertexStore[vertex_ind], newRay, split_var);
            newValOnSplit = dotProductWithOptTableau(tmpPoint, split_row_ind,
                solnInfo);
            notOnSplit =
                (!isVal(newValOnSplit, std::floor(split_val), param.getEPS()))
                    && (!isVal(newValOnSplit, std::ceil(split_val), param.getEPS()));
            if (hplane.row >= 0) {
              newValOnHplane = dotProductWithOptTableau(tmpPoint, hplane.row,
                  solnInfo);
              notOnHplane = (!isVal(newValOnHplane, colLB, param.getEPS()))
                  && (!isVal(newValOnHplane, colUB, param.getEPS()));
            }
            param.setRAYEPS(prevRAYEPS);
          }

          if (notOnSplit) {
            const bool reallyNotOnSplit = 
                (!isVal(newValOnSplit, std::floor(split_val), param.getAWAY()))
                    && (!isVal(newValOnSplit, std::ceil(split_val), param.getAWAY()));
            if (reallyNotOnSplit) {
              error_msg(errstr,
                  "Point %s not on bd of split %d (var %d), created by activating var %d (hplane %d). newVal: %s, splitVal: %s.\n",
                  tmpPoint.vectorString().c_str(), split_ind, split_var,
                  hplane.var, hplane_ind, stringValue(newValOnSplit, "%e").c_str(),
                stringValue(split_val).c_str());
              writeErrorToII(errstr, GlobalVariables::log_file);
              exit(1);
            } else {
              warning_msg(errstr,
                  "Point %s not on bd of split %d (var %d), created by activating var %d (hplane %d). newVal: %s, splitVal: %s. Small enough violation that we only send a warning.\n",
                  tmpPoint.vectorString().c_str(), split_ind, split_var,
                  hplane.var, hplane_ind, stringValue(newValOnSplit, "%e").c_str(),
                stringValue(split_val).c_str());
            }
          }
          if (notOnHplane) {
            const bool reallyNotOnHplane = (!isVal(newValOnHplane, colLB, param.getAWAY()))
              && (!isVal(newValOnHplane, colUB, param.getAWAY()));
            if (reallyNotOnHplane) {
              error_msg(errstr,
                  "Basic var %d corresponding to hplane %d activated on split %d yielding pt %s is not on the activated hyperplane. newVal: %s, colLB: %s, colUB: %s.\n",
                  hplane.var, hplane_ind, split_ind,
                  tmpPoint.vectorString().c_str(),
                  stringValue(newValOnHplane, "%e").c_str(), stringValue(colLB).c_str(),
                  stringValue(colUB).c_str());
              writeErrorToII(errstr, GlobalVariables::log_file);
              exit(1);
            } else {
              warning_msg(errstr,
                  "Basic var %d corresponding to hplane %d activated on split %d yielding pt %s is not on the activated hyperplane. newVal: %s, colLB: %s, colUB: %s. Small enough violation that we only send a warning.\n",
                  hplane.var, hplane_ind, split_ind,
                  tmpPoint.vectorString().c_str(),
                  stringValue(newValOnHplane, "%e").c_str(), stringValue(colLB).c_str(),
                  stringValue(colUB).c_str());
            }
          }

          // Ensure that no ``future'' hyperplane cuts the point
          bool isCut = false;
          for (int h_ind = hplane_ind_in_split;
              h_ind < (int) interPtsAndRays.allHplaneToAct.size(); h_ind++) {
            const int curr_hplane_ind = interPtsAndRays.allHplaneToAct[h_ind];
            // The hyperplane needs to have cut the current ray
            const bool hplaneCutsThisRay =
                (find_val(ray_to_cut,
                    hplaneStore[curr_hplane_ind].rayToBeCutByHplane[split_ind])
                    >= 0);
            if (hplaneCutsThisRay
                && pointIsCutByHplane(solver, solnInfo, tmpPoint, 1.0,
                    hplaneStore[curr_hplane_ind])) {
              isCut = true;
              break;
            }
          }

          if (!isCut) {
            numNewPointsAdded++;
            const bool pointIsFinal = pointInP(solver, solnInfo,
                tmpPoint.getVectorNumElements(), tmpPoint.getVectorIndices(),
                tmpPoint.getVectorElements());

            // Get depth wrt SIC
            const double SICactivity = getSICActivity(split_var, vertexStore,
                rayStore, solver, solnInfo, tmpPoint);
            double curr_SIC_depth =
                (NBSIC) ? (SICactivity - NBSIC->rhs()) / NBSICnorm : 0.0;
            if (lessThanVal(curr_SIC_depth, 0.0, param.getDIFFEPS())) {
              error_msg(errstr,
                  "Depth %e is negative! Point cannot be above the SIC.\n",
                  curr_SIC_depth);
              writeErrorToII(errstr, GlobalVariables::log_file);
              exit(1);
            } else if (lessThanVal(curr_SIC_depth, 0.0, param.getEPS())) {
              warning_msg(errstr,
                  "Depth %e is negative! Point cannot be above the SIC. Small enough violation that only sending warning.\n",
                  curr_SIC_depth);
            }
            if (isZero(curr_SIC_depth, param.getDIFFEPS())) {
              curr_SIC_depth = +0.0;
            } 

            // Get depth wrt objective
            const double objActivity = objCut.getActivity(tmpPoint);
            double curr_obj_depth = (objActivity - objCut.rhs()) / objNorm;
            if (isZero(curr_obj_depth)) {
              curr_obj_depth = +0.0;
            }
            if ((objActivity - objCut.rhs()) < this->min_nb_obj_val) {
              this->min_nb_obj_val = curr_obj_depth;
            }
            // interPtsAndRays.objDepth.push_back(curr_obj_depth); // this caused the objective depth to be improperly calculated

            interPtsAndRays.addPointOrRayToIntersectionInfo(tmpPoint, 1.0,
                ray_to_cut, newRayColIndices[r], hplane_ind,
                hplane_ind_in_split, vertex_ind, -1, curr_SIC_depth,
                curr_obj_depth, pointIsFinal);

            sumAddedDepth += curr_SIC_depth;
          }
        }
      } /* if new ray intersects split bd, it creates a new point */
      else if (!CHECK_FOR_DUP_RAYS
          || !duplicateRay(interPtsAndRays, interPtsAndRays.RHS, newRay)) { // else, add the ray
        // Ensure that no ``future'' hyperplane cuts the point
        bool isCut = false;
        for (int h_ind = hplane_ind_in_split;
            h_ind < (int) interPtsAndRays.allHplaneToAct.size(); h_ind++) {
          const int curr_hplane_ind = interPtsAndRays.allHplaneToAct[h_ind];
          // The hyperplane needs to have cut the current ray
          const bool hplaneCutsThisRay = (find_val(ray_to_cut,
              hplaneStore[curr_hplane_ind].rayToBeCutByHplane[split_ind]) >= 0);
          if (hplaneCutsThisRay
              && pointIsCutByHplane(solver, solnInfo, newRay, 0.0,
                  hplaneStore[curr_hplane_ind])) {
            isCut = true;
            break;
          }
        }
        if (!isCut) {
          const bool rayIsFinal = isRayFinal(solver, solnInfo,
              newRay.getVectorNumElements(), newRay.getVectorIndices(),
              newRay.getVectorElements());
          interPtsAndRays.addPointOrRayToIntersectionInfo(newRay, 0.0,
              ray_to_cut, newRayColIndices[r], hplane_ind, hplane_ind_in_split,
              vertex_ind, -1, 0.0, 0.0, rayIsFinal);
        }
      } /* if new ray does not intersect split bd, add the ray */
    } /* end iterating over new rays to be created */
  } /* end iterating over rays intersected by hplane in this activation */

  if (numNewPointsAdded > 0) {
    avgAddedDepth += sumAddedDepth / numNewPointsAdded;
    totalNumPointsAdded += numNewPointsAdded;
  }
} /* algorithm1 (PHA1) */

/**
 * ALGORITHM 2
 *
// * Input (from paper, NOT to this function):
// * convex rational polyhedron P defined by a set of hyperplanes H;
// * set I \subseteq [n] denoting the integer variables;
// * objective coefficient vector c;
// * classes of objectives O;
// * number non-tilted hyperplanes k_h;
// * fractionality threshold \eps_{frac}.
 */
void CglPHA::algorithm2(AdvCuts &structPHA, const SolutionInfo& solnInfo,
    const AdvCuts &structSICs, const AdvCuts &NBSICs,
    Stats& phaTimeStats, std::vector<int>& num_gics_per_split, AdvCut& NBObj,
    const double objNorm,
    std::vector<std::vector<double> >& distToSplitAlongRay,
    std::vector<Hplane>& hplaneStore, std::vector<Vertex>& vertexStore,
    std::vector<Ray>& rayStore, std::vector<int>& sortedRaysIndex,
//    std::vector<int>& allRayIndices,
    std::vector<int>& allRaysCutAcrossSplits,
    std::vector<std::vector<int> >& raysCutForSplit) {
  HplaneScoringFn option = HplaneScoringFn::final;

  const int num_rays_to_be_cut = param.getParamVal(ParamIndices::NUM_RAYS_CUT_PARAM_IND);
  const int num_alg3_rounds =
        param.getParamVal(ParamIndices::NUM_ALG2_ROUNDS_PARAM_IND) < 0 ?
            0 : std::ceil(std::log2(num_rays_to_be_cut + 1));
  const int num_alg2_rounds = std::abs(param.getParamVal(ParamIndices::NUM_ALG2_ROUNDS_PARAM_IND));
  const int num_total_act_rounds = num_alg3_rounds + num_alg2_rounds;

  for (int alg2_act_ind = 0; alg2_act_ind < num_alg2_rounds; alg2_act_ind++) {
    const bool isTimeLimitReached = reachedTimeLimit(phaTimeStats, "PHA_TOTAL", param.getTIMELIMIT());
    if (reachedCutLimit(0, num_total_gics) || isTimeLimitReached) {
      break;
    }
    const int int_option = static_cast<int>(option);
    printf(
        "\n## Start act %d/%d (``alg2'' round %d/%d): adds extra hyperplane per split via PHA1(option:%d,%s), %s. ##\n",
        num_alg3_rounds + alg2_act_ind + 1, num_total_act_rounds,
        alg2_act_ind + 1, num_total_act_rounds,int_option,
        HplaneScoringFnName[int_option].c_str(),
        HplaneScoringFnDescription[int_option].c_str());
    PHAHelper(0, solnInfo.numNB, solnInfo, structSICs, NBSICs, NBObj, objNorm,
        distToSplitAlongRay, hplaneStore, vertexStore,
        rayStore, sortedRaysIndex, /*allRayIndices,*/ allRaysCutAcrossSplits,
        raysCutForSplit, interPtsAndRays, num_alg3_rounds + alg2_act_ind,
        num_total_act_rounds, true, option, structPHA, num_gics_per_split,
        num_total_gics, phaTimeStats);
  } /* algorithm 2 activations (one hyperplane per split per round) */
} /* algorithm2 */

/**
 * ALGORITHM 3
 * TARGETED TILTING
 * Description: Activations that focus on cutting one ray at a time
 *
// * Input (from paper, NOT to this function):
// * polyhedron P defined by a set of hyperplanes H;
// * a vertex \bar{x} of P;
// * indices of fractional integer variables \sigma;
// * hyperplane selection criterion SC;
// * classes of objectives O;
// * set R_A \subseteq \bar{R} for the rays of \bar{C} to cut.
 */
void CglPHA::algorithm3(AdvCuts &structPHA, const SolutionInfo& solnInfo,
    const AdvCuts &structSICs, const AdvCuts &NBSICs,
    Stats& phaTimeStats, std::vector<int>& num_gics_per_split, AdvCut& NBObj,
    const double objNorm,
    std::vector<std::vector<double> >& distToSplitAlongRay,
    std::vector<Hplane>& hplaneStore, std::vector<Vertex>& vertexStore,
    std::vector<Ray>& rayStore, std::vector<int>& sortedRaysIndex,
//    std::vector<int>& allRayIndices,
    std::vector<int>& allRaysCutAcrossSplits,
    std::vector<std::vector<int> >& raysCutForSplit) {
  const int num_rays_to_be_cut = param.getParamVal(ParamIndices::NUM_RAYS_CUT_PARAM_IND);
  const int num_alg3_rounds = std::ceil(std::log2(num_rays_to_be_cut + 1)); // see below
  std::vector<int> numRaysCutInFirstActRound;

  // Set up number rays cuts in each of the first activation rounds
  // We double the number of rays cut each time (we generate cuts after each round)
  numRaysCutInFirstActRound.resize(num_alg3_rounds, 0);
  numRaysCutInFirstActRound[0] = 1;
  if (num_alg3_rounds > 1) {
    for (int i = 1; i < num_alg3_rounds - 1; i++) {
      numRaysCutInFirstActRound[i] = (2 * numRaysCutInFirstActRound[i - 1]);
    }
    numRaysCutInFirstActRound[num_alg3_rounds - 1] = num_rays_to_be_cut
        - std::pow(2, num_alg3_rounds - 1) + 1;
  }

  const int num_alg2_rounds = std::abs(param.getParamVal(ParamIndices::NUM_ALG2_ROUNDS_PARAM_IND));
  const int num_total_act_rounds = num_alg3_rounds + num_alg2_rounds;

  const int int_option = param.getParamVal(ParamIndices::HPLANE_SCORING_FN_PARAM_IND);
  HplaneScoringFn option = static_cast<HplaneScoringFn>(int_option);
  int curr_num_rays_cut = 0;
  for (int alg3_act_ind = 0; alg3_act_ind < num_alg3_rounds; alg3_act_ind++) {
    const bool isTimeLimitReached = reachedTimeLimit(phaTimeStats, "PHA_TOTAL", param.getTIMELIMIT());
    if (reachedCutLimit(0, num_total_gics) || isTimeLimitReached) {
      break;
    }
    printf(
        "\n## Start act %d/%d (``alg3'' activation round %d/%d) of PHA1(option:%d,%s), %s. ##\n",
        alg3_act_ind + 1, num_total_act_rounds, alg3_act_ind + 1, num_alg3_rounds, int_option,
        HplaneScoringFnName[int_option].c_str(),
        HplaneScoringFnDescription[int_option].c_str());
    PHAHelper(curr_num_rays_cut, numRaysCutInFirstActRound[alg3_act_ind],
        solnInfo, structSICs, NBSICs, NBObj, objNorm, distToSplitAlongRay,
        hplaneStore, vertexStore, rayStore, sortedRaysIndex,
//        allRayIndices,
        allRaysCutAcrossSplits, raysCutForSplit, interPtsAndRays,
        alg3_act_ind, num_total_act_rounds, false, option, structPHA,
        num_gics_per_split, num_total_gics, phaTimeStats);
    curr_num_rays_cut += numRaysCutInFirstActRound[alg3_act_ind];
  } /* algorithm 3 activations */
} /* algorithm3 (targeted tilting) */

void CglPHA::runAnalyses(const SolutionInfo* const solnInfoPtr,
    std::vector<Hplane>& hplaneStore, std::vector<int>& num_SIC_points,
    std::vector<int>& num_SIC_final_points, 
    std::vector<int>& num_SIC_final_rays, 
    std::vector<int>& num_SIC_rays,
    const AdvCuts& structSICs, const AdvCuts& NBSICs, 
		const AdvCuts& structPHA) {
  const int numSICs = structSICs.sizeCuts();
  const int numGICs = structPHA.sizeCuts();
  const int numSplits = solnInfoPtr->numFeasSplits;

#ifdef SHOULD_CALC_ACTIVE_FINAL_POINTS
  int numActiveSICs = 0, numActiveGICs = 0;
  {
  // Per request, calculate average number of tight final points across cuts and splits
  // We calculate two things
  // (1) for each cut (SICs and GICs), number of final points & rays that the cut is tight on (across all splits)
  //      --> we report the avg number of tight final points/rays, averaging over the number of cuts and splits
  // (2) for each _active_ cut (SICs and GICs), number of final points & rays that the cut is tight on (across all splits)
  PHASolverInterface* PHASolver = dynamic_cast<PHASolverInterface*>(solver->clone()); // this will fail if we use anything but OsiClp
  setMessageHandler(PHASolver);
  PHASolver->applyCuts(structSICs);
  PHASolver->applyCuts(structPHA);
  PHASolver->isProvenOptimal() ? PHASolver->resolve() : PHASolver->initialSolve();
  checkSolverOptimality(PHASolver, false);

  // Check activity of SICs and GICs
  std::vector<bool> SIC_is_active(structSICs.sizeCuts());
  for (int cut_ind = 0; cut_ind < structSICs.sizeCuts(); cut_ind++) {
    const CoinPackedVector currCut = structSICs.cuts[cut_ind].row();
    const double activity = dotProduct(currCut.getNumElements(),
        currCut.getIndices(), currCut.getElements(),
        PHASolver->getColSolution());
    if (isVal(activity, structSICs[cut_ind].rhs())) {
      numActiveSICs++;
      SIC_is_active[cut_ind] = true;
    } else {
      SIC_is_active[cut_ind] = false;
    }
  }

  std::vector<bool> GIC_is_active(structPHA.sizeCuts());
  for (int cut_ind = 0; cut_ind < structPHA.sizeCuts(); cut_ind++) {
    const CoinPackedVector currCut = structPHA.cuts[cut_ind].row();
    const double activity = dotProduct(currCut.getNumElements(),
        currCut.getIndices(), currCut.getElements(),
        PHASolver->getColSolution());
    if (isVal(activity, structPHA[cut_ind].rhs())) {
      numActiveGICs++;
      GIC_is_active[cut_ind] = true;
    } else {
      GIC_is_active[cut_ind] = false;
    }
  }

  //assert(structPHA.sizeCuts() == (int) GlobalVariables::cutCoeffs.size());
  //assert(structPHA.sizeCuts() == (int) GlobalVariables::cutRHS.size());
  GlobalVariables::numTightFinalPointsGICs.resize(numSplits, 0);
  GlobalVariables::numTightFinalRaysGICs.resize(numSplits, 0);
  GlobalVariables::numTightFinalPointsSICs.resize(numSplits, 0);
  GlobalVariables::numTightFinalRaysSICs.resize(numSplits, 0);
  GlobalVariables::numTightFinalPointsActiveGICs.resize(numSplits, 0);
  GlobalVariables::numTightFinalRaysActiveGICs.resize(numSplits, 0);
  GlobalVariables::numTightFinalPointsActiveSICs.resize(numSplits, 0);
  GlobalVariables::numTightFinalRaysActiveSICs.resize(numSplits, 0);
  int totalNumFinalPoints = 0, totalNumFinalRays = 0;
  for (int s = 0; s < numSplits; s++) {
    for (int p = 0; p < interPtsAndRays[s].getNumRows(); p++) {
      if (!interPtsAndRays[s].finalFlag[p]) {
        continue;
      }
      bool isRay = isZero(interPtsAndRays[s].RHS[p]);
      if (!isRay) {
        totalNumFinalPoints++;
      } else {
        totalNumFinalRays++;
      }
      const CoinShallowPackedVector vec = interPtsAndRays[s].getVector(p);

      // Calculate tight points for GICs
      for (int cut_ind = 0; cut_ind < structPHA.sizeCuts(); cut_ind++) {
        const double row_activity = dotProduct(vec.getNumElements(),
            vec.getIndices(), vec.getElements(),
            GlobalVariables::cutCoeffsNB[cut_ind].data());
        if (isVal(row_activity, GlobalVariables::cutRHSNB[cut_ind])) {
          if (!isRay) {
            GlobalVariables::numTightFinalPointsGICs[s]++;
            GlobalVariables::totalNumTightFinalPointsGICs++;
          } else {
            GlobalVariables::numTightFinalRaysGICs[s]++;
            GlobalVariables::totalNumTightFinalRaysGICs++;
          }

          // Now increment the relevant vectors if the cut is also active at the new optimum
          if (GIC_is_active[cut_ind]) {
            if (!isRay) {
              GlobalVariables::numTightFinalPointsActiveGICs[s]++;
              GlobalVariables::totalNumTightFinalPointsActiveGICs++;
            } else {
              GlobalVariables::numTightFinalRaysActiveGICs[s]++;
              GlobalVariables::totalNumTightFinalRaysActiveGICs++;
            }
          }
        } /* check if GIC is tight for this point/ray */
      } /* iterate over GICs */

      // Calculate tight points for SICs
      for (int cut_ind = 0; cut_ind < (int) NBSICs.cuts.size(); cut_ind++) {
        const CoinPackedVector currCut = NBSICs.cuts[cut_ind].row();
        const double row_activity = dotProduct(vec.getNumElements(),
            vec.getIndices(), vec.getElements(), currCut.getNumElements(),
            currCut.getIndices(), currCut.getElements());
        if (isVal(row_activity, NBSICs.cuts[cut_ind].rhs())) {
          if (!isRay) {
            GlobalVariables::numTightFinalPointsSICs[s]++;
            GlobalVariables::totalNumTightFinalPointsSICs++;
          } else {
            GlobalVariables::numTightFinalRaysSICs[s]++;
            GlobalVariables::totalNumTightFinalRaysSICs++;
          }

          // Now increment the relevant vectors if the cut is also active at the new optimum
          if (SIC_is_active[cut_ind]) {
            if (!isRay) {
              GlobalVariables::numTightFinalPointsActiveSICs[s]++;
              GlobalVariables::totalNumTightFinalPointsActiveSICs++;
            } else {
              GlobalVariables::numTightFinalRaysActiveSICs[s]++;
              GlobalVariables::totalNumTightFinalRaysActiveSICs++;
            }
          }
        } /* check if SIC is tight for this point/ray */
      } /* iterate over SICs */
    } /* iterate over points/rays */
  } /* iterate over splits */

  if (PHASolver) {
    delete PHASolver;
  }

#if 0
  { // DEBUG DEBUG DEBUG
    int totalNumSICPoints = 0, totalNumSICRays = 0;
    int totalNumPoints = 0, totalNumRays = 0;
    int tmpTotalNumFinalPoints = 0, tmpTotalNumFinalRays = 0;
    int totalNumSICFinalPoints = 0, totalNumSICFinalRays = 0;
    for (int s = 0; s < numSplits; s++) {
      totalNumSICPoints += num_SIC_points[s]++;
      totalNumSICRays += num_SIC_rays[s]++;
      totalNumPoints += interPtsAndRays[s].numPoints;
      totalNumRays += interPtsAndRays[s].numRays;
      tmpTotalNumFinalPoints += interPtsAndRays[s].numFinalPoints;
      tmpTotalNumFinalRays += interPtsAndRays[s].numFinalRays;
      totalNumSICFinalPoints += num_SIC_final_points[s]++;
      totalNumSICFinalRays += num_SIC_final_rays[s]++;
    }
    assert(totalNumFinalPoints == tmpTotalNumFinalPoints);
    assert(totalNumFinalRays == tmpTotalNumFinalRays);
    assert(totalNumSICFinalPoints <= tmpTotalNumFinalPoints);
    assert(totalNumSICFinalRays <= tmpTotalNumFinalRays);
    printf("Num splits: %d\tNum vars: %d\tNum SICs: %d\tNum GICs: %d\n",
        numSplits, solnInfoPtr->numCols, (int) NBSICs.cuts.size(),
        structPHA.sizeCuts());
    printf(
        "SIC points: %d\tSIC final points: %d\tSIC rays: %d\tSIC final rays: %d\n",
        totalNumSICPoints, totalNumSICFinalPoints, totalNumSICRays,
        totalNumSICFinalRays);
    printf("Total points: %d\tTotal rays: %d\n", totalNumPoints, totalNumRays);
    printf("Total final points: %d\tAvg tight on SICs: %.2f\tAvg tight on GICs: %.2f\n",
        totalNumFinalPoints,
        (double) totalNumTightFinalPointsSICs / (numSICs * numSplits),
        (double) totalNumTightFinalPointsGICs / (numGICs * numSplits));
    printf("Total final rays: %d\tAvg tight on SICs: %.2f\tAvg tight on GICs: %.2f\n",
        totalNumFinalRays,
        (double) totalNumTightFinalRaysSICs / (numSICs * numSplits),
        (double) totalNumTightFinalRaysGICs / (numGICs * numSplits));
    printf("Total active SICs: %d\tAvg tight on active SICs: %.2f\n",
        numActiveSICs,
        (double) totalNumTightFinalPointsActiveSICs / (numActiveSICs * numSplits));
    printf("Total active GICs: %d\tAvg tight on active GICs: %.2f\n",
        numActiveGICs,
        (double) totalNumTightFinalPointsActiveGICs / (numActiveGICs * numSplits));
    exit(1);
  } // DEBUG DEBUG DEBUG
#endif
  } /* big ugly container so that we can keep some variables local */
#endif

#ifdef TRACE
  printf("\n## Analyzing last point-ray collection. ##\n");
#endif
  GlobalVariables::timeStats.start_timer(COLLECT_STATS_TIME);
  std::vector<double> minSICDepth(solnInfoPtr->numFeasSplits), maxSICDepth(
      solnInfoPtr->numFeasSplits), avgSICDepth(solnInfoPtr->numFeasSplits);
  std::vector<double> minObjDepth(solnInfoPtr->numFeasSplits), maxObjDepth(
      solnInfoPtr->numFeasSplits), avgObjDepth(solnInfoPtr->numFeasSplits);
  for (int split_ind = 0; split_ind < (int) interPtsAndRays.size();
      split_ind++) {
    minSICDepth[split_ind] = interPtsAndRays[split_ind].minSICDepth;
    maxSICDepth[split_ind] = interPtsAndRays[split_ind].maxSICDepth;
    avgSICDepth[split_ind] = interPtsAndRays[split_ind].totalSICDepth
        / interPtsAndRays[split_ind].numPoints;
    minObjDepth[split_ind] = interPtsAndRays[split_ind].minObjDepth;
    maxObjDepth[split_ind] = interPtsAndRays[split_ind].maxObjDepth;
    avgObjDepth[split_ind] = interPtsAndRays[split_ind].totalObjDepth
        / interPtsAndRays[split_ind].numPoints;
  }

#if defined TRACE && defined SHOULD_WRITE_INTPOINTS
  printf("\n## Save information about the points and rays. ##\n");
  char filename[300];
  snprintf(filename, sizeof(filename) / sizeof(char),
      "%s-PointsAndRays.csv", GlobalVariables::out_f_name_stub.c_str());
  FILE *ptsRaysOut = fopen(filename, "w");
  writeInterPtsAndRaysInJSpace(ptsRaysOut, *solnInfoPtr,
      interPtsAndRays, hplaneStore, num_SIC_points,
      num_SIC_final_points, num_SIC_rays);
  fclose(ptsRaysOut);
#endif

  // Write GIC point info to II
  int totalNumPoints = 0, totalNumFinalPoints = 0;
  int totalNumRays = 0, totalNumFinalRays = 0;
  std::vector<int> num_points_per_split(solnInfoPtr->numFeasSplits);
  std::vector<int> num_final_points_per_split(solnInfoPtr->numFeasSplits);
  std::vector<int> num_rays_per_split(solnInfoPtr->numFeasSplits);
  std::vector<int> num_final_rays_per_split(solnInfoPtr->numFeasSplits);
  for (int s = 0; s < solnInfoPtr->numFeasSplits; s++) {
    num_points_per_split[s] = interPtsAndRays[s].numPoints;
    num_final_points_per_split[s] = interPtsAndRays[s].numFinalPoints;
    num_rays_per_split[s] = interPtsAndRays[s].numRays;
    num_final_rays_per_split[s] = interPtsAndRays[s].numFinalRays;
    totalNumPoints += num_points_per_split[s];
    totalNumFinalPoints += num_final_points_per_split[s];
    totalNumRays += num_rays_per_split[s];
    totalNumFinalRays += num_final_rays_per_split[s];
  }

  const double avg_num_SIC_final_points = computeAverage(num_SIC_final_points);
  const double avg_num_SIC_final_rays = computeAverage(num_SIC_final_rays);
  const double avg_num_par_rays = computeAverage(num_SIC_rays);
#ifdef SHOULD_CALC_ACTIVE_FINAL_POINTS
  const double avgNumTightFinalPointsSICs =
      (numSICs > 0) ? (double) GlobalVariables::totalNumTightFinalPointsSICs / (numSICs * numSplits) : 0.;
  const double avgNumTightFinalPointsGICs =
      (numGICs > 0) ? (double) GlobalVariables::totalNumTightFinalPointsGICs / (numGICs * numSplits) : 0.;
  const double avgNumTightFinalRaysSICs =
      (numSICs > 0) ? (double) GlobalVariables::totalNumTightFinalRaysSICs / (numSICs * numSplits) : 0.;
  const double avgNumTightFinalRaysGICs =
      (numGICs > 0) ? (double) GlobalVariables::totalNumTightFinalRaysGICs / (numGICs * numSplits) : 0.;
  const double avgNumTightFinalPointsActiveSICs =
      (numActiveSICs > 0) ?
          (double) GlobalVariables::totalNumTightFinalPointsActiveSICs
              / (numActiveSICs * numSplits) :
          0.;
  const double avgNumTightFinalPointsActiveGICs =
      (numActiveGICs > 0) ?
          (double) GlobalVariables::totalNumTightFinalPointsActiveGICs
              / (numActiveGICs * numSplits) :
          0.;
  const double avgNumTightFinalRaysActiveSICs =
      (numActiveSICs > 0) ?
          (double) GlobalVariables::totalNumTightFinalRaysActiveSICs
              / (numActiveSICs * numSplits) :
          0.;
  const double avgNumTightFinalRaysActiveGICs =
      (numActiveGICs > 0) ?
          (double) GlobalVariables::totalNumTightFinalRaysActiveGICs
              / (numActiveGICs * numSplits) :
          0.;
#else
  const double avgNumTightFinalPointsSICs = 0.;
  const double avgNumTightFinalPointsGICs = 0.;
  const double avgNumTightFinalRaysSICs = 0.;
  const double avgNumTightFinalRaysGICs = 0.;
  const double avgNumTightFinalPointsActiveSICs = 0.;
  const double avgNumTightFinalPointsActiveGICs = 0.;
  const double avgNumTightFinalRaysActiveSICs = 0.;
  const double avgNumTightFinalRaysActiveGICs = 0.;
#endif
  writeGICPointInfoToII(param.getParamVal(ParamIndices::HPLANE_SCORING_FN_PARAM_IND),
      param.getParamVal(ParamIndices::NUM_ALG2_ROUNDS_PARAM_IND),
      solnInfoPtr->numNB, avg_num_par_rays, avg_num_SIC_final_points, avg_num_SIC_final_rays,
      param.getParamVal(ParamIndices::NUM_RAYS_CUT_PARAM_IND),
      param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND),
      num_points_per_split, num_final_points_per_split, num_rays_per_split, num_final_rays_per_split,
      totalNumPoints, totalNumFinalPoints, totalNumRays, totalNumFinalRays,
      avgNumTightFinalPointsSICs, avgNumTightFinalPointsGICs,
      avgNumTightFinalRaysSICs, avgNumTightFinalRaysGICs,
      avgNumTightFinalPointsActiveSICs, avgNumTightFinalPointsActiveGICs,
      avgNumTightFinalRaysActiveSICs, avgNumTightFinalRaysActiveGICs,
      avgSICDepth, minSICDepth, maxSICDepth, avgObjDepth, minObjDepth, maxObjDepth,
      GlobalVariables::log_file);
  GlobalVariables::timeStats.end_timer(COLLECT_STATS_TIME);
} /* runAnalyses */

/**
 * PHA Helper
 */
void CglPHA::PHAHelper(const int indexOfFirstRayToCut, const int numRaysToBeCut,
    const SolutionInfo& solnInfo, const AdvCuts &structSICs,
    const AdvCuts &NBSICs, const AdvCut& NBObj, const double objNorm,
    const std::vector<std::vector<double> >& distToSplitAlongRay,
    //std::vector<std::vector<bool> >& pointIsFinal,
    std::vector<Hplane>& hplaneStore, std::vector<Vertex>& vertexStore,
    std::vector<Ray>& rayStore, const std::vector<int>& sortIndex,
//    std::vector<int>& allRayIndices,
    std::vector<int> &allRaysCutAcrossSplits,
    std::vector<std::vector<int> > &raysCutForSplit,
    std::vector<IntersectionInfo>& interPtsAndRays, const int act_ind,
    const int total_num_act, const bool isAlg2,
    const HplaneScoringFn option, AdvCuts& structPHA,
    std::vector<int>& num_gics_per_split, int& num_total_pha, Stats& phaTimeStats) {
#ifdef TRACE
  const int int_option = (int) option;
#endif

  // Choose hyperplanes and update vertices
  GlobalVariables::timeStats.start_timer(CHOOSE_HPLANES_TIME);

  /***********************************************************************************
   * Alg 2 step 7 or Alg 3 steps 9 to 11
   ***********************************************************************************/
  // Either one hyperplane is chosen per each split for this activation,
  // or one ray is cut (at a time) and a hyperplane is chosen for each ray
  if (isAlg2) {
#ifdef TRACE
    printf(
        "Act %d/%d: PHA1(option:%d,%s): Choosing extra hplane for all splits.\n",
        act_ind + 1, total_num_act, int_option,
        HplaneScoringFnName[int_option].c_str());
#endif
    chooseHplaneBySplit(//allRayIndices,
        0, numRaysToBeCut, sortIndex, castsolver, solnInfo, hplaneStore,
        vertexStore, rayStore, NBSICs, distToSplitAlongRay,
        allRaysCutAcrossSplits, raysCutForSplit, interPtsAndRays, act_ind,
        (option == HplaneScoringFn::neighbor),
        (option == HplaneScoringFn::depth), (option == HplaneScoringFn::final),
        (option == HplaneScoringFn::mostremoved));
  } /* choose hyperplane by split? (when in alg2) */
  else {
    for (int ray_ind = indexOfFirstRayToCut; ray_ind < indexOfFirstRayToCut + numRaysToBeCut; ray_ind++) {
#ifdef TRACE
      const int nonBasicVarIndex = solnInfo.nonBasicVarIndex[sortIndex[ray_ind]];
      const std::string nb_var_name =
          (nonBasicVarIndex < solnInfo.numCols) ?
              solver->getColName(nonBasicVarIndex) :
              solver->getRowName(nonBasicVarIndex - solnInfo.numCols);
      printf(
          "Act %d/%d: PHA1(option:%d,%s): Choosing hplane for ray %d/%d (index %d, variable %d, and name %s).\n",
          act_ind + 1, total_num_act, int_option,
          HplaneScoringFnName[int_option].c_str(),
          ray_ind - indexOfFirstRayToCut + 1, numRaysToBeCut,
          sortIndex[ray_ind], nonBasicVarIndex, nb_var_name.c_str());
#endif

      Ray rayOut;
      const int rayIndices[1] = { sortIndex[ray_ind] };
      const double rayValues[1] = { 1.0 };
      rayOut.setVector(1, rayIndices, rayValues);

      // Find hyperplane that intersects ray
      // either (0) first, (1) yielding highest avg depth of resulting points, or (2) yielding most final points
      // Last option allows you to also calculate number of points removed and use that as a secondary criterion
      // Returns hplane.var = -1 if no hyperplane is intersected by the ray before bd S for all splits
      chooseHplaneByRay(//allRayIndices,
          0, rayOut, sortIndex[ray_ind],
          castsolver, solnInfo, hplaneStore, vertexStore, rayStore, NBSICs,
          distToSplitAlongRay, allRaysCutAcrossSplits, raysCutForSplit,
          interPtsAndRays, act_ind, (option == HplaneScoringFn::neighbor),
          (option == HplaneScoringFn::depth),
          (option == HplaneScoringFn::final), false);
    } /* end iterating over rays */
  } /* choose hyperplane by ray? (when in alg3) */

  GlobalVariables::timeStats.end_timer(CHOOSE_HPLANES_TIME);

#ifdef SHOULD_WRITE_HPLANES
  printf(
      "\n## Save information about which hyperplanes will be activated. ##\n");
  char filename[300];
  snprintf(filename, sizeof(filename) / sizeof(char), "%s-hplanesAct%dOpt%d.csv",
      GlobalVariables::out_f_name_stub.c_str(), act_ind,
      static_cast<int>(option));
  FILE *hplaneOut = fopen(filename, "w");
  printActivatedHplaneDetails(hplaneOut, solnInfo, hplaneStore,
      interPtsAndRays);
  fclose(hplaneOut);
#endif

  /***********************************************************************************
   * Call Algorithm 1 (Alg 2 step 8 or Alg 3 step 12)
   ***********************************************************************************/
  // Create new intersection points and rays
  GlobalVariables::timeStats.start_timer(GEN_INTERSECTION_POINTS_TIME);
#ifdef TRACE
  printf("\n## Act %d/%d: Perform activations. ##\n", act_ind + 1,
      total_num_act);
#endif
  performActivations(solnInfo, vertexStore, rayStore, hplaneStore,
      interPtsAndRays, NBSICs, NBObj, objNorm, act_ind);
  GlobalVariables::timeStats.end_timer(GEN_INTERSECTION_POINTS_TIME);

  /***********************************************************************************
   * Generate cuts (Alg 2 step 9 or Alg 3 step 14)
   ***********************************************************************************/
  if (!GlobalConstants::NO_GEN_CUTS) {
    // Generate cuts after activation
#ifdef TRACE
    printf(
        "\n## Act %d/%d, PHA1(option:%d,%s): starting to try objectives (total so far: %d) ##\n",
        act_ind + 1, total_num_act, int_option,
        HplaneScoringFnName[int_option].c_str(), num_total_pha);
#endif
    tryObjectives(structPHA, num_gics_per_split, num_total_pha, num_obj_tried,
        solnInfo, /*solver->getNumCols(),*/ interPtsAndRays, /*rayStore,*/
        vertexStore, hplaneStore, structSICs, phaTimeStats);
    printf("\n## Finishing act %d/%d with %d cuts. ##\n", act_ind + 1, total_num_act, num_total_pha);
  }
} /* PHAHelper */

/**
 * Given a particular ray, we will intersect it with a hyperplane based on a given selection criterion
 * This is the method used for Algorithm 3 (step 9)
 *
 * The selection criteria are:
 * (1) selectByDist: choose the hplane that cuts the given ray first
 * (2) selectByDepth: choose the hplane for which the resulting points have greatest average depth
 * (3) selectByFinalPoints: choose the hplane that yields the most number of final points
 */
bool CglPHA::chooseHplaneByRay(//const std::vector<int>& allRayIndices,
    const int originatingVertexIndex, const Ray& rayOut, const int ray_ind,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
    std::vector<Hplane>& hplaneStore, std::vector<Vertex>& vertexStore,
    std::vector<Ray>& rayStore, const AdvCuts& NBSIC,
    const std::vector<std::vector<double> >& distToSplitAlongRay,
    std::vector<int>& allRaysCutAcrossSplits,
    std::vector<std::vector<int> >& raysCutForSplit,
    std::vector<IntersectionInfo>& interPtsAndRays, const int act_ind,
    const bool selectByDist, const bool selectByDepth,
    const bool selectByFinalPoints, const bool secondarySelectByMostRemoved) {
  // Initialize values, in case no hplane is found at all
  const int numSplitsToCheck = solnInfo.numFeasSplits; //splitsAffByHplane.size();
  Hplane tmpHplane(HplaneVarFlag::none, HplaneRowFlag::none,
      HplaneBoundFlag::none, 1); // var, row, ubflag, numsplits
  std::vector<Hplane> hplane(numSplitsToCheck);

  // For each of the splits we care about (splitsAffByHplane) look at:
  // (1) min distance to the hyperplane
  double minDistOverall = solver->getInfinity();
  std::vector<double> distToHplaneChosenForSplit(numSplitsToCheck, -1.0);
  // (2) the average depth of the points that would be created if we activated each hplane
  std::vector<double> addedDepthByHplaneChosenForSplit(numSplitsToCheck, -1.0);
  // (3) total number of final points that would be created
  std::vector<int> numFinalPointsByHplaneChosenForSplit(numSplitsToCheck, -1);
  // (4) number of points removed by the hyperplane
  int numPointsCut = 0; // not used by this function

  std::vector<Vertex> vertexCreatedForSplit(numSplitsToCheck);

  // Check the basic variable in each row to see if it can be activated
  if (ACTIVATE_BASIC_BOUNDS) {
    for (int row_ind = 0; row_ind < solver->getNumRows(); row_ind++) {
      if (solver->getRowSense()[row_ind] == 'E') {
        continue; // Do not activate equality rows
      }
      tmpHplane.var = solnInfo.varBasicInRow[row_ind];
      const bool isStructVar = (tmpHplane.var < solver->getNumCols());
      tmpHplane.row = row_ind;

      const double varVal = getVarVal(solver, tmpHplane.var);
      const double rayVal = solnInfo.raysOfC1[ray_ind][row_ind]; // rayMult * -1 * currRay[r];

      tmpHplane.bound = 0.0;
      if (isZero(rayVal, param.getRAYEPS())) {
        continue; // ray does not intersect this hyperplane
      } else if (lessThanVal(rayVal, 0.0, param.getRAYEPS())) {
        tmpHplane.bound = getVarLB(solver, tmpHplane.var);
        tmpHplane.ubflag = static_cast<int>(HplaneBoundFlag::toLB);
      } else if (greaterThanVal(rayVal, 0.0, param.getRAYEPS())) {
        tmpHplane.bound = getVarUB(solver, tmpHplane.var);
        tmpHplane.ubflag = static_cast<int>(HplaneBoundFlag::toUB);
      }

      // If this bound is actually at infinity; we will never intersect it
      if (isInfinity(std::abs(tmpHplane.bound), solver->getInfinity())) {
        continue;
      }

      if (!isStructVar) {
        tmpHplane.ubflag = static_cast<int>(HplaneBoundFlag::toLBSlack);
      }

      // Otherwise, find the distance to the hyperplane, and ensure it is actually intersected
      const double tmpDistToHplane = (tmpHplane.bound - varVal) / rayVal;
      if (isInfinity(tmpDistToHplane, solver->getInfinity())) {
        continue; // ray does not intersect this hyperplane
      }

      // Set packed vector corresponding to this vertex
      Vertex vertCreated;
      calcNewVectorCoordinates(vertCreated, tmpDistToHplane,
          vertexStore[originatingVertexIndex], rayOut);

      std::vector<double> addedDepthForThisHplane(numSplitsToCheck, 0.0);
      std::vector<int> numFinalPointsThisHplane(numSplitsToCheck, 0);

      for (int split_ind = 0; split_ind < numSplitsToCheck; split_ind++) {
        const int split_var = solnInfo.feasSplitVar[split_ind];
        if (split_var == tmpHplane.var) {
          continue; // Skip if the split var bound is the one being activated
        }

        // Check distance to split
        const double distToSplit =
            (originatingVertexIndex > 0) ?
                distanceToSplit(vertexStore[originatingVertexIndex], rayOut,
                    split_var, solver, solnInfo) :
                distToSplitAlongRay[split_ind][ray_ind];

        // Check split value
        const double split_row = solnInfo.rowOfVar[split_var];
        const double splitRayVal = solnInfo.raysOfC1[ray_ind][split_row];
        const double vertexSplitVal = solnInfo.a0[split_row]
            + splitRayVal * tmpDistToHplane;
        const bool isVertexOnFloor = !greaterThanVal(vertexSplitVal,
            std::floor(solnInfo.a0[split_row]), param.getAWAY());
        const bool isVertexOnCeil = !lessThanVal(vertexSplitVal,
            std::ceil(solnInfo.a0[split_row]), param.getAWAY());
        const bool isVertexOnSplit = isVertexOnFloor || isVertexOnCeil;

        if (chooseHplaneHelper(split_ind, distToSplit, isVertexOnSplit, ray_ind,
            hplane, tmpHplane, tmpDistToHplane, minDistOverall,
            distToHplaneChosenForSplit[split_ind],
            addedDepthByHplaneChosenForSplit[split_ind],
            numFinalPointsByHplaneChosenForSplit[split_ind], numPointsCut,
            vertCreated, addedDepthForThisHplane[split_ind],
            numFinalPointsThisHplane[split_ind], numPointsCut, //allRayIndices,
            solver, solnInfo, interPtsAndRays[split_ind], vertexStore, rayStore,
            NBSIC, raysCutForSplit, selectByDist, selectByDepth,
            selectByFinalPoints, secondarySelectByMostRemoved, false, true)) {
          vertexCreatedForSplit[split_ind] = vertCreated;
        }
      } /* end looping over splits */
    } /* end looping over rows */
  } /* checking if should activate basic bounds */

  // Check nonbasic bound
  const int rayVarStructIndex = solnInfo.nonBasicVarIndex[ray_ind];
  if (ACTIVATE_NB_BOUNDS && rayVarStructIndex < solver->getNumCols()) {
    // distToHplane is (UB - LB) / 1.0
    const double LB = getVarLB(solver, rayVarStructIndex);
    const double UB = getVarUB(solver, rayVarStructIndex);
    if (!isNegInfinity(LB, solver->getInfinity())
        && !isInfinity(UB, solver->getInfinity())) {
      tmpHplane.var = rayVarStructIndex;
      tmpHplane.row = static_cast<int>(HplaneRowFlag::NB);;
      tmpHplane.ubflag = static_cast<int>(HplaneBoundFlag::toUBNB);
      tmpHplane.bound =
          true ? UB - LB : UB;

      const double tmpDistToHplane = UB - LB;
      // Set packed vector corresponding to this vertex
      Vertex vertCreated;
      calcNewVectorCoordinates(vertCreated, tmpDistToHplane,
          vertexStore[originatingVertexIndex], rayOut);

//      int numFinalPointsThisHplane = 0;
//      double avgDepthForThisHplane = 0.0;
      std::vector<double> addedDepthForThisHplane(numSplitsToCheck, 0.0);
      //      double avgDepthForThisHplane = 0.0;
      std::vector<int> numFinalPointsThisHplane(numSplitsToCheck, 0);

//      bool hplaneAffectsASplit = false;
//      tmpSplitsAffByHplane.clear();
      for (int split_ind = 0; split_ind < numSplitsToCheck; split_ind++) {
        const double distToSplit =
            (originatingVertexIndex > 0) ?
                distanceToSplit(vertexStore[originatingVertexIndex], rayOut,
                    solnInfo.feasSplitVar[split_ind], solver, solnInfo) :
                distToSplitAlongRay[split_ind][ray_ind];
        // Check split value
        const int split_var = solnInfo.feasSplitVar[split_ind];
        const double split_row = solnInfo.rowOfVar[split_var];
        const double splitRayVal = solnInfo.raysOfC1[ray_ind][split_row];
        const double vertexSplitVal = solnInfo.a0[split_row]
            + splitRayVal * tmpDistToHplane;
        const bool isVertexOnFloor = !greaterThanVal(vertexSplitVal,
            std::floor(solnInfo.a0[split_row]), param.getAWAY());
        const bool isVertexOnCeil = !lessThanVal(vertexSplitVal,
            std::ceil(solnInfo.a0[split_row]), param.getAWAY());
        const bool isVertexOnSplit = isVertexOnFloor || isVertexOnCeil;

        if (chooseHplaneHelper(split_ind, distToSplit, isVertexOnSplit, ray_ind,
            hplane, tmpHplane, tmpDistToHplane, minDistOverall,
            distToHplaneChosenForSplit[split_ind],
            addedDepthByHplaneChosenForSplit[split_ind],
            numFinalPointsByHplaneChosenForSplit[split_ind], numPointsCut,
            vertCreated, addedDepthForThisHplane[split_ind],
            numFinalPointsThisHplane[split_ind], numPointsCut, //allRayIndices,
            solver, solnInfo, interPtsAndRays[split_ind], vertexStore, rayStore,
            NBSIC, raysCutForSplit, selectByDist, selectByDepth,
            selectByFinalPoints, secondarySelectByMostRemoved, true, true)) {
          vertexCreatedForSplit[split_ind] = vertCreated;
        }
      }
    }
  } /* check nb bound */

  if (isInfinity(minDistOverall, solver->getInfinity())) {
//    hplane.structvar = static_cast<int>(HplaneVarFlag::none);
//    hplane.row = static_cast<int>(HplaneRowFlag::none);
//    hplane.ubflag = static_cast<int>(HplaneBoundFlag::none);
    return false;
  }

  // If we did find a hyperplane
  // Recall var >= 0 if found, -1 if not found
  for (int split_ind = 0; split_ind < numSplitsToCheck; split_ind++) {
    // Set packed vector corresponding to this vertex
    if (hplane[split_ind].var >= 0) {
      updateHplaneStore(interPtsAndRays, vertexStore, 0, rayStore,
          allRaysCutAcrossSplits, raysCutForSplit, hplaneStore, solver,
          solnInfo, vertexCreatedForSplit[split_ind], hplane[split_ind],
          split_ind, ray_ind, act_ind);
    }
  }
  return true;
} /* chooseHplaneByRay */

/**
 * Choose a hyperplane to cut a number of rays all at once
 * This is the method used for Algorithm 2 (step 9)
 *
 * The selection criteria are:
 * (1) selectByDist: choose the hplane that cuts any uncut ray first
 * (2) selectByDepth: choose the hplane for which the resulting points have greatest average depth
 * (3) selectByFinalPoints: choose the hplane that yields the most number of final points
 */
bool CglPHA::chooseHplaneBySplit(//const std::vector<int>& allRayIndices,
    const int originatingVertexIndex, const int numRaysToCut,
    const std::vector<int>& sortIndex,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
    std::vector<Hplane>& hplaneStore, std::vector<Vertex>& vertexStore,
    std::vector<Ray>& rayStore, const AdvCuts& NBSIC,
    const std::vector<std::vector<double> >& distToSplitAlongRay,
    std::vector<int>& allRaysCutAcrossSplits,
    std::vector<std::vector<int> >& raysCutForSplit,
    std::vector<IntersectionInfo>& interPtsAndRays, const int act_ind,
    const bool selectByDist, const bool selectByDepth,
    const bool selectByFinalPoints, const bool secondarySelectByMostRemoved) {
  // Initialize values, in case no hplane is found at all
  const int numSplitsToCheck = solnInfo.numFeasSplits; //splitsAffByHplane.size();
  Hplane tmpHplane(HplaneVarFlag::none, HplaneRowFlag::none,
      HplaneBoundFlag::none, numSplitsToCheck); // var, row, ubflag, numsplits
  std::vector<Hplane> hplane(numSplitsToCheck);
  for (int split_ind = 0; split_ind < numSplitsToCheck; split_ind++) {
    hplane[split_ind].resizeVectors(1);
  }

  // For each of the splits we care about (splitsAffByHplane) look at:
  // (1) min distance to the hyperplane
  double minAvgDistOverall = solver->getInfinity();
  std::vector<double> avgDistToHplaneChosenForSplit(numSplitsToCheck, -1.0);
  // (2) the average depth of the points that would be created if we activated each hplane
  std::vector<double> addedDepthByHplaneChosenForSplit(numSplitsToCheck, -1.0);
  // (3) total number of final points that would be created
  std::vector<int> numFinalPointsByHplaneChosenForSplit(numSplitsToCheck, -1);
  // (4) number of points removed by the hyperplane
  std::vector<int> numPointsCutByHplaneChosenForSplit(numSplitsToCheck, -1);

  // Check the basic variable in each row to see if it can be activated
  if (ACTIVATE_BASIC_BOUNDS) {
    for (int row_ind = 0; row_ind < solver->getNumRows(); row_ind++) {
      tmpHplane.var = solnInfo.varBasicInRow[row_ind];
      const bool isStructVar = (tmpHplane.var < solver->getNumCols());
      tmpHplane.row = row_ind;
      if (solver->getRowSense()[row_ind] == 'E') {
        continue; // Do not activate equality rows
      }
      const int numBoundsToCheck = 1 + isStructVar;

      for (int b = 0; b < numBoundsToCheck; b++) {
        std::vector<double> addedDepthForThisHplane(numSplitsToCheck, 0.0);
        std::vector<int> numFinalPointsThisHplane(numSplitsToCheck, 0);
        std::vector<int> numRaysIntersectedThisHplane(numSplitsToCheck, 0);
        std::vector<double> totalDistToThisHplane(numSplitsToCheck, 0.0);
        std::vector<int> numPointsCutByThisHplane(numSplitsToCheck, -1);

        if (b == 0) {
          tmpHplane.bound = getVarLB(solver, tmpHplane.var);
          tmpHplane.ubflag = static_cast<int>(HplaneBoundFlag::toLB);
        } else {
          tmpHplane.bound = getVarUB(solver, tmpHplane.var);
          tmpHplane.ubflag = static_cast<int>(HplaneBoundFlag::toUB);
        }

        // If this bound is actually at infinity; we will never intersect it
        if (isInfinity(std::abs(tmpHplane.bound), solver->getInfinity())) {
          continue;
        }

        if (!isStructVar) {
          tmpHplane.ubflag = static_cast<int>(HplaneBoundFlag::toLBSlack);
        }

        const double varVal = getVarVal(solver, tmpHplane.var);

        for (int tmp_ray_ind = 0; tmp_ray_ind < numRaysToCut; tmp_ray_ind++) {
          const int ray_ind = sortIndex[tmp_ray_ind];
          const double rayVal = solnInfo.raysOfC1[ray_ind][row_ind]; // rayMult * -1 * currRay[r];

          // Check that ray intersects this hyperplane
          if (isZero(rayVal, param.getRAYEPS())) {
            continue;
          }
          if (b == 1 && lessThanVal(rayVal, 0.0, param.getRAYEPS())) {
            continue;
          }
          if (b == 0 && greaterThanVal(rayVal, 0.0, param.getRAYEPS())) {
            continue;
          }

          // Otherwise, find the distance to the hyperplane, and ensure it is actually intersected
          const double tmpDistToHplane = (tmpHplane.bound - varVal) / rayVal;
          if (isInfinity(tmpDistToHplane, solver->getInfinity())) {
            continue; // ray does not intersect this hyperplane
          }

          // Set packed vector corresponding to this vertex
          Vertex vertCreated;
          Ray rayOut;
          const int rayIndices[1] = { ray_ind };
          const double rayValues[1] = { 1.0 };
          rayOut.setVector(1, rayIndices, rayValues);
          calcNewVectorCoordinates(vertCreated, tmpDistToHplane,
              vertexStore[originatingVertexIndex], rayOut);

          for (int split_ind = 0; split_ind < numSplitsToCheck; split_ind++) {
            const int split_var = solnInfo.feasSplitVar[split_ind];
            if (split_var == tmpHplane.var) {
              continue; // Skip if the split var bound is the one being activated
            }

            // If hypeprlane has been activated before, check that it does not already cut the current ray
            const int oldHplaneForSplitOverall = hasHplaneBeenActivated(
                hplaneStore, interPtsAndRays[split_ind].allHplaneToAct,
                tmpHplane, true);
            if (oldHplaneForSplitOverall >= 0) {
              const bool rayCutByCurrHplane =
                  (find_val(ray_ind,
                      hplaneStore[oldHplaneForSplitOverall].rayToBeCutByHplane[split_ind])
                      >= 0);
              if (rayCutByCurrHplane) {
                continue; // Skip if this hplane has already been activated on this ray
              }
            }

            // Check distance to split
            const double distToSplit =
                (originatingVertexIndex > 0) ?
                    distanceToSplit(vertexStore[originatingVertexIndex], rayOut,
                        split_var, solver, solnInfo) :
                    distToSplitAlongRay[split_ind][ray_ind];

            // Check split value
            const double split_row = solnInfo.rowOfVar[split_var];
            const double splitRayVal = solnInfo.raysOfC1[ray_ind][split_row];
            const double vertexSplitVal = solnInfo.a0[split_row]
                + splitRayVal * tmpDistToHplane;
            const bool isVertexOnFloor = !greaterThanVal(vertexSplitVal,
                std::floor(solnInfo.a0[split_row]), param.getAWAY());
            const bool isVertexOnCeil = !lessThanVal(vertexSplitVal,
                std::ceil(solnInfo.a0[split_row]), param.getAWAY());
            const bool isVertexOnSplit = isVertexOnFloor || isVertexOnCeil;

            // Decide whether to use this hplane for this split
            if (chooseHplaneHelper(split_ind, distToSplit, isVertexOnSplit,
                ray_ind, hplane, tmpHplane, tmpDistToHplane, minAvgDistOverall,
                avgDistToHplaneChosenForSplit[split_ind],
                addedDepthByHplaneChosenForSplit[split_ind],
                numFinalPointsByHplaneChosenForSplit[split_ind],
                numPointsCutByHplaneChosenForSplit[split_ind], vertCreated,
                addedDepthForThisHplane[split_ind],
                numFinalPointsThisHplane[split_ind],
                numPointsCutByThisHplane[split_ind], //allRayIndices,
                solver, solnInfo, interPtsAndRays[split_ind], vertexStore,
                rayStore, NBSIC, raysCutForSplit, selectByDist, selectByDepth,
                selectByFinalPoints, secondarySelectByMostRemoved, false, false)) {
              // If this ray is cut by the hplane, count the distance
              numRaysIntersectedThisHplane[split_ind]++;
              totalDistToThisHplane[split_ind] += tmpDistToHplane;
              tmpHplane.rayToBeCutByHplane[split_ind].push_back(ray_ind);
            } /* end hyperplane update for the specific split */
          } /* end iterating over splits to check if chosen hplane affects ray on each split */
        } /* end iterating over rays */

        // Check if this hyperplane is better than the previous one chosen for each split
        for (int split_ind = 0; split_ind < numSplitsToCheck; split_ind++) {
          if ((int) tmpHplane.rayToBeCutByHplane[split_ind].size() == 0) {
            continue;
          }
          const double avgDistToThisHplane = totalDistToThisHplane[split_ind]
              / numRaysIntersectedThisHplane[split_ind];
          const bool hplaneIsCloser = (selectByDist
              && (lessThanVal(avgDistToHplaneChosenForSplit[split_ind], 0.0)
                  || lessThanVal(avgDistToThisHplane, avgDistToHplaneChosenForSplit[split_ind])));
          const bool higherDepth = (selectByDepth
              && greaterThanVal(addedDepthForThisHplane[split_ind],
                  addedDepthByHplaneChosenForSplit[split_ind]));
          const bool moreFinalPoints = (selectByFinalPoints
              && (numFinalPointsThisHplane[split_ind]
                  > numFinalPointsByHplaneChosenForSplit[split_ind]));
          const bool cutMorePoints = (secondarySelectByMostRemoved
              && (numPointsCutByThisHplane[split_ind]
                  > numPointsCutByHplaneChosenForSplit[split_ind]));
          const bool tieBreaking = tieBreakingForChooseHplane(selectByDist,
              selectByDepth, selectByFinalPoints,
              avgDistToHplaneChosenForSplit[split_ind],
              addedDepthByHplaneChosenForSplit[split_ind],
              numFinalPointsByHplaneChosenForSplit[split_ind],
              avgDistToThisHplane, addedDepthForThisHplane[split_ind],
              numFinalPointsThisHplane[split_ind]);

          if ((hplaneIsCloser || higherDepth || moreFinalPoints || cutMorePoints
              || tieBreaking)) {
#ifdef TRACE
            char temp[500];
            if (hplaneIsCloser) {
              snprintf(temp, sizeof(temp) / sizeof(char), "Closer, %f < %f.",
                  avgDistToThisHplane,
                  avgDistToHplaneChosenForSplit[split_ind]);
            } else if (higherDepth) {
              snprintf(temp, sizeof(temp) / sizeof(char),
                  "Higher depth, %f > %f.", addedDepthForThisHplane[split_ind],
                  addedDepthByHplaneChosenForSplit[split_ind]);
            } else if (moreFinalPoints) {
              snprintf(temp, sizeof(temp) / sizeof(char),
                  "More final points, %d > %d.",
                  numFinalPointsThisHplane[split_ind],
                  numFinalPointsByHplaneChosenForSplit[split_ind]);
            } else if (cutMorePoints) {
              snprintf(temp, sizeof(temp) / sizeof(char),
                  "More cut points, %d > %d.",
                  numPointsCutByThisHplane[split_ind], numPointsCutByHplaneChosenForSplit[split_ind]);
            } else if (tieBreaking) {
              snprintf(temp, sizeof(temp) / sizeof(char),
                  "Tiebreak win. Dist: %f (vs %f). Depth: %f (vs %f). Final points: %d (vs %d). Cut points: %d (vs %d).",
                  avgDistToThisHplane, avgDistToHplaneChosenForSplit[split_ind],
                  addedDepthForThisHplane[split_ind],
                  addedDepthByHplaneChosenForSplit[split_ind],
                  numFinalPointsThisHplane[split_ind],
                  numFinalPointsByHplaneChosenForSplit[split_ind],
                  numPointsCutByThisHplane[split_ind],
                  numPointsCutByHplaneChosenForSplit[split_ind]);
            }
            printf("\tSplit %d. Found better hyperplane. %s\n", split_ind, temp);
            printf("\tOld hyperplane: ");
            hplane[split_ind].printAll();
            printf("\tNew hyperplane: ");
            tmpHplane.printAll();
#endif
            if (lessThanVal(avgDistToThisHplane, minAvgDistOverall)) {
              minAvgDistOverall = avgDistToThisHplane;
            }
            avgDistToHplaneChosenForSplit[split_ind] = avgDistToThisHplane;
            addedDepthByHplaneChosenForSplit[split_ind] =
                addedDepthForThisHplane[split_ind];
            numFinalPointsByHplaneChosenForSplit[split_ind] =
                numFinalPointsThisHplane[split_ind];
            numPointsCutByHplaneChosenForSplit[split_ind] =
                numPointsCutByThisHplane[split_ind];

            hplane[split_ind].var = tmpHplane.var;
            hplane[split_ind].row = tmpHplane.row;
            hplane[split_ind].ubflag = tmpHplane.ubflag;
            hplane[split_ind].bound = tmpHplane.bound;
            if (!(tmpHplane.var < solver->getNumCols())) {
              hplane[split_ind].ubflag =
                  static_cast<int>(HplaneBoundFlag::toLBSlack);
            } else {
              // If this variable was originally NB and at its upper-bound,
              // then we had complemented it to be at its lower-bound
              // Hence, if it is at its upper-bound now, this is the same as its lower-bound
              // in the complemented nonbasic space
              if (solnInfo.isNBUBVar(tmpHplane.var)) {
                hplane[split_ind].ubflag *= -1 * tmpHplane.ubflag + 1; // Switch between 0 and 1
              } // TODO Double check this
            }
            hplane[split_ind].rayToBeCutByHplane[0] =
                tmpHplane.rayToBeCutByHplane[split_ind];
          } /* end hyperplane update for the specific split */
          tmpHplane.rayToBeCutByHplane[split_ind].clear();
        } /* end looping over splits to update hyperplanes */
      } /* end looping over bounds to check */
    } /* end looping over rows */
  } /* checking if should activate basic bounds */

  // Check nonbasic bounds
  if (ACTIVATE_NB_BOUNDS) {
    for (int tmp_ray_ind = 0; tmp_ray_ind < numRaysToCut; tmp_ray_ind++) {
      const int ray_ind = sortIndex[tmp_ray_ind];
      const int rayVarStructIndex = solnInfo.nonBasicVarIndex[ray_ind];
      if (rayVarStructIndex < solver->getNumCols()) {
        Ray rayOut;
        const int rayIndices[1] = { ray_ind };
        const double rayValues[1] = { 1.0 };
        rayOut.setVector(1, rayIndices, rayValues);

        // distToHplane is (UB - LB) / 1.0
        const double LB = getVarLB(solver, rayVarStructIndex);
        const double UB = getVarUB(solver, rayVarStructIndex);
        if (!isNegInfinity(LB, solver->getInfinity())
            && !isInfinity(UB, solver->getInfinity())) {
          tmpHplane.var = rayVarStructIndex;
          tmpHplane.row = static_cast<int>(HplaneRowFlag::NB);
          tmpHplane.ubflag = static_cast<int>(HplaneBoundFlag::toUBNB);
          tmpHplane.bound =
//              param.getParamVal(ParamIndices::NB_SPACE_PARAM_IND)
              true ?
                  UB - LB : UB;

          const double tmpDistToHplane = UB - LB;
          // Set packed vector corresponding to this vertex
          Vertex vertCreated;
          calcNewVectorCoordinates(vertCreated, tmpDistToHplane,
              vertexStore[originatingVertexIndex], rayOut);

          //      int numFinalPointsThisHplane = 0;
          //      double avgDepthForThisHplane = 0.0;
          std::vector<double> addedDepthForThisHplane(numSplitsToCheck, 0.0);
          //      double avgDepthForThisHplane = 0.0;
          std::vector<int> numFinalPointsThisHplane(numSplitsToCheck, 0);
          std::vector<int> numPointsCutThisHplane(numSplitsToCheck, 0);

          //      bool hplaneAffectsASplit = false;
          //      tmpSplitsAffByHplane.clear();
          for (int split_ind = 0; split_ind < numSplitsToCheck; split_ind++) {
            const int oldHplaneForSplitOverall = hasHplaneBeenActivated(
                hplaneStore, interPtsAndRays[split_ind].allHplaneToAct,
                tmpHplane, true);
            if (oldHplaneForSplitOverall >= 0) {
              const bool rayCutByCurrHplane =
                  (find_val(ray_ind,
                      hplaneStore[oldHplaneForSplitOverall].rayToBeCutByHplane[split_ind])
                      >= 0);
              if (rayCutByCurrHplane) {
                continue; // Skip if this hplane has already been activated on this ray
              }
            }
            const double distToSplit =
                (originatingVertexIndex > 0) ?
                    distanceToSplit(vertexStore[originatingVertexIndex], rayOut,
                        solnInfo.feasSplitVar[split_ind], solver, solnInfo) :
                    distToSplitAlongRay[split_ind][ray_ind];

            // Check split value
            const int split_var = solnInfo.feasSplitVar[split_ind];
            const double split_row = solnInfo.rowOfVar[split_var];
            const double splitRayVal = solnInfo.raysOfC1[ray_ind][split_row];
            const double vertexSplitVal = solnInfo.a0[split_row]
                + splitRayVal * tmpDistToHplane;
            const bool isVertexOnFloor = !greaterThanVal(vertexSplitVal,
                std::floor(solnInfo.a0[split_row]), param.getAWAY());
            const bool isVertexOnCeil = !lessThanVal(vertexSplitVal,
                std::ceil(solnInfo.a0[split_row]), param.getAWAY());
            const bool isVertexOnSplit = isVertexOnFloor || isVertexOnCeil;

            if (chooseHplaneHelper(split_ind, distToSplit, isVertexOnSplit,
                ray_ind, hplane, tmpHplane, tmpDistToHplane, minAvgDistOverall,
                avgDistToHplaneChosenForSplit[split_ind],
                addedDepthByHplaneChosenForSplit[split_ind],
                numFinalPointsByHplaneChosenForSplit[split_ind],
                numPointsCutByHplaneChosenForSplit[split_ind], vertCreated,
                addedDepthForThisHplane[split_ind],
                numFinalPointsThisHplane[split_ind],
                numPointsCutThisHplane[split_ind], //allRayIndices,
                solver, solnInfo, interPtsAndRays[split_ind], vertexStore,
                rayStore, NBSIC, raysCutForSplit, selectByDist, selectByDepth,
                selectByFinalPoints, secondarySelectByMostRemoved, true, true)) {
              // Update rays cut; other items updated within chooseHplaneHelper since selecting by "ray"
              hplane[split_ind].rayToBeCutByHplane[0].assign(rayIndices,
                  rayIndices + 1);
            } /* end hyperplane update for the specific split */
          } /* end looping over splits */
        } /* ensure ray is bounded from both sides */
      } /* ensure the ray is not corresponding to a slack variable */
    } /* looping over rays to activate nonbasic bounds */
  } /* checking if should activate nonbasic bounds */

  if (isInfinity(minAvgDistOverall, solver->getInfinity())) { // did we not find a hplane?
    return false;
  }

  // If we did find a hyperplane
  // Recall var >= 0 if found, -1 if not found
  for (int split_ind = 0; split_ind < numSplitsToCheck; split_ind++) {
    // Set packed vector corresponding to this vertex
    if (hplane[split_ind].var >= 0) {
      const int numRaysToBeCutByThisHplane =
          hplane[split_ind].rayToBeCutByHplane[0].size();
      for (int r = 0; r < numRaysToBeCutByThisHplane; r++) {
        const int ray_ind = hplane[split_ind].rayToBeCutByHplane[0][r];
        const int row_ind = hplane[split_ind].row;
        const double tmpDistToHplane = (hplane[split_ind].bound
            - getVarVal(solver, hplane[split_ind].var))
            / ((row_ind < 0) ? 1. : solnInfo.raysOfC1[ray_ind][row_ind]);
        Vertex vertCreated;
        Ray rayOut;
        const int rayIndices[1] = { ray_ind };
        const double rayValues[1] = { 1.0 };
        rayOut.setVector(1, rayIndices, rayValues);
        calcNewVectorCoordinates(vertCreated, tmpDistToHplane,
            vertexStore[originatingVertexIndex], rayOut);
        updateHplaneStore(interPtsAndRays, vertexStore, 0, rayStore,
            allRaysCutAcrossSplits, raysCutForSplit, hplaneStore, solver,
            solnInfo, vertCreated, hplane[split_ind], split_ind,
            ray_ind, act_ind, true);
      } /* loop over rays cut by the hyperplane for this split */
    }
  }
  return true;
} /* chooseHplaneBySplit */

/**
 * chooseHplaneHelper
 *
 * If "chooseByRay" is selected, then the decision of whether to replace the current
 * hplane cutting the ray is made within this method;
 * otherwise, the method computes the relevant stats,
 * and the decision is made within the chooseHplaneBySplit method
 */
bool CglPHA::chooseHplaneHelper(const int split_ind, const double distToSplit,
    const bool isVertexOnSplit, const int ray_ind, std::vector<Hplane>& hplane,
    const Hplane& tmpHplane, const double tmpDistToHplane,
    double& minDistOverall, double& distToHplaneChosenForSplit,
    double& addedDepthByHplaneChosenForSplit,
    int& numFinalPointsByHplaneChosenForSplit,
    int& numPointsCutByHplaneChosenForSplit,
    const Vertex& vertCreated, double& addedDepthForThisHplane,
    int& numFinalPointsThisHplane, int& numPointsCutThisHplane,
    /*const std::vector<int>& allRayIndices,*/
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
    const IntersectionInfo& interPtsAndRays, const std::vector<Vertex>& vertexStore,
    const std::vector<Ray>& rayStore, const AdvCuts& NBSIC,
    const std::vector<std::vector<int> >& raysCutForSplit, const bool selectByDist,
    const bool selectByDepth, const bool selectByFinalPoints,
    const bool calcPointsRemoved, const bool isNBBound,
    const bool chooseByRay) {
  if (solnInfo.feasSplitVar[split_ind] == tmpHplane.var) {
    return false; // Skip if the split var bound is the one being activated
  }

  if (!CUT_PARALLEL_RAYS && isInfinity(distToSplit, solver->getInfinity())) {
    return false;
  }

  // Check whether the point hits the split first
  if (!lessThanVal(tmpDistToHplane, distToSplit) || isVertexOnSplit) {
    return false;
  }

    // The hplane is intersected before bd S
  IntersectionInfo* tmpInterPtsAndRays;
  if (!selectByDist) { // Currently skip this for selectByDist
    tmpInterPtsAndRays = new IntersectionInfo;
    tmpInterPtsAndRays->reverseOrdering();
    tmpInterPtsAndRays->setDimensions(0, solnInfo.numNB);
    fakeActivateHplane(numFinalPointsThisHplane, selectByFinalPoints,
        addedDepthForThisHplane, selectByDepth, /*tmpDistToHplane,*/ solver,
        solnInfo, split_ind, vertexStore, vertCreated, rayStore, ray_ind,
        tmpHplane, (NBSIC.size() > 0) ? &(NBSIC.cuts[split_ind]) : NULL, /*allRayIndices,*/
        *tmpInterPtsAndRays, interPtsAndRays, raysCutForSplit[split_ind]);
  }

  if (calcPointsRemoved) {
    std::vector<int> tmpPointsToRemove;
    pointsFromRayViolatingHplane(numPointsCutByHplaneChosenForSplit,
        tmpPointsToRemove, ray_ind, /*vertexStore, rayStore,*/ tmpHplane,
        interPtsAndRays, /*split_ind,*/ solver, solnInfo);
  }

  if (!selectByDist) {
    if (tmpInterPtsAndRays) {
      delete tmpInterPtsAndRays;
    }
  }

  // Exit early if we are not coming from chooseHplaneByRay (from algorithm 3)
  // This is because we need to activate the hyperplane on all the rays to get the full picture
  if (!chooseByRay) {
    return true;
  }

  // Check various selection criteria
  const bool hplaneIsCloser = (selectByDist
      && (lessThanVal(distToHplaneChosenForSplit, 0.0)
          || lessThanVal(tmpDistToHplane, distToHplaneChosenForSplit)));
  const bool higherDepth = (selectByDepth
      && greaterThanVal(addedDepthForThisHplane,
          addedDepthByHplaneChosenForSplit));
  const bool moreFinalPoints = (selectByFinalPoints
      && (numFinalPointsThisHplane > numFinalPointsByHplaneChosenForSplit));
  const bool cutMorePoints = (calcPointsRemoved
      && (numPointsCutThisHplane > numPointsCutByHplaneChosenForSplit));

  // Check tie breaking rules
  const bool tieBreaking = tieBreakingForChooseHplane(selectByDist,
      selectByDepth, selectByFinalPoints, distToHplaneChosenForSplit,
      addedDepthByHplaneChosenForSplit, numFinalPointsByHplaneChosenForSplit,
      tmpDistToHplane, addedDepthForThisHplane, numFinalPointsThisHplane);

  if ((hplaneIsCloser || higherDepth || moreFinalPoints || cutMorePoints
      || tieBreaking)) {
#ifdef TRACE
    char temp[500];
    if (hplaneIsCloser) {
      snprintf(temp, sizeof(temp) / sizeof(char), "Closer, %f < %f.",
          tmpDistToHplane, distToHplaneChosenForSplit);
    } else if (higherDepth) {
      snprintf(temp, sizeof(temp) / sizeof(char), "Higher depth, %f > %f.",
          addedDepthForThisHplane, addedDepthByHplaneChosenForSplit);
    } else if (moreFinalPoints) {
      snprintf(temp, sizeof(temp) / sizeof(char), "More final points, %d > %d.",
          numFinalPointsThisHplane, numFinalPointsByHplaneChosenForSplit);
    } else if (cutMorePoints) {
      snprintf(temp, sizeof(temp) / sizeof(char), "More cut points, %d > %d.",
          numPointsCutThisHplane, numPointsCutByHplaneChosenForSplit);
    } else if (tieBreaking) {
      snprintf(temp, sizeof(temp) / sizeof(char),
          "Tiebreak win. Dist: %f (vs %f). Depth: %f (vs %f). Final points: %d (vs %d). Cut points: %d (vs %d).",
          tmpDistToHplane, distToHplaneChosenForSplit, addedDepthForThisHplane,
          addedDepthByHplaneChosenForSplit, numFinalPointsThisHplane,
          numFinalPointsByHplaneChosenForSplit, numPointsCutThisHplane,
          numPointsCutByHplaneChosenForSplit);
    }
    printf("\tSplit %d. Found better hyperplane. %s\n", split_ind, temp);
    printf("\tOld hyperplane: ");
    hplane[split_ind].printAll();
    printf("\tNew hyperplane: ");
    tmpHplane.printAll();
#endif
    if (lessThanVal(tmpDistToHplane, minDistOverall)) {
      minDistOverall = tmpDistToHplane;
    }
    distToHplaneChosenForSplit = tmpDistToHplane;
    addedDepthByHplaneChosenForSplit = addedDepthForThisHplane;
    numFinalPointsByHplaneChosenForSplit = numFinalPointsThisHplane;
    numPointsCutByHplaneChosenForSplit = numPointsCutThisHplane;

    hplane[split_ind].var = tmpHplane.var;
    hplane[split_ind].row = tmpHplane.row;
    hplane[split_ind].ubflag = tmpHplane.ubflag;
    hplane[split_ind].bound = tmpHplane.bound;
    if (!isNBBound) {
      if (!(tmpHplane.var < solver->getNumCols())) {
        hplane[split_ind].ubflag = static_cast<int>(HplaneBoundFlag::toLBSlack);
      } else {
        // If this variable was originally NB and at its upper-bound,
        // then we had complemented it to be at its lower-bound
        // Hence, if it is at its upper-bound now, this is the same as its lower-bound
        // in the complemented nonbasic space
        if ( //param.getParamVal(ParamIndices::NB_SPACE_PARAM_IND) &&
            solnInfo.isNBUBVar(tmpHplane.var)) {
          hplane[split_ind].ubflag *= -1 * tmpHplane.ubflag + 1; // Switch between 0 and 1
        } // TODO Double check this
      }
    }
    return true;
  } else {
    return false;
  }
} /* chooseHplaneHelper */

void CglPHA::fakeActivateHplane(int& numFinalPointsAdded,
    const bool calcFinalPoints, double& avgAddedDepth, const bool calcDepth,
    /*const double distToHplane,*/ const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const int split_ind,
    const std::vector<Vertex>& vertexStore, const Vertex& vertCreated,
    const std::vector<Ray>& rayStore, const int ray_ind, const Hplane& hplane,
    const AdvCut* NBSIC, /*const std::vector<int>& allRayIndices,*/
    IntersectionInfo& interPtsAndRays,
    const IntersectionInfo& oldInterPtsAndRays,
    const std::vector<int>& raysPrevCut) {
//  int totalNumPointsAdded = 0;
//  int totalNumPointsRemoved = 0;
//  double avgRemovedDepth = 0.0;

  const std::vector<int> raysToBeCutByHplane(1, ray_ind);
  std::vector<int> vertOfAllRaysIntByHplane, allRaysIntByHplane;

  // Calculate totalSOS for the SIC for this split
  CoinPackedVector SICVec;
  if (calcDepth && NBSIC) {
    SICVec = NBSIC->row();
  }
  const double norm = (calcDepth) ? std::sqrt(SICVec.normSquare()) : 0;

  // Get all rays intersected by this hplane
  addRaysIntByHplaneBeforeBdS(vertOfAllRaysIntByHplane, allRaysIntByHplane,
      solver, solnInfo, hplane, vertexStore, 0, rayStore, split_ind);

  // We we want to track the average depth of the removed intersection points
//  double sumRemovedDepth = 0.0;
  double sumAddedDepth = 0.0;
//  int numNewPointsRemoved = 0;
  int numNewPointsAdded = 0;

  // Rays cut: ray_ind + all rays from allRaysIntByHplane that have previously been cut
  std::vector<bool> raysCutFlag(solnInfo.numNB, false);
  std::vector<bool> rayIntersectedFlag(solnInfo.numNB, false);

  raysCutFlag[ray_ind] = true;
  const int numRaysIntByHplane = allRaysIntByHplane.size();
  for (int r = 0; r < numRaysIntByHplane; r++) {
    rayIntersectedFlag[allRaysIntByHplane[r]] = true;
  }
  for (int r = 0; r < (int) raysPrevCut.size(); r++) {
    if (rayIntersectedFlag[raysPrevCut[r]])
      raysCutFlag[raysPrevCut[r]] = true;
  }

  // Remove all intersection points cut off by this hyperplane
  const int split_var = solnInfo.feasSplitVar[split_ind];
//  if (calcDepth) {
////    numNewPointsRemoved +=
//    removePtsAndRaysCutByHplane(solver, solnInfo, raysCutFlag, split_var,
//        hplane, vertexStore, rayStore, interPtsAndRays, &NBSIC, norm,
//        sumRemovedDepth, true);
//  }

  // Calculate average depth of the removed points
//  if (numNewPointsRemoved > 0) {
//    avgRemovedDepth += sumRemovedDepth / numNewPointsRemoved;
//    totalNumPointsRemoved += numNewPointsRemoved;
//  }

  /***********************************************************************************
   * Alg 1, step 4 (find indices of new edges emanating from activated hyperplane)
   ***********************************************************************************/
  // Create the set J \setminus K
  std::vector<int> newRayColIndices;
  storeNewRayNBIndices(solnInfo.numNB, newRayColIndices, allRaysIntByHplane,
      raysCutFlag);
  const int numNewRays = newRayColIndices.size();

  //For ray intersected by hyperplane H generate new rays corresponding to each element in J \setminus K
  //This is of course available in closed form from the tableau
  for (int r = 0; r < numNewRays; r++) {
    Ray newRay;
    Point tmpPoint;
    calcNewRayCoordinatesRank1(newRay, solnInfo, ray_ind, newRayColIndices[r],
        hplane, CoinMin(param.getRAYEPS(), solnInfo.EPS));
    bool ptFlagToAdd = calcIntersectionPointWithSplit(tmpPoint, solver,
        solnInfo, vertCreated, newRay, split_var);

    if (ptFlagToAdd) {
      if (!duplicatePoint(interPtsAndRays, interPtsAndRays.RHS, tmpPoint)) {
        numNewPointsAdded++;
        interPtsAndRays.addPointOrRayToIntersectionInfo(tmpPoint, 1.0, ray_ind,
            newRayColIndices[r], -1, -1, -1, -1, -1., -1., false);

        // Check whether it is a final point
        if (calcFinalPoints && pointInP(solver, solnInfo, tmpPoint)) {
          if (!duplicatePoint(oldInterPtsAndRays, oldInterPtsAndRays.RHS,
              tmpPoint)) {
            numFinalPointsAdded++;
          }
        }

        if (calcDepth) {
          // Calculate depth of added point
          const double SIC_activity = getSICActivity(split_var, vertexStore,
              rayStore, solver, solnInfo, tmpPoint);
          double curr_SIC_depth =
              (NBSIC) ? (SIC_activity - NBSIC->rhs()) / norm : 0.0;
          if (lessThanVal(curr_SIC_depth, 0.0, param.getDIFFEPS())) {
            error_msg(errstr,
                "Depth %e is negative! Point cannot be above the SIC.\n",
                curr_SIC_depth);
            writeErrorToII(errstr, GlobalVariables::log_file);
            exit(1);
          } else if (lessThanVal(curr_SIC_depth, 0.0, param.getEPS())) {
            warning_msg(errstr,
                "Depth %e is negative! Point cannot be above the SIC. Small enough violation that only sending warning.\n",
                curr_SIC_depth);
          }
          if (isZero(curr_SIC_depth, param.getDIFFEPS())) {
            curr_SIC_depth = +0.0;
          } 

          sumAddedDepth += curr_SIC_depth;
        }
      } /* check that not duplicate point */
    } /* if new ray intersects split bd, it creates a new point */
  } /* end iterating over new rays to be created */
  if (numNewPointsAdded > 0) {
    avgAddedDepth += sumAddedDepth / numNewPointsAdded;
  }
} /* fakeActivateHplane */

/**
 * Finds points created from rayInd that are cut by this hyperplane
 */
void CglPHA::pointsFromRayViolatingHplane(int& numPointsCutByHplane,
    std::vector<int>& pointsToRemove,
    const int rayInd,
//    const std::vector<Vertex>& vertexStore, const std::vector<Ray>& rayStore,
    const Hplane& hplane, const IntersectionInfo& interPtsAndRays,
//    const int splitsAffByHplane,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo) {
  pointsToRemove.clear();
  int numViolPoints = 0;
  for (int row_ind = 0; row_ind < interPtsAndRays.getNumRows(); row_ind++) {
    const int ray_ind =
        (interPtsAndRays.cutRay[row_ind] >= 0) ?
            interPtsAndRays.cutRay[row_ind] :
            interPtsAndRays.newRay[row_ind];
    if (ray_ind != rayInd) {
      continue;
    }

    // Ensure this is a point, not a ray
    if (isZero(interPtsAndRays.RHS[row_ind])) {
      continue;
    }

    const CoinShallowPackedVector tmpVector = interPtsAndRays.getVector(
        row_ind);
    const int numElmts = tmpVector.getNumElements();
    const int* elmtIndices = tmpVector.getIndices();
    const double* elmtVals = tmpVector.getElements();

    const bool NBFlag = (hplane.row == static_cast<int>(HplaneRowFlag::NB));

    // The NBFlag prevents us from using a -1 row index
    const double xVal = !NBFlag * solnInfo.a0[!NBFlag * hplane.row];
    double NBRowActivity = 0.0;
    for (int k = 0; k < numElmts; k++) {
      if (!NBFlag) {
        NBRowActivity += solnInfo.raysOfC1[elmtIndices[k]][hplane.row]
            * elmtVals[k];
      } else if (elmtIndices[k] == hplane.var) {
        NBRowActivity += elmtVals[k];
      }
    }
    const double newVal = xVal + NBRowActivity;
    double colLB, colUB;

    if (hplane.var < solver->getNumCols()) {
      if (!NBFlag) {
        colLB = getVarLB(solver, hplane.var);
        colUB = getVarUB(solver, hplane.var);
      } else {
        colLB = 0.0;
        colUB = getVarUB(solver, hplane.var) - getVarLB(solver, hplane.var);
      }
    } else {
      colLB = 0.0;
      colUB = solver->getInfinity();
    }

    if (lessThanVal(newVal, colLB) || greaterThanVal(newVal, colUB)) {
      numViolPoints++;
      pointsToRemove.push_back(row_ind);
      if (numPointsCutByHplane == -1) {
        numPointsCutByHplane++;
      }
      numPointsCutByHplane++;
    }
  }
} /* pointsFromRayViolatingHplane */

/***********************************************************************/
bool CglPHA::addInterPtFromUncutRayOfC1(IntersectionInfo& interPtsAndRays,
    double& distToSplitAlongRay, const SolutionInfo& solnInfo,
    const int split_ind, const std::vector<Ray>& rayStore, const int ray_ind,
    const std::vector<Vertex>& vertexStore, const bool recalc_dist) {
  const char* rowSense = solver->getRowSense();

  if (rayStore[ray_ind].cutFlag[split_ind]) {
    return false;
  }

  const int nb_var = solnInfo.nonBasicVarIndex[ray_ind];
  if ((nb_var >= solnInfo.numCols)
      && (rowSense[nb_var - solnInfo.numCols] == 'E')) {
    return false;
  }

  Point tmpPoint;
  bool rayIntSplitFlag = calcIntersectionPointWithSplit(tmpPoint,
      distToSplitAlongRay, castsolver, solnInfo, vertexStore[0],
      rayStore[ray_ind], solnInfo.feasSplitVar[split_ind], recalc_dist);
  if (!rayIntSplitFlag) {
//    // Make sure this is not a duplicate
//    if (!CHECK_FOR_DUP_RAYS
//        || !duplicateRay(interPtsAndRays, interPtsAndRays.RHS,
//            rayStore[ray_ind])) {
    const bool currRayIsFinal = isRayFinal(castsolver, solnInfo, ray_ind);
    interPtsAndRays.addPointOrRayToIntersectionInfo(rayStore[ray_ind], 0.0, -1,
        ray_ind, -1, -1, 0, -1, 0.0, 0.0, currRayIsFinal);
//    }
  } else {
//    if (!duplicatePoint(interPtsAndRays, interPtsAndRays.RHS, tmpPoint)) {
    // TODO Check if the ray value is equal to the distance to split; should be the same
    const bool currPointIsFinal = pointInP(castsolver, solnInfo, tmpPoint);
    interPtsAndRays.addPointOrRayToIntersectionInfo(tmpPoint, 1.0, -1, ray_ind,
        -1, -1, 0, -1, 0.0, 0.0, currPointIsFinal);
//    }
  }
  return true;
} /* addInterPtFromUncutRayOfC1 */

/***********************************************************************/
void CglPHA::performActivations(const SolutionInfo& solnInfo,
    const std::vector<Vertex>& vertexStore, std::vector<Ray>& rayStore,
    const std::vector<Hplane>& hplaneStore,
    std::vector<IntersectionInfo>& interPtsAndRays, const AdvCuts& NBSIC,
    const AdvCut& objCut, const double objNorm, const int act_ind) {
  for (int split_ind = 0; split_ind < solnInfo.numFeasSplits; split_ind++) {
    const int numHplanesActForSplit =
        interPtsAndRays[split_ind].hplaneIndexByAct[act_ind].size();

    // Calculate norm for the SIC for this split
    const AdvCut* currNBSIC = (NBSIC.size() > 0) ? NBSIC.getCutPointer(split_ind) : NULL;
    const double NBSICNorm = (currNBSIC) ? (currNBSIC->row()).twoNorm() : 0.0;

    for (int h = 0; h < numHplanesActForSplit; h++) {
      const int hplane_ind =
          interPtsAndRays[split_ind].hplaneIndexByAct[act_ind][h];
      const int hplane_ind_in_split =
          interPtsAndRays[split_ind].indexOfHplaneInOverallOrderForSplit[act_ind][h];
#ifdef TRACE
      printf(
          "Act %d, split %d (var %d), hplane %d/%d: hplane var %d, hplane row index %d.\n",
          act_ind + 1, split_ind, solnInfo.feasSplitVar[split_ind], h + 1,
          numHplanesActForSplit, hplaneStore[hplane_ind].var,
          hplaneStore[hplane_ind].row);
#endif
      algorithm1(interPtsAndRays[split_ind].numPointsRemoved[act_ind],
          interPtsAndRays[split_ind].avgRemovedDepth[act_ind],
          interPtsAndRays[split_ind].numPointsAdded[act_ind],
          interPtsAndRays[split_ind].avgAddedDepth[act_ind],
          interPtsAndRays[split_ind], split_ind, hplane_ind,
          hplane_ind_in_split, hplaneStore, hplaneStore[hplane_ind], rayStore,
          vertexStore, castsolver, solnInfo, currNBSIC, NBSICNorm, objCut,
          objNorm, act_ind);
    } /* activate hyperplanes for this split */

    // Average average depth
    if (numHplanesActForSplit > 0) {
      interPtsAndRays[split_ind].avgRemovedDepth[act_ind] /=
          numHplanesActForSplit;
      interPtsAndRays[split_ind].avgAddedDepth[act_ind] /=
          numHplanesActForSplit;
    }
  } /* iterate over splits */
} /* performActivations */

/***********************************************************************/
void CglPHA::tryObjectives(AdvCuts& structPHA,
    std::vector<int>& num_cuts_per_cgs, int& num_total_gics, int& num_obj_tried,
    const SolutionInfo& solnInfo, //const int dim,
    const std::vector<IntersectionInfo>& interPtsAndRays,
//    const std::vector<Ray>& rayStore,
    const std::vector<Vertex>& vertexStore,
    const std::vector<Hplane>& hplaneStore, const AdvCuts& structSICs,
    Stats& phaTimeStats) {
  GlobalVariables::timeStats.start_timer(GEN_MPHA_TIME);

  // Set up solvers with points and rays found
  std::vector<PHASolverInterface*> cutSolver(solnInfo.numFeasSplits);

  // This is in the * complemented * nonbasic space,
  // so the origin is the solution to the LP relaxation
  // Moreover, every cut will have right-hand side 1,
  // since we want to separate the origin
  
  // Set up scaling
//  const bool useScale = false && !isInfinity(std::abs(this->min_nb_obj_val)) && !lessThanVal(this->min_nb_obj_val, 1.);
//  const double beta = (useScale) ? this->min_nb_obj_val : 1.0;
#ifdef TRACE
//  printf("Beta for cut generation set to %.3f (min_nb_obj_val = %.3f).\n", beta, this->min_nb_obj_val);
#endif
  std::vector<bool> isCutSolverPrimalFeasible(solnInfo.numFeasSplits);
  for (int s = 0; s < solnInfo.numFeasSplits; s++) {
    if (reachedCutLimit(num_cuts_per_cgs[s], num_total_gics)) {
      isCutSolverPrimalFeasible[s] = false;
      continue;
    }

    cutSolver[s] = new PHASolverInterface;
    isCutSolverPrimalFeasible[s] = setupCutSolver(cutSolver[s],
        interPtsAndRays[s], s,
        solver->getColName(solnInfo.feasSplitVar[s]));
  }

  bool isTimeLimitReached = reachedTimeLimit(phaTimeStats, "PHA_TOTAL", param.getTIMELIMIT());
  if (param.getParamVal(ParamIndices::USE_SPLIT_SHARE_PARAM_IND) > 0
      && !isTimeLimitReached) {
    // Try split sharing idea
#ifdef TRACE
    printf("\n## Try split sharing idea. ##\n");
#endif
    splitSharingCutGeneration(cutSolver, structPHA, num_cuts_per_cgs,
        num_total_gics, num_obj_tried, interPtsAndRays, castsolver, solnInfo,
        vertexStore, hplaneStore, structSICs,
        isCutSolverPrimalFeasible);
  }

  for (int s = 0; s < solnInfo.numFeasSplits; s++) {
    if (reachedTimeLimit(phaTimeStats, "PHA_TOTAL", param.getTIMELIMIT())) {
      break;
    }
    if (reachedCutLimit(num_cuts_per_cgs[s], num_total_gics)) {
      isCutSolverPrimalFeasible[s] = false;
      continue;
    }
    if (isCutSolverPrimalFeasible[s]) {
      tryObjectivesForCgs(cutSolver[s], structPHA, num_cuts_per_cgs[s],
          num_total_gics, num_obj_tried, solnInfo, //solver->getNumCols(),
          s, solver->getColName(solnInfo.feasSplitVar[s]), //interPtsAndRays[s],
          vertexStore, structSICs);
    }
  }

  isTimeLimitReached = reachedTimeLimit(phaTimeStats, "PHA_TOTAL", param.getTIMELIMIT());
  if (param.getParamVal(ParamIndices::USE_SPLIT_SHARE_PARAM_IND) < 0
      && !isTimeLimitReached) {
    // Try split sharing idea
#ifdef TRACE
    printf("\n## Try split sharing idea. ##\n");
#endif
    splitSharingCutGeneration(cutSolver, structPHA, num_cuts_per_cgs,
        num_total_gics, num_obj_tried, interPtsAndRays, castsolver, solnInfo,
        vertexStore, hplaneStore, structSICs,
        isCutSolverPrimalFeasible);
  }

  // Clear things
  for (int s = 0; s < solnInfo.numFeasSplits; s++) {
    if (cutSolver[s]) {
      delete cutSolver[s];
    }
  }
  GlobalVariables::timeStats.end_timer(GEN_MPHA_TIME);
} /* tryObjectives */

/***********************************************************************/
void CglPHA::tryObjectivesForCgs(PHASolverInterface* cutSolver,
    AdvCuts& structPHA, int& num_cuts_per_cgs, int& num_cuts_generated,
    int& num_obj_tried, const SolutionInfo& probData,
//    const int dim, const bool inNBSpace,
    const int cgs_ind, const std::string& cgsName,
//    const IntersectionInfo& interPtsAndRays,
    const std::vector<Vertex>& vertexStore, const AdvCuts& structSICs) {
  const double beta = 1.0;
  tightOnPointsRaysCutGeneration(cutSolver, beta, structPHA, num_cuts_per_cgs,
      num_cuts_generated, num_obj_tried,
      dynamic_cast<PHASolverInterface*>(solver), probData, cgs_ind,
      cgsName, structSICs);

  // This is in the * complemented * nonbasic space,
  // so the origin is the solution to the LP relaxation
  // Moreover, every cut will have right-hand side 1,
  // since we want to separate the origin
  genCutsFromPointsRaysNB(cutSolver, beta, structPHA, num_cuts_per_cgs,
      num_cuts_generated, num_obj_tried, castsolver, probData, cgs_ind, cgsName,
      structSICs);
  if (param.paramVal[USE_CUT_VERT_HEUR_PARAM_IND] == 1) {
    cutVerticesCutGeneration(cutSolver, beta, structPHA, num_cuts_per_cgs,
        num_cuts_generated, num_obj_tried, castsolver, probData,
        vertexStore, cgs_ind, cgsName, structSICs);
  }
} /* tryObjectivesForCgs */

bool CglPHA::tieBreakingForChooseHplane(const bool selectByDist,
    const bool selectByDepth, const bool selectByFinalPoints,
    const double distToHplaneChosenForSplit,
    const double addedDepthByHplaneChosenForSplit,
    const int numFinalPointsByHplaneChosenForSplit,
    const double tmpDistToHplane, const double& addedDepthForThisHplane,
    const int numFinalPointsThisHplane) const {
  bool tieBreaking = false;
  if (selectByDist && isVal(tmpDistToHplane, distToHplaneChosenForSplit)) {
    tieBreaking = greaterThanVal(addedDepthForThisHplane,
        addedDepthByHplaneChosenForSplit)
        || (numFinalPointsThisHplane > numFinalPointsByHplaneChosenForSplit);
  } else if (selectByDepth
      && isVal(addedDepthForThisHplane, addedDepthByHplaneChosenForSplit)) {
    tieBreaking = lessThanVal(tmpDistToHplane, distToHplaneChosenForSplit)
        || (numFinalPointsThisHplane > numFinalPointsByHplaneChosenForSplit);
  } else if (selectByFinalPoints
      && (numFinalPointsThisHplane == numFinalPointsByHplaneChosenForSplit)) {
    tieBreaking = lessThanVal(tmpDistToHplane, distToHplaneChosenForSplit)
        || greaterThanVal(addedDepthForThisHplane,
            addedDepthByHplaneChosenForSplit);
  }
  return tieBreaking;
} /* tieBreakingForChooseHplane */

void CglPHA::storeNewRayNBIndices(int numNonBasicCols,
    std::vector<int> &newRayColIndices,
    const std::vector<int> &allRaysCutByHplane,
    std::vector<bool> &raysCutFlag) const {
  //Clear array and fill all raysCutFlag elements to false
  newRayColIndices.clear();
  std::fill(raysCutFlag.begin(), raysCutFlag.end(), 0);

  for (unsigned j = 0; j < allRaysCutByHplane.size(); j++)
    raysCutFlag[allRaysCutByHplane[j]] = true;

  for (int j = 0; j < numNonBasicCols; j++)
    if (raysCutFlag[j] == false)
      newRayColIndices.push_back(j);
} /* storeNewRayNBIndices */

/***********************************************************************/
void CglPHA::generateSICs(AdvCuts& NBSICs, AdvCuts& structSICs, SolutionInfo& solnInfo) {
  // Only feasible ones are returned
  CglSIC SICGen(param);
  SICGen.generateCuts(*solver, NBSICs, structSICs, solnInfo);
#ifdef TRACE
  printf("***** SICs generated: %d.\n", structSICs.sizeCuts());
#endif
}

/***********************************************************************/
void CglPHA::setParam(const CglGICParam &source) {
  param = source;
} /* setParam */

/*********************************************************************/
// Returns true if needs optimal basis to do cuts
bool CglPHA::needsOptimalBasis() const {
  return true;
}
