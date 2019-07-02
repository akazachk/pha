//============================================================================
// Name        : CutHelper.cpp
// Author      : akazachk
// Copyright   : Your copyright notice
// Description : Helper functions of cuts
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#include <cmath> // preferred over math.h
#include <numeric> // inner_prod
#include "GlobalConstants.hpp"
#include "Utility.hpp"
#include "Output.hpp"
#include "CutHelper.hpp"
#include "CoinTime.hpp"
#include "CglSIC.hpp"

/***********************************************************************/
/***********************************************************************/
/** Methods related to generating cuts from non-basic space */
/***********************************************************************/
/***********************************************************************/

/**
 * Generates cuts based on points and rays in IntersectionInfo, which are in the *non-basic* space
 * Cuts are generated in the structural space however
 * @return Number of cuts generated
 */
int genCutsFromPointsRaysNB(PHASolverInterface* const cutSolver,
    const double beta, AdvCuts & cuts, int& num_cuts_this_cgs,
    int& num_cuts_total, int& num_obj_tried,
    const PHASolverInterface* const solver, const SolutionInfo & probData,
    const int cgsIndex, const std::string& cgsName, const AdvCuts& structSICs) {
  if (reachedCutLimit(num_cuts_this_cgs, num_cuts_total)) {
    return 0;
  }
  std::vector<double> obj(cutSolver->getNumCols(), 0.0);

  // Generate cuts from this solver
  const int init_num_cuts = (int) cuts.size();

  if (param.paramVal[NUM_CUTS_ITER_BILINEAR_PARAM_IND] >= 1) {
    // Try the bilinear optimization
    GlobalVariables::timeStats.start_timer(
        CutHeuristicsName[CutHeuristics::ITER_BILINEAR_CUT_HEUR] + "_TIME");
#ifdef TRACE
    printf(
        "\n## Using iterative procedure to get deepest cut post SICs (cgs %d/%d with name %s). ##\n",
        cgsIndex + 1, param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND),
        cgsName.c_str());
#endif

    iterateDeepestCutPostMSIC(cuts, num_cuts_this_cgs, num_cuts_total,
        num_obj_tried, cutSolver, solver, probData, beta, cgsIndex, cgsName,
        structSICs, true);
    GlobalVariables::timeStats.end_timer(
        CutHeuristicsName[CutHeuristics::ITER_BILINEAR_CUT_HEUR] + "_TIME");
  }

  if (reachedCutLimit(num_cuts_this_cgs, num_cuts_total)) {
    return (int) cuts.size() - init_num_cuts;;
  }

  if (param.paramVal[ParamIndices::USE_UNIT_VECTORS_HEUR_PARAM_IND] >= 1) {
    // Try each of the axis directions (corresponds to cutting away SIC points/rays in nb space)
    GlobalVariables::timeStats.start_timer(
        CutHeuristicsName[CutHeuristics::UNIT_VECTORS_CUT_HEUR] + "_TIME");
#ifdef TRACE
    printf(
        "\n## Try finding deepest cut with respect to each of the non-basic axis directions (cgs %d/%d with name %s). ##\n",
        cgsIndex + 1, param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND),
        cgsName.c_str());
#endif
    int nb_dir = 0;
    while (!reachedCutLimit(num_cuts_this_cgs, num_cuts_total)
        && (nb_dir < cutSolver->getNumCols())) { // right now we are using either all or none
      obj[nb_dir] = 1.0;
      cutSolver->setObjective(obj.data());
      num_obj_tried++;
      genCutsFromCutSolver(cuts, cutSolver, solver, probData, beta, cgsIndex,
          cgsName, structSICs, num_cuts_this_cgs, num_cuts_total, true,
          CutHeuristics::UNIT_VECTORS_CUT_HEUR);
      obj[nb_dir] = 0.0;
      nb_dir++;
    }
    GlobalVariables::timeStats.end_timer(
        CutHeuristicsName[CutHeuristics::UNIT_VECTORS_CUT_HEUR] + "_TIME");
  } /* try axis directions (corresponds to cutting away SIC points/rays in nb space) */

  return (int) cuts.size() - init_num_cuts;
} /* genCutsFromPointsRaysNB */

/**
 * Set the cut solver objective based on option, in the *non-basic* space
 * Options are:
 *  1. max_eff
 */
void setCutSolverObjectiveFromVertexNB(const int vert_ind,
    PHASolverInterface* const cutSolver,
    const std::vector<Vertex>& vertexStore) {
  // We set the objective function to be the one to maximize "effectiveness"
  // That is, we want to cut off objVertex by as much as possible
  // To do this we *minimize* alpha * objVertex
  // If < beta, then alpha x >= beta valid and cuts off objVertex
  // However, objVertex is just the origin, which will not lead to any meaningful cut.
  // We instead optimize in the direction (1,1,...,1), which may actually be feasible,
  // but certainly for some t, t * (1,1,...,1) is cut by, say, the SIC.
  cutSolver->setObjSense(1.0);

  const int numElem = vertexStore[vert_ind].getNumElements();
  if (numElem == 0) {
    std::vector<double> obj(cutSolver->getNumCols(), 1.0);
    cutSolver->setObjective(obj.data());
  } else {
    std::vector<double> obj(cutSolver->getNumCols(), 0.0);
    setFullVecFromPackedVec(obj, vertexStore[vert_ind]);
    cutSolver->setObjective(obj.data());
  }
} /* setCutSolverObjectiveFromVertexNB */

/**
 * Set up solver from points and rays provided.
 * The rows are \alpha (p or r) \ge \beta = 1 (or 0 or -1)
 * PHA version
 */
bool setupCutSolver(PHASolverInterface* const cutSolver,
    const IntersectionInfo& interPtsAndRays,
    const int cgsIndex, const std::string& cgsName) {
  const double scale = 1.;

  // Build the matrix as well as col lower/upper bounds
  const int num_cols = interPtsAndRays.getNumCols();
  std::vector<double> colLB(num_cols, -1 * cutSolver->getInfinity());
  std::vector<double> colUB(num_cols, cutSolver->getInfinity());
  std::vector<double> rowLB;
  rowLB.reserve(num_cols + 1);
  int numElementsEstimate = 0; // estimate max size of sparse coefficient matrix, for reserve

  std::vector<int> rowOfPoint, rowOfRay;
  std::vector<double> objDepthPoints, objDepthRays;

  for (int r = 0; r < interPtsAndRays.getNumRows(); r++) {
    const int num_vec_el = interPtsAndRays.getVectorSize(r);
    const double vec_rhs =
        (interPtsAndRays.RHS[r] != 0.) ?
            interPtsAndRays.RHS[r] * scale : interPtsAndRays.RHS[r];
    // Check if it is a bound
    if (num_vec_el == 1) {
      // row is coeff * var >= rhs
      const int start_ind = interPtsAndRays.getVectorFirst(r);
      const double coeff = interPtsAndRays.getElements()[start_ind];
      const int var = interPtsAndRays.getIndices()[start_ind];
      const double old_lb = colLB[var];
      const double old_ub = colUB[var];
      if (isZero(vec_rhs)) {
        if (lessThanVal(coeff, 0.)) {
          if (old_ub > 0.) {
            colUB[var] = 0.;
          }
        } else if (greaterThanVal(coeff, 0.)) {
          if (old_lb < 0.) {
            colLB[var] = 0.;
          }
        } else {
          // Do nothing, because this is a 0 >= 0 row (shouldn't happen... but who knows)
        }
      } /* if zero rhs */
      else {
        if (isVal(coeff, 0.)) {
          // Do nothing, because this is 0 >= something
          if (lessThanVal(vec_rhs, 0.)) {
            error_msg(errorstring,
                "We have an infeasible row: %.2f * x%d >= %.2f.\n", coeff, var,
                vec_rhs);
            writeErrorToII(errorstring, GlobalVariables::log_file);
            exit(1);
          }
        } else {
          const double new_bound = vec_rhs / coeff;
          if (lessThanVal(coeff, 0.)) {
            if (old_ub > new_bound) {
              colUB[var] = new_bound;
            }
          } else {
            if (old_lb < new_bound) {
              colLB[var] = new_bound;
            }
          }
        }
      } /* else it is a non-zero rhs */
    } /* if num_vec_el == 1, it is a bound */

    // When rhs non-zero, we input it anyway, to keep for objective function purposes
    if (num_vec_el > 1 || !isZero(vec_rhs)) {
      numElementsEstimate += num_vec_el;
      if (!isZero(vec_rhs)) {
        rowOfPoint.push_back(r);
        objDepthPoints.push_back(interPtsAndRays.objDepth[r]);
      } else {
        rowOfRay.push_back(r);
        objDepthRays.push_back(interPtsAndRays.objDepth[r]);
      }
    } /* we add it for processing */
  } /* iterate over the rows in the matrix */

  // Sort by orthogonality
  const int numPointsToProcess = objDepthPoints.size();
  std::vector<int> sortIndexPoints(numPointsToProcess);
  for (int i = 0; i < numPointsToProcess; i++) {
    sortIndexPoints[i] = i;
  }
  sort(sortIndexPoints.begin(), sortIndexPoints.end(),
      index_cmp_asc<const double*>(objDepthPoints.data())); // ascending order

  const int numRaysToProcess = objDepthRays.size();
  std::vector<int> sortIndexRays(numRaysToProcess);
  for (int i = 0; i < numRaysToProcess; i++) {
    sortIndexRays[i] = i;
  }
  sort(sortIndexRays.begin(), sortIndexRays.end(),
      index_cmp_asc<const double*>(objDepthRays.data())); // ascending order

  // Now input into the matrix
  CoinPackedMatrix mx; // create a new col-ordered matrix
  mx.reverseOrdering(); // make it row-ordered
  mx.setDimensions(0, num_cols);
  mx.reserve(numPointsToProcess + numRaysToProcess, numElementsEstimate, false);
  //    (numPointsToProcess + numRaysToProcess) * num_cols, false);
  int numPoints = 0, numRays = 0;

  // First the points
  std::vector<int> rowInMxOfPoint;
//    ortho.clear();
//    ortho.reserve(numPointsToProcess + numRaysToProcess);
  for (int p = 0; p < numPointsToProcess; p++) {
//      if (isDuplicatePoint[p]) {
//        continue;
//      }
    const int ind = sortIndexPoints[p];
//      const int t = disjTermOfPoint[ind];
    const int r = rowOfPoint[ind];
    const double rhs = interPtsAndRays.RHS[r] * scale;
    mx.appendRow(interPtsAndRays.getVectorSize(r),
        interPtsAndRays.getIndices() + interPtsAndRays.getVectorFirst(r),
        interPtsAndRays.getElements() + interPtsAndRays.getVectorFirst(r));
    rowLB.push_back(rhs);
    numPoints++;
//    ortho.push_back(objDepthPoints[ind]);
  } /* input points */

  // Then the rays
  for (int p = 0; p < numRaysToProcess; p++) {
//    if (isDuplicateRay[p]) {
//      continue;
//    }
    const int ind = sortIndexRays[p];
//    const int t = disjTermOfRay[ind];
    const int r = rowOfRay[ind];
    const double rhs = interPtsAndRays.RHS[r];
    mx.appendRow(interPtsAndRays.getVectorSize(r),
        interPtsAndRays.getIndices() + interPtsAndRays.getVectorFirst(r),
        interPtsAndRays.getElements()
            + interPtsAndRays.getVectorFirst(r));
    rowLB.push_back(rhs);
    numRays++;
//    ortho.push_back(objDepthRays[ind]);
  } /* input rays */

  // Clean matrix
  mx.cleanMatrix();

//  // Filter out columns that are empty
//  const bool shouldDelete = true;
//  if (shouldDelete) {
//    std::vector<int> deleteColIndex;
//    const int* length = mx.countOrthoLength();
//    for (int col_ind = 0; col_ind < mx.getNumCols(); col_ind++) {
//      if (length[col_ind] == 0) {
//        const bool ignore_lb = isNegInfinity(colLB[col_ind]) || isZero(colLB[col_ind]);
//        const bool ignore_ub = isInfinity(colUB[col_ind]) || isZero(colUB[col_ind]);
//        if (!ignore_lb || !ignore_ub) {
//          nonZeroColIndex.push_back(col_ind);
//        } else {
//          deleteColIndex.push_back(col_ind);
//        }
//      } else {
//        nonZeroColIndex.push_back(col_ind);
//      }
//    }
//    if (length) {
//      delete[] length;
//      length = NULL;
//    }
//    // Delete the columns from mx, as well as from colLB and colUB
//    mx.deleteCols(deleteColIndex.size(), deleteColIndex.data());
//    for (int pos = deleteColIndex.size() - 1; pos >= 0; pos--) {
//      colLB.erase(colLB.begin() + deleteColIndex[pos]);
//      colUB.erase(colUB.begin() + deleteColIndex[pos]);
//    }
//  }

  // Load problem into cutSolver
  // Defaults (in COIN-OR):
  // colLB: 0 **** This should be -inf (except when changed above)
  // colUB: inf
  // obj coeff: 0
  // rowSense: >=
  // rhs: 0 **** We want to change this to beta for the points
  // rowRange: 0
  cutSolver->loadProblem(mx, colLB.data(), colUB.data(), NULL, NULL,
      rowLB.data(), NULL);
  cutSolver->disableFactorization();

  setMessageHandler(cutSolver);
  cutSolver->getModelPtr()->setMaximumIterations(param.getCUTSOLVER_MAX_ITER());
  cutSolver->getModelPtr()->setNumberIterations(0);
  cutSolver->getModelPtr()->setMaximumSeconds(param.getCUTSOLVER_TIMELIMIT());
  cutSolver->setHintParam(OsiDoPresolveInInitial,
      (param.getParamVal(ParamIndices::CUT_PRESOLVE_PARAM_IND) > 0) ?
          true : false);
  cutSolver->setHintParam(OsiDoPresolveInResolve,
      (param.getParamVal(ParamIndices::CUT_PRESOLVE_PARAM_IND) > 0) ?
          true : false);

//  // Defaults:
//  // colLB: 0 **** This should be -inf
//  // colUB: inf
//  // obj coeff: 0
//  // rowSense: >=
//  // rhs: 0 **** We want to change this to beta for the points
//  // rowRange: 0
//  cutSolver->loadProblem(interPtsAndRays, NULL, NULL, NULL, NULL,
//      interPtsAndRays.RHS.data(), NULL);
//
//  for (int j = 0; j < cutSolver->getNumCols(); j++) {
//    cutSolver->setColLower(j, -1.0 * cutSolver->getInfinity());
//  }
//
//  cutSolver->disableFactorization();

  // Check that the cutSolver is feasible for the zero objective function
  // If the number of cuts generated so far is zero, we should try a dummy objective
  GlobalVariables::timeStats.start_timer(
      CutHeuristicsName[CutHeuristics::DUMMY_OBJ_CUT_HEUR] + "_TIME");
#ifdef TRACE
  printf(
      "\n## Check feasibility with all zeroes objective (cgs %d/%d with name %s). ##\n",
      cgsIndex + 1, param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND),
      cgsName.c_str());
#endif
//  std::vector<double> obj(cutSolver->getNumCols(), 0.0);
//  cutSolver->setObjective(obj.data());
  cutSolver->initialSolve();
  GlobalVariables::timeStats.end_timer(
      CutHeuristicsName[CutHeuristics::DUMMY_OBJ_CUT_HEUR] + "_TIME");

  // Check that cutSolver is not primal infeasible
  // This can happen if 0 \in cone(\rayset), for instance, in the non-basic space.
  // We might also get that the cut LP is "bad" and does not solve quickly
  // in which case we abort with numerical issues as the error.
  bool retval;
  if (cutSolver->isProvenOptimal()) {
    retval = true;
  } else if (cutSolver->isProvenPrimalInfeasible()
      && !cutSolver->isAbandoned()) {
    warning_msg(errorstring,
        "PRLP (cgs %d/%d) is primal infeasible; giving up on this point-ray collection.\n",
        cgsIndex + 1, param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND));
    GlobalVariables::numCutSolverFails[CutSolverFails::PRIMAL_CUTSOLVER_NO_OBJ_FAIL_IND]++;
    retval = false;
  } else if (cutSolver->isIterationLimitReached()
      || cutSolver->getModelPtr()->hitMaximumIterations()) {
    warning_msg(errorstring,
        "PRLP (cgs %d/%d) is having numerical issues; giving up on this point-ray collection.\n",
        cgsIndex + 1, param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND));
    GlobalVariables::numCutSolverFails[CutSolverFails::NUMERICAL_ISSUES_NO_OBJ_FAIL_IND]++;
    retval = false;
  } else if (cutSolver->isProvenDualInfeasible() && !cutSolver->isAbandoned()) {
    // This should never happen, but we are going to chock it up to numerical problems
    warning_msg(errorstring,
        "PRLP (cgs %d/%d) is having numerical issues; giving up on this point-ray collection.\n",
        cgsIndex + 1, param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND));
    GlobalVariables::numCutSolverFails[CutSolverFails::NUMERICAL_ISSUES_NO_OBJ_FAIL_IND]++;
    retval = false;
  } else {
    // Something else happened; still bad...
    warning_msg(errorstring,
        "PRLP (cgs %d/%d) is having numerical issues; giving up on this point-ray collection.\n",
        cgsIndex + 1, param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND));
    GlobalVariables::numCutSolverFails[CutSolverFails::NUMERICAL_ISSUES_NO_OBJ_FAIL_IND]++;
    retval = false;
  }

  return retval;
} /* setupCutSolver */

/**
 * Generate a cut, if possible, from the solver
 * Rhs is given as beta
 *
 * return -1 * (fail index + 1) if it is something that means to discontinue this solver
 */
int genCutsFromCutSolver(AdvCuts & cuts,
    PHASolverInterface* const cutSolver,
    const PHASolverInterface* const origSolver,
    const SolutionInfo & origProbData, const double beta, const int cgsIndex,
    const std::string& cgsName, const AdvCuts& structSICs,
    int& num_cuts_this_cgs, int& num_cuts_total, const bool inNBSpace,
    const CutHeuristics cutHeur, const bool tryExtraHard) {
//  GlobalVariables::numObjFromHeur[cutHeur]++;
  if (reachedCutLimit(num_cuts_this_cgs, num_cuts_total)) {
    return exitGenCutsFromCutSolver(cutSolver,
        -1 * (CutSolverFails::CUT_LIMIT_FAIL_IND + 1));
  }

  // First we save the old solution; this is just to have a quick way to check that the same solution was found
  std::vector<double> old_solution;
  if (cuts.size() > 0 && cutSolver->isProvenOptimal()) {
    old_solution.assign(cutSolver->getColSolution(),
        cutSolver->getColSolution() + cutSolver->getNumCols());
  }

  // Resolve
  const double timeLimit = (tryExtraHard) ?  -1. : param.getCUTSOLVER_TIMELIMIT(); // maybe go unlimited time this round
  cutSolver->getModelPtr()->setNumberIterations(0);
  cutSolver->getModelPtr()->setMaximumSeconds(param.getCUTSOLVER_TIMELIMIT());
  cutSolver->resolve();

  if (!checkSolverOptimality(cutSolver, false, timeLimit)) {
  //    if (cutSolver->getModelPtr()->hitMaximumIterations()) {
  //      // Sometimes Clp gets stuck
  //      cutSolver->getModelPtr()->setMaximumSeconds(timeLimit);
  //      cutSolver->initialSolve();
  //      checkSolverOptimality(cutSolver, false, timeLimit, true);
  //    }
      if (cutSolver->isIterationLimitReached()) {
        return exitGenCutsFromCutSolver(cutSolver,
            -1 * (CutSolverFails::ITERATION_CUTSOLVER_FAIL_IND + 1));
      }
      if (cutSolver->getModelPtr()->hitMaximumIterations()) {
        return exitGenCutsFromCutSolver(cutSolver,
            -1 * (CutSolverFails::TIMELIMIT_CUTSOLVER_FAIL_IND + 1));
      }
      if (cutSolver->isAbandoned()) {
        return exitGenCutsFromCutSolver(cutSolver,
            -1 * (CutSolverFails::ABANDONED_CUTSOLVER_FAIL_IND + 1));
      }
    }

  // If "nan" objective, then something went terribly wrong; try again
  if (std::isnan(cutSolver->getObjValue())) {
#ifdef TRACE
    warning_msg(warnstr,
        "Encountered NaN objective value on cgs %d/%d.\n",
        cgsIndex + 1, param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND));
#endif
    cutSolver->getModelPtr()->setNumberIterations(0);
    cutSolver->getModelPtr()->setMaximumSeconds(timeLimit);
    cutSolver->initialSolve();
    checkSolverOptimality(cutSolver, false, timeLimit);
  }

  // Too expensive
  //  if (!cutSolver->isProvenOptimal() && !cutSolver->isProvenDualInfeasible()) {
  //    // Try resolving since this seems to help sometimes
  //    cutSolver->getModelPtr()->setMaximumSeconds(timeLimit);
  //    cutSolver->initialSolve();
  //    checkSolverOptimality(cutSolver, false, timeLimit, true);
  //  }

  if (cutSolver->isProvenOptimal()
      && isNegInfinity(cutSolver->getObjValue(), param.getINFINITY())) {
    // We essentially have a ray, though it is not being reported as such
    // Try resolving because this fixes the issue sometimes
    cutSolver->getModelPtr()->setMaximumSeconds(timeLimit);
    cutSolver->initialSolve();
    checkSolverOptimality(cutSolver, false, timeLimit);
  }

  if (cutSolver->isProvenDualInfeasible()) {
    //    // Is it really true?
    //    cutSolver->getModelPtr()->setMaximumSeconds(param.getCUTSOLVER_TIMELIMIT());
    //    cutSolver->resolve();
    //
    //    // This may stop it from being dual infeasible
    //    if (!cutSolver->isProvenOptimal()
    //        && !cutSolver->isProvenDualInfeasible()) {
    //      // Sometimes Clp gets stuck
    //      cutSolver->getModelPtr()->setMaximumSeconds(param.getCUTSOLVER_TIMELIMIT());
    //      cutSolver->initialSolve();
    //    }
    //
    // Homogeneous cuts never cut the LP optimum in NB space so no point trying rays
    if (cutSolver->isProvenDualInfeasible() && inNBSpace) {
      return exitGenCutsFromCutSolver(cutSolver,
          -1 * (CutSolverFails::DUAL_CUTSOLVER_FAIL_IND + 1));
    }
  }

  int num_cuts_generated = 0;
  const int num_rays_to_get = 1; // Actually this is the only possibility currently

  // If it is proven optimal, we can simply obtain the solution
  // If it is dual infeasible, then it is unbounded (or infeasible, but it would say so)
  // in which case we need to get the ray and use that as the solution
  if (cutSolver->isProvenOptimal()) {
    // Try to generate the cut (but first, test against previous solution)
    if (!old_solution.empty()) {
      double maxDiff = 0.;
      for (int ind = 0; ind < cutSolver->getNumCols(); ind++) {
        const double diff = std::abs(cutSolver->getColSolution()[ind] - old_solution[ind]);
        if (diff > maxDiff) {
          maxDiff = diff;
        }
      }

      if (isZero(maxDiff, origProbData.EPS)) {
        return exitGenCutsFromCutSolver(cutSolver,
            -1 * (CutSolverFails::DUPLICATE_GIC_CUTSOLVER_FAIL_IND + 1));
      }
    } /* check whether it is the same solution as before */

    const int return_code = genCutsHelper(cuts, cutSolver, origSolver,
        origProbData, beta, cgsIndex, cgsName, structSICs, num_cuts_this_cgs,
        num_cuts_total, inNBSpace, cutHeur);
    if (return_code > 0) {
      num_cuts_generated += return_code;
    } else {
      return exitGenCutsFromCutSolver(cutSolver, return_code);
    }
  } else if (cutSolver->isProvenDualInfeasible() && !inNBSpace) {
    // For the time being this seems to not be working well
    // For bell3a for instance, when cutting 64 rays, there is a cut that appears
    // both as a primal ray and as a feasible solution, so it has both >= 0 and >= -1
    // constant sides, and it causes infeasibility of the original problem when added
    return exitGenCutsFromCutSolver(cutSolver,
        -1 * (CutSolverFails::DUAL_CUTSOLVER_FAIL_IND + 1));

    AdvCut currCut(false, 1);
    // We need to get primal rays
    // In case of dual infeasibility, there should be at least one
    std::vector<double*> primalRay(num_rays_to_get);
    for (int i = 0; i < num_rays_to_get; i++) {
      primalRay[i] = new double[cutSolver->getNumCols()]; // TODO valgrind says this isn't getting freed for some reason?
    }
    primalRay = cutSolver->getPrimalRays(num_rays_to_get); // There may be more than one, which all lead to cuts

    const int scaleStatus = scalePrimalRay(cutSolver->getNumCols(), &primalRay[0]);

    if (scaleStatus >= 0) {
      double obj_val;
      if (inNBSpace) {
        AdvCut currCutNB(true, 1);
        currCutNB.cgsIndex = cgsIndex;
        currCutNB.cgsName = cgsName;
        currCutNB.splitVarIndex[0] = cgsIndex;
        currCutNB.setOsiRowCut(cutSolver->getNumCols(), primalRay[0], 0.0,
            origProbData.EPS); // This is an extreme ray, so rhs is 0
        obj_val = currCutNB.row().dotProduct(
            cutSolver->getObjCoefficients());
        AdvCut::convertCutFromJSpaceToStructSpace(currCut, origSolver,
            origProbData.nonBasicVarIndex, currCutNB, true, origProbData.EPS);
      } else {
        obj_val = std::inner_product(origSolver->getObjCoefficients(),
            origSolver->getObjCoefficients()
                + origSolver->getNumCols(), primalRay[0], 0.0);
        currCut.cgsIndex = cgsIndex;
        currCut.cgsName = cgsName;
        currCut.splitVarIndex[0] = cgsIndex;
        currCut.setOsiRowCut(cutSolver->getNumCols(), primalRay[0], 0.0, origProbData.EPS); // This is an extreme ray, so rhs is 0
      }

      const bool cuttingSolutionFlag = lessThanVal(obj_val, 0.0);
      const bool duplicateSICFlag = structSICs.isDuplicateCut(currCut);
      const bool duplicateGICFlag = cuts.isDuplicateCut(currCut);
      const bool toAdd = cuttingSolutionFlag && !duplicateSICFlag
          && !duplicateGICFlag
          && !reachedCutLimit(num_cuts_this_cgs, num_cuts_total);
      if (toAdd) {
        currCut.cutHeur = cutHeur;
        GlobalVariables::numCutsFromHeur[cutHeur]++;
        cuts.addNonDuplicateCut(currCut, false);
        num_cuts_generated++;
        num_cuts_this_cgs++;
        num_cuts_total++;
      } else {
        GlobalVariables::numCutSolverFails[CutSolverFails::NON_CUTTING_CUTSOLVER_FAIL_IND] +=
            !cuttingSolutionFlag;
        GlobalVariables::numCutSolverFails[CutSolverFails::DUPLICATE_SIC_CUTSOLVER_FAIL_IND] +=
            duplicateSICFlag;
        GlobalVariables::numCutSolverFails[CutSolverFails::DUPLICATE_GIC_CUTSOLVER_FAIL_IND] +=
            duplicateGICFlag;
      }
      for (int i = 0; i < num_rays_to_get; i++) {
        delete[] primalRay[i];
      }
      primalRay.clear();
    } else {
      return exitGenCutsFromCutSolver(cutSolver,
          -1 * (CutSolverFails::DUAL_CUTSOLVER_FAIL_IND + 1));
    }
//    primalRay.clear();
//    std::vector<double*>().swap(primalRay);
  }
  return exitGenCutsFromCutSolver(cutSolver, num_cuts_generated);
} /* genCutsFromCutSolver */

int genCutsHelper(AdvCuts & cuts, PHASolverInterface* const cutSolver,
    const PHASolverInterface* const origSolver,
    const SolutionInfo & origProbData, const double beta, const int cgsIndex,
    const std::string& cgsName, const AdvCuts& structSICs,
    int& num_cuts_this_cgs, int& num_cuts_total, const bool inNBSpace,
    const CutHeuristics cutHeur) {
  // We might get here without cutSolver being optimal
  if (!cutSolver->isProvenOptimal()) {
    return 0; // errors are tabulated in exitGenCutsFromCutSolver
  } /* may have run into issues */
  if (reachedCutLimit(num_cuts_this_cgs, num_cuts_total)) {
    return 0; // will tabulate cut limit fails later
  }

  int num_cuts_generated = 0;
  AdvCut currCut(false, 1);
  currCut.cgsIndex = cgsIndex;
  currCut.cgsName = cgsName;
  currCut.splitVarIndex[0] = cgsIndex;
  if (inNBSpace) {
    currCut.setCutFromJSpaceCoefficients(cutSolver->getColSolution(), beta,
        origSolver, origProbData.nonBasicVarIndex, true, origProbData.EPS);
  } else {
    currCut.setOsiRowCut(cutSolver->getNumCols(), cutSolver->getColSolution(),
        beta, origProbData.EPS);
  }

  // Clean
  const int clean_error = currCut.cleanCut(origSolver, origProbData/*, beta*/);
  if (clean_error != 0) {
//    GlobalVariables::numCutSolverFails[clean_error]++;
//    return 0;
    return clean_error;
  }

  // Ensure the cut is not a duplicate, and some other checks
  const bool generateAnywayIfDuplicateSIC = true; // in principle, when this is true, we do not need to check ortho with SICs; in our tests, we might anyway if we want to get the relevant stats
  const bool checkIfDuplicateSIC = false;
  const double currCutNorm = currCut.getTwoNorm();
  const double violation = currCut.violated(origSolver->getColSolution());
  int duplicateSICIndex = -1, duplicateGICIndex = -1;
  int minOrthogonalityIndex = -1;
  double orthogonalityWithSICs = 1., orthogonalityWithGICs = 1.;
  if (checkIfDuplicateSIC) {
    structSICs.howDuplicate(currCut, duplicateSICIndex, minOrthogonalityIndex,
        orthogonalityWithSICs, origProbData.EPS);
  }
  const bool duplicateSICFlag = (duplicateSICIndex >= 0);
  const bool orthogonalitySICFailFlag = !duplicateSICFlag
      && (orthogonalityWithSICs < param.getMINORTHOGONALITY());
  if (generateAnywayIfDuplicateSIC
      || (!duplicateSICFlag && !orthogonalitySICFailFlag)) {
    cuts.howDuplicate(currCut, duplicateGICIndex, minOrthogonalityIndex,
        orthogonalityWithGICs, origProbData.EPS);
  }
  const bool duplicateGICFlag = (duplicateGICIndex >= 0);
  bool orthogonalityGICFailFlag = !duplicateGICFlag
      && (orthogonalityWithGICs < param.getMINORTHOGONALITY());
  if (orthogonalityGICFailFlag) {
    // If rhs is better or sparser than existing cut, replace the cut
    const int cut_num_el =
        cuts.cuts[minOrthogonalityIndex].row().getNumElements();
    //    const int* fpc_indices =
    //        fpcs.cuts[minOrthogonalityIndex].row().getIndices();
    //    const double* fpc_vals =
    //        fpcs.cuts[minOrthogonalityIndex].row().getElements();
    const bool isSparser = currCut.row().getNumElements() < cut_num_el;
    const double old_effectiveness =
        cuts.cuts[minOrthogonalityIndex].effectiveness();
    const bool isMoreEffective = (violation / currCutNorm) > old_effectiveness;

    if (isSparser || isMoreEffective) {
      // Adjust statistics
      GlobalVariables::numCutsFromHeur[cuts.cuts[minOrthogonalityIndex].cutHeur]--;
      GlobalVariables::numCutsFromHeur[cutHeur]++;

      // Replace old cut
      currCut.cutHeur = cutHeur;
      currCut.setEffectiveness(violation / currCutNorm);
      cuts.cuts[minOrthogonalityIndex] = currCut;

      orthogonalityGICFailFlag = true;

#ifdef SHOULD_CALC_ACTIVE_FINAL_POINTS
      std::vector<double> tmpvec(cutSolver->getColSolution(), cutSolver->getColSolution() + cutSolver->getNumCols());
      GlobalVariables::cutCoeffsNB[minOrthogonalityIndex] = tmpvec; // modified from code below 2018/11/03
//      GlobalVariables::cutCoeffsNB[minOrthogonalityIndex].insert(
//          GlobalVariables::cutCoeffsNB[minOrthogonalityIndex].begin(),
//          cutSolver->getColSolution(),
//          cutSolver->getColSolution() + cutSolver->getNumCols());
      GlobalVariables::cutRHSNB[minOrthogonalityIndex] = beta;
#endif
    }
  } /* orthogonalityGICFail (replace cut possibly) */
  const bool toAdd = (generateAnywayIfDuplicateSIC
      || (!duplicateSICFlag && !orthogonalitySICFailFlag)) && !duplicateGICFlag
      && !orthogonalityGICFailFlag;
  if (toAdd) {
    currCut.cutHeur = cutHeur;
    currCut.setEffectiveness(violation / currCutNorm);
    GlobalVariables::numCutsFromHeur[cutHeur]++;
    cuts.addNonDuplicateCut(currCut, false);
    num_cuts_generated++;
    num_cuts_this_cgs++;
    num_cuts_total++;

#ifdef SHOULD_CALC_ACTIVE_FINAL_POINTS
      std::vector<double> tmpvec(cutSolver->getColSolution(), cutSolver->getColSolution() + cutSolver->getNumCols());
      GlobalVariables::cutCoeffsNB.push_back(tmpvec);
      GlobalVariables::cutRHSNB.push_back(beta);
#endif
  } else {
    // Note that these are all mutually exclusive IF we are checking for duplicates
    if (generateAnywayIfDuplicateSIC) {
      GlobalVariables::numCutSolverFails[CutSolverFails::DUPLICATE_SIC_CUTSOLVER_FAIL_IND] +=
          duplicateSICFlag;
      GlobalVariables::numCutSolverFails[CutSolverFails::ORTHOGONALITY_CUTSOLVER_FAIL_IND] +=
          orthogonalitySICFailFlag;
    } else if (duplicateSICFlag) {
      return -1 * (CutSolverFails::DUPLICATE_SIC_CUTSOLVER_FAIL_IND + 1);
    } else if (orthogonalitySICFailFlag) {
      return -1 * (CutSolverFails::ORTHOGONALITY_CUTSOLVER_FAIL_IND + 1);
    }
    if (duplicateGICFlag) {
      return -1 * (CutSolverFails::DUPLICATE_GIC_CUTSOLVER_FAIL_IND + 1);
    }
    if (orthogonalityGICFailFlag) {
      return -1 * (CutSolverFails::ORTHOGONALITY_CUTSOLVER_FAIL_IND + 1);
    }
  }
  return num_cuts_generated;
} /* genCutsHelper */

int exitGenCutsFromCutSolver(PHASolverInterface* const cutSolver,
    const int return_status) {
  if (return_status < 0) {
    const int error_index = -1 * (return_status + 1);
    GlobalVariables::numCutSolverFails[error_index]++;
  }
  // Reset timer and number iterations
  cutSolver->getModelPtr()->setNumberIterations(0);
  cutSolver->getModelPtr()->setMaximumSeconds(param.getCUTSOLVER_TIMELIMIT());
  return return_status;
} /* exitGenCutsFromCutSolver */

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
    const bool inNBSpace) {
  const int init_num_cuts = (int) cuts.size();

  GlobalVariables::timeStats.start_timer(CutHeuristicsName[CutHeuristics::TIGHT_POINTS_CUT_HEUR] + "_TIME");
#ifdef TRACE
  printf(
      "\n## Try finding a cut tight on each intersection point or ray (cgs %d/%d with name %s): %d total. ##\n",
      cgsIndex + 1, param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND),
      cgsName.c_str(), cutSolver->getNumRows());
#endif
  const CoinPackedMatrix* mat = cutSolver->getMatrixByRow();
  const int paramNumRowsToTry =
      param.paramVal[ParamIndices::USE_TIGHT_POINTS_HEUR_PARAM_IND];
  const int numRowsToTry =
      (mat->getNumRows() < paramNumRowsToTry) ?
          mat->getNumRows() : paramNumRowsToTry;
  for (int ind = 0; ind < numRowsToTry; ind++) {
    if (reachedCutLimit(num_cuts_this_cgs, num_cuts_total)) {
      break; // cut limit reached
    }
    const int row_ind = mat->getNumRows() - 1 - ind;
    std::vector<double> obj(cutSolver->getNumCols(), 0.0);
    setFullVecFromPackedVec(obj, mat->getVector(row_ind));
    cutSolver->setObjective(obj.data());
//    cutSolver->initialSolve();
//    cutSolver->resolve();
    num_obj_tried++;
    const int return_code = genCutsFromCutSolver(cuts, cutSolver, solver,
        probData, beta, cgsIndex, cgsName, structSICs, num_cuts_this_cgs,
        num_cuts_total, inNBSpace, CutHeuristics::TIGHT_POINTS_CUT_HEUR);
    if (return_code == -1 * (CutSolverFails::PRIMAL_CUTSOLVER_FAIL_IND + 1)) {
      break; // abandon, may be primal infeasible
    }
  }
  GlobalVariables::timeStats.end_timer(CutHeuristicsName[CutHeuristics::TIGHT_POINTS_CUT_HEUR] + "_TIME");

  return (int) cuts.size() - init_num_cuts;;
} /* tightOnPointsAndRaysGeneration */

/**
 * @brief Try cutting off each of the vertices generated in the intermediate stages
 */
int cutVerticesCutGeneration(PHASolverInterface* const cutSolver,
    const double beta, AdvCuts & cuts, int& num_cuts_this_cgs,
    int& num_cuts_total, int& num_obj_tried,
    const PHASolverInterface* const solver, const SolutionInfo & probData,
    const std::vector<Vertex>& vertexStore, const int cgsIndex,
    const std::string& cgsName, const AdvCuts& structSICs) {
  const int init_num_cuts = (int) cuts.size();

  GlobalVariables::timeStats.start_timer(
      CutHeuristicsName[CutHeuristics::CUT_VERTICES_CUT_HEUR] + "_TIME");
#ifdef DEBUG
  printf(
      "\n## Try cutting each of the generated vertices by as much as possible (cgs %d/%d with name %s). ##\n",
      cgsIndex + 1, param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND),
      cgsName.c_str());
#endif
  for (int vert_ind = 0; vert_ind < (int) vertexStore.size(); vert_ind++) {
    if (reachedCutLimit(num_cuts_this_cgs, num_cuts_total)) {
      break; // cut limit reached
    }
    setCutSolverObjectiveFromVertexNB(vert_ind, cutSolver, vertexStore);
    num_obj_tried++;
    const int return_code = genCutsFromCutSolver(cuts, cutSolver, solver,
        probData, beta, cgsIndex, cgsName, structSICs, num_cuts_this_cgs,
        num_cuts_total, true, CutHeuristics::CUT_VERTICES_CUT_HEUR);
    if (return_code == -1 * (CutSolverFails::PRIMAL_CUTSOLVER_FAIL_IND + 1)) {
      break; // abandon, may be primal infeasible
    }
  }
  GlobalVariables::timeStats.end_timer(
      CutHeuristicsName[CutHeuristics::CUT_VERTICES_CUT_HEUR] + "_TIME");

  return (int) cuts.size() - init_num_cuts;;
} /* cutVerticesCutGeneration */

/**
 * For each ray, look at set of intersection points (across all splits) generated by this ray
 * Call this set P(r), corresponding to ray r
 * Let p(r,s) be the point generated by ray r for split s
 * Let d(r,s) = distance along ray r to split s
 * Then for any split s, we can try to cut points p(r,s') for all s' : d(r,s') < d(r,s)
 */
void splitSharingCutGeneration(std::vector<PHASolverInterface*>& cutSolver,
    AdvCuts & cuts, std::vector<int>& num_cuts_per_cgs, int& num_cuts_total,
    int& num_obj_tried, const std::vector<IntersectionInfo>& interPtsAndRays,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
    const std::vector<Vertex>& vertexStore,
    const std::vector<Hplane>& hplaneStore, const AdvCuts& structSICs,
    const std::vector<bool>& isCutSolverPrimalFeasible) {
  GlobalVariables::timeStats.start_timer(
      CutHeuristicsName[CutHeuristics::SPLIT_SHARE_CUT_HEUR] + "_TIME");
  // point analysis must have been done prior to entering this function
  // (this fills the objective bound vector and identifies points)
  const double beta = 1.0; // We are in the complemented nb space
  for (int s = 0; s < solnInfo.numFeasSplits; s++) {
    if (reachedCutLimit(num_cuts_per_cgs[s], num_cuts_total)) {
      continue;
    }
    std::vector<double> obj(interPtsAndRays[s].getNumCols(), 0.0);

    // Sort points based on obj cost
//    std::vector<int> sortIndex(interPtsAndRays[s].pointRowIndex);
    std::vector<int> sortIndex(interPtsAndRays[s].getNumRows());
    for (int i = 0; i < interPtsAndRays[s].getNumRows(); i++) {
      sortIndex[i] = i;
    }

//    std::vector<double> objDepth(interPtsAndRays[s].objDepth);
    std::sort(sortIndex.begin(), sortIndex.end(),
        index_cmp_asc<const std::vector<double>&>(interPtsAndRays[s].objDepth));

//    const int num_points = interPtsAndRays[s].pointRowIndex.size();
    int num_points_tried = 0;
    for (int ind = 0;
        (ind < interPtsAndRays[s].getNumRows())
            && (num_points_tried
                < std::abs(param.paramVal[ParamIndices::USE_SPLIT_SHARE_PARAM_IND]));
        ind++) {
      if (reachedCutLimit(num_cuts_per_cgs[s], num_cuts_total)) {
        break;
      }
      const int row_ind = sortIndex[ind];
      if (isZero(interPtsAndRays[s].RHS[row_ind])) {
        continue;
      }
      num_points_tried++;
      const int deep_split_index = findDeepestSplit(interPtsAndRays[s],
                row_ind, s, solver, solnInfo,
                vertexStore, hplaneStore);
      if ((deep_split_index == -1) || !isCutSolverPrimalFeasible[s]
          || reachedCutLimit(num_cuts_per_cgs[deep_split_index],
              num_cuts_total)) {
        continue;
      }

      // Cut this point using the deepest split
      setFullVecFromPackedVec(obj, interPtsAndRays[s].getVector(row_ind));

      // Set the objective for the deepest split solver
      cutSolver[deep_split_index]->setObjective(obj.data());
//      cutSolver[deep_split_index]->initialSolve();
//      cutSolver[deep_split_index]->resolve();
      num_obj_tried++;
      genCutsFromCutSolver(cuts, cutSolver[deep_split_index], solver, solnInfo,
          beta, solnInfo.feasSplitVar[deep_split_index],
          solver->getColName(solnInfo.feasSplitVar[deep_split_index]),
          structSICs, num_cuts_per_cgs[deep_split_index], num_cuts_total, true,
          CutHeuristics::SPLIT_SHARE_CUT_HEUR);

      // Clear the objective vector back to zeroes
      clearFullVecFromPackedVec(obj, interPtsAndRays[s].getVector(row_ind));
    }
  }
  GlobalVariables::timeStats.end_timer(
      CutHeuristicsName[CutHeuristics::SPLIT_SHARE_CUT_HEUR] + "_TIME");
} /* splitSharing */

/**
 * @brief For the ray that generated the given point, find the deepest split the ray intersects
 */
int findDeepestSplit(const IntersectionInfo& interPtsAndRays, const int pt_ind,
    const int split_ind, const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo,
//    const std::vector<Ray>& rayStore,
    const std::vector<Vertex>& vertexStore,
    const std::vector<Hplane>& hplaneStore) {
  const int hplane_ind =
      (interPtsAndRays.hplane[pt_ind] >= 0) ?
          interPtsAndRays.hplane[pt_ind] : 0;

  Ray newRay;
  calcNewRayCoordinatesRank1(newRay, /*solver,*/ solnInfo,
        interPtsAndRays.cutRay[pt_ind], interPtsAndRays.newRay[pt_ind], /*rayStore,*/
        hplaneStore[hplane_ind], param.getRAYEPS());

  // Calculate deepest split
  int max_split_ind = -1;
  double maxDistToSplit = -1.;
  double distToThisSplit = -1.;
  for (int s = 0; s < solnInfo.numFeasSplits; s++) {
    const double distToSplit = distanceToSplit(
        vertexStore[interPtsAndRays.vertex[pt_ind]], newRay,
        solnInfo.feasSplitVar[s], solver, solnInfo);
    if (!isInfinity((distToSplit))
        && greaterThanVal(distToSplit, maxDistToSplit)) {
      maxDistToSplit = distToSplit;
      max_split_ind = s;
    }
    if (s == split_ind) {
      distToThisSplit = distToSplit;
    }
  }

  if (lessThanVal(distToThisSplit, maxDistToSplit))
    return max_split_ind;
  else
    return -1;
} /* findDeepestSplit */

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
    const bool inNBSpace) {
  if (reachedCutLimit(num_cuts_this_cgs, num_cuts_total)
      || (param.paramVal[NUM_CUTS_ITER_BILINEAR_PARAM_IND] == 0)) {
    return 0;
  }

  bool stationary_point_flag = false;

  // Get solution after adding all structSICs
  PHASolverInterface* MSICSolver = dynamic_cast<PHASolverInterface*>(origSolver->clone());
  setMessageHandler(MSICSolver);
  if (isInfinity(applyCuts(MSICSolver, structSICs))) {
    error_msg(error_str,
      "Adding SICs causes infeasible linear program. This should not happen!\n");
    exit(1);
  }

  std::vector<double> NBSICSolution;
  if (inNBSpace) {
    compNBCoor(NBSICSolution, MSICSolver->getColSolution(), MSICSolver,
        origSolver, origSolnInfo);
    cutSolver->setObjective(NBSICSolution.data());
  } else {
    cutSolver->setObjective(MSICSolver->getColSolution());
  }

  num_obj_tried++;
  int num_cuts_gen = genCutsFromCutSolver(cuts, cutSolver, origSolver,
      origSolnInfo, beta, cgsIndex, cgsName, structSICs, num_cuts_this_cgs,
      num_cuts_total, true, CutHeuristics::ITER_BILINEAR_CUT_HEUR);

  if (num_cuts_gen < 0) {
    if (MSICSolver) {
      delete MSICSolver;
    }
    return 0;
  }

  // We may need to get primal rays
  // In case of dual infeasibility, there should be at least one
  // During this procedure, since we are using dual rays, do not use presolve
  cutSolver->setHintParam(OsiDoPresolveInInitial, false);
  cutSolver->setHintParam(OsiDoPresolveInResolve, false);

  // During this procedure, since we are using dual rays, do not use presolve
  MSICSolver->setHintParam(OsiDoPresolveInInitial, false);
  MSICSolver->setHintParam(OsiDoPresolveInResolve, false);


//  if (isNegInfinity(cutSolver->getObjValue(), param.getINFINITY())) {
//    // We essentially have a ray, though it is not being reported as such
//    // Try resolving because this fixes the issue sometimes
//    cutSolver->resolve();
//
//    // Sometimes it simply takes too long to solve optimally
//    if (cutSolver->isIterationLimitReached() || cutSolver->getModelPtr()->hitMaximumIterations()) {
////      GlobalVariables::numCutSolverFails[CutSolverFails::ITERATION_CUTSOLVER_FAIL_IND]++;
//      return exitFromIterateDeepestCutPostMSIC(num_cuts_gen,
//          cutSolver, MSICSolver, NULL);
//    }
//    if (!cutSolver->isProvenDualInfeasible()) {
//      warning_msg(warnstring,
//          "Scaling issues in cutSolver for cgs %d/%d with name %s when trying to iterate deepest cut post SICs. Abandoning this function.\n",
//          cgsIndex + 1, param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND),
//          cgsName.c_str());
//      // We don't exit out because maybe this is still salvageable
//      //GlobalVariables::numCutSolverFails[CutSolverFails::DUAL_CUTSOLVER_FAIL_IND]++;
//      return exitFromIterateDeepestCutPostMSIC(num_cuts_gen,
//          cutSolver, MSICSolver, NULL);
//    }
//  }
//
//  PHASolverInterface* iterativeSICSolver = dynamic_cast<PHASolverInterface*>(MSICSolver->clone());
//  setMessageHandler(iterativeSICSolver);
//  // During this procedure, since we are using dual rays, do not use presolve
//  iterativeSICSolver->setHintParam(OsiDoPresolveInInitial, false);
//  iterativeSICSolver->setHintParam(OsiDoPresolveInResolve, false);

  std::vector<double> SICSolution(MSICSolver->getColSolution(),
      MSICSolver->getColSolution() + MSICSolver->getNumCols());
  std::vector<double> structSpaceCut;
  for (int j = 1;
      (j < param.paramVal[NUM_CUTS_ITER_BILINEAR_PARAM_IND])
          && !reachedCutLimit(num_cuts_this_cgs, num_cuts_total); j++) {
    for (int i = 0; i < param.getITER_PER_CUT_BILINEAR(); i++) {
      if (cutSolver->isProvenOptimal()) {
        // Recall that iterativeSICSolver is in the structural space,
        // while cutSolver is in the NB space potentially
        if (inNBSpace) {
          AdvCut::convertCutCoeffFromJSpaceToStructSpaceNonPacked(
              structSpaceCut, origSolver,
              origSolnInfo.nonBasicVarIndex,
              cutSolver->getColSolution(), 1.0, true);
          MSICSolver->setObjective(structSpaceCut.data());
        } else {
          MSICSolver->setObjective(cutSolver->getColSolution());
        }
      } else if (cutSolver->isProvenDualInfeasible()) {
        //cutSolverPrimalRay.resize(1);
        //cutSolverPrimalRay[0] = new double[cutSolver->getNumCols()]; // valgrind says this isn't getting freed for some reason
        //cutSolverPrimalRay = cutSolver->getPrimalRays(1); // There may be more than one, which all lead to cuts
        std::vector<double*> cutSolverPrimalRay = cutSolver->getPrimalRays(1); // There may be more than one, which all lead to cuts
        
        // It has happened (e.g., instance rlp2)
        if (cutSolverPrimalRay[0] == NULL) {
          cutSolver->initialSolve();
          if (cutSolver->isProvenDualInfeasible()) {
            cutSolverPrimalRay = cutSolver->getPrimalRays(1);
          }
          if (!cutSolver->isProvenDualInfeasible() || !cutSolverPrimalRay[0]) {
            warning_msg(warnstring,
                "Issues obtaining ray when cut solver is unbounded on cgs %d/%d with name %s when trying to iterate deepest cut post SICs to get cut %d in iteration %d. Abandoning this function.\n",
                cgsIndex + 1,
                param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND),
                cgsName.c_str(), j, i);
          return exitFromIterateDeepestCutPostMSIC(num_cuts_gen,
              cutSolver, MSICSolver);
          }
        }

        scalePrimalRay(cutSolver->getNumCols(), &cutSolverPrimalRay[0]);
        if (inNBSpace) {
          AdvCut::convertCutCoeffFromJSpaceToStructSpaceNonPacked(
              structSpaceCut, origSolver,
              origSolnInfo.nonBasicVarIndex,
              cutSolverPrimalRay[0], 1.0, true);
          MSICSolver->setObjective(structSpaceCut.data());
        } else {
          MSICSolver->setObjective(cutSolverPrimalRay[0]);
        }

        // Free memory
        if (cutSolverPrimalRay.size() > 0 && cutSolverPrimalRay[0]) {
          delete[] cutSolverPrimalRay[0];
        }
        cutSolverPrimalRay.clear();
      }

      MSICSolver->initialSolve();
      checkSolverOptimality(MSICSolver, false, param.getCUTSOLVER_TIMELIMIT());
//      MSICSolver->resolve(); // Because in pk1, using fpcs, we get a time in which the iterSICSolver is dual infeasible but does not say it is proven so
      if (!MSICSolver->isProvenOptimal()
          && !MSICSolver->isProvenDualInfeasible()) {
        error_msg(errorstring,
            "cgs %d/%d with name %s: Iterative SIC solver not proven optimal or dual infeasible (j: %d, iter: %d).\n",
            cgsIndex + 1, param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND),
            cgsName.c_str(), j, i);
        writeErrorToII(errorstring, GlobalVariables::log_file);
        exit(1);
      }

      if (isNegInfinity(MSICSolver->getObjValue(), param.getINFINITY())) {
        // We essentially have a ray, though it is not being reported as such
        // Try resolving because this fixes the issue sometimes
        MSICSolver->resolve();
        checkSolverOptimality(MSICSolver, false, param.getCUTSOLVER_TIMELIMIT());

        if (!MSICSolver->isProvenDualInfeasible()) {
          warning_msg(warnstring,
              "Scaling issues in iterativeSICSolver for split %d when trying to iterate deepest cut post SICs to get cut %d in iteration %d. Abandoning this function.\n",
              cgsIndex, j, i);
          return exitFromIterateDeepestCutPostMSIC(num_cuts_gen,
              cutSolver, MSICSolver);
        }
      }

      if (MSICSolver->isProvenOptimal()) {
        // Check if we have reached a stationary point
        if (vectorsEqual(MSICSolver->getNumCols(), SICSolution.data(),
            MSICSolver->getColSolution())) {
          stationary_point_flag = true;
          break;
        }

        // Otherwise, reset SICSolution to be current solution and set objective
        SICSolution.assign(MSICSolver->getColSolution(),
            MSICSolver->getColSolution() + MSICSolver->getNumCols());
      } else if (MSICSolver->isProvenDualInfeasible()) {
        // Obtain the dual ray
        std::vector<double*> SICSolverPrimalRay = MSICSolver->getPrimalRays(1);
        
        // It has happened (e.g., instance rlp2)
        if (SICSolverPrimalRay[0] == NULL) {
          MSICSolver->initialSolve();
          if (MSICSolver->isProvenDualInfeasible()) {
            SICSolverPrimalRay = MSICSolver->getPrimalRays(1);
          }
          if (!MSICSolver->isProvenDualInfeasible()
              || !SICSolverPrimalRay[0]) {
            warning_msg(warnstring,
                "Issues obtaining ray when iterative SIC solver is unbounded on cgs %d/%d with name %s when trying to iterate deepest cut post SICs to get cut %d in iteration %d. Abandoning this function.\n",
                cgsIndex + 1,
                param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND),
                cgsName.c_str(), j, i);
            return exitFromIterateDeepestCutPostMSIC(num_cuts_gen, cutSolver,
                MSICSolver);
          }
        }
        scalePrimalRay(MSICSolver->getNumCols(), &SICSolverPrimalRay[0]);

        // Check if we have reached a stationary point
        if (vectorsEqual(MSICSolver->getNumCols(),
            SICSolution.data(), SICSolverPrimalRay[0])) {
          stationary_point_flag = true;
          // Free memory
          if (SICSolverPrimalRay.size() > 0 && SICSolverPrimalRay[0]) {
            delete[] SICSolverPrimalRay[0];
          }
          SICSolverPrimalRay.clear();
          break;
        }

        // Otherwise, reset SICSolution to be current solution and set objective
        SICSolution.assign(SICSolverPrimalRay[0],
            SICSolverPrimalRay[0] + MSICSolver->getNumCols());
        
        // Free memory
        if (SICSolverPrimalRay.size() > 0 && SICSolverPrimalRay[0]) {
          delete[] SICSolverPrimalRay[0];
        }
        SICSolverPrimalRay.clear();
      }

      if (inNBSpace) {
        compNBCoor(NBSICSolution, SICSolution.data(), MSICSolver,
            origSolver, origSolnInfo);
        cutSolver->setObjective(NBSICSolution.data());
      } else {
        cutSolver->setObjective(SICSolution.data());
      }

      // Resolve PRLP without generating a cut
      cutSolver->getModelPtr()->setNumberIterations(0);
      cutSolver->getModelPtr()->setMaximumSeconds(param.getCUTSOLVER_TIMELIMIT());
      cutSolver->resolve();
      checkSolverOptimality(cutSolver, false, param.getCUTSOLVER_TIMELIMIT());
      // Sometimes it simply takes too long to solve optimally
      if (cutSolver->isIterationLimitReached() || cutSolver->getModelPtr()->hitMaximumIterations()) {
        return exitFromIterateDeepestCutPostMSIC(num_cuts_gen, cutSolver,
            MSICSolver);
      }
      if (!cutSolver->isProvenOptimal()
          && !cutSolver->isProvenDualInfeasible()) {
        error_msg(errorstring,
            "cgs %d/%d with name %s: cutSolver issue in iter bilinear: not proven optimal and not dual infeasible (unbounded) (j: %d, iter: %d).\n",
            cgsIndex + 1, param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND),
            cgsName.c_str(), j, i);
        writeErrorToII(errorstring, GlobalVariables::log_file);
        exit(1);
      }

      if (isNegInfinity(cutSolver->getObjValue(), param.getINFINITY())) {
        // We essentially have a ray, though it is not being reported as such
        // Try resolving because this fixes the issue sometimes
        cutSolver->resolve();
        checkSolverOptimality(cutSolver, false, param.getCUTSOLVER_TIMELIMIT());

        // Sometimes it simply takes too long to solve optimally
        if (cutSolver->isIterationLimitReached()
            || cutSolver->getModelPtr()->hitMaximumIterations()) {
//          GlobalVariables::numCutSolverFails[CutSolverFails::ITERATION_CUTSOLVER_FAIL_IND]++;
          return exitFromIterateDeepestCutPostMSIC(num_cuts_gen, cutSolver,
              MSICSolver);
        }
        if (!cutSolver->isProvenDualInfeasible()) {
          warning_msg(warnstring,
              "Scaling issues for cgs %d/%d with name %s when trying to iterate deepest cut post SICs to get cut %d in iteration %d. Abandoning this function.\n",
              cgsIndex + 1, param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND),
              cgsName.c_str(), j, i);
          return exitFromIterateDeepestCutPostMSIC(num_cuts_gen, cutSolver,
              MSICSolver);
        }
      }
    } /* do this many iterations per obj */
    num_obj_tried++;
    const int return_code = genCutsFromCutSolver(cuts, cutSolver, origSolver,
        origSolnInfo, beta, cgsIndex, cgsName, structSICs, num_cuts_this_cgs,
        num_cuts_total, true, CutHeuristics::ITER_BILINEAR_CUT_HEUR);
    if (return_code >= 0) {
      num_cuts_gen += return_code;
    } else {
      return exitFromIterateDeepestCutPostMSIC(num_cuts_gen, cutSolver,
          MSICSolver);
    }
    if (stationary_point_flag)
      break;
  } /* loop at most num_cuts_iter_bilinear times */

  return exitFromIterateDeepestCutPostMSIC(num_cuts_gen, cutSolver, MSICSolver);
} /* iterateDeepestCutPostMSIC */

int exitFromIterateDeepestCutPostMSIC(const int num_cuts_gen,
    PHASolverInterface* const cutSolver,
    PHASolverInterface* const MSICSolver) {
  // Delete solvers
  if (MSICSolver) {
    delete MSICSolver;
  }

  cutSolver->getModelPtr()->setMaximumIterations(param.getCUTSOLVER_MAX_ITER());
  cutSolver->getModelPtr()->setMaximumSeconds(param.getCUTSOLVER_TIMELIMIT());

  cutSolver->setHintParam(OsiDoPresolveInInitial,
      (param.getParamVal(ParamIndices::CUT_PRESOLVE_PARAM_IND) > 0) ?
          true : false);
  cutSolver->setHintParam(OsiDoPresolveInResolve,
      (param.getParamVal(ParamIndices::CUT_PRESOLVE_PARAM_IND) > 0) ?
          true : false);

  return num_cuts_gen;
} /* exitFromIterateDeepestCutPostMSIC */

/**
 * Scales ray generated.
 * If max is not too big, then we do not bother rescaling
 * If min / norm too small, then we do not do that scaling
 *
 * @return status: -1 is bad ray, don't add it; 0 is no scaling; 1 is okay scaling
 */
int scalePrimalRay(int num_cols, double** rayToScale,
    const double SENSIBLE_MAX) {
  double totalSOS = 0.0, minAbsElem = 0.0, maxAbsElem = 0.0;
  for (int i = 0; i < num_cols; i++) {
    double curr = (*rayToScale)[i];
    if (isZero(curr, param.getEPS())) {
      curr = +0.0;
      (*rayToScale)[i] = +0.0;
      continue;
    }

    const double absCurr = std::abs(curr);
    totalSOS += absCurr * absCurr;

    if ((minAbsElem > absCurr) || isZero(minAbsElem)) {
      minAbsElem = absCurr;
    }
    if (maxAbsElem < absCurr) {
      maxAbsElem = absCurr;
    }
  }

  // If we already have a sensible max or if the smallest element is really small...
  // Should we maybe check that it is numerically unstable?
  // If it has elements that are too small, we simply reject it
  // We can check this by ensuring the ratio minAbsElem / maxAbsElem is not too small
  if (lessThanVal(maxAbsElem, SENSIBLE_MAX)) {
    return 0;
  }
  if (isZero(maxAbsElem) || isZero(minAbsElem, 1000 * param.getEPS())
      || isZero(minAbsElem / maxAbsElem)) {
    return -1;
  }

  double norm = std::sqrt(totalSOS);
  double divFactor;

  double newAbsMin = minAbsElem / norm; // Do not want this to be too small
  if (isZero(newAbsMin)) {
    divFactor = minAbsElem / (100 * param.getEPS());
  } else {
    divFactor = norm;
  }

  for (int i = 0; i < num_cols; i++) {
    (*rayToScale)[i] = (*rayToScale)[i] / divFactor;
  }
  return 1;
} /* scalePrimalRay */

/**
 * Compute score as in SCIP from effectiveness and orthogonalities
 */
void computeScore(AdvCuts& cs, const int pos) {
  cs.score[pos] = cs.cuts[pos].effectiveness() + cs.orthogonality[pos];
} /* computeScore */

/**
 * As in SCIP
 */
void updateOrthogonalities(AdvCuts& cs, const OsiRowCut& cut,
    const double min_cut_orthogononality,
    const std::vector<int>& roundWhenCutApplied) {
  for (int pos = 0; pos < cs.sizeCuts(); pos++) {
    if (roundWhenCutApplied[pos] >= 0) {
      continue;
    }
    // Update orthogonality
    const double this_ortho = getOrthogonality(cut.row(), cs.cuts[pos].row());
    if (this_ortho < cs.orthogonality[pos]) {
      if (this_ortho < min_cut_orthogononality) {
        /* cut is too parallel: release the row and delete the cut */
        //SCIPdebugMessage("    -> deleting parallel cut <%s> after adding <%s> (pos=%d, len=%d, orthogonality=%g, score=%g)\n",
        //SCIProwGetName(sepastore->cuts[pos]), SCIProwGetName(cut), pos, SCIProwGetNNonz(cut), thisortho, sepastore->scores[pos]);
        //SCIP_CALL( sepastoreDelCut(sepastore, blkmem, set, eventqueue, eventfilter, lp, pos) );
        cs.orthogonality[pos] = -1;
        computeScore(cs, pos); // we actually only penalize it a lot
      } else {
        /* recalculate score */
        cs.orthogonality[pos] = this_ortho;
        computeScore(cs, pos);
      }
    }
  }
} /* updateOrthogonalities */

/**
 * As in SCIP, returns the position of the best non-forced cut in the cuts array
 * Currently not forcing any cuts
 * */
int getBestCutIndex(const AdvCuts& cs,
    const std::vector<int>& roundWhenCutApplied) {
  int bestpos = -1;
  double bestscore = -1;
  for (int pos = 0; pos < cs.sizeCuts(); pos++) {
    // Check if cut has not been applied and is current best cut
    if ((roundWhenCutApplied[pos] < 0) && (cs.score[pos] > bestscore)) {
      bestscore = cs.score[pos];
      bestpos = pos;
    }
  }
  return bestpos;
} /* getBestCutIndex */

/**
 * Generic way to apply a set of cuts to a solver
 * @return Solver value after adding the cuts, if they were added successfully
 */
double applyCuts(PHASolverInterface* solver, const OsiCuts& cs,
    const int startCutIndex, const int numCutsOverload) {
  if (solver && !solver->isProvenOptimal()) {
    solver->resolve();
    checkSolverOptimality(solver, false);
  }
  if (solver == NULL || !solver->isProvenOptimal()) {
    return solver->getInfinity();
  }

  const int numCuts = (numCutsOverload >= 0) ? numCutsOverload : cs.sizeCuts() - startCutIndex;
  const double compare_obj = solver->getObjValue();

  // If the number of cuts is 0, we simply get the value without these cuts applied
  if (numCuts > 0) {
#ifdef TRACE
    OsiSolverInterface::ApplyCutsReturnCode code;
    int num_applied = 0;
#endif
    if ((startCutIndex == 0) && (numCuts == cs.sizeCuts())) {
#ifdef TRACE
      code = solver->applyCuts(cs);
      num_applied = code.getNumApplied();
#else
      solver->applyCuts(cs);
#endif
    } else {
#ifdef TRACE
      const int old_num_rows = solver->getNumRows();
#endif
      OsiRowCut* cuts;
      cuts = new OsiRowCut[numCuts];

      for (int i = startCutIndex; i < startCutIndex + numCuts; i++) {
        cuts[i - startCutIndex] = cs.rowCut(i);
      }
      solver->applyRowCuts(numCuts, cuts);
      delete[] cuts;
#ifdef TRACE
      num_applied = (solver->getNumRows() - old_num_rows);
#endif
    }

#ifdef TRACE
    printf("Cuts applied: %d (of the %d possible).\n",
        num_applied, numCuts);
#endif
  } /* check if num cuts > 0 */

  const bool useResolve = solver->isProvenOptimal();
  if (useResolve) {
    solver->resolve();
  } else {
    solver->initialSolve();
  }

  // At this point it is useful to check for integer feasibility quickly
	if (!solver->isProvenOptimal()) {
		// Sometimes need to resolve
		solver->resolve();
		if (!solver->isProvenOptimal()) {
			error_msg(errstr,
					"Not in subspace and solver not optimal after adding cuts, so we have an integer infeasible problem.\n");
			if (param.getEXIT_ON_INFEAS()) {
				writeErrorToII(errstr, GlobalVariables::log_file);
				exit(1);
			}
		}
	}

  // Sometimes there are slight inaccuracies at the end that we can get rid of
  if (solver->isProvenOptimal() && solver->getObjValue() < compare_obj) {
    solver->resolve();
    checkSolverOptimality(solver, false);
	}

  if (solver->isProvenOptimal()) {
    return solver->getObjValue();
  } else {
    return solver->getInfinity();  // Minimization problem, so this means infeas
  }
} /* applyCuts */

/**
 * Generic way to apply a set of cuts to a solver
 * @return Solver value after adding the cuts, if they were added successfully
 */
double applyCuts(PHASolverInterface* solver, const OsiCuts& cs,
    const int numCutsOverload) {
  const int numCuts = (numCutsOverload >= 0) ? numCutsOverload : cs.sizeCuts();
  return applyCuts(solver, cs, 0, numCuts);
} /* applyCuts */

/**
 * Generic way to apply a set of cuts to a solver
 * @return Solver value after adding the cuts, if they were added successfully
 */
double applyCutsInRounds(PHASolverInterface* solver, AdvCuts& cs,
    std::vector<int>& numCutsByRound, std::vector<double>& boundByRound,
    std::vector<double>& basisCond, const int numCutsToAddPerRound,
    std::vector<int>* const oldRoundWhenCutApplied,
    bool populateRoundWhenCutApplied, 
    const int maxRounds, double compare_obj) {
  if (solver && !solver->isProvenOptimal()) {
    solver->resolve();
    checkSolverOptimality(solver, false);
  }
  if (solver == NULL || !solver->isProvenOptimal()) {
    return solver->getInfinity();
  }
  
  const double initObjValue = solver->getObjValue();
  if (isNegInfinity(compare_obj)) {
    compare_obj = initObjValue;
  }

  const int numCuts = cs.sizeCuts();
  const int numRounds =
      std::ceil((double) numCuts / numCutsToAddPerRound) <= maxRounds ?
          std::ceil((double) numCuts / numCutsToAddPerRound) : maxRounds + 1;
  populateRoundWhenCutApplied = populateRoundWhenCutApplied
      && (oldRoundWhenCutApplied != NULL);
  const bool updateScores = (oldRoundWhenCutApplied == NULL)
      || populateRoundWhenCutApplied;
  if (populateRoundWhenCutApplied) {
    numCutsByRound.resize(numRounds);
  }
  basisCond.resize(numRounds);
  boundByRound.resize(numRounds);

#ifdef TRACE
    int numOrigRows = solver->getNumRows();
#endif

  OsiRowCut* cuts = NULL;
//  std::vector<int> sortedCutIndex(numCuts);
  if (numCuts > 0) {
    cuts = new OsiRowCut[numCuts];

//    if (numRounds > 1) { // Sort and compute scores with initialized orthogonalities
    if (updateScores) {
      cs.score.clear();
      cs.score.resize(numCuts);
      cs.orthogonality.clear();
      cs.orthogonality.resize(numCuts, 1.); // initialized orthogonality
    }
    for (int pos = 0; pos < numCuts; pos++) {
      if (updateScores) {
        computeScore(cs, pos);
      }
//        sortedCutIndex[pos] = pos;
    }
//      std::sort(sortedCutIndex.begin(), sortedCutIndex.end(),
//          index_cmp<const std::vector<double>&>(cs.score));
//    }

    for (int i = 0; i < numCuts; i++) {
//      cuts[i] = cs.rowCut(sortedCutIndex[i]);
      cuts[i] = cs.rowCut(i);
    }
  }

  // If the number of cuts is 0, we simply get the value without these cuts applied
  //const std::vector<OsiRowCut> cutVec = cs.cuts;
  std::vector<int> roundWhenCutApplied(numCuts, -1);
  for (int round_ind = 0; round_ind < numRounds; round_ind++) {
    // Re-start scores
    if (updateScores) {
      for (int pos = 0; pos < numCuts; pos++) {
        if (roundWhenCutApplied[pos] < 0) {
          cs.orthogonality[pos] = 1.; // re-start orthogonality
          computeScore(cs, pos);
        }
      }
    }
    const int start_ind = round_ind * numCutsToAddPerRound;
    const int numCutsToAddInThisRound =
        (round_ind < numRounds - 1) ?
            numCutsToAddPerRound : numCuts - start_ind;

//    if (numRounds > 1) {
    if (updateScores) {
      if (round_ind < numRounds - 1) {
        //std::vector<int> indicesAppliedThisRound;
        //indicesAppliedThisRound.reserve(numCutsToAddInThisRound);
        for (int i = 0; i < numCutsToAddInThisRound; i++) {
          const int pos = getBestCutIndex(cs, roundWhenCutApplied);
          if (pos >= 0) {
            solver->applyRowCuts(1, &(cuts[pos]));
            roundWhenCutApplied[pos] = round_ind;
            numCutsByRound[round_ind]++;
            if (updateScores) {
              updateOrthogonalities(cs, cs.cuts[pos], 0, roundWhenCutApplied);
            }
            //indicesAppliedThisRound.push_back(pos);
          } else {
            break;
          }
        }

        // 2018-09-16: update effectiveness wrt new opt; need to un-update at end
        for (int pos = 0; pos < (int) cs.size(); pos++) {
          if (roundWhenCutApplied[pos] >= 0) {
            continue;
          }
          //const int pos = indicesAppliedThisRound[i];
          const double violation = cs.cuts[pos].violated(solver->getColSolution());
          const double currCutNorm = cs.cuts[pos].getTwoNorm();
          cs.cuts[pos].setEffectiveness(violation / currCutNorm);
        }
      } else {
        for (int pos = 0; pos < numCuts; pos++) {
          if (roundWhenCutApplied[pos] < 0) {
            solver->applyRowCuts(1, &(cuts[pos]));
            roundWhenCutApplied[pos] = round_ind;
            numCutsByRound[round_ind]++;
            if (numCutsByRound[round_ind] == numCutsToAddInThisRound) {
              break;
            }
          } else {
            // Reset effectiveness
            cs.cuts[pos].setEffectiveness(cuts[pos].effectiveness());
          }
        }
      }
    } else {
      int numCutsAdded = 0;
      for (int pos = 0; pos < numCuts; pos++) {
        if ((*oldRoundWhenCutApplied)[pos] == round_ind) {
          solver->applyRowCuts(1, &(cuts[pos]));
          numCutsAdded++;
          if (numCutsAdded == numCutsToAddInThisRound) {
            break;
          }
        }
      }
    }
//    } else {
//      solver->applyRowCuts(numCutsToAddInThisRound, &(cuts[start_ind]));
//    }

    const bool useResolve = solver->isProvenOptimal();
    if (useResolve) {
      solver->resolve();
    } else {
      solver->initialSolve();
    }
    checkSolverOptimality(solver, false);

    // At this point it is useful to check for integer feasibility quickly
    if (!checkSolverOptimality(solver, false)) {
      error_msg(errstr,
          "Solver not optimal after adding cuts, so we have an integer infeasible problem. Exit status: %d.\n", solver->getModelPtr()->status());
      if (param.getEXIT_ON_INFEAS()) {
        if (numCuts && cuts) {
          delete[] cuts;
        }
        writeErrorToII(errstr, GlobalVariables::log_file);
        exit(1);
      }
    }
     
    // Sometimes there are slight inaccuracies at the end that we can get rid of
    if (solver->isProvenOptimal() && solver->getObjValue() < compare_obj) {
      solver->resolve();
      checkSolverOptimality(solver, false);
    }

    if (solver->isProvenOptimal()) {
      boundByRound[round_ind] = solver->getObjValue();
      //basisCond[round_ind] = compute_condition_number_norm2(solver);
    } else {
      for (int round_ind2 = round_ind; round_ind2 < numRounds; round_ind2++) {
        boundByRound[round_ind2] = solver->getInfinity(); // minimization problem, so this means infeas
        basisCond[round_ind2] = -1;
      }
      if (numCuts && cuts) {
        delete[] cuts;
      }
      return (solver->getInfinity()); // minimization problem, so this means infeas
    }
  } /* end iterating over rounds */

  if (numCuts && cuts) {
    delete[] cuts;
  }
  if (populateRoundWhenCutApplied) {
    for (int pos = 0; pos < numCuts; pos++) {
      (*oldRoundWhenCutApplied)[pos] = roundWhenCutApplied[pos];
    }
  }

#ifdef TRACE
  printf("Cuts applied: %d (of the %d possible).\n",
      solver->getNumRows() - numOrigRows, numCuts);
#endif
  return (solver->getObjValue());
} /* applyCutsInRounds */

/**
 * Get bound from applying
 * 1. SICs
 * 2. PHA
 * 3. SICs + PHA
 */
void compareCuts(const int num_obj_tried,
    const OsiSolverInterface* const solver, AdvCuts& structSICs,
    AdvCuts& structPHAs, const int numCutsToAddPerRound, std::string& output) {
#ifdef TRACE
  printf("\n## Comparing cuts. ##\n");
#endif
  const int numSICs = structSICs.sizeCuts();
  const int numPHA = structPHAs.sizeCuts();
  double sic_opt = 0.0, pha_opt = 0.0;
  double sic_pha_opt = 0.0;
  std::vector<int> numCutsByRoundSIC, numCutsByRoundPHA;
  std::vector<double> boundByRoundSIC, boundByRoundPHA;
  std::vector<double> boundByRoundPHAPlus;
  std::vector<double> basisCondByRoundSIC, basisCondByRoundPHA;
  std::vector<double> basisCondByRoundPHAPlus;
  int max_rounds = 0;
  double minSICOrthogonalityInit = -1.;
  double maxSICOrthogonalityInit = -1.;
  double avgSICOrthogonalityInit = -1.;
  double minPHAOrthogonalityInit = -1.;
  double maxPHAOrthogonalityInit = -1.;
  double avgPHAOrthogonalityInit = -1.;
  double minSICOrthogonalityFinal = -1.;
  double maxSICOrthogonalityFinal = -1.;
  double avgSICOrthogonalityFinal = -1.;
  double minPHAOrthogonalityFinal = -1.;
  double maxPHAOrthogonalityFinal = -1.;
  double avgPHAOrthogonalityFinal = -1.;

  // SICS
#ifdef TRACE
  if (numSICs > 0) {
    printf("\nGetting SICs bound.\n");
  } else {
    printf("\nNo SICs to apply.\n");
  }
#endif
  GlobalVariables::timeStats.start_timer(APPLY_MSICS_TIME);
  PHASolverInterface* SICSolver;
  SICSolver = dynamic_cast<PHASolverInterface*>(solver->clone());
  setMessageHandler(SICSolver);
  sic_opt = applyCuts(SICSolver, structSICs);
  GlobalVariables::timeStats.end_timer(APPLY_MSICS_TIME);

  // PHA
  GlobalVariables::timeStats.start_timer(APPLY_MPHA_TIME);
  if (numPHA > 0) {
#ifdef TRACE
    printf("\nGetting PHA bound.\n");
#endif
    PHASolverInterface* PHASolver;
    PHASolver = dynamic_cast<PHASolverInterface*>(solver->clone());
    setMessageHandler(PHASolver);
    std::vector<int> roundWhenCutAppliedPHA(numPHA);
    pha_opt = applyCutsInRounds(PHASolver, structPHAs, numCutsByRoundPHA,
        boundByRoundPHA, basisCondByRoundPHA, numCutsToAddPerRound,
        &(roundWhenCutAppliedPHA), true);

    // Get max orthogonality for round 1
    getOrthogonalityStatsForRound(structPHAs,
        minPHAOrthogonalityInit, maxPHAOrthogonalityInit,
        avgPHAOrthogonalityInit, &roundWhenCutAppliedPHA, 0);
    // Get max orthogonality for all PHA cuts
    getOrthogonalityStatsForRound(structPHAs,
        minPHAOrthogonalityFinal, maxPHAOrthogonalityFinal,
        avgPHAOrthogonalityFinal);

    if (PHASolver) {
      delete PHASolver;
    }
    
    if (numSICs > 0) {
#ifdef TRACE
      printf("\nGetting SIC+PHA bound.\n");
#endif
    }
    sic_pha_opt = applyCutsInRounds(SICSolver, structPHAs,
        numCutsByRoundPHA, boundByRoundPHAPlus, basisCondByRoundPHAPlus,
        numCutsToAddPerRound, &(roundWhenCutAppliedPHA), false);
    if (max_rounds < (int) boundByRoundPHAPlus.size()) {
      max_rounds = boundByRoundPHAPlus.size();
    }

    //applyCuts(SICSolver, structPHAs);
  } else {
    sic_pha_opt = sic_opt;
  }
  GlobalVariables::timeStats.end_timer(APPLY_MPHA_TIME);

  // Get number rounds of SICs needed to meet bound from GICs+SICs
#ifdef TRACE
  printf("\nGetting number rounds of SICs req'd to get bound.\n");
#endif
  int total_num_sics = 0;
  double final_sic_bound = 0.;
  const int min_sic_rounds = 0;
  int num_sic_rounds = 0;
  if (numSICs > 0) {
    PHASolverInterface* roundSICSolver = dynamic_cast<PHASolverInterface*>(solver->clone());
    setMessageHandler(roundSICSolver);
    double curr_sic_opt = applyCuts(roundSICSolver, structSICs);
    num_sic_rounds++;
    boundByRoundSIC.push_back(curr_sic_opt);
    //basisCondByRoundSIC.push_back(compute_condition_number_norm2(roundSICSolver));
    numCutsByRoundSIC.push_back(numSICs);
    getOrthogonalityStatsForRound(structSICs, minSICOrthogonalityInit,
        maxSICOrthogonalityInit, avgSICOrthogonalityInit);
    total_num_sics += numSICs;
    OsiCuts allSICs = structSICs;
    while (num_sic_rounds < min_sic_rounds || 
        (lessThanVal(curr_sic_opt, sic_pha_opt)
          && total_num_sics < numPHA)) {
#ifdef TRACE
      printf("\n##### Starting round %d of SIC generation #####\n", num_sic_rounds+1);
#endif
      OsiCuts currStructSICs;
      CglGICParam tmpParam(param);
      tmpParam.setParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND, 0); // will cause number cuts added per round to not be the same
      CglSIC SICGen(tmpParam);
      SICGen.generateCuts(*roundSICSolver, currStructSICs);

      // Check if no SICs generated
      if (currStructSICs.sizeCuts() == 0) {
        break;
      }

      num_sic_rounds++;
      curr_sic_opt = applyCuts(roundSICSolver, currStructSICs);
      boundByRoundSIC.push_back(curr_sic_opt);
      //basisCondByRoundSIC.push_back(compute_condition_number_norm2(roundSICSolver));
      numCutsByRoundSIC.push_back(currStructSICs.sizeCuts());
      total_num_sics += currStructSICs.sizeCuts();
      allSICs.insert(currStructSICs);

      // Other stopping conditions:
      // Bound does not improve at all after one round
      if (isVal(curr_sic_opt, boundByRoundSIC[num_sic_rounds - 2])) {
        break;
      }

      // Bound does not significantly improve after five rounds
      if (num_sic_rounds > 4) {
        const double delta = curr_sic_opt - boundByRoundSIC[num_sic_rounds - 4];
        if (!greaterThanVal(delta, 1e-3)) {
          break;
        }
      }
    } /* one round of Gomory cuts */
    final_sic_bound = roundSICSolver->getObjValue();
    getOrthogonalityStatsForRound(allSICs, minSICOrthogonalityFinal,
        maxSICOrthogonalityFinal, avgSICOrthogonalityFinal);
    if (roundSICSolver) {
      delete roundSICSolver;
    }
    if (max_rounds < num_sic_rounds) {
      max_rounds = boundByRoundSIC.size();
    }
  } /* add Gomory cuts in rounds */
      
  // After adding all cuts, count active cuts
//#ifdef TRACE
//  printf("\nCounting number of active cuts.\n");
//#endif
  std::vector<int> active_cut_heur_sicsolver(CutHeuristics::NUM_CUT_HEUR, 0);
  int num_active_pha_sicsolver = 0;

  if (!isInfinity(pha_opt)) {
    // Count number of active SICs
    for (int cut_ind = 0; cut_ind < structSICs.sizeCuts(); cut_ind++) {
			OsiRowCut* currCut = structSICs.rowCutPtr(cut_ind);
      const double activity = currCut->row().dotProduct(
          SICSolver->getColSolution());
      if (isVal(activity, currCut->rhs())) {
        active_cut_heur_sicsolver[CutHeuristics::SIC_CUT_GEN]++;
      }
    }

    // Count number of active PHA
    for (int cut_ind = 0; cut_ind < structPHAs.sizeCuts(); cut_ind++) {
      if (numSICs > 0) {
        const double activity = structPHAs[cut_ind].row().dotProduct(
            SICSolver->getColSolution());
        if (isVal(activity, structPHAs[cut_ind].rhs())) {
          num_active_pha_sicsolver++;
          active_cut_heur_sicsolver[structPHAs[cut_ind].cutHeur]++;
        }
      }
    }
  }

  const double lp_opt = solver->getObjValue();
  char tmpstring[300];
  sprintf(tmpstring, "\n## Results from adding cuts: ##\n");
  output += tmpstring;
  const int NAME_WIDTH = 25;
  const int NUM_DIGITS_BEFORE_DEC = 7;
  const int NUM_DIGITS_AFTER_DEC = 7;
  const double ip_opt = GlobalVariables::bestObjValue;
  const bool calcGapClosed = !isInfinity(std::abs(ip_opt));
  const int NUM_DIGITS_BEFORE_DEC_GAP_CLOSED = calcGapClosed ? 3 : 0;
  const int NUM_DIGITS_AFTER_DEC_GAP_CLOSED = calcGapClosed ? 3 : 0;
  char gapstring[300];
  sprintf(gapstring, " ");
  sprintf(tmpstring, "%-*.*s% -*.*f\n", NAME_WIDTH, NAME_WIDTH,
      "Original solver: ", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC, lp_opt);
  output += tmpstring;
  if (calcGapClosed) {
    sprintf(tmpstring, "%-*.*s% -*.*f\n", NAME_WIDTH, NAME_WIDTH,
        "IP: ", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC,
        ip_opt);
    output += tmpstring;
  }
  if (numSICs > 0 && !isInfinity(sic_opt)) {
    if (calcGapClosed) {
      snprintf(gapstring, sizeof(gapstring) / sizeof(char), " [%*.*f%%] ",
          NUM_DIGITS_BEFORE_DEC_GAP_CLOSED, NUM_DIGITS_AFTER_DEC_GAP_CLOSED,
          100. * (sic_opt - lp_opt) / (ip_opt - lp_opt));
    }
    sprintf(tmpstring, "%-*.*s% -*.*f%s(%d cuts)\n", NAME_WIDTH, NAME_WIDTH,
        "SICs: ", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC, sic_opt,
        gapstring, numSICs);
    output += tmpstring;
  }
  if (numPHA > 0 && !isInfinity(pha_opt)) {
    if (calcGapClosed) {
      snprintf(gapstring, sizeof(gapstring) / sizeof(char), " [%*.*f%%] ",
          NUM_DIGITS_BEFORE_DEC_GAP_CLOSED, NUM_DIGITS_AFTER_DEC_GAP_CLOSED,
          100. * (pha_opt - lp_opt) / (ip_opt - lp_opt));
    }
    sprintf(tmpstring, "%-*.*s% -*.*f%s(%d cuts)\n", NAME_WIDTH, NAME_WIDTH,
        "PHA: ", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC, pha_opt,
        gapstring, numPHA);
    output += tmpstring;
  }
//  if (isInfinity(pha_opt)) {
//    sprintf(tmpstring, "%-*.*s+inf (%d cuts)\n", NAME_WIDTH, NAME_WIDTH,
//        "*** All but SICs: ", numPHA);
//    output += tmpstring;
//  } else {
//    sprintf(tmpstring, "%-*.*s% -*.*f (%d cuts)\n", NAME_WIDTH, NAME_WIDTH,
//        "*** All but SICs: ", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC,
//        pha_opt, numPHA);
//    output += tmpstring;
//  }

  if (numSICs > 0) {
    if (isInfinity(sic_pha_opt)) {
      sprintf(tmpstring, "%-*.*s+inf (%d cuts", NAME_WIDTH, NAME_WIDTH,
          "*** All cuts with SICs:: ", numSICs + numPHA);
      output += tmpstring;
    } else {
      if (calcGapClosed) {
        snprintf(gapstring, sizeof(gapstring) / sizeof(char), " [%*.*f%%] ",
            NUM_DIGITS_BEFORE_DEC_GAP_CLOSED, NUM_DIGITS_AFTER_DEC_GAP_CLOSED,
            100. * (sic_pha_opt - lp_opt) / (ip_opt - lp_opt));
      }
      sprintf(tmpstring, "%-*.*s% -*.*f%s(%d cuts", NAME_WIDTH, NAME_WIDTH,
          "*** All cuts with SICs: ", NUM_DIGITS_BEFORE_DEC,
          NUM_DIGITS_AFTER_DEC, sic_pha_opt, gapstring, numSICs + numPHA);
      output += tmpstring;
    }
    if (numSICs > 0) {
      sprintf(tmpstring, ", %d / %d active SICs", active_cut_heur_sicsolver[CutHeuristics::SIC_CUT_GEN], numSICs);
      output += tmpstring;
    }
    if (numPHA > 0) {
      sprintf(tmpstring, ", %d / %d active PHA", num_active_pha_sicsolver, numPHA);
      output += tmpstring;
    }
    sprintf(tmpstring, ")\n");
    output += tmpstring;
  }

  // BOUND BY ROUND
  if (max_rounds > 0) {
    printf("\n## Bound by round ##\n");
    const int COL_WIDTH = 10
        + sizeof(stringValue(solver->getObjValue()).c_str()) / sizeof(char);
    printf("%-*.*s", 6, 6, "Round");
    if (numSICs > 0) {
      printf("%-*.*s", COL_WIDTH, COL_WIDTH, "SIC");
    }
    if (numPHA > 0) {
      printf("%-*.*s", COL_WIDTH, COL_WIDTH, "PHA");
    }
    if (numPHA > 0 && numSICs > 0) {
      printf("%-*.*s", COL_WIDTH, COL_WIDTH, "PHA+");
    }
    printf("\n");
    for (int round_ind = 0; round_ind < max_rounds; round_ind++) {
      printf("%-*.*s", 6, 6, stringValue(round_ind).c_str());
      if (round_ind < (int) boundByRoundSIC.size()) {
        const std::string tmpstring = stringValue(boundByRoundSIC[round_ind],
            "%.2f") + "(" + stringValue(numCutsByRoundSIC[round_ind]) + ")";
        printf("%-*.*s", COL_WIDTH, COL_WIDTH, tmpstring.c_str());
      } else if (numSICs > 0) {
        printf("%-*.*s", COL_WIDTH, COL_WIDTH, "");
      }
      if (round_ind < (int) boundByRoundPHA.size()) {
        printf("%-*.*s", COL_WIDTH, COL_WIDTH,
            stringValue(boundByRoundPHA[round_ind], "%.2f").c_str());
      } else if (numPHA > 0) {
        printf("%-*.*s", COL_WIDTH, COL_WIDTH, "");
      }
      if (round_ind < (int) boundByRoundPHAPlus.size()) {
        printf("%-*.*s", COL_WIDTH, COL_WIDTH,
            stringValue(boundByRoundPHAPlus[round_ind], "%.2f").c_str());
      } else if (numPHA > 0 && numSICs > 0) {
        printf("%-*.*s", COL_WIDTH, COL_WIDTH, "");
      }
      printf("\n");
    }
  }

  // Also print if we improved
  // Recall everything has been standardized to be a minimization problem
  // Thus, larger objective means we have closed more gap
  printf("\n## Improvements: ##\n");
  bool PHAOverSICs = greaterThanVal(pha_opt, sic_opt);
  bool allOverSICs = greaterThanVal(sic_pha_opt, sic_opt);
  if (numSICs == 0) {
    printf("No SICs were generated.\n");
  } else if (PHAOverSICs || allOverSICs) {
    if (PHAOverSICs) {
      printf("New cuts > SICs by %f.\n", pha_opt - sic_opt);
    }
    if (allOverSICs) {
      printf("All cuts (with SICs) > SICs by %f.\n",
          sic_pha_opt - sic_opt);
    }
  } else {
    printf("No improvements over SICs.\n");
  }

  if ((param.getParamVal(ParamIndices::PHA_PARAM_IND) > 0)) {
    printf("\n## Cut summary: ##\n");
    // Report how many objectives were tried
    printf(
        "Generated %d cuts from %d objectives tried "
            "(on the %d cut-generating sets point-ray collections).\n",
        numPHA,
        num_obj_tried
            - GlobalVariables::numCutSolverFails[CutSolverFails::PRIMAL_CUTSOLVER_FAIL_IND],
        param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND));
  }

  // Also print which heuristic had active cuts
  printf("\n## Active cuts from which heuristics: ##\n");
  for (int heur_ind = 0; heur_ind < CutHeuristics::NUM_CUT_HEUR; heur_ind++) {
    printf("%s:", CutHeuristicsName[heur_ind].c_str());
    if (numSICs > 0) {
      printf(" [%d]", active_cut_heur_sicsolver[heur_ind]);
    }
    if (numSICs == 0) {
      printf(" 0");
    }
    printf(" (of %d) \n", GlobalVariables::numCutsFromHeur[heur_ind]);
  }

  // Save all this to ii
  writeCutHeurInfoToII(numSICs + numPHA,
      num_obj_tried, active_cut_heur_sicsolver, 
      GlobalVariables::log_file);
  writeRootBoundToII(solver->getObjValue(), numSICs, sic_opt, 
      numPHA, pha_opt, sic_pha_opt, 
      active_cut_heur_sicsolver[CutHeuristics::SIC_CUT_GEN],
      num_active_pha_sicsolver,
      GlobalVariables::log_file);

  // Save the round info to ii as well
  writeRoundInfoToII(sic_pha_opt, total_num_sics, num_sic_rounds,
      final_sic_bound, boundByRoundPHAPlus,
//      compute_condition_number_norm2(solver),
//      basisCondByRoundSIC, basisCondByRoundPHA, basisCondByRoundPHAPlus,
      minSICOrthogonalityInit, maxSICOrthogonalityInit, avgSICOrthogonalityInit, 
      minPHAOrthogonalityInit, maxPHAOrthogonalityInit, avgPHAOrthogonalityInit,
      minSICOrthogonalityFinal, maxSICOrthogonalityFinal, avgSICOrthogonalityFinal, 
      minPHAOrthogonalityFinal, maxPHAOrthogonalityFinal, avgPHAOrthogonalityFinal, 
      GlobalVariables::log_file);

  // Free
  if (SICSolver) {
    delete SICSolver;
  }
} /* compareCuts */
