//============================================================================
// Name        : Utility.cpp
// Author      : akazachk
// Version     : 2018.11
// Description : Used for functions that are useful but not particular to any
//               file/class/etc.
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "GlobalConstants.hpp"

#include "Utility.hpp"
#include <sys/stat.h>

/** Overload solve from hot start because of issues */
bool solveFromHotStart(OsiSolverInterface* const solver, const int col,
    const bool isChangedUB, const double origBound, const double newBound) {
  solver->solveFromHotStart();
  if (solver->isIterationLimitReached()) {
    // This sometimes happens, e.g., with arki001
    solver->unmarkHotStart();
    //solver->disableFactorization();
    solver->resolve();
    if (isChangedUB) {
      solver->setColUpper(col, origBound);
    } else {
      solver->setColLower(col, origBound);
    }
    solver->resolve();
    //solver->enableFactorization();
    solver->markHotStart();
    if (isChangedUB) {
      solver->setColUpper(col, newBound);
    } else {
      solver->setColLower(col, newBound);
    }
    solver->solveFromHotStart();
  }
  return (solver->isProvenOptimal());
} /* solveFromHotStart overload */

/**
 * Calculate strong branching lower bound
 */
double calculateStrongBranchingLB(double& worstOptVal, double& bestOptVal,
    OsiSolverInterface* solver, const SolutionInfo& probData) {
  const double* colSolution = solver->getColSolution();
  double strong_branching_lb = solver->getInfinity();
  worstOptVal = solver->getInfinity();
  bestOptVal = -1 * solver->getInfinity();

  try {
		OsiClpSolverInterface* clpsolver = dynamic_cast<OsiClpSolverInterface*>(solver);
		const int hot_start_iter_limit = std::numeric_limits<int>::max();
		clpsolver->setIntParam(OsiMaxNumIterationHotStart, hot_start_iter_limit);
		clpsolver->setSpecialOptions(16); // use standard strong branching rather than clp's
  } catch (std::exception& e) {
    // It's okay, we can continue
  }

  solver->enableFactorization();
  solver->markHotStart();
  //for (int col = 0; col < solver->getNumCols(); col++) {
  //int cols[2] = {877,924};
  //for ( auto col : cols ) {
	for (auto col : probData.fractionalCore) {
		const double origLB = solver->getColLower()[col];
		const double origUB = solver->getColUpper()[col];
		double optValDown, optValUp;
		bool downBranchFeasible = true, upBranchFeasible = true;

		// Check down branch
		solver->setColUpper(col, std::floor(colSolution[col]));
		solveFromHotStart(solver, col, true, origUB, std::floor(colSolution[col]));
		if (solver->isProvenOptimal()) {
			optValDown = solver->getObjValue();
		} else if (solver->isProvenPrimalInfeasible()) {
			downBranchFeasible = false;
		} else {
			// Something strange happened
			error_msg(errorstring,
					"Down branch is neither optimal nor primal infeasible on variable %d (value %e).\n",
					col, colSolution[col]);
			writeErrorToII(errorstring, GlobalVariables::log_file);
			exit(1);
		}
		solver->setColUpper(col, origUB);

		// Return to previous state
		solver->solveFromHotStart();

		// Check up branch
		solver->setColLower(col, std::ceil(colSolution[col]));
		solveFromHotStart(solver, col, true, origUB, std::ceil(colSolution[col]));
		if (solver->isProvenOptimal()) {
			optValUp = solver->getObjValue();
		} else if (solver->isProvenPrimalInfeasible()) {
			upBranchFeasible = false;
		} else {
			// Something strange happened
			error_msg(errorstring,
					"Up branch is neither optimal nor primal infeasible on variable %d (value %e).\n",
					col, colSolution[col]);
			writeErrorToII(errorstring, GlobalVariables::log_file);
			exit(1);
		}
		solver->setColLower(col, origLB);

		// Return to original state
		solver->solveFromHotStart();

		// Check if some side of the split is infeasible
		if (!downBranchFeasible || !upBranchFeasible) {
			if (!downBranchFeasible && !upBranchFeasible) {
				// Infeasible problem
				error_msg(errorstring,
						"Infeasible problem due to integer variable %d (value %e).\n",
						col, colSolution[col]);
				writeErrorToII(errorstring, GlobalVariables::log_file);
				exit(1);
			}
		} else {
			// Update strong branching value, if appropriate
			const double lowerOptVal = (optValDown < optValUp) ? optValDown : optValUp;
			if (lowerOptVal < strong_branching_lb) {
				strong_branching_lb = lowerOptVal;
			}
			if (lowerOptVal > bestOptVal) {
			  bestOptVal = lowerOptVal;
			}
		}
	} /* get strong branching lb */
  solver->unmarkHotStart();
  solver->disableFactorization();

  worstOptVal = strong_branching_lb;
  return strong_branching_lb;
} /* calculateStrongBranchingLB */

/**
 * @brief Checks whether a solver is optimal and resolves a few times if there are some infeasibilities; doRedCostFixing is disabled no matter what
 */
bool checkSolverOptimality(OsiSolverInterface* const solver,
    const bool exitOnDualInfeas, const double timeLimit,
    const int maxNumResolves) {
  OsiClpSolverInterface* clpsolver = NULL;
  try {
    clpsolver = dynamic_cast<OsiClpSolverInterface*>(solver);
  } catch (std::exception& e) {
    // Disregard
  }

  // If it is supposedly proven primal infeasible, might be good to do a resolve first
  // This is, e.g., something that arises with miplib2003/pp08aCUTS_presolved -32
  int resolve_count = 0;
  if (clpsolver->isProvenPrimalInfeasible()) {
    clpsolver->getModelPtr()->setNumberIterations(0);
    clpsolver->getModelPtr()->setMaximumSeconds(timeLimit);
    clpsolver->resolve();
    resolve_count++;
  }

  // First clean
  if (resolve_count < maxNumResolves && clpsolver) {
    int status = clpsolver->getModelPtr()->secondaryStatus();
    bool is_cleaned = false;
    if (status == 2 || !isZero(clpsolver->getModelPtr()->sumPrimalInfeasibilities())) {
      clpsolver->getModelPtr()->cleanup(1);
      is_cleaned = true;
    } else if (status == 3 || status == 9 || !isZero(clpsolver->getModelPtr()->sumDualInfeasibilities())) {
      clpsolver->getModelPtr()->cleanup(2);
      is_cleaned = true;
    } else if (status == 4) {
      clpsolver->getModelPtr()->cleanup(1);
      is_cleaned = true;
    }
    if (is_cleaned) {
      //      clpsolver->getModelPtr()->setNumberIterations(0);
      //      clpsolver->getModelPtr()->setMaximumSeconds(timeLimit);
      //      clpsolver->resolve();
      resolve_count++;
//      return checkSolverOptimality(solver, exitOnDualInfeas, timeLimit,
//          doRedCostFixing, maxNumResolves - resolve_count);
      return checkSolverOptimality(solver, exitOnDualInfeas, timeLimit, 0);
    }

    // Do some resolves, but first save whether the initial status is dual infeasible
    // The reason is that for neos15 w/str=2, we were getting a dual infeasible problem
    // turn into a primal infeasible problem after resolves
    // (probably the strengthening causes numerical issues)
    const bool oldStatusDualInfeasible = (clpsolver->isProvenDualInfeasible());
    const bool oldStatusOptimal = (clpsolver->isProvenOptimal());
    bool resolve = true;
    double infeas = std::numeric_limits<double>::max();
    while (resolve && resolve_count < maxNumResolves) {
      const double curr_infeas =
          clpsolver->getModelPtr()->sumPrimalInfeasibilities()
              + clpsolver->getModelPtr()->sumDualInfeasibilities();
      resolve = (resolve_count == 0)
          || (!isZero(curr_infeas) && (curr_infeas < infeas));
      if (resolve) {
        clpsolver->getModelPtr()->setNumberIterations(0);
        clpsolver->getModelPtr()->setMaximumSeconds(timeLimit);
        clpsolver->resolve();
        resolve_count++;
      }
      infeas = curr_infeas;
    }

    // If dual infeas -> primal infeas, do a dual resolve, which should fix the issue
    if (oldStatusDualInfeasible && clpsolver->isProvenPrimalInfeasible()) {
      clpsolver->getModelPtr()->dual(2,0); // just do values pass
      if (!clpsolver->isProvenDualInfeasible() && !clpsolver->isProvenOptimal()) {
        clpsolver->getModelPtr()->setProblemStatus(2);
        resolve_count = maxNumResolves; // stop resolving
      }
    }

    if (oldStatusOptimal && clpsolver->isProvenPrimalInfeasible()) {
      clpsolver->enableFactorization();
      clpsolver->getModelPtr()->setNumberIterations(0);
      clpsolver->getModelPtr()->setMaximumSeconds(timeLimit);
      clpsolver->resolve();
      resolve_count++;
      clpsolver->disableFactorization();
    }

    // Clean once more if needed and possible
    status = clpsolver->getModelPtr()->secondaryStatus();
    is_cleaned = false;
    if (status == 2
        || !isZero(clpsolver->getModelPtr()->sumPrimalInfeasibilities())) {
      clpsolver->getModelPtr()->cleanup(1);
      is_cleaned = true;
    } else if (status == 3 || status == 9
        || !isZero(clpsolver->getModelPtr()->sumDualInfeasibilities())) {
      clpsolver->getModelPtr()->cleanup(2);
      is_cleaned = true;
    } else if (status == 4) {
      clpsolver->getModelPtr()->cleanup(1);
      is_cleaned = true;
    }
    if (is_cleaned) {
      resolve_count++;
      return checkSolverOptimality(solver, exitOnDualInfeas, timeLimit, 0);
    }
  }

  if (solver->isProvenPrimalInfeasible()) {
    return false;
  } else if (!(solver->isProvenOptimal())) {
    // Sometimes need to resolve once more to get the correct status
    if (resolve_count < maxNumResolves) {
      solver->resolve();
    }
    if (solver->isProvenPrimalInfeasible()) {
      return false;
    } else if (solver->isProvenDualInfeasible()) {
      if (exitOnDualInfeas) {
        error_msg(errstr,
            "Solver is dual infeasible. Check why this happened!\n");
        writeErrorToII(errstr, GlobalVariables::log_file);
        exit(1);
      } else {
        return false;
      }
    } else if (!(solver->isProvenOptimal())) {
      return false;
    }
  }

  return true;
} /* checkSolverOptimality */

bool enableFactorization(OsiSolverInterface* const solver, const int resolveFlag) {
  solver->enableFactorization();

  if (resolveFlag > 0) {
    try {
      OsiClpSolverInterface* clpsolver =
          dynamic_cast<OsiClpSolverInterface*>(solver);

      // After enabling factorization, things sometimes look different and it is worth resolving
      // This is motivated by seymour-disj-10; after adding one round of GMICs, row 5077 had negative row price
      if ((clpsolver->getModelPtr()->sumPrimalInfeasibilities() > param.getEPS())
          || (clpsolver->getModelPtr()->sumDualInfeasibilities()
              > param.getEPS())) {
        if (resolveFlag == 1) {
          clpsolver->initialSolve(); // hopefully this fixes things rather than breaks things; actually it breaks things, e.g., for rd2 SICs for coral/neos17, the variable basic in row 648 changes from 164 to 23
        } else {
          clpsolver->resolve();
        }
        clpsolver->disableFactorization();
        clpsolver->enableFactorization(); // this will hopefully solve the neos17 problem
      }
    } catch (std::exception& e) {
      // Disregard
    }
  }

  return solver->isProvenOptimal();
} /* enableFactorization */

/**
 * @brief Parses int from string using strtol
 */
bool parseInt(const char *str, int &val) {
  char *temp;
  bool rc = true;
  errno = 0;
  long tmpval = strtol(str, &temp, 0);

  if (temp == str || *temp != '\0' ||
      ((tmpval == std::numeric_limits<long>::min() || tmpval == std::numeric_limits<long>::max()) && errno == ERANGE))
    rc = false;

  val = static_cast<int>(tmpval);

  return rc;
}

/**
 * @brief Parses double from string using strtod
 */
bool parseDouble(const char *str, double &val) {
  char *temp;
  bool rc = true;
  errno = 0;
  val = strtod(str, &temp);

  if (temp == str || *temp != '\0' ||
      ((val == std::numeric_limits<double>::min()
          || val == std::numeric_limits<double>::max()) && errno == ERANGE))
    rc = false;

  return rc;
}

/**
 * @brief Universal way to check whether we reached the limit for the number of cuts for each split
 * This allows us to change between restricting number of cuts per split and total number of cuts easily
 */
bool reachedCutLimit(const int num_cuts_this_cgs, const int num_cuts_total) {
  // The cut limit is either across all cut-generating sets
  // or it is divided among the cut-generating sets (either as the limit / num_cgs, or as a fixed number per cgs)
  // If CUT_LIMIT = 0 => no cut limit
  // If CUT_LIMIT > 0 => absolute cut limit
  // If CUT_LIMIT < 0 => cut limit per cgs
  if (param.paramVal[CUT_LIMIT_PARAM_IND] == 0) {
    return false;
  } else if (param.paramVal[CUT_LIMIT_PARAM_IND] > 0) {
    return (num_cuts_total >= param.paramVal[CUT_LIMIT_PARAM_IND]);
  } else {
    const int num_splits = param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND);
    return (num_cuts_this_cgs >= -1. * param.paramVal[CUT_LIMIT_PARAM_IND] / num_splits);
  }
} /* reachedCutLimit */

/**
 * @brief Universal way to check whether we reached the limit for the time for some subroutine
 */
bool reachedTimeLimit(const Stats& timer, const std::string& timeName,
    const double max_time) {
  return (timer.get_total_time(timeName) > max_time);
} /* reachedTimeLimit */

/**
 * Two vectors are parallel iff u.v/|u|*|v| = 1
 */
double getOrthogonality(const CoinPackedVector& vec1,
    const CoinPackedVector& vec2) {
  const double scalarprod = sparseDotProduct(vec1, vec2);
  const double parallelism = std::abs(scalarprod)
      / (vec1.twoNorm() * vec2.twoNorm());
  return 1. - parallelism;
}

void getOrthogonalityStatsForRound(const OsiCuts& cs, double& minOrtho,
    double& maxOrtho, double& avgOrtho,
    const std::vector<int>* const roundWhenCutApplied, const int this_round) {
  const bool check_round = roundWhenCutApplied != NULL && this_round >= 0;
  int num_checked = 0;
  minOrtho = 2.;
  maxOrtho = -1.;
  avgOrtho = 0.;
  for (int pos = 0; pos < cs.sizeCuts(); pos++) {
    if (check_round && (*roundWhenCutApplied)[pos] != this_round) {
      continue;
    }
    for (int pos2 = pos + 1; pos2 < cs.sizeCuts(); pos2++) {
      if (check_round && (*roundWhenCutApplied)[pos2] != this_round) {
        continue;
      }
      const double this_ortho = getOrthogonality(cs.rowCutPtr(pos)->row(),
          cs.rowCutPtr(pos2)->row());
      if (this_ortho < minOrtho) {
        minOrtho = this_ortho;
      }
      if (this_ortho > maxOrtho) {
        maxOrtho = this_ortho;
      }
      avgOrtho += this_ortho;
      num_checked++;
    }
  }
  avgOrtho /= num_checked;
}

double getMinOrthogonality(const CoinPackedVector& vec, const OsiCuts& cs) {
  double min_ortho = 2.;
  for (int pos = 0; pos < cs.sizeCuts(); pos++) {
    const double this_ortho = getOrthogonality(cs.rowCutPtr(pos)->row(), vec);
    if (this_ortho < min_ortho) {
      min_ortho = this_ortho;
    }
  }
  return min_ortho;
}

/**
 * Compare two double arrays for equality
 * No effort is made to check that both have the same number of elements
 */
bool vectorsEqual(const int n, const double* vec1, const double* vec2) {
  for (int i = 0; i < n; i++) {
    if (!isZero(vec1[i] - vec2[i])) {
      return false;
    }
  }
  return true;
}

/********************************************************************************/
/**
 * @brief Decide if two rows are the same.
 * @return
 *   0: seem same in coeff and rhs,
 *   +/-1: seem same in coeff, but diff in rhs (-1: cut1 better, +1: cut2 better),
 *   2: seem different in coeff and rhs.
 *
 * The value of eps will be used to determine whether a cut coefficient is zero or not.
 * We are roughly checking that whether cut1 = ratio * cut2 for some ratio.
 * Let ratio = rhs1 / rhs2. (If both are non-zero. Typically both are = 1.)
 * We say the cuts are "truly different" if |coeff^1_j - ratio * coeff^2_j| >= diffeps.
 * However, might be that the cut coefficients are scaled versions, but the rhs values are different.
 * So we also compute ratio_coeff = coeff^1_j / coeff^2_j for the first j where both are non-zero.
 * If |coeff^1_j - ratio * coeff^2_j| >= diffeps for some j,
 * but |coeff^1_j - ratio_coeff * coeff^2_j| < diffeps for all j,
 * then the return code will be -1 or 1 (essentially, one cut will dominate the other).
 *
 * By the way, this is all essentially checking orthogonality...
 */
int isRowDifferent(const CoinPackedVector& cut1Vec, const double cut1rhs,
    const CoinPackedVector& cut2Vec, const double cut2rhs, const double eps) {
  const int numElem1 = cut1Vec.getNumElements();
  const int* index1 = cut1Vec.getIndices();
  const double* value1 = cut1Vec.getElements();

  const int numElem2 = cut2Vec.getNumElements();
  const int* index2 = cut2Vec.getIndices();
  const double* value2 = cut2Vec.getElements();

  double ratio = -1.0;
  double ratio_coeff = -1.0;
  double maxDiff = 0.;
  double maxDiffCoeff = 0.;
  double maxDiffZero = 0.;

  // First, we try to get the ratio from RHS
  if (isZero(cut1rhs, eps) || isZero(cut2rhs, eps)) {
    // If one but not both are 0, then the cuts are clearly different
    // But they may only differ in the rhs
    if (!isZero(cut1rhs, eps)) {
      maxDiff = std::abs(cut1rhs);
    } else if (!isZero(cut2rhs, eps)) {
      maxDiff = std::abs(cut2rhs);
    }
    // Otherwise, they're both near zero,
    // and we're not going to be able to get
    // any kind of ratio from here
  } else {
    ratio = cut1rhs / cut2rhs;
  }

  // Now check if this ratio (if we found one) holds for
  // all the coefficients. In particular, we check that
  // std::abs(cut1.coeff - ratio * cut2.coeff) < eps
  int ind2 = 0;
  for (int ind1 = 0; ind1 < numElem1; ind1++) {
    // If this entry is zero, then we cannot have a violation
    // (unless the corresponding entry of cut2 is non-zero,
    // but this will be detected later)
    //if (isZero(value1[ind1], eps)) {
    //  continue;
    //}

    // Find matching indices
    while (ind2 < numElem2 && index1[ind1] > index2[ind2]) {
      if (!isZero(value2[ind2], eps)) {
        return 2;
      }
      if (value2[ind2] > maxDiffZero) {
        maxDiffZero = std::abs(value2[ind2]);
      }
      ind2++;
    }

    // If we reached end of the second cut,
    // or index1[ind1] does not exist in the second cut,
    // capture maxDiffZero, or exit early if something "really different" detected
    if (ind2 >= numElem2 || index2[ind2] > index1[ind1]) {
      if (!isZero(value1[ind1], eps)) {
        return 2;
      }
      if (value1[ind1] > maxDiffZero) {
        maxDiffZero = std::abs(value1[ind1]);
      }
      continue;
    }

    // Either one or both of the coefficients may be zero,
    // in which case both better be zero.
    if (isZero(value1[ind1], eps) || isZero(value2[ind2], eps)) {
      if (!isZero(value1[ind1], eps) || !isZero(value2[ind2], eps)) {
        return 2; // truly different
      }
      /*
       if (isZero(value1[ind1], eps) && isZero(value2[ind2], eps)) {
       // Both are zero, nothing to see here
       } else if (isZero(value2[ind2], eps)) {
       if (value1[ind1] > maxDiffZero) {
       maxDiffZero = value1[ind1];
       }
       } else {
       if (value2[ind2] > maxDiffZero) {
       maxDiffZero = value2[ind1];
       }
       }
       */
      ind2++;
      continue;
    }

    // Alternatively, we make sure inequality is not scaled version of other
    // Set scale factor if it has not been set before
    if (lessThanVal(ratio, 0.0)) {
      ratio = value1[ind1] / value2[ind2];
      if (lessThanVal(ratio, 0.0)) {
        // The scale is negative
        // Potentially one cut is the negative of the other,
        // but they are thus still distinct
        // (all cuts are ax >= b)
        return 2;
      }
    } else {
      const double diff = std::abs(value1[ind1] - ratio * value2[ind2]);
      if (diff > maxDiff) {
        maxDiff = diff;
      }
    }

    // Also check the ratio_coeff in case one cut dominates the other
    if (lessThanVal(ratio_coeff, 0.0)) {
      ratio_coeff = value1[ind1] / value2[ind2];
      if (lessThanVal(ratio_coeff, 0.0)) {
        // The scale is negative
        // Potentially one cut is the negative of the other,
        // but they are thus still distinct
        // (all cuts are ax >= b)
        return 2;
      }
    } else {
      const double diff = std::abs(value1[ind1] - ratio_coeff * value2[ind2]);
      if (diff > maxDiffCoeff) {
        maxDiffCoeff = diff;
      }
    }

    if (!isZero(maxDiff, eps) && !isZero(maxDiffCoeff, eps)) {
      return 2;
    }

    ind2++;
  } /* iterate over cut1 coefficients */

  // Check any remaining cut2 indices that were not matched
  while (ind2 < numElem2) {
    if (!isZero(value2[ind2], eps)) {
      return 2;
    }
    if (value2[ind2] > maxDiffZero) {
      maxDiffZero = std::abs(value2[ind2]);
    }
    ind2++;
  }

  // Are the cuts scaled versions of one another?
  if (isZero(maxDiff, eps)) {
    return 0; // they seem the same
  }

  // Are the cut coefficients scaled versions of one another?
  if (isZero(maxDiffCoeff, eps)) {
    if (lessThanVal(ratio_coeff, 0., eps)) {
      return 2; // cuts face opposite directions
    }

    // Either one or the other rhs will be better
    const double rhs_diff = cut1rhs - ratio_coeff * cut2rhs;
    if (isZero(rhs_diff, eps)) {
      return 0; // they seem the same, not sure how ratio missed it
    } else if (greaterThanVal(rhs_diff, 0., eps)) {
      return -1;
    } else {
      return 1;
    }
  }

  return 0;
} /* isRowDifferent */

/**
 * @brief Set full std::vector coordinates from packed vector
 *
 * No effort is made to reset the indices that are not present in packed vector
 * so these should be initialized as desired prior to calling this function
 *
 * No check is performed for indices existing,
 * so this should be ensured before calling this function
 */
void setFullVecFromPackedVec(std::vector<double>& fullVec,
    const CoinPackedVector& vec) {
  const int numElem = vec.getNumElements();
  const int* vertIndex = vec.getIndices();
  const double* vertElem = vec.getElements();
  for (int el = 0; el < numElem; el++) {
    const int curr_ind = vertIndex[el];
    fullVec[curr_ind] = vertElem[el];
  }
}

/**
 * @brief Any indices that appear in vec will be put to zero in fullVec
 *
 * No effort is made to reset the indices that are not present in packed std::vector
 * so these should be initialized as desired prior to calling this function
 *
 * No check is performed for indices existing,
 * so this should be ensured before calling this function
 */
void clearFullVecFromPackedVec(std::vector<double>& fullVec,
    const CoinPackedVector& vec) {
  const int numElem = vec.getNumElements();
  const int* vertIndex = vec.getIndices();
  for (int el = 0; el < numElem; el++) {
    const int curr_ind = vertIndex[el];
    fullVec[curr_ind] = +0.0;
  }
}

void copySolver(PHASolverInterface* const solverCopy,
    const PHASolverInterface* const solverOld) {
  solverCopy->enableFactorization();
  // This is to prevent a different basis from being chosen for solverCopy, which sometimes happens.
  int* cstat = new int[solverOld->getNumCols()];
  int* rstat = new int[solverOld->getNumRows()];
  solverOld->getBasisStatus(cstat, rstat);
  solverCopy->setBasisStatus(cstat, rstat);
  delete[] cstat;
  delete[] rstat;
  solverCopy->disableFactorization();
  //  solverCopy->enableSimplexInterface(true); // This caused problems with non-negative variables getting negative values in the solution.
  if (!solverCopy->basisIsAvailable())
    solverCopy->resolve();
  solverCopy->enableFactorization(); // Seemingly is broken by the resolve or setBasisStatus
}

/**
 * @brief returns lb of col (recall complementing of slacks)
 */
double getVarLB(const OsiSolverInterface* const solver, const int col) {
//    const bool in_subspace_or_compl_space) {
  if (col < solver->getNumCols()) {
    return solver->getColLower()[col];
  } else {
    return 0.0;
  }
}

/**
 * @brief returns ub of col (recall complementing of slacks)
 */
double getVarUB(const OsiSolverInterface* const solver, const int col) {
//    const bool in_subspace_or_compl_space) {
  if (col < solver->getNumCols()) {
    return solver->getColUpper()[col];
  } else {
    return solver->getInfinity();
  }
}

/**
 * @brief returns value of col (recall complementing of slacks)
 */
double getVarVal(const PHASolverInterface* const solver, const int col) {
  if (col < solver->getNumCols()) {
    return solver->getColSolution()[col];
  } else {
    const int row = col - solver->getNumCols();
    // Absolute value since we assume all slacks are non-negative
    return std::abs(solver->getRowActivity()[row] - solver->getRightHandSide()[row]);
//    return (solver->getRightHandSide()[row] - solver->getRowActivity()[row]);
  }
}


/**
 * @brief returns value of col (recall complementing of slacks)
 */
double getCNBVarVal(const PHASolverInterface* const solver, const int col) {
  return getCNBVarVal(solver, col, isNonBasicUBVar(solver, col));
}

/**
 * @brief returns value of col (recall complementing of slacks)
 */
double getCNBVarVal(const PHASolverInterface* const solver, const int col, const bool complement) {
  double value;
  if (col < solver->getNumCols()) {
    value = solver->getColSolution()[col];
    const double LB = getVarLB(solver, col);
    const double UB = getVarUB(solver, col);

    if (complement) {
      value = UB - value;
    } else {
      value = value - LB;
    }
  } else {
    const int row = col - solver->getNumCols();
    // Absolute value since we assume all slacks are non-negative
    value = (solver->getRightHandSide()[row] - solver->getRowActivity()[row]);
    if (complement) {
      value = std::abs(value);
    }
  }
  return value;
}

/**
 * The variables that lead to rays are non-basic, non-fixed structural variables
 * and non-basic slacks not coming from equality constraints.
 */
bool isRayVar(const OsiClpSolverInterface* const solver, const int var) {
  if (var < solver->getNumCols())
    return (!isBasicCol(solver, var));
  else
    return (!isBasicSlack(solver, var - solver->getNumCols()));
}

/**
 * Returns the index in raysCutByHplane of the first ray that is cut by this hyperplane in this activation
 */
int getStartIndex(const Hplane& hplane, const int split_ind, const int act_ind) {
  return (act_ind > 0) ? hplane.numRaysCutByHplaneInclPrevAct[split_ind][act_ind - 1] : 0;
}

/**
 * Returns the index in raysCutByHplane of the last ray that is cut by this hyperplane in this activation
 */
int getEndIndex(const Hplane& hplane, const int split_ind, const int act_ind) {
  return hplane.numRaysCutByHplaneInclPrevAct[split_ind][act_ind];
}

bool duplicatePoint(const CoinPackedMatrix& interPtsAndRays,
    const std::vector<double>& rhs, const Point& newPoint) {
  const int numPtsAndRays = interPtsAndRays.getNumRows();
  for (int r = 0; r < numPtsAndRays; r++) {
    if (isZero(rhs[r]))
      continue;

    const CoinShallowPackedVector vec = interPtsAndRays.getVector(r);
    if (newPoint == vec) // This invokes the operator from Point.hpp
      return true;
  }
  return false;
}

bool duplicateRay(const CoinPackedMatrix& interPtsAndRays,
    const std::vector<double>& rhs, const Ray& newRay) {
  const int numPtsAndRays = interPtsAndRays.getNumRows();
  for (int r = 0; r < numPtsAndRays; r++) {
    if (!isZero(rhs[r]))
      continue;

    const CoinShallowPackedVector vec = interPtsAndRays.getVector(r);
    if (newRay == vec) // This invokes the operator from Ray.hpp
      return true;
  }
  return false;
}

int findDuplicateRay(const CoinPackedMatrix& interPtsAndRays,
    const std::vector<double>& rhs, const Ray& newRay) {
  const int numPtsAndRays = interPtsAndRays.getNumRows();
  for (int r = 0; r < numPtsAndRays; r++) {
    if (!isZero(rhs[r]))
      continue;

    const CoinShallowPackedVector vec = interPtsAndRays.getVector(r);
    if (newRay == vec) // This invokes the operator from Ray.hpp
      return r;
  }
  return -1;
}

/**
 * Determine if ray was previously cut by this hplane
 */
bool rayCutByHplaneForSplit(const int oldHplaneIndex,
    const std::vector<Hplane>& hplaneStore, const int ray_ind,
    const std::vector<Ray>& rayStore, const int split_ind) {
  if (oldHplaneIndex < 0) {
    return false;
  }

  const int numRaysCutForSplit =
      hplaneStore[oldHplaneIndex].rayToBeCutByHplane[split_ind].size();
  for (int r = 0; r < numRaysCutForSplit; r++) {
    if (rayStore[ray_ind]
        == rayStore[hplaneStore[oldHplaneIndex].rayToBeCutByHplane[split_ind][r]]) {
      return true;
    }
  }
  return false;
}

/*
 * Method checks if a hyperplane has been previous activated.
 * Only checks indices specified
 * @return the index if yes and -1 otherwise
 */
int hasHplaneBeenActivated(const std::vector<Hplane>& hplanesPrevAct,
    const std::vector<int> actHplaneIndex, const Hplane& tmpHplane,
    const bool returnIndexWIActivatedVec) {
  for (int i = 0; i < (int) actHplaneIndex.size(); i++) {
    const int hplane_ind = actHplaneIndex[i];
    const int prev_hplane_var = hplanesPrevAct[hplane_ind].var;
    const int prev_hplane_ubflag = hplanesPrevAct[hplane_ind].ubflag;
    if ((prev_hplane_var == tmpHplane.var)
        && (prev_hplane_ubflag == tmpHplane.ubflag)) {
      return returnIndexWIActivatedVec ? i : hplane_ind;
    }
  }
  return -1;
}

/*
 * Method checks if a hyperplane has been previous activated.
 * @return the index if yes and -1 otherwise
 */
int hasHplaneBeenActivated(const std::vector<Hplane>& hplanesPrevAct,
    const Hplane& tmpHplane, const int start_ind) {
  for (int hplane_ind = start_ind; hplane_ind < (int) hplanesPrevAct.size();
      hplane_ind++) {
    if (hplanesPrevAct[hplane_ind] == tmpHplane) {
      return hplane_ind;
    }
  }
  return -1;
}

/**
 * Method checks if a vertex has been previously created
 * @return the index if yes and -1 otherwise
 */
int hasVertexBeenActivated(const std::vector<Vertex>& verticesPrevAct,
    const Vertex& tmpVertex) {
  for (int vert_ind = 0; vert_ind < (int) verticesPrevAct.size();
      vert_ind++) {
    if (tmpVertex == verticesPrevAct[vert_ind])
      return vert_ind;
  }
  return -1;
}

/**
 *
 */
bool calcIntersectionPointWithSplit(Point& intPt, double& dist,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
    const Vertex& vertex, const Ray& ray, const int split_var,
    const bool recalc_dist) {
  if (recalc_dist) {
    dist = distanceToSplit(vertex, ray, split_var, solver, solnInfo);
  }
  if (isInfinity(dist, solver->getInfinity())) {
    return false;
  }

  calcNewVectorCoordinates(intPt, dist, vertex, ray, CoinMin(param.getEPS(), solnInfo.EPS));
  return true;
}

/**
 *
 */
bool calcIntersectionPointWithSplit(Point& intPt,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
    const Vertex& vertex, const Ray& ray, const int split_var) {
  const double dist = distanceToSplit(vertex, ray, split_var, solver,
      solnInfo);
  if (isInfinity(dist, solver->getInfinity())) {
    return false;
  }

  calcNewVectorCoordinates(intPt, dist, vertex, ray, CoinMin(param.getEPS(), solnInfo.EPS));
  return true;
}

void calcNewRayCoordinatesRank1(CoinPackedVector& rayIn, //const PointCutsSolverInterface* const solver,
    const SolutionInfo& solnInfo, const int rayOutIndex,
    const int newRayIndex, //const std::vector<Ray>& rayStore,
    const Hplane& hplane, const double RAYEPS) {
//    const int hplane_ind, const std::vector<Hplane>& hplaneStore) {
  std::vector<int> rayInIndex;
  std::vector<double> rayInVal;

  // First component
  rayInIndex.push_back(newRayIndex);
  rayInVal.push_back(1.0);

  // Second component (ignored if activating nb bound)
  // Also ignored if rayOutIndex is -1, meaning no ray was cut, so hplane is some dummy hplane
  if (rayOutIndex != -1 && hplane.row >= 0) {
    const double rayDirn = -1.0 * solnInfo.raysOfC1[newRayIndex][hplane.row]
        / solnInfo.raysOfC1[rayOutIndex][hplane.row];
    if (!isZero(rayDirn, RAYEPS)) {
      rayInIndex.push_back(rayOutIndex);
      rayInVal.push_back(rayDirn);
    }
  }

  // Add and sort
  rayIn.setVector(rayInIndex.size(), rayInIndex.data(), rayInVal.data());
  rayIn.sortIncrIndex();
}

/**
 * @brief Reimplemented from AdvCuts class because of precision issues with taking inverses
 *
 * @param split_var :: split
 * @param vec :: must be sorted in increasing order (in terms of indices) and be in NB space
 */
double getSICActivity(const int split_var,
    const std::vector<Vertex>& vertexStore,
    const std::vector<Ray>& rayStore, const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const CoinPackedVector& vec) {
  const double EPS = CoinMin(param.getEPS(), solnInfo.EPS);
  const int numVecElem = vec.getNumElements();
  if (numVecElem == 0)
    return +0.0;

  const int* vecIndex = vec.getIndices();
  const double* vecElem = vec.getElements();

  double activity = +0.0;

  for (int vec_el = 0; vec_el < numVecElem; vec_el++) {
    const int vec_ind = vecIndex[vec_el];
    const double val = vecElem[vec_el] 
        / distanceToSplit(vertexStore[0], rayStore[vec_ind],
            split_var, solver, solnInfo);
    if (!isZero(vecElem[vec_el], EPS) || !isZero(val, EPS)) {
      activity += val;
    }
  }

  return activity;
}

/**
 * Finds coordinates of new vertex or point obtained by
 * proceeding along rayOut from vertex a distance of dist
 * Depends crucially on both being sorted in increasing order
 */
void calcNewVectorCoordinates(CoinPackedVector& newVector, const double dist,
    const CoinPackedVector& vertex, const CoinPackedVector& rayOut, const double EPS) {
  const int numVertexEls = vertex.getNumElements();
  const int* vertexIndex = vertex.getIndices();
  const double* vertexVal = vertex.getElements();

  const int numRayEls = rayOut.getNumElements();
  const int* rayIndex = rayOut.getIndices();
  const double* rayVal = rayOut.getElements();

  std::vector<int> newVectorIndex;
  std::vector<double> newVectorVal;
  newVectorIndex.reserve(numVertexEls + numRayEls);
  newVectorVal.reserve(numVertexEls + numRayEls);

  int vert_el = 0;
  for (int ray_el = 0; ray_el < numRayEls; ray_el++) {
    const int ray_ind = rayIndex[ray_el];
    const double ray_val = rayVal[ray_el];

    while ((vert_el < numVertexEls) && (vertexIndex[vert_el] < ray_ind)) {
      // The new vertex has not changed along this coordinate
      if (!isZero(vertexVal[vert_el], EPS)) {
        newVectorIndex.push_back(vertexIndex[vert_el]);
        newVectorVal.push_back(vertexVal[vert_el]);
      }

      vert_el++;
    }

    double init_vert_val = 0.0;
    if ((vert_el < numVertexEls) && (vertexIndex[vert_el] == ray_ind)) {
      init_vert_val = vertexVal[vert_el];
      vert_el++;
    }

    const double new_val = init_vert_val + dist * ray_val;
    if (!isZero(new_val, EPS)) {
      newVectorIndex.push_back(ray_ind);
      newVectorVal.push_back(new_val);
    }
  }

  // Add any remaining vertex elements
  while (vert_el < numVertexEls) {
    // The new vertex has not changed along this coordinate
    newVectorIndex.push_back(vertexIndex[vert_el]);
    newVectorVal.push_back(vertexVal[vert_el]);

    vert_el++;
  }

  newVector.setVector(newVectorIndex.size(), newVectorIndex.data(),
      newVectorVal.data());
} /* calcNewVectorCoordinates */

double distanceToSplit(const CoinPackedVector& vertex,
    const CoinPackedVector& rayOut, const int split_var,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo) {
  double dist = solver->getInfinity();
  const int split_row = solnInfo.rowOfVar[split_var];

  Hplane hplane(split_var, split_row, static_cast<int>(HplaneBoundFlag::none));

  const double raySplitDirn = hplane.dotProductWithHplane(rayOut, solnInfo,
      true);

  // Which side of the hyperplane are we trying to hit?
  const double splitVal = solnInfo.a0[hplane.row];

  double boundBasicSpace;
  const double EPS = param.getRAYEPS(); //CoinMin(param.getRAYEPS(), solnInfo.EPS);
  if (isZero(raySplitDirn, EPS)) {
    // Ray does not intersect the bd of the split
    return dist; // infinity
  } else if (greaterThanVal(raySplitDirn, 0.0, EPS)) {
    hplane.ubflag = static_cast<int>(HplaneBoundFlag::toUB);
    boundBasicSpace = std::ceil(splitVal);
  } else {
    hplane.ubflag = static_cast<int>(HplaneBoundFlag::toLB);
    boundBasicSpace = std::floor(splitVal);
  }

  const double vDotH = hplane.dotProductWithHplane(vertex, solnInfo, false);
  const double xVal = splitVal + vDotH;

  // Ensure the vertex is not violating the hyperplane
  // Will not happen when we intersect each ray with the first hplane
  if (hplane.ubflag <= 0) {
    if (lessThanVal(xVal, boundBasicSpace)) {
      return dist;
    }
  } else {
    if (greaterThanVal(xVal, boundBasicSpace)) {
      return dist;
    }
  }

  dist = (boundBasicSpace - xVal) / (raySplitDirn);

//  if (isZero(dist)) {
//    dist = 0.0;
//  }

  if (lessThanVal(dist, 0.0)) {
    // Negative distance means the ray is in the opposite direction
    return solver->getInfinity();
  } else {
    return dist;
  }
} /* distanceToSplit */

/**
 * Used to calculate value of new point on split
 */
double dotProductWithOptTableau(const CoinPackedVector& vec, const int row_ind,
    const SolutionInfo& solnInfo) {
  const int numElmts = vec.getNumElements();
  const int* elmtIndices = vec.getIndices();
  const double* elmtVals = vec.getElements();

  double NBRowActivity = 0.0;
  for (int k = 0; k < numElmts; k++) {
    //      if (!isZero(solnInfo.raysOfC1[elmtIndices[k]][row_ind], param.getRAYEPS())
    //          && !isZero(elmtVals[k], param.getEPS())) {
    NBRowActivity += solnInfo.raysOfC1[elmtIndices[k]][row_ind]
        * elmtVals[k];
    //      }
  }
  return (solnInfo.a0[row_ind] + NBRowActivity);
} /* dotProductWithOptTableau */

/**
 * Checks that the point p does not violate any hyperplane or bound
 */
bool pointInP(const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const IntersectionInfo& interPtsAndRays,
    const int pt_ind) {
  if (isZero(interPtsAndRays.RHS[pt_ind]))
    return false;

  const CoinShallowPackedVector tmpVector = interPtsAndRays.getVector(pt_ind);
  const int numElmts = tmpVector.getNumElements();
  const int* elmtIndices = tmpVector.getIndices();
  const double* elmtVals = tmpVector.getElements();

  return pointInP(solver, solnInfo, numElmts, elmtIndices, elmtVals);
} /* pointInP */

bool pointInP(const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const int numElmts, const int* elmtIndices,
    const double* elmtVals) {
  // We check to make sure each of the rays stays within its bounds
  for (int k = 0; k < numElmts; k++) {
    if (lessThanVal(elmtVals[k], 0., param.getDIFFEPS())) {
      error_msg(errstr,
          "Point has negative component corresponding to index %d: %e.\n",
          elmtIndices[k], elmtVals[k]);
      writeErrorToII(errstr, GlobalVariables::log_file);
      exit(1);
    } else if (lessThanVal(elmtVals[k], 0., param.getEPS())) {
      warning_msg(errstr,
          "Point has negative component corresponding to index %d: %e. Small enough violation that only sending warning.\n",
          elmtIndices[k], elmtVals[k]);
    }

    // Make sure we have not proceeded past the bound for each ray
    // Actually, this could happen using Selva's PHA method
    double colLB = 0.0, colUB = solver->getInfinity();
    int currVar = solnInfo.nonBasicVarIndex[elmtIndices[k]];
    if (currVar < solnInfo.numCols) {
      colLB = getVarLB(solver, currVar);
      colUB = getVarUB(solver, currVar);
    }
    if (greaterThanVal(elmtVals[k], colUB - colLB, param.getEPS())) {
      // We have a violation
      // This is possible, since a ray stemming from a new vertex may go past its bound
      // However, if it is the cut ray, then this should not happen
      // (unless using Selva's PHA procedure
      // --- might intersect this ray by a hplane activated later and lying "past the bound")
      return false;
    }
  }

  // For each basic variable, compute its new value, given how far
  // traveled for each of the ray directions in comp NB space
  // Check if any bounds are violated
//  const int split_var = solnInfo.feasSplitVar[split_ind];
//  const int split_row = solnInfo.rowOfVar[split_var];

  for (int row_ind = 0; row_ind < solnInfo.numRows; row_ind++) {
    const int varBasicInRow = solnInfo.varBasicInRow[row_ind];
    const double xVal = solnInfo.a0[row_ind];
    double NBRowActivity = 0.0;
    for (int k = 0; k < numElmts; k++) {
//      if (!isZero(solnInfo.raysOfC1[elmtIndices[k]][row_ind], param.getRAYEPS())
//          && !isZero(elmtVals[k], param.getEPS())) {
      NBRowActivity += solnInfo.raysOfC1[elmtIndices[k]][row_ind]
          * elmtVals[k];
//      }
    }
    double newVal = xVal + NBRowActivity;

//    // We check whether the point lies on the boundary of the split
//    if (row_ind == split_row) {
//      // Value should be either floor or ceiling of xVal
//      const double floorVal = std::floor(xVal);
//      const double ceilVal = std::ceil(xVal);
//
//      if (!isVal(newVal, floorVal) && !isVal(newVal, ceilVal)) {
//        error_msg(errstr,
//            "Point %d not on bd of split %d (var %d), created by activating hplane %d. newVal: %s, floor: %s, ceil: %s.\n",
//            pt_ind, split_ind, varBasicInRow,
//            interPtsAndRays.hplane[pt_ind],
//            stringValue(newVal).c_str(),
//            stringValue(floorVal).c_str(),
//            stringValue(ceilVal).c_str());
//        writeErrorToII(errstr, GlobalVariables::inst_info_out);
//        exit(1);
//      }
//    }

    const double colLB = getVarLB(solver, varBasicInRow);
    const double colUB = getVarUB(solver, varBasicInRow);

    // What we really want to ensure is that the variable corresponding
    // to the hplane we activated (if there was one) stays within its bounds
    if (lessThanVal(newVal, colLB, param.getEPS())
        || greaterThanVal(newVal, colUB, param.getEPS())) {
//      if (hplane_row_ind <= -1)
//        return false;
//
//      if (row_ind == hplane_row_ind) {
//        error_msg(errstr,
//            "Basic var %d corresponding to hplane %d activated on split %d yielding pt %d is outside its bounds. newVal: %s, colLB: %s, colUB: %s.\n",
//            varBasicInRow, interPtsAndRays.hplane[pt_ind],
//            split_ind, pt_ind, stringValue(newVal).c_str(),
//            stringValue(colLB).c_str(), stringValue(colUB).c_str());
//        writeErrorToII(errstr, GlobalVariables::inst_info_out);
//        exit(1);
//      }

      return false;
    }
  }

  return true;
} /* pointInP */

/**
 * We assume ray starts at origin (LP optimum)
 */
bool isRayFinal(const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const int ray_ind) {
  // Check whether the ray hits its own other bound
  const double eps = CoinMin(param.getEPS(), solnInfo.EPS);
  const int rayVarStructIndex = solnInfo.nonBasicVarIndex[ray_ind];
  if (rayVarStructIndex < solver->getNumCols() && 
      !isNegInfinity(getVarLB(solver, rayVarStructIndex), solver->getInfinity()) && 
      !isInfinity(getVarUB(solver, rayVarStructIndex), solver->getInfinity())) {
    return false;
  }

  // Now check the basic variables
  bool rayHitsSomething = false;
  for (int row_ind = 0; row_ind < solnInfo.numRows; row_ind++) {
    if (solver->getRowSense()[row_ind] == 'E') {
      continue; // Do not activate equality rows
    }
    const int var = solnInfo.varBasicInRow[row_ind];

    const double varVal = getVarVal(solver, var);
    const double rayVal = solnInfo.raysOfC1[ray_ind][row_ind];
    double bound = 0.;
    if (isZero(rayVal, eps)) {
      continue; // ray does not intersect this hyperplane
    } else if (lessThanVal(rayVal, 0.0, eps)) {
      bound = getVarLB(solver, var);
    } else if (greaterThanVal(rayVal, 0.0, eps)) {
      bound = getVarUB(solver, var);
    }

    // If this bound is actually at infinity; we will never intersect it
    if (isInfinity(std::abs(bound), solver->getInfinity())) {
      continue;
    }

    // Otherwise, find the distance to the hyperplane, and check if it is actually intersected
    const double tmpDistToHplane = (bound - varVal) / rayVal;
    if (isInfinity(tmpDistToHplane, solver->getInfinity())) {
      continue; // ray does not intersect this hyperplane
    } else {
      rayHitsSomething = true;
      break;
    }
  } /* iterate over basic variables */

  return !rayHitsSomething;
} /* isRayFinal */

/**
 * Checks that the ray r does not violate any hyperplane or bound
 */
bool isRayFinal(const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const IntersectionInfo& interPtsAndRays,
    const int ray_ind) {
  if (!isZero(interPtsAndRays.RHS[ray_ind]))
    return false;

  const CoinShallowPackedVector tmpVector = interPtsAndRays.getVector(ray_ind);
  const int numElmts = tmpVector.getNumElements();
  const int* elmtIndices = tmpVector.getIndices();
  const double* elmtVals = tmpVector.getElements();

  return isRayFinal(solver, solnInfo, numElmts, elmtIndices, elmtVals);
} /* isRayFinal */

/**
 *
 */
bool isRayFinal(const PHASolverInterface* const solver, 
    const SolutionInfo& solnInfo, const int n,
    const int* indices, const double* values) {
  const double eps = CoinMin(param.getEPS(), solnInfo.EPS);

  // Check nonbasic bounds
  // If any of the directions is negative, we know new ray hits lower bound on that nb dirn
  // Else we check ub and lb both existing as before for structural variables
  for (int i = 0; i < n; i++) {
    if (lessThanVal(values[i], 0., eps)) {
      return false;
    }
    const int ray_ind = indices[i];
    const int rayVarStructIndex = solnInfo.nonBasicVarIndex[ray_ind];
    if (rayVarStructIndex < solver->getNumCols() &&
        !isNegInfinity(getVarLB(solver, rayVarStructIndex), solver->getInfinity()) &&
        !isInfinity(getVarUB(solver, rayVarStructIndex), solver->getInfinity())) {
      return false;
    }
  }

  // Now check the basic variables
  bool rayHitsSomething = false;
  for (int row_ind = 0; row_ind < solnInfo.numRows; row_ind++) {
    if (solver->getRowSense()[row_ind] == 'E') {
      continue; // Do not activate equality rows
    }
    const int var = solnInfo.varBasicInRow[row_ind];
    const double varVal = getVarVal(solver, var);
    double rayVal = 0.;
    for (int i = 0; i < n; i++) {
      const int ray_ind = indices[i];
      rayVal += solnInfo.raysOfC1[ray_ind][row_ind];
    }

    double bound = 0.;
    if (isZero(rayVal, eps)) {
      continue; // ray does not intersect this hyperplane
    } else if (lessThanVal(rayVal, 0.0, eps)) {
      bound = getVarLB(solver, var);
    } else if (greaterThanVal(rayVal, 0.0, eps)) {
      bound = getVarUB(solver, var);
    }

    // If this bound is actually at infinity; we will never intersect it
    if (isInfinity(std::abs(bound), solver->getInfinity())) {
      continue;
    }

    // Otherwise, find the distance to the hyperplane, and check if it is actually intersected
    const double tmpDistToHplane = (bound - varVal) / rayVal;
    if (isInfinity(tmpDistToHplane, solver->getInfinity())) {
      continue; // ray does not intersect this hyperplane
    } else {
      rayHitsSomething = true;
      break;
    }
  } /* iterate over basic variables */

  return !rayHitsSomething;
} /* isRayFinal */

/*********************/
/* GENERAL FUNCTIONS */
/*********************/

/***********************************************************************/
/**
 * @brief Get environment variable by key.
 * @param key :: Key to look up.
 */
std::string get_env_var(std::string const & key) {
  char * val;
  val = getenv(key.c_str());
  std::string retval = "";
  if (val != NULL) {
    retval = val;
  }
  return retval;
}

/***********************************************************************/
/**
 * @brief Find value in a std::vector and return where it is in the std::vector.
 * @param key :: Key to look up.
 */
int find_val(const int &key, const std::vector<int> &vec) {
  for (int k = 0; k < (int) vec.size(); k++) {
    if (vec[k] == key)
      return k;
  }
  return -1;
}

/***********************************************************************/
/**
 * @brief Find value in a std::vector and return where it is in the std::vector.
 * @param key :: Key to look up.
 */
int find_val(const double &key, const std::vector<double> &vec) {
  for (int k = 0; k < (int) vec.size(); k++) {
    if (vec[k] == key)
      return k;
  }
  return -1;
}

/***********************************************************************/
/**
 * @brief Find value in a std::vector and return where it is in the std::vector.
 * @param key :: Key to look up.
 */
int find_val(const std::string &key, const std::vector<std::string> &vec) {
  for (int k = 0; k < (int) vec.size(); k++) {
    if (vec[k].compare(key) == 0)
      return k;
  }
  return -1;
}

/****************************************************************************/
/**
 * @brief Check if a given file exists in the system.
 *
 * @returns True if exists, False otherwise
 */
bool fexists(const char *filename) {
//  std::ifstream ifile(filename);
//  return ifile; // This worked before g++5
  struct stat buffer;
  return (stat (filename, &buffer) == 0);
} /* fexists */

/********************************************************************************/
/**
 *
 */
double computeAverage(const std::vector<double> &vec) {
  return computeAverage(vec, 0, (int) vec.size() - 1);
}

/********************************************************************************/
/**
 *
 */
double computeAverage(const std::vector<int> &vec) {
  return computeAverage(vec, 0, (int) vec.size() - 1);
}

/********************************************************************************/
/**
 *
 */
double computeAverage(const std::vector<double> &vec, const int start_index,
    const int end_index) {
//  if (start_index < 0 || start_index > vecsize || end_index < 0
//      || end_index > vecsize || start_index > end_index)
//    return 0.0;

  double sum = 0.0;
  int n = end_index - start_index + 1;
  for (int i = start_index; i <= end_index; i++) {
    sum += vec[i];
  }

  return sum / n;
}

/********************************************************************************/
/**
 *
 */
double computeAverage(const std::vector<int> &vec, const int start_index,
    const int end_index) {
//  if (start_index < 0 || start_index > vecsize || end_index < 0
//      || end_index > vecsize || start_index > end_index)
//    return 0.0;

  double sum = 0.0;
  int n = end_index - start_index + 1;
  for (int i = start_index; i <= end_index; i++) {
    sum += vec[i];
  }

  return sum / n;
}

/********************************************************************************/
/**
 *
 */
double computeAverage(const std::vector<double> &vec,
    const std::vector<int> &index) {
  int n = (int) index.size();
//  if (n == 0)
//    return 0.0;

  double sum = 0.0;
  for (int i = 0; i < (int) index.size(); i++) {
    sum += vec[index[i]];
  }

  return sum / n;
}

/********************************************************************************/
/**
 *
 */
double computeAverage(const std::vector<int> &vec,
    const std::vector<int> &index) {
  int n = (int) index.size();
//  if (n == 0)
//    return 0.0;

  double sum = 0.0;
  for (int i = 0; i < (int) index.size(); i++) {
    sum += vec[index[i]];
  }

  return sum / n;
}

/********************************************************************************/
/**
 *
 */
double computeAverageDiff(const std::vector<double> &vec1,
    const std::vector<double> &vec2) {
  return computeAverageDiff(vec1, vec2, 0, (int) vec1.size() - 1);
}

/********************************************************************************/
/**
 *
 */
double computeAverageDiff(const std::vector<int> &vec1,
    const std::vector<int> &vec2) {
  return computeAverageDiff(vec1, vec2, 0, (int) vec1.size() - 1);
}

/********************************************************************************/
/**
 *
 */
double computeAverageDiff(const std::vector<double> &vec1,
    const std::vector<double> &vec2, const int start_index,
    const int end_index) {
//  if (start_index < 0 || start_index > vecsize || end_index < 0
//      || end_index > vecsize || start_index > end_index)
//    return 0.0;

  double sum = 0.0;
  int n = end_index - start_index + 1;
  for (int i = start_index; i <= end_index; i++) {
    sum += vec1[i] - vec2[i];
  }

  return sum / n;
}

/********************************************************************************/
/**
 *
 */
double computeAverageDiff(const std::vector<int> &vec1,
    const std::vector<int> &vec2, const int start_index,
    const int end_index) {
//  if (start_index < 0 || start_index > vecsize || end_index < 0
//      || end_index > vecsize || start_index > end_index)
//    return 0.0;

  double sum = 0.0;
  int n = end_index - start_index + 1;
  for (int i = start_index; i <= end_index; i++) {
    sum += vec1[i] - vec2[i];
  }

  return sum / n;
}

/********************************************************************************/
/**
 *
 */
double computeAverageDiff(const std::vector<double> &vec1,
    const std::vector<double> &vec2, const std::vector<int> &index) {
  int n = (int) index.size();
//  if (n == 0)
//    return 0.0;

  double sum = 0.0;
  for (int i = 0; i < (int) index.size(); i++) {
    sum += vec1[index[i]] - vec2[index[i]];
  }

  return sum / n;
}

/********************************************************************************/
/**
 *
 */
double computeAverageDiff(const std::vector<int> &vec1,
    const std::vector<int> &vec2, const std::vector<int> &index) {
  int n = (int) index.size();
//  if (n==0)
//    return 0.0;

  double sum = 0.0;
  for (int i = 0; i < (int) index.size(); i++) {
    sum += vec1[index[i]] - vec2[index[i]];
  }

  return sum / n;
}

/********************************************************************************/
/**
 *
 */
std::vector<double>::const_iterator computeMin(const std::vector<double> &vec) {
  std::vector<double>::const_iterator minit = vec.begin();
  for (std::vector<double>::const_iterator it = vec.begin(); it != vec.end();
      ++it) {
    if (*it < *minit) {
      minit = it;
    }
  }

  return minit;
}

/********************************************************************************/
/**
 *
 */
std::vector<int>::const_iterator computeMin(const std::vector<int> &vec) {
  std::vector<int>::const_iterator minit = vec.begin();
  for (std::vector<int>::const_iterator it = vec.begin(); it != vec.end();
      ++it) {
    if (*it < *minit) {
      minit = it;
    }
  }

  return minit;
}

double computeMinAbs(const int num_el, const double* vec) {
  int min_ind = 0;
  if (num_el == 0) {
    return std::numeric_limits<double>::lowest();
  }
  for (int i = 1; i < num_el; ++i) {
    if (std::abs(vec[i]) < std::abs(vec[min_ind])) {
      min_ind = i;
    }
  }

  return std::abs(vec[min_ind]);
}

double computeMaxAbs(const int num_el, const double* vec) {
  int max_ind = 0;
  if (num_el == 0) {
    return std::numeric_limits<double>::max();
  }
  for (int i = 1; i < num_el; ++i) {
    if (std::abs(vec[i]) > std::abs(vec[max_ind])) {
      max_ind = i;
    }
  }

  return std::abs(vec[max_ind]);
}

/********************************************************************************/
/**
 *
 */
std::vector<double>::const_iterator computeMin(const std::vector<double> &vec,
    const std::vector<int> &index) {
  std::vector<double>::const_iterator minit = vec.begin() + index[0];
  for (std::vector<int>::const_iterator it = index.begin(); it != index.end();
      ++it) {
    if (*(vec.begin() + *it) < *minit) {
      minit = vec.begin() + *it;
    }
  }

  return minit;
}

/********************************************************************************/
/**
 *
 */
std::vector<int>::const_iterator computeMin(const std::vector<int> &vec,
    const std::vector<int> &index) {
  std::vector<int>::const_iterator minit = vec.begin() + index[0];
  for (std::vector<int>::const_iterator it = index.begin(); it != index.end();
      ++it) {
    if (*(vec.begin() + *it) < *minit) {
      minit = vec.begin() + *it;
    }
  }

  return minit;
}

/********************************************************************************/
/**
 *
 */
std::vector<double>::const_iterator computeMax(const std::vector<double> &vec) {
  std::vector<double>::const_iterator minit = vec.begin();
  for (std::vector<double>::const_iterator it = vec.begin(); it != vec.end();
      ++it) {
    if (*it > *minit) {
      minit = it;
    }
  }

  return minit;
}

/********************************************************************************/
/**
 *
 */
std::vector<int>::const_iterator computeMax(const std::vector<int> &vec) {
  std::vector<int>::const_iterator minit = vec.begin();
  for (std::vector<int>::const_iterator it = vec.begin(); it != vec.end();
      ++it) {
    if (*it > *minit) {
      minit = it;
    }
  }

  return minit;
}

/********************************************************************************/
/**
 *
 */
std::vector<double>::const_iterator computeMax(const std::vector<double> &vec,
    const std::vector<int> &index) {
  std::vector<double>::const_iterator minit = vec.begin() + index[0];
  for (std::vector<int>::const_iterator it = index.begin(); it != index.end();
      ++it) {
    if (*(vec.begin() + *it) > *minit) {
      minit = vec.begin() + *it;
    }
  }

  return minit;
}

/********************************************************************************/
/**
 *
 */
std::vector<int>::const_iterator computeMax(const std::vector<int> &vec,
    const std::vector<int> &index) {
  std::vector<int>::const_iterator minit = vec.begin() + index[0];
  for (std::vector<int>::const_iterator it = index.begin(); it != index.end();
      ++it) {
    if (*(vec.begin() + *it) > *minit) {
      minit = vec.begin() + *it;
    }
  }

  return minit;
}

/********************************************************************************/
/**
 * @brief Computes minimum difference between the two std::vectors.
 * This is not commutative. We compute min_{i} (vec1[i] - vec2[i]).
 *
 */
double computeMinDiff(const std::vector<double> &vec1,
    const std::vector<double> &vec2) {
  double currmin = vec1[0] - vec2[0];

  std::vector<double>::const_iterator it1 = vec1.begin();
  std::vector<double>::const_iterator it2 = vec2.begin();
  while (it1 != vec1.end() && it2 != vec2.end()) {
    double tmp = *it1 - *it2;
    if (tmp < currmin) {
      currmin = tmp;
    }
  }

  return currmin;
}

/********************************************************************************/
/**
 * @brief Computes minimum difference between the two std::vectors.
 * This is not commutative. We compute min_{i} (vec1[i] - vec2[i]).
 *
 */
int computeMinDiff(const std::vector<int> &vec1, const std::vector<int> &vec2) {
  int currmin = vec1[0] - vec2[0];

  std::vector<int>::const_iterator it1 = vec1.begin();
  std::vector<int>::const_iterator it2 = vec2.begin();
  while (it1 != vec1.end() && it2 != vec2.end()) {
    int tmp = *it1 - *it2;
    if (tmp < currmin) {
      currmin = tmp;
    }
  }

  return currmin;
}

/********************************************************************************/
/**
 * @brief Computes minimum difference between the two std::vectors, only considering
 * those indices provided.
 * This is not commutative. We compute min_{i \in index} (vec1[i] - vec2[i]).
 *
 */
double computeMinDiff(const std::vector<double> &vec1,
    const std::vector<double> &vec2, const std::vector<int> &index) {
  double currmin = vec1[index[0]] - vec2[index[0]];

  for (int i = 0; i < (int) index.size(); i++) {
    double tmp = vec1[index[i]] - vec2[index[i]];
    if (tmp < currmin)
      currmin = tmp;
  }

  return currmin;
}

/********************************************************************************/
/**
 * @brief Computes minimum difference between the two std::vectors, only considering
 * those indices provided.
 * This is not commutative. We compute min_{i \in index} (vec1[i] - vec2[i]).
 *
 */
int computeMinDiff(const std::vector<int> &vec1, const std::vector<int> &vec2,
    const std::vector<int> &index) {
  int currmin = vec1[index[0]] - vec2[index[0]];

  for (int i = 0; i < (int) index.size(); i++) {
    int tmp = vec1[index[i]] - vec2[index[i]];
    if (tmp < currmin)
      currmin = tmp;
  }

  return currmin;
}

/********************************************************************************/
/**
 * @brief Computes maximum difference between the two std::vectors, only considering
 * those indices provided.
 * This is not commutative. We compute max_{i \in index} (vec1[i] - vec2[i]).
 *
 */
double computeMaxDiff(const std::vector<double> &vec1,
    const std::vector<double> &vec2) {
  double currmax = vec1[0] - vec2[0];

  for (int i = 0; i < (int) vec1.size(); i++) {
    double tmp = vec1[i] - vec2[i];
    if (tmp > currmax)
      currmax = tmp;
  }

  return currmax;
}

/********************************************************************************/
/**
 * @brief Computes maximum difference between the two std::vectors, only considering
 * those indices provided.
 * This is not commutative. We compute max_{i \in index} (vec1[i] - vec2[i]).
 *
 */
int computeMaxDiff(const std::vector<int> &vec1, const std::vector<int> &vec2) {
  int currmax = vec1[0] - vec2[0];

  for (int i = 0; i < (int) vec1.size(); i++) {
    int tmp = vec1[i] - vec2[i];
    if (tmp > currmax)
      currmax = tmp;
  }

  return currmax;
}

/********************************************************************************/
/**
 * @brief Computes maximum difference between the two std::vectors, only considering
 * those indices provided.
 * This is not commutative. We compute max_{i \in index} (vec1[i] - vec2[i]).
 *
 */
double computeMaxDiff(const std::vector<double> &vec1,
    const std::vector<double> &vec2, const std::vector<int> &index) {
  double currmax = vec1[index[0]] - vec2[index[0]];

  for (int i = 0; i < (int) index.size(); i++) {
    double tmp = vec1[index[i]] - vec2[index[i]];
    if (tmp > currmax)
      currmax = tmp;
  }

  return currmax;
}

/********************************************************************************/
/**
 * @brief Computes maximum difference between the two std::vectors, only considering
 * those indices provided.
 * This is not commutative. We compute max_{i \in index} (vec1[i] - vec2[i]).
 *
 */
int computeMaxDiff(const std::vector<int> &vec1, const std::vector<int> &vec2,
    const std::vector<int> &index) {
  int currmax = vec1[index[0]] - vec2[index[0]];

  for (int i = 0; i < (int) index.size(); i++) {
    int tmp = vec1[index[i]] - vec2[index[i]];
    if (tmp > currmax)
      currmax = tmp;
  }

  return currmax;
}

/********************/
/* CONVERSION FUNCS */
/********************/

/**
 * Create coordinates of the point from the std::vector structSolution,
 * using the cobasis defined in solnInfo
 */
void compNBCoor(std::vector<double>& NBSolution, const double* structSolution,
    const PHASolverInterface* const tmpSolver,
    const PHASolverInterface* const origSolver,
    const SolutionInfo& origSolnInfo) {
  NBSolution.resize(origSolnInfo.numNB);

  // All coefficients are going to be non-negative, since we work in the complemented NB space
  for (int i = 0; i < origSolnInfo.numNB; i++) {
    int currVar = origSolnInfo.nonBasicVarIndex[i];
    double tmpVal, origVal;
    if (currVar < origSolnInfo.numCols) {
      tmpVal = structSolution[currVar];
      origVal = origSolver->getColSolution()[currVar];
    } else {
      int currRow = currVar - origSolnInfo.numCols;
      tmpVal = tmpSolver->getRightHandSide()[currRow]
          - tmpSolver->getRowActivity()[currRow];
      origVal = 0.0;
    }
    // ***
    // For *struct* var:
    // If the variable is lower-bounded, then we need to shift it up for the comp space
    // That is, xComp = x - lb
    // If the variable is upper-bounded, we need to get the value in the complemented space
    // That is, xComp = ub - x
    // What to do if the variable is fixed NB or free NB?
    // In the former case, thinking of it as a lower-bounded variable, it will simply be zero
    // in any solution, including in the original one at v, in which case we don't need to
    // worry about it for the packed ray.
    // In the latter case, we may need to do something more clever, as in wrt original val
    // of variable in the NB space at v.
    // Actually, value should not change from that in v...?
    // TODO Figure out what to do for free vars
    // ****
    // For *slack* var:
    // In the original NB space at v, all slacks were at zero.
    // We simply look at b[row] - activity[row]
    // For >= rows, slack is non-positive
    // For <= rows, slack is non-negative
    // We take absolute value to deal with this
    const double newVal = std::abs(tmpVal - origVal);
    if (!isZero(newVal))
      NBSolution[i] = newVal;
    else
      NBSolution[i] = 0.0;
  }
//  CoinPackedVector::setVector((int) packedIndex.size(), packedIndex.data(),
//      packedVal.data(), false);
} /* compNBCoor */

/********************/
/* STRING FUNCTIONS */
/********************/

/***********************************************************************/
/**
 * Parses std::string line by the set of delimiters given.
 * (Only works well with one for now.)
 * Parses std::string based on delimiters provided and puts output into std::vector.
 */
void stringparse(std::string &line, std::string delimiters,
    std::vector<std::string> &strings) {
  if (line.empty())
    return;

  std::string curr;
//  cout << "Here." << " Reading: " << line << ".\n";
  // Get all characters to be used as delimiters
//  cout << "Number of delimiters is " << delimiter.length() << std::endl;
  int numdelims = delimiters.length();
//  numdelims = strlen(delims);
  const char *delims = delimiters.c_str();

  // Make istringstream from line
  for (int i = 0; i < numdelims; i++) {
    std::istringstream tmp;
    tmp.str(line);
    while (getline(tmp, curr, delims[i])) {
//      printf("%s\n", curr.c_str());
//      cout << curr << std::endl;
      strings.push_back(curr);
    }
  }
}

/************************************************************************/
// Singular Value Decomposition (from Lapack);
// needed to compute condition number in norm 2
extern "C" int dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a,
    int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt,
    double *work, int *lwork, int *info);

/**********************************************************/
double compute_condition_number_norm2(const OsiSolverInterface* const solver) {
  double condition = 0.0;
  if (!solver->basisIsAvailable()) {
    return condition;
  }
  int m = solver->getNumRows();
  int n = solver->getNumCols();
  int* rstat = new int[m];
  int* cstat = new int[n];
  solver->getBasisStatus(cstat, rstat);
  double* basis = new double[m * m];
  memset(basis, 0, m * m * sizeof(double));
  const CoinPackedMatrix* byCol = solver->getMatrixByCol();
  const int * rowPos = byCol->getIndices();
  const CoinBigIndex * columnStart = byCol->getVectorStarts();
  const int * columnLength = byCol->getVectorLengths();
  const double * columnElements = byCol->getElements();
  int bindex = 0;
  // fill the structural columns of the basis
  for (int i = 0; i < n; ++i) {
    if (cstat[i] == 1) {
      for (int j = columnStart[i]; j < columnStart[i] + columnLength[i];
          ++j) {
        basis[bindex * m + rowPos[j]] = columnElements[j];
      }
      bindex++;
    }
  }
  // fill the artificial columns of the basis
  for (int i = 0; i < m; ++i) {
    if (rstat[i] == 1) {
      basis[bindex * m + i] = 1.0;
      bindex++;
    }
  }
  // prepare data for SVD; get the optimal size of the work array
  int info = 0;
  char job = 'N';
  double* sv = new double[m];
  double sizeWork = 0.0;
  // we first call the SVD routine with lwork = -1; this way, the function
  // returns the optimal size of the work array that is needed for the
  // computation (the value is written in sizeWork)
  int lwork = -1;
  dgesvd_(&job, &job, &m, &m, basis, &m, sv, NULL, &m, NULL, &m, &sizeWork,
      &lwork, &info);
  if (info != 0) {
    printf(
        "### WARNING: can not obtain optimal size of work array for SVD\n");
    return condition;
  }
  // now we know the optimal size, so we allocate work array and do SVD
  lwork = (int) sizeWork;
  double* work = new double[lwork];
  dgesvd_(&job, &job, &m, &m, basis, &m, sv, NULL, &m, NULL, &m, work, &lwork,
      &info);
  if (info != 0) {
    printf("### WARNING: can not compute condition number\n");
    return condition;
  }
  condition = sv[0] / sv[m - 1];
  delete[] rstat;
  delete[] cstat;
  delete[] basis;
  delete[] sv;
  delete[] work;
  return condition;
} /* compute_condition_number_norm2 */

double compute_condition_number_norm1(OsiSolverInterface* solver) {
  double condition = 0.0;
  if (!solver->basisIsAvailable()) {
    return condition;
  }
  int m = solver->getNumRows();
  int n = solver->getNumCols();
  int* rstat = new int[m];
  int* cstat = new int[n];
  solver->getBasisStatus(cstat, rstat);
  const CoinPackedMatrix* byCol = solver->getMatrixByCol();
  const CoinBigIndex * columnStart = byCol->getVectorStarts();
  const int * columnLength = byCol->getVectorLengths();
  const double * columnElements = byCol->getElements();
  double norm1B = 0.0;
  double currnorm;
  // norm of structural columns
  for (int i = 0; i < n; ++i) {
    if (cstat[i] == 1) {
      currnorm = 0.0;
      for (int j = columnStart[i]; j < columnStart[i] + columnLength[i];
          ++j) {
        currnorm += std::abs(columnElements[j]);
      }
      if (currnorm > norm1B) {
        norm1B = currnorm;
      }
    }
  }
  // norm of artificial columns
  if (norm1B < 1.0) {
    norm1B = 1.0;
  }
  solver->enableSimplexInterface(true);
  solver->enableFactorization();
  double norm1Binv = 0.0;
  double* binvcol = new double[m];
  for (int i = 0; i < m; ++i) {
    solver->getBInvCol(i, binvcol);
    currnorm = 0.0;
    for (int j = 0; j < m; ++j) {
      currnorm += std::abs(binvcol[j]);
    }
    if (currnorm > norm1Binv) {
      norm1Binv = currnorm;
    }
  }
  solver->disableFactorization();
  solver->disableSimplexInterface();
  condition = norm1B * norm1Binv;
  delete[] rstat;
  delete[] cstat;
  delete[] binvcol;
  return condition;
} /* compute_condition_number_norm1 */

// compute the condition number in the given norm (1 or 2)
double compute_condition_number(OsiSolverInterface* solver, const int norm) {
  if (norm < 1 || norm > 2) {
    printf("### WARNING: can not compute condition number in norm %d\n",
        norm);
    return 0.0;
  }
  if (norm == 1) {
    return compute_condition_number_norm1(solver);
  }
  return compute_condition_number_norm2(solver);
} /* compute_condition_number */

std::string toLowerString(const std::string str) {
  std::string tmp = str;
  for (int i = 0; i < (int) tmp.length(); i++) {
    tmp[i] = tolower(tmp[i]);
  }
  return tmp;
}

std::string toUpperString(const std::string str) {
  std::string tmp = str;
  for (int i = 0; i < (int) tmp.length(); i++) {
    tmp[i] = toupper(tmp[i]);
  }
  return tmp;
}

std::string toLowerStringNoUnderscore(const std::string str) {
  std::string tmp = str;
  for (int i = 0; i < (int) tmp.length(); i++) {
    if (tmp[i] != '_') {
      tmp[i] = tolower(tmp[i]);
    }
  }
  return tmp;
}

bool remove_underscore(char c) {
  return (c == '_');
}

std::string toUpperStringNoUnderscore(const std::string str) {
  std::string tmp = str;
  std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
  tmp.erase(std::remove_if(tmp.begin(), tmp.end(), remove_underscore), tmp.end());
  return tmp;
}

std::vector<std::string> lowerCaseStringVector(const std::vector<std::string>& strVec) {
  std::vector<std::string> tmp = strVec;
  for (int i = 0; i < (int) strVec.size(); i++) {
    std::transform(tmp[i].begin(), tmp[i].end(), tmp[i].begin(), ::tolower);
  }
  return tmp;
}

/** We assume it is comma separated */
double getObjValueFromFile(const char* filename, const std::string& this_inst_name) {
  if (!filename) {
    return std::numeric_limits<double>::lowest();
  }

  std::ifstream infile(filename);
  if (infile.is_open()) {
    std::string line;
    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      if (line.empty()) {
        continue;
      }
      std::string inst_name;
      if (!(std::getline(iss, inst_name, ','))) {
        warning_msg(warnstring,
            "Could not read instance name. String is %s.\n",
            line.c_str());
        continue;
      }
      if (inst_name == this_inst_name) {
        try {
          std::string token;
          if (!(std::getline(iss, token, ','))) {
            throw;
          }
          return std::stod(token);
        } catch (std::exception& e) {
          warning_msg(warnstring,
              "Could not read optimal value. String is %s.\n",
              line.c_str());
          continue;
        }
      }
    }
    infile.close();
  } else {
    // If we were not able to open the file, throw an error
    error_msg(errorstring, "Not able to open obj file.\n");
    writeErrorToII(errorstring, GlobalVariables::log_file);
    exit(1);
  }

  return std::numeric_limits<double>::lowest();
} /* getObjValueFromFile */

/** We assume it is tab separated */
void getOptSolFromFile(const char* filename, std::vector<double>& opt_sol) {
  if (!filename) {
    return;
  }

  std::ifstream infile(filename);
  if (infile.is_open()) {
    std::string line;
    std::getline(infile, line); // discard the first line
    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      if (line.empty()) {
        continue;
      }
      if (line == "ENDATA") {
        break;
      }

      std::string name, val;
      iss >> std::skipws >> name >> val;
      opt_sol.push_back( std::stod(val) );
      //std::cout << line.c_str() << '\n'; // DEBUG
      //std::cout << "Name: " << name.c_str() << std::endl; // DEBUG
      //std::cout << "Val: " << val.c_str() << std::endl; // DEBUG
      //std::cout << opt_sol.size() << '\n'; // DEBUG
    } /* run through the file */
    infile.close();
  } else {
    // If we were not able to open the file, throw an error
    error_msg(errorstring, "Not able to open solution file.\n");
    writeErrorToII(errorstring, GlobalVariables::log_file);
    exit(1);
  }

  return;
} /* getOptSolFromFile */

/**********************************************************/
/** from Giacomo Nannicini */
double dotProduct(const double* a, const double* b, int dimension) {
  // Accumulate the dot product here
  double accumulator = 0.0;
  // Compensation term for low-order bits lost
  double compensation = 0.0;
  // Temporary number for storing current term
  double currterm = 0.0;
  // Temporary number for storing new accumulator
  double nextacc = 0.0;
  for (int i = 0; i < dimension; ++i) {
    // Compute current term of the dot product adding compensation
    currterm = (a[i] * b[i]) - compensation;
    // This is the value of the sum
    nextacc = accumulator + currterm;
    // Recover what we just lost adding currterm to accumulator
    compensation = (nextacc - accumulator) - currterm;
    // Now save new value of the accumulator
    accumulator = nextacc;
  }
  return accumulator;
} /* dotProduct (dense x dense) */

double dotProduct(int sizea, const int* indexa, const double* a,
    const double* b) {
  // Accumulate the dot product here
  double accumulator = 0.0;
  // Compensation term for low-order bits lost
  double compensation = 0.0;
  // Temporary number for storing current term
  double currterm = 0.0;
  // Temporary number for storing new accumulator
  double nextacc = 0.0;
  for (int i = 0; i < sizea; ++i) {
    // Compute current term of the dot product adding compensation
    currterm = (a[i] * b[indexa[i]]) - compensation;
    // This is the value of the sum
    nextacc = accumulator + currterm;
    // Recover what we just lost adding currterm to accumulator
    compensation = (nextacc - accumulator) - currterm;
    // Now save new value of the accumulator
    accumulator = nextacc;
  }
  return accumulator;
} /* dotProduct (sparse x dense) */

double dotProduct(int sizea, const int* indexa, const double* a, int sizeb,
    const int* indexb, const double* b) {
  // Accumulate the dot product here
  double accumulator = 0.0;
  // Compensation term for low-order bits lost
  double compensation = 0.0;
  // Temporary number for storing current term
  double currterm = 0.0;
  // Temporary number for storing new accumulator
  double nextacc = 0.0;
  // Current position in vectors a and b
  int posa = 0, posb = 0;
  while (posa < sizea && posb < sizeb) {
    // If it is the same component, compute dot product
    if (indexa[posa] == indexb[posb]) {
      // Compute current term of the dot product adding compensation
      currterm = (a[posa] * b[posb]) - compensation;
      // This is the value of the sum
      nextacc = accumulator + currterm;
      // Recover what we just lost adding currterm to accumulator
      compensation = (nextacc - accumulator) - currterm;
      // Now save new value of the accumulator
      accumulator = nextacc;
      // Increment both position indices
      ++posa;
      ++posb;
    } else if (indexa[posa] < indexb[posb]) {
      // Increment only smaller position index
      ++posa;
    } else if (indexa[posa] > indexb[posb]) {
      // Increment only smaller position index
      ++posb;
    }
  }
  return accumulator;
} /* dotProduct (sparse x sparse) */
