//============================================================================
// Name        : SolutionInfo.cpp
// Author      : Aleksandr M. Kazachkov (edited from Selvaprabu Nadarajah)
// Version     : 0.2013.06.16
// Copyright   : Your copyright notice
// Description : Information about a problem
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#include <cmath> // sqrt
#include <numeric> // inner_product
#include "SolutionInfo.hpp"
#include "Utility.hpp"

/**
 * @brief brief: if true, generate a small solution information file; specifically, it does not collect raysOfC1.
 */
SolutionInfo::SolutionInfo(OsiSolverInterface* const solver,
    const int max_frac_core, const bool brief) {
  initialize(solver, max_frac_core, brief);
} /* constructor */

void SolutionInfo::copyOurStuff(const SolutionInfo* const rhs) {
  if (rhs) {
    // primitives
    SPACE = rhs->SPACE;
//    basis = rhs->basis;
    LPopt = rhs->LPopt;
    objNorm = rhs->objNorm;
    minAbsRayVal = rhs->minAbsRayVal;
    minAbsRayValNoEPS = rhs->minAbsRayValNoEPS;
    minAbsCoeff = rhs->minAbsCoeff;
    maxAbsCoeff = rhs->maxAbsCoeff;
    minReferenceValue = rhs->minReferenceValue;
    maxReferenceValue = rhs->maxReferenceValue;
    EPS = rhs->EPS;
    numCols = rhs->numCols;
    numRows = rhs->numRows;
    numANonzero = rhs->numANonzero;
    numNBOrig = rhs->numNBOrig;
    numNBSlack = rhs->numNBSlack;
    numNB = rhs->numNB;
    numBasicCols = rhs->numBasicCols;
    numIntegerNonBinary = rhs->numIntegerNonBinary;
    numBinary = rhs->numBinary;
    numContinuous = rhs->numContinuous;
    numFixed = rhs->numFixed;
    numEqRows = rhs->numEqRows;
    numIneqRows = rhs->numIneqRows;
    numBoundRows = rhs->numBoundRows;
    numAssignRows = rhs->numAssignRows;
    numLB = rhs->numLB;
    numUB = rhs->numUB;
    numPrimalDegeneratePivots = rhs->numPrimalDegeneratePivots;
    numDualDegeneratePivots = rhs->numDualDegeneratePivots;
    numSplits = rhs->numSplits;
    numFeasSplits = rhs->numFeasSplits;
    num0feas = rhs->num0feas;
    num1feas = rhs->num1feas;
    num2feas = rhs->num2feas;
    numFloorInfeas = rhs->numFloorInfeas;
    numCeilInfeas = rhs->numCeilInfeas;
    // vectors
    nonBasicReducedCost = rhs->nonBasicReducedCost;
    raysOfC1 = rhs->raysOfC1;
    basicOrigVarIndex = rhs->basicOrigVarIndex;
    basicSlackVarIndex = rhs->basicSlackVarIndex;
    nonBasicVarIndex = rhs->nonBasicVarIndex;
    varBasicInRow = rhs->varBasicInRow;
    rowOfVar = rhs->rowOfVar;
    cstat = rhs->cstat;
    rstat = rhs->rstat;
    fractionalCore = rhs->fractionalCore;
    a0 = rhs->a0;
    splitFeasFlag = rhs->splitFeasFlag;
    fracCoreFloorFeasFlag = rhs->fracCoreFloorFeasFlag;
    fracCoreCeilFeasFlag = rhs->fracCoreCeilFeasFlag;
    lopsidedFlag = rhs->lopsidedFlag;
    feasSplitFracCoreIndex = rhs->feasSplitFracCoreIndex;
    feasSplitVar = rhs->feasSplitVar;
    infeasSplitVar = rhs->infeasSplitVar;
    lopsidedSplitVar = rhs->lopsidedSplitVar;
    count_intersecting_rays = rhs->count_intersecting_rays;
    count_intersecting_rays_floor = rhs->count_intersecting_rays_floor;
    count_intersecting_rays_ceil = rhs->count_intersecting_rays_ceil;
    rayDensitySortIndex = rhs->rayDensitySortIndex;
    rayDensity = rhs->rayDensity;
  } else {
    SPACE = 20;
    LPopt = 0.;
    objNorm = 0.;
    minAbsRayVal = GlobalVariables::param.getINFINITY();
    minAbsRayValNoEPS = GlobalVariables::param.getINFINITY();
    minAbsCoeff = 0.;
    maxAbsCoeff = 0.;
    minReferenceValue = 1;
    maxReferenceValue = 1;
    EPS = GlobalVariables::param.getEPS();
    numCols = 0;
    numRows = 0;
    numANonzero = 0;
    numNBOrig = 0;
    numNBSlack = 0;
    numNB = 0;
    numBasicCols = 0;
    numIntegerNonBinary = 0;
    numBinary = 0;
    numContinuous = 0;
    numFixed = 0;
    numEqRows = 0;
    numIneqRows = 0;
    numBoundRows = 0;
    numAssignRows = 0;
    numLB = 0;
    numUB = 0;
    numPrimalDegeneratePivots = 0;
    numDualDegeneratePivots = 0;
    numSplits = 0;
    numFeasSplits = 0;
    num0feas = 0;
    num1feas = 0;
    num2feas = 0;
    numFloorInfeas = 0;
    numCeilInfeas = 0;
  }
} /* copyOurStuff */

void SolutionInfo::initialize(OsiSolverInterface* const solver,
    const int max_frac_core, const bool brief) {
  copyOurStuff(NULL);
  OsiClpSolverInterface* clpsolver = NULL;
  try {
    clpsolver = dynamic_cast<OsiClpSolverInterface*>(solver);
  } catch (std::exception& e) {
    error_msg(errorstring,
        "There are several places I assume OsiClpSolverInterface, such as getting the basis status quickly. Edit those before proceeding.\n");
    writeErrorToII(errorstring, GlobalVariables::log_file);
    exit(1);
  }
  checkSolverOptimality(solver, false);
  enableFactorization(solver);

  if (!solver->isProvenOptimal()) {
    warning_msg(warnstr, "Solver to be used is not proven optimal.\n");
  }

  minAbsCoeff = computeMinAbs(solver->getMatrixByCol()->getNumElements(),
      solver->getMatrixByCol()->getElements());
  minReferenceValue = 1;
  maxReferenceValue = 1;
  minReferenceValue =
      (minAbsCoeff > 0) ?
          CoinMin(minReferenceValue, minAbsCoeff) : minReferenceValue;
  maxAbsCoeff = computeMaxAbs(solver->getMatrixByCol()->getNumElements(),
      solver->getMatrixByCol()->getElements());
  maxReferenceValue = CoinMax(maxReferenceValue, maxAbsCoeff);

  // Objective value
  LPopt = solver->getObjValue();
//  minReferenceValue =
//      (LPopt != 0) ?
//          CoinMin(minReferenceValue, std::abs(LPopt)) : minReferenceValue;
//  maxReferenceValue = CoinMax(maxReferenceValue, std::abs(LPopt));
  objNorm = std::sqrt(
      std::inner_product(solver->getObjCoefficients(),
          solver->getObjCoefficients() + solver->getNumCols(),
          solver->getObjCoefficients(), 0.0));

  // Information about data
  numCols = solver->getNumCols();
  numRows = solver->getNumRows();
  numANonzero = solver->getNumElements();

  // Get basic / non-basic information
//  basis = clpsolver->getModelPtr()->getBasis();
  cstat.resize(numCols);
  rstat.resize(numRows);
  solver->getBasisStatus(&cstat[0], &rstat[0]);
  getVariableAndRowInfo(solver);
  EPS = CoinMin(EPS, minReferenceValue / maxReferenceValue);

  // Initialize fracCoreFloorFeas and fracCoreCeilFeas
  numSplits = (int) fractionalCore.size();
  std::vector<double> splitCost(numSplits);
  for (int i = 0; i < numSplits; i++) {
    splitCost[i] = solver->getObjCoefficients()[fractionalCore[i]];
  }
  trimFractionalCore(splitCost, max_frac_core, true);
  numFeasSplits = numSplits; // may be adjusted from CglSIC

  feasSplitVar = fractionalCore;
  splitFeasFlag.resize(numSplits, true);
  feasSplitFracCoreIndex.resize(numSplits);
  for (int i = 0; i < numFeasSplits; i++) {
    feasSplitFracCoreIndex[i] = i;
  }
  fracCoreFloorFeasFlag.resize(numSplits, false);
  fracCoreCeilFeasFlag.resize(numSplits, false);
  lopsidedFlag.resize(numSplits, false);

  checkFracCoreFeas(solver);

  // Save rays of C1, where the coefficients of inv(B) * A
  // are sometimes negated because we have
  // \bar x = inv(B) * b - inv(B) * A * x_N
  int tempIndex = 0;
  const char* rowSense = solver->getRowSense();
  if (!brief) {
    raysOfC1.resize(numNB);
//    rayMultiplier.resize(numNB, 1.);
  }
  nonBasicReducedCost.resize(numNB);
  for (int j = 0; j < numNB; j++) {
    const int NBVar = nonBasicVarIndex[j];
    if (!brief) {
      raysOfC1[j].resize(numRows);
      solver->getBInvACol(NBVar, &raysOfC1[j][0]);

      // We are looking at inv(B)*A, jth col
      // Call it \bar a^j
      // Moving it to the right side we get \bar a^0 - \bar a^j
      // So call ray j r^j := -\bar a^j
      // (EXCEPT for the ub vars, because for those, we actually
      // proceed along the negative direction).
      // For >= slacks, Clp has them as non-positive variables
      // (enforcing +1 coefficient when adding slacks, but this
      // destroys the property that only cuts with rhs 1 exist in nb space),
      // whereas we treat them as non-NEGATIVE variables.
      // Thus, the basis inverse entry is already the negation of what we'd see in
      // the tableau, so we do not need to do anything.
      if (isNonBasicLBVar(clpsolver, NBVar)) { // This checks if it is at the lower-bound
        for (int k = 0; k < numRows; k++) {
          raysOfC1[j][k] *= -1.0;
        }
      } else {
//        rayMultiplier[j] = -1.;
      }
    } /* not brief */

    if (NBVar < numCols) {
      if (isNonBasicUBCol(clpsolver, NBVar))
        nonBasicReducedCost[j] = -1.0 * solver->getReducedCost()[NBVar];
      else
        nonBasicReducedCost[j] = solver->getReducedCost()[NBVar];
    } else {
      tempIndex = NBVar - numCols;
      if (isNonBasicUBSlack(clpsolver, tempIndex))
        nonBasicReducedCost[j] = solver->getRowPrice()[tempIndex];
      else
        nonBasicReducedCost[j] = -1.0 * solver->getRowPrice()[tempIndex];
    }
    if (lessThanVal(nonBasicReducedCost[j], 0.)) {
      if (lessThanVal(nonBasicReducedCost[j], 0., 1e-3)) {
        error_msg(errorstring,
            "Non-basic reduced cost should be >= 0 in the complemented non-basic space. " "However, for nb var %d (real var %d), it is %e.\n",
            j, NBVar, nonBasicReducedCost[j]);
        writeErrorToII(errorstring, GlobalVariables::log_file);
        exit(1);
      } else {
        warning_msg(warnstring,
            "Non-basic reduced cost should be >= 0 in the complemented non-basic space. " "However, for nb var %d (real var %d), it is %f. Small enough error that we only send warning.\n",
            j, NBVar, nonBasicReducedCost[j]);
        GlobalVariables::numCutSolverFails[CutSolverFails::NUMERICAL_ISSUES_WARNING_NO_OBJ_FAIL_IND]++;
      }
    }

    if (isZero(nonBasicReducedCost[j], EPS)) {
      numDualDegeneratePivots++;
    }
  } /* loop over non-basic vars to collect rays and reduced costs */

  if (!brief) {
    // For the rows that are >=, if the variable basic in the row is slack,
    // recall that clp has that as a non-positive variable,
    // but we want it to be non-negative.
    // Thus, we multiply the row through by minus one.
    // That is, the row is s = b_i + r x, but the "real" situation is e = -b_i - r x.
    // Hence, we multiply through by -1.
    for (int r = 0; r < numRows; r++) {
      if (varBasicInRow[r] >= numCols) {
        if (rowSense[varBasicInRow[r] - numCols] == 'G') {
          for (int j = 0; j < numNB; j++) {
            raysOfC1[j][r] *= -1.0;
          }
        }
      }
    }

    // TODO This could cause problems if we are changing a value that is supposed to be non-zero to zero
    // Conversely, we could have problems if a zero value appears as non-zero, so both ways may have issues
    // Because of numerical issues, go through all rows and explicitly make zero all near-zero elements
    for (int j = 0; j < numNB; j++) {
      for (int r = 0; r < numRows; r++) {
        if (!isZero(raysOfC1[j][r], 0.0)
            && lessThanVal(std::abs(raysOfC1[j][r]), minAbsRayValNoEPS, 0.0)) {
          minAbsRayValNoEPS = std::abs(raysOfC1[j][r]);
        }
        if (isZero(raysOfC1[j][r], param.getRAYEPS())) { // TODO should maybe check against maxReferenceValue
          raysOfC1[j][r] = 0.0;
        }
        if (!isZero(raysOfC1[j][r], 0.0)
            && lessThanVal(std::abs(raysOfC1[j][r]), minAbsRayVal, 0.0)) {
          minAbsRayVal = std::abs(raysOfC1[j][r]);
        }
      }
    }
  } /* check if !brief */

//  if (solverFactorizationStatus == 0) {
    solver->disableFactorization();
//  }
} /* initialize */

/**
 * @brief This should probably do some freeing.
 */
SolutionInfo::~SolutionInfo() {
//  if (basis) {
//    delete basis;
//  }
} /* destructor */

/**
 * @brief Counts the number of rays intersecting the bd of each split
 */
void SolutionInfo::countIntersectingRays(const std::vector<int>& varIndex) {
  if (raysOfC1.empty()) {
    return;
  }

  count_intersecting_rays.resize(varIndex.size(), 0);
  count_intersecting_rays_floor.resize(varIndex.size(), 0);
  count_intersecting_rays_ceil.resize(varIndex.size(), 0);
  for (int split = 0; split < (int) varIndex.size(); split++) {
    int currSplit = varIndex[split];
    int currSplitRow = rowOfVar[currSplit];
//    double pi0 = floor(primalSoln[currSplit]);
    for (int j = 0; j < numNB; j++) {
//      int currCol = nonBasicVarIndex[j];
      double rayGrad = raysOfC1[j][currSplitRow];
      //Compute slack ratio
//      double splitSR;
      if (greaterThanVal(rayGrad, 0.0)) {
//        splitSR = (pi0 + 1 - primalSoln[currSplit]) / (rayGrad);
        count_intersecting_rays[split]++;
        count_intersecting_rays_ceil[split]++;
      } else {
        if (lessThanVal(rayGrad, 0.0)) {
//          splitSR = (pi0 - primalSoln[currSplit]) / (rayGrad);
          count_intersecting_rays[split]++;
          count_intersecting_rays_floor[split]++;
        } else {
//          splitSR = param.getINFINITY();
        }
      }
    }
  }
}

void SolutionInfo::trimFractionalCore(const std::vector<double>& splitCost,
    const int max_frac_core, const bool shouldSort) {
  std::vector<int> sortIndex;
  if (shouldSort) {
    sortIndex.resize(numSplits);
    for (int i = 0; i < numSplits; i++) {
      sortIndex[i] = i;
    }
    std::sort(sortIndex.begin(), sortIndex.end(), index_cmp<const double*>(&splitCost[0]));
  }
  std::vector<int> tmpFracCore;
  tmpFracCore.reserve(numSplits);
  for (int s = 0; s < numSplits; s++) {
    if ((max_frac_core > 0)
        && ((int) tmpFracCore.size() >= max_frac_core)) {
      break;
    }
    if (shouldSort) {
      tmpFracCore.push_back(fractionalCore[sortIndex[s]]);
    } else {
      tmpFracCore.push_back(fractionalCore[s]);
    }
  }
//  tmpFracCore.clear(); tmpFracCore.push_back(4);
  numSplits = tmpFracCore.size();
  fractionalCore = tmpFracCore;
}

void SolutionInfo::sortFractionalCore(const std::vector<double>& splitCost) {
  std::vector<int> sortIndex;
  sortIndex.resize(numSplits);
  for (int i = 0; i < numSplits; i++) {
    sortIndex[i] = i;
  }
  sort(sortIndex.begin(), sortIndex.end(),
      index_cmp<const double*>(&splitCost[0]));
  std::vector<int> tmpFracCore;
  tmpFracCore.reserve(numSplits);
  for (int s = 0; s < numSplits; s++) {
      tmpFracCore.push_back(fractionalCore[sortIndex[s]]);
  }
  fractionalCore = tmpFracCore;
}

void SolutionInfo::checkFracCoreFeas(const OsiSolverInterface* const tmpsolver,
    const bool working_in_subspace) {
  if (!working_in_subspace) {
    return;
  }

  // Now we go through each split, and force each split to be
  // at its lower bound, and at its upper bound. We count how many
  // of those are feasible.
  num0feas = 0;
  num1feas = 0;
  num2feas = 0;
  numFloorInfeas = 0;
  numCeilInfeas = 0;

  PHASolverInterface* solver = dynamic_cast<PHASolverInterface*>(tmpsolver->clone());
  setMessageHandler(solver);
  int numrows = solver->getNumRows();
  const double* colSolution = tmpsolver->getColSolution();

  for (int i = 0; i < numSplits; i++) {
    const int split_var = fractionalCore[i];

    // If we are working in the subspace, then check feasibility
//    bool toAdd = false;
    if (!working_in_subspace) {
      //      toAdd = true;
    } else {
      int numvalid = 0;
      double currval = colSolution[split_var];
      std::vector<int> cols;
      std::vector<double> vals;
      double rowlb, rowub;

      cols.push_back(split_var);
      vals.push_back(1);
      rowlb = floor(currval);
      rowub = floor(currval);

      solver->addRow(1, cols.data(), vals.data(), rowlb, rowub);
      solver->initialSolve();

      if (solver->isProvenOptimal()) {
        numvalid++;
        fracCoreFloorFeasFlag[split_var] = true;
      } else {
        numFloorInfeas++;
      }

      std::vector<int> rowsToDelete(1);
      rowsToDelete[0] = numrows;
      solver->deleteRows(1, rowsToDelete.data());

      rowlb = ceil(currval);
      rowub = ceil(currval);

      solver->addRow(1, cols.data(), vals.data(), rowlb, rowub);
      solver->initialSolve();

      if (solver->isProvenOptimal()) {
        numvalid++;
        fracCoreCeilFeasFlag[i] = true;
      } else {
        numCeilInfeas++;
      }

      solver->deleteRows(1, rowsToDelete.data());

      if (numvalid == 0) {
        num0feas++;
        //        toAdd = USE_EMPTY;
      } else if (numvalid == 1) {
        num1feas++;
        lopsidedFlag[i] = true;
        lopsidedSplitVar.push_back(split_var);
        //        toAdd = YESLOP;
      } else {
        num2feas++;
        //        toAdd = true;
      }
    }
  }

  // Free memory
  if (solver) {
    delete solver;
  }
}

/**
 * @brief Prints brief solution info.
 */
void SolutionInfo::printBriefSolutionInfo(
    const OsiSolverInterface* const solver, const std::string outputFileName) {
  std::string ext = "-info.alex";
//#ifdef TRACE
//  std::cout
//    << "\n############### Starting printBriefSolutionInfo routine to generate "
//    << outputFileName + ext << ". ###############" << std::endl;
//#endif

  FILE *fptr = fopen((outputFileName + ext).c_str(), "w");
  if (fptr == NULL) {
    error_msg(errorstring,
        "While printing brief solution info for %s encountered error.\n",
        outputFileName.c_str());
    writeErrorToII(errorstring, GlobalVariables::log_file);
    exit(1);
  }

  fprintf(fptr, "Solution information (brief):\n\n");
  fprintf(fptr, "%-*s%-*d\n", SPACE, "ROWS:", SPACE, solver->getNumRows());
  fprintf(fptr, "%-*s%-*d\n", SPACE, "COLUMNS:", SPACE, solver->getNumCols());
  fprintf(fptr, "%-*s%-*d\n", SPACE, "INT NON-BIN:", SPACE,
      numIntegerNonBinary);
  fprintf(fptr, "%-*s%-*d\n", SPACE, "BINARY:", SPACE, numBinary);
  fprintf(fptr, "%-*s%-*d\n", SPACE, "NUM EQ ROWS:", SPACE, numEqRows);
  fprintf(fptr, "%-*s%-*d\n", SPACE, "NUM INEQ ROWS:", SPACE, numIneqRows);
  fprintf(fptr, "%-*s%-*d\n", SPACE, "NUM BOUND ROWS:", SPACE, numBoundRows);
  fprintf(fptr, "%-*s%-*d\n", SPACE, "NUM ASSIGN ROWS:", SPACE,
      numAssignRows);
  fprintf(fptr, "%-*s%-*f\n", SPACE, "LP OPT:", SPACE, solver->getObjValue());
  fprintf(fptr, "%-*s%-*d\n", SPACE, "FRAC CORE:", SPACE,
      (int) fractionalCore.size());
  //  fprintf(fptr, "%-*s%-*d\n", SPACE, "ORIG BASICS:", SPACE,
  //      (int) basicOrigVarIndex.size());
  fprintf(fptr, "%-*s%-*d\n", SPACE, "ORIG NON-BASICS:", SPACE,
      (int) nonBasicOrigVarIndex.size());
  //  fprintf(fptr, "%-*s%-*d\n", SPACE, "SLACK BASICS:", SPACE,
  //      (int) basicSlackVarIndex.size());
  fprintf(fptr, "%-*s%-*d\n", SPACE, "SLACK NON-BASICS:", SPACE,
      (int) nonBasicSlackVarIndex.size());

  fclose(fptr);

#ifdef TRACE
  std::cout << "\n## Brief solution written to " << outputFileName + ext << ". ##"
    << std::endl;
#endif
}

void SolutionInfo::getVariableAndRowInfo(const OsiSolverInterface* const solver/*, const bool brief*/) {
  const double* rhs = solver->getRightHandSide();
  const double* rowActivity = solver->getRowActivity();
  const double* tempPrimalSoln = solver->getColSolution();
  //  const double* tempDualSoln = solver->getRowPrice();
  //  const double* tempOrigRedCosts = solver->getReducedCost();

  varBasicInRow.resize(numRows);
  rowOfVar.clear();
  //rowOfSlack.clear();
  rowOfVar.resize(numCols + numRows, -1);
  //rowOfSlack.resize(numRows, -1);
  solver->getBasics(&varBasicInRow[0]); // Get which variable is basic in each row

  // Clear data structures
  basicOrigVarIndex.clear();
  basicSlackVarIndex.clear();
  nonBasicVarIndex.clear();
  nonBasicOrigVarIndex.clear();
  nonBasicSlackVarIndex.clear();
//  nonBasicFixedVarIndex.clear();
  a0.clear();

  // Allocate space
  basicOrigVarIndex.reserve(numCols);
  basicSlackVarIndex.reserve(numRows);
  nonBasicVarIndex.reserve(numCols);
  nonBasicOrigVarIndex.reserve(numCols);
  nonBasicSlackVarIndex.reserve(numRows);
  a0.resize(numRows);

  // Get the basic and non-basic original variable info
  // Note that fixed non-basic variables might not need to be added to the non-basic index vector
  // since they do not correspond to any rays...
  // but right now we typically assume we get n rays in some parts of the code
  for (int col = 0; col < numCols; col++) {
    // Count how many non-inf lower and upper bounds
    double colLB, colUB;
    colLB = solver->getColLower()[col];
    colUB = solver->getColUpper()[col];
    if (colLB > -1 * solver->getInfinity() + param.getEPS()) {
      numLB++;
      minReferenceValue =
          (colLB != 0.) ?
              CoinMin(minReferenceValue, std::abs(colLB)) : minReferenceValue;
      maxReferenceValue = CoinMax(maxReferenceValue, std::abs(colLB));
    }

    if (colUB < solver->getInfinity() - param.getEPS()) {
      numUB++;
      minReferenceValue =
          (colUB != 0.) ?
              CoinMin(minReferenceValue, std::abs(colUB)) : minReferenceValue;
      maxReferenceValue = CoinMax(maxReferenceValue, std::abs(colUB));
    }

    if (solver->isBinary(col)) {
      numBinary++;
    } else if (solver->isContinuous(col)) {
      numContinuous++;
    } else if (solver->isIntegerNonBinary(col)) {
      numIntegerNonBinary++;
    }

    if (isBasicVar(col)) {
      basicOrigVarIndex.push_back(col);
      if (solver->isInteger(col)) {
        if (!isVal(tempPrimalSoln[col], std::floor(tempPrimalSoln[col]),
              param.getAWAY())
            && !isVal(tempPrimalSoln[col],
              std::ceil(tempPrimalSoln[col]), param.getAWAY())) {
          fractionalCore.push_back(col);
        }
      }

      // If it is at its lower or upper bound, it is primal degenerate
      if (isVal(solver->getColSolution()[col], colLB)
          || isVal(solver->getColSolution()[col], colUB)) {
        numPrimalDegeneratePivots++;
      }
    } else {
      // Recall that rowOfVar stores -1 - nb index for nb variables
      // The -1 is to prevent the conflict of the 0th nb var and var basic in row 0
      rowOfVar[col] = -1 - (int) nonBasicVarIndex.size();
      nonBasicOrigVarIndex.push_back(col);
      numNBOrig++;
      nonBasicVarIndex.push_back(col);
      
      // If it is a fixed variable, we should be able to handle it more efficiently...
      if (isVal(colLB, colUB)) {
//        nonBasicFixedVarIndex.push_back(col);
        numFixed++;
      }
    }
  } /* iterate over columns */

  // Here only store the non-basic slacks corresponding to inequalities since
  // the equality slack columns don't correspond to any rays
  const CoinPackedMatrix* mx = solver->getMatrixByRow();
  const char* rowSense = solver->getRowSense();
  for (int row = 0; row < numRows; row++) {
    const double absrhs = std::abs(solver->getRightHandSide()[row]);
    minReferenceValue =
        (absrhs > 0) ? CoinMin(minReferenceValue, absrhs) : minReferenceValue;
    maxReferenceValue = CoinMax(maxReferenceValue, absrhs);
    //    slackVar.push_back(rhs[row] - rowActivity[row]);
    //    if (rowSense[row] == 'G')
    //      slackVar[row] *= -1.0;
    //    dualSoln.push_back(tempDualSoln[row]);
    //
    //if (basicSlack(row)) { // rstat[row] == 0 should also be basic? Seemingly not from ClpSimplex model status definitions.
    if (isBasicVar(row + numCols)) {
      basicSlackVarIndex.push_back(row + numCols);

      // If the slack in the row is zero, then this is a primal degenerate pivot
      if (isVal(solver->getRowActivity()[row],
          solver->getRightHandSide()[row])) {
        numPrimalDegeneratePivots++;
      }
    } else {
      rowOfVar[row + numCols] = -1 - (int) nonBasicVarIndex.size();
//      if (rowSense[row] != 'E') {
        nonBasicSlackVarIndex.push_back(row + numCols);
//      }
      numNBSlack++;
      nonBasicVarIndex.push_back(row + numCols);
    //  }
    }
    if (rowSense[row] == 'E') {
      //    nonBasicEqualitySlackVarIndex.push_back(row + numcols);
      numEqRows++;
    }

    // Also decide how many numEqRows, numIneqRows there are
    if (mx->getVectorSize(row) == 1) {
      if (rowSense[row] == 'E') {
        numAssignRows++;
      } else {
        numBoundRows++;
      }
    }

    const int var = varBasicInRow[row];
    const int currrow = var - numCols;
    if (isBasicVar(row + numCols) && currrow != row) {
      // Elsewhere we are using that each slack is basic in its own row,
      // so if this is not true, we will have to adjust
      // We use this, for example, in PCut, to calculate the right order
      // for the packed NB rays in genCornerNB
      error_msg(errstr,
          "Basic variable in row %d is variable %d, but it should be the slack on this row (%d).\n",
          row, var, row + numCols);
      writeErrorToII(errstr, GlobalVariables::log_file);
      exit(1);
    }
    if (var < numCols) {
      rowOfVar[var] = row;
      a0[row] = (tempPrimalSoln[var]);
    } else {
      rowOfVar[var] = row;
      //rowOfSlack[currrow] = row;
      if (rowSense[row] == 'G') {
        a0[row] = (rowActivity[row] - rhs[row]);
      } else {
        a0[row] = (rhs[row] - rowActivity[row]);
      }
    }
    // 2017-08-23: The minReferenceValue set from the solution is going to be far too small typically, with epsilons everywhere
//    minReferenceValue =
//        (a0[row] != 0) ?
//            CoinMin(minReferenceValue, std::abs(a0[row])) : minReferenceValue;
//    maxReferenceValue = CoinMax(maxReferenceValue, std::abs(a0[row]));
  }
  numIneqRows = numRows - numEqRows;
  numBasicCols = (int) basicOrigVarIndex.size()
    + (int) basicSlackVarIndex.size();
  numNB = (int) nonBasicVarIndex.size();
} /* getVariableAndRowInfo */

bool SolutionInfo::isBasicVar(const int var) const {
  if (var < numCols) {
    return isBasicCol(this->cstat, var);
  } else {
    return isBasicSlack(this->rstat, var - numCols);
  }
}

bool SolutionInfo::isNBUBVar(const int var) const {
  if (var < numCols) {
    return isNonBasicUBCol(this->cstat, var);
  } else {
    return isNonBasicUBSlack(this->rstat, var - numCols);
  }
}

bool SolutionInfo::isNBLBVar(const int var) const {
  if (var < numCols) {
    return isNonBasicLBCol(this->cstat, var);
  } else {
    return isNonBasicLBSlack(this->rstat, var - numCols);
  }
}
