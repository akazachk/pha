///============================================================================
// Name        : SolutionInfo.hpp
// Author      : Aleksandr M. Kazachkov (edited from Selvaprabu Nadarajah)
// Version     : 0.2013.06.16
// Copyright   : Your copyright notice
// Description : Information about a problem
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================
#pragma once

#include <vector>
#include "typedefs.hpp"
#include "GlobalConstants.hpp"

class SolutionInfo {
public:
  /** INFORMATION THAT WE WILL KEEP **/
//  CoinWarmStartBasis basis;
  int SPACE;
  std::vector<double> nonBasicReducedCost;
  std::vector<std::vector<double> > raysOfC1;
  double LPopt;
  double objNorm;
  double minAbsRayValNoEPS, minAbsRayVal;
  double minAbsCoeff, maxAbsCoeff;
  double minReferenceValue, maxReferenceValue;
  double EPS; // this will be used to decide what epsilon is used for zeroing out coefficients and such

  /*
   * basicOrigVarIndex              ::  For the kth basic orig variable, what is its original index
   * nonBasicOrigVarIndex           ::  For the kth non-basic orig var, what is it its original index
   * basicSlackVarIndex             ::  For the kth basic slack var, what is it its original index
   * nonBasicSlackVarIndex          ::  For the kth non-basic slack [not artificial] var, the original index
   * nonBasicEqualitySlackVarIndex  ::  For the kth non-basic artificial var, the original index
   * nonBasicVarIndex               ::  The non-basic variables and their indices (first original then slack [not artificial])
   * fractionalCore                 ::  The set of original integer variables with non-integer current vals
   * rowOfVar                       ::  The row of each original variable (-(1 + nb index) if not basic)
   * varBasicInRow                  ::  For each row, which variable is basic in it
   * basicOrigVarRowIndex           ::  For each basic orig variable, what row is it basic in
   * basicSlackVarRowIndex          ::  For each basic slack var, what row is it basic in
   * cstat                          ::  For each variable, whether it is free (0), basic (1), nb at upper-bound (2), or nb at lb (3)
   * rstat                          ::
   * structVarWithUBs               ::
   * mapNonBasicStructToJ           ::
   */
  std::vector<int> basicOrigVarIndex, basicSlackVarIndex;
  std::vector<int> nonBasicVarIndex, nonBasicOrigVarIndex,
      nonBasicSlackVarIndex, nonBasicFixedVarIndex;
  std::vector<int> varBasicInRow; //, mapStructVarIndexToNBVarIndex,
  // Note, what we store in rowOfVar for nb vars will always be negative,
  // but as a result, will not exactly be its nb index
  // (the true nb index will -1 * (rowOfVar[i] + 1))
  std::vector<int> rowOfVar; // This will keep the row that var is basic in, or (-1 - its non-basic index)

  std::vector<int> cstat, rstat;
  std::vector<int> fractionalCore;

  std::vector<double> a0;
  int numCols, numRows, numANonzero;
  int numNBOrig, numNBSlack, numNB;
  int numBasicCols, numIntegerNonBinary, numBinary, numContinuous, numFixed; 
  int numEqRows, numIneqRows, numBoundRows, numAssignRows, numLB, numUB;

  int numPrimalDegeneratePivots, numDualDegeneratePivots;

  // Things added to ensure the fractional core we chose is feasible
  int numSplits, numFeasSplits;
  int num0feas, num1feas, num2feas, numFloorInfeas, numCeilInfeas;

  // This is a vector of booleans, one per split, to check if that is feasible.
  // I.e. if i \in feasSplitSub, splitFeasSub[i] = true. Else, false.
  std::vector<bool> splitFeasFlag;
  std::vector<bool> fracCoreFloorFeasFlag, fracCoreCeilFeasFlag, lopsidedFlag;

  // We need this so that we don't try to compile stats on the infeasible splits.
  // infeasSplit is simply the complement of the set of feasSplitVar indices.
  std::vector<int> feasSplitFracCoreIndex;
  std::vector<int> feasSplitVar, infeasSplitVar, lopsidedSplitVar;

  // For counting the number of intersecting rays
  std::vector<int> count_intersecting_rays, count_intersecting_rays_floor,
      count_intersecting_rays_ceil;

  // Number of the desired (subspace split) rows each ray intersects
  std::vector<int> rayDensitySortIndex;
  std::vector<int> rayDensity;
  
  void copyOurStuff(const SolutionInfo* const rhs);

  /** CONSTRUCTORS **/
  SolutionInfo() {
    copyOurStuff(NULL);
  }
  ;

  SolutionInfo(OsiSolverInterface* const solver, const bool brief = true) :
      SolutionInfo(solver, solver->getNumCols(), brief) {
  }
  ;

  SolutionInfo(OsiSolverInterface* const solver, const int max_frac_core,
      const bool brief);
  void initialize(OsiSolverInterface* const solver, const int max_frac_core,
      const bool brief = true);

  ~SolutionInfo();
  
  /** OTHER METHODS **/
  /**
   * @brief Extracts information about basic and non-basic variables,
   *        as well as how many inequality, non-bound rows we have.
   *        Also identify the fractional core.
   *        Also get number integer / binary / continuous
   */
  void getVariableAndRowInfo(const OsiSolverInterface* const solver/*, const bool brief*/);

  /**
   * @brief Counts the number of rays intersecting the bd of each split
   */
  void countIntersectingRays(const std::vector<int> &index);

  void trimFractionalCore(const std::vector<double>& splitCost,
      const int max_frac_core, const bool shouldSort = true);

  void sortFractionalCore(const std::vector<double>& splitCost);

  void checkFracCoreFeas(const OsiSolverInterface* const tmpsolver,
      const bool working_in_subspace = false);
  
  bool isBasicVar(const int var) const;
  bool isNBUBVar(const int var) const;
  bool isNBLBVar(const int var) const;

  /** OUTPUT **/
  void printBriefSolutionInfo(const OsiSolverInterface* const solver,
      const std::string outputFileName);

  /** INLINE FUNCTIONS **/
  inline const std::vector<int>& getCstat() const {
    return cstat;
  }
  inline const std::vector<int>& getRstat() const {
    return rstat;
  }
  inline const std::vector<int>& getFractionalCore() const {
    return fractionalCore;
  }
  
  /*
  inline bool basicCol(int col) const {
    return isBasicCol(cstat, col);
  }

  inline bool NBUBCol(int col) const {
    return isNonBasicUBCol(cstat, col);
  }

  inline bool NBLBCol(int col) const {
    return isNonBasicLBCol(cstat, col);
  }

  inline bool basicSlack(int row) const {
    return isBasicSlack(rstat, row);
  }

  inline bool NBUBSlack(int row) const {
    return isNonBasicUBSlack(rstat, row);
  }

  inline bool NBLBSlack(int row) const {
    return isNonBasicLBSlack(rstat, row);
  }
  */
  
  inline int getVarNBIndex(const int var) const {
    if ((var < numCols + numRows) && (rowOfVar[var] <= -1)) {
      return -1 * (1 + rowOfVar[var]);
    } else {
      return -1;
    }
  }
};
