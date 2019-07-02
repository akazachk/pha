//============================================================================
// Name        : Hplane.cpp
// Author      : Aleksandr M. Kazachkov
// Version     : 0.2014.05.23
// Copyright   : Your copyright notice
// Description : Hyperplane class
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#include "Hplane.hpp"
#include "GlobalConstants.hpp"
#include "Utility.hpp"

Hplane::Hplane(const int tmpVar, const int tmpRow, const int tmpUBFlag,
    const double tmpBound, const int numSplits) {
  var = tmpVar;
  row = tmpRow;
  ubflag = tmpUBFlag;
  bound = tmpBound;

  resizeVectors(numSplits);
}

Hplane::~Hplane() {

}

void Hplane::resizeVectors(const int numSplits) {
  if (numSplits > 0) {
    vertOfRayToBeCutByHplane.resize(numSplits);
    rayToBeCutByHplane.resize(numSplits);
    vertCreatedByRayToBeCutByHplane.resize(numSplits);
    vertOfAllRaysIntByHplane.resize(numSplits);
    allRaysIntByHplane.resize(numSplits);
    numRaysCutByHplaneInclPrevAct.resize(numSplits);
//    for (int s = 0; s < numSplits; s++) {
//      const int old_size_vec = numRaysCutByHplaneInclPrevAct[s].size();
//      const int prevTotal =
//          (old_size_vec > 0) ?
//              numRaysCutByHplaneInclPrevAct[s][old_size_vec - 1] :
//              0;
//      numRaysCutByHplaneInclPrevAct[s].resize(num_hplane_act, prevTotal);
//    }
  }
}

/**
 * Based on UB flag, set bound of the hplane
 */
void Hplane::setHplaneBoundFromUBFlag(
    const PHASolverInterface* const solver) {
  // Which side of the hyperplane are we trying to hit?
  if (ubflag == static_cast<int>(HplaneBoundFlag::toLBNB)) {
    bound = 0.0;
  } else if (ubflag == static_cast<int>(HplaneBoundFlag::toLBSlack)) {
    bound = 0.0;
  } else if (ubflag == static_cast<int>(HplaneBoundFlag::toLB)) {
    bound = getVarLB(solver, var);
  } else if (ubflag == static_cast<int>(HplaneBoundFlag::toUB)) {
    bound = getVarUB(solver, var);
  } else if (ubflag == static_cast<int>(HplaneBoundFlag::toUBNB)) {
    const double UB = getVarUB(solver, var);
    const double LB = getVarLB(solver, var);
    if (!isInfinity(UB, solver->getInfinity())
        && !isNegInfinity(LB, solver->getInfinity())) {
      bound = UB - LB;
    }
  }
}

/**
 * distanceToHplane 0
 * Given current value varVal and the value of the ray rayVal for this hplane,
 * calculate the distance to the current hyperplane
 */
double Hplane::distanceToHplane(const double varVal, const double rayVal,
    const PHASolverInterface* const solver) const {
  // Ensure bound is finite
  if (isInfinity(std::abs(bound), solver->getInfinity())
      || isZero(rayVal, param.getRAYEPS())) {
    return solver->getInfinity();
  }

  // Which side of the hyperplane are we trying to hit?
  // This is important to check because if bound and varVal are equal,
  // we will get a distance of zero even if the ray is in the wrong direction
  if (greaterThanVal(rayVal, 0.0, param.getRAYEPS())
      && !isVal(bound, getVarUB(solver, var))) {
    return solver->getInfinity();
  } else if (lessThanVal(rayVal, 0.0, param.getRAYEPS())
      && !isVal(bound, getVarLB(solver, var))) {
    return solver->getInfinity();
  }

  // Find the distance to the hyperplane, and ensure it is actually intersected
  const double dist = (bound - varVal) / rayVal;

  if (isInfinity(dist, solver->getInfinity()) || lessThanVal(dist, 0.0)) {
    return solver->getInfinity(); // ray does not intersect this hyperplane
  } else if (isZero(dist)) {
    return 0.0;
  } else {
    return dist;
  }
} /* distanceToHplane 0 */

/**
 * distanceToHplane 1
 * Assume starting at optimal vertex
 * The ray is given by the variable (in the structural space)
 */
double Hplane::distanceToHplane(const int structRayOutVar,
    const PHASolverInterface* const solver) const {
  if (var < solver->getNumCols()) {
    return distanceToHplane(solver->getColSolution()[var],
        structRayOutVar, solver);
  } else {
    const int currRow = var - solver->getNumCols();
    return distanceToHplane(
        std::abs(
            solver->getRowActivity()[currRow]
                - solver->getRightHandSide()[currRow]), structRayOutVar, solver);
  }
} /* distanceToHplane 1 */

/**
 * distanceToHplane 2
 * Given current value varVal of the hyperplane's variable
 * and the ray given by the variable (in the struct space),
 * calculate the distance to the current hyperplane
 */
double Hplane::distanceToHplane(const double varVal, const int structRayOutVar,
     const PHASolverInterface* const solver) const {
  // Ensure bound is finite
  if (isInfinity(std::abs(bound), solver->getInfinity())) {
    return solver->getInfinity();
  }

  // Get ray element
  double rayVal;
  const int rayMult = isNonBasicUBVar(solver, structRayOutVar) ? -1 : 1;
  if (row >= 0) {
//    const int solverFactorizationStatus = solver->canDoSimplexInterface();
//    if (solverFactorizationStatus == 0) {
    solver->enableFactorization();
//    }
    std::vector<double> BInvA(solver->getNumRows());
    solver->getBInvACol(structRayOutVar, &BInvA[0]);
//    if (solverFactorizationStatus == 0) {
    solver->disableFactorization();
//    }
    rayVal = rayMult * -1 * BInvA[row];
  } else if (var == structRayOutVar) {
    rayVal = rayMult;
  } else {
    return solver->getInfinity(); // Ray and hyperplane do not intersect
  }

  return distanceToHplane(varVal, rayVal, solver);
} /* distanceToHplane 2 */

/**
 * distanceToHplane 3
 * Assume starting at optimal vertex
 * The ray is given by the CoinPackedVector rayOut
 * This could be in the structural space or the non-basic space
 * If it is in the non-basic space, then the set of non-basic variables needs to be provided
 */
double Hplane::distanceToHplane(const CoinPackedVector& rayOut,
    const PHASolverInterface* const solver, const bool inNBSpace,
    const std::vector<int>* nonBasicVarIndex) const {
  if (var < solver->getNumCols()) {
    return distanceToHplane(solver->getColSolution()[var], rayOut, solver,
        inNBSpace, nonBasicVarIndex);
  } else {
    const int currRow = var - solver->getNumCols();
    return distanceToHplane(
        std::abs(
            solver->getRowActivity()[currRow]
                - solver->getRightHandSide()[currRow]), rayOut, solver,
        inNBSpace, nonBasicVarIndex);
  }
} /* distanceToHplane 3 */

/**
 * distanceToHplane 4
 * Starting at given vertex, which could be provided in the struct or non-basic space
 * The ray is given by the CoinPackedVector rayOut
 * This also could be in either the struct or non-basic space
 */
double Hplane::distanceToHplane(const CoinPackedVector &vertex,
    const CoinPackedVector& rayOut, const PHASolverInterface* const solver,
    const bool inNBSpace, const std::vector<int>* nonBasicVarIndex) const {
  // Ensure bound is finite
  if (isInfinity(std::abs(bound), solver->getInfinity())) {
    return solver->getInfinity();
  }

  double varVal = 0.0;
  if (inNBSpace) {
    varVal = getVarVal(solver, var)
        + dotProductWithHplane(vertex, solver, false, inNBSpace,
            nonBasicVarIndex);
  } else { // Else vectors are in the structural space
    if (var < solver->getNumCols()) {
      varVal = vertex[var];
    } else {
      // We want to know the value of the slack on this hyperplane
      const int curr_row = var - solver->getNumCols();
      varVal = std::abs(
          getRowActivity(vertex, solver, false)
              - solver->getRightHandSide()[curr_row]);
    }
  }
  return distanceToHplane(varVal, rayOut, solver, inNBSpace, nonBasicVarIndex);
} /* distanceToHplane 4 */

/**
 * distanceToHplane 4 WITH SOLUTION INFO
 * Starting at given vertex, which could be provided in the struct or non-basic space
 * The ray is given by the CoinPackedVector rayOut
 * This also could be in either the struct or non-basic space
 */
double Hplane::distanceToHplane(const CoinPackedVector &vertex,
    const CoinPackedVector& rayOut, const PHASolverInterface* const solver,
    const bool inNBSpace, const SolutionInfo& solnInfo) const {
  // Ensure bound is finite
  if (isInfinity(std::abs(bound), solver->getInfinity())) {
    return solver->getInfinity();
  }

  double varVal = 0.0;
  if (inNBSpace) {
    varVal = getVarVal(solver, var)
        + dotProductWithHplane(vertex, solnInfo, false);
  } else { // Else vectors are in the structural space
    if (var < solver->getNumCols()) {
      varVal = vertex[var];
    } else {
      // We want to know the value of the slack on this hyperplane
      const int curr_row = var - solver->getNumCols();
      varVal = std::abs(
          getRowActivity(vertex, solver, false)
              - solver->getRightHandSide()[curr_row]);
    }
  }

  return distanceToHplane(varVal, rayOut, solver, /*inNBSpace,*/ solnInfo); // distToHplane 5.2
} /* distanceToHplane 4 */

/**
 * distanceToHplane 5.1
 * Given current value varVal of the hyperplane's variable
 * and the ray given by the CoinPackedVector rayOut,
 * which could be in either the struct or non-basic space,
 * calculate the distance to the current hyperplane
 * If it is in the non-basic space, then the set of non-basic variables needs to be provided
 */
double Hplane::distanceToHplane(const double varVal,
    const CoinPackedVector& rayOut, const PHASolverInterface* const solver,
    const bool inNBSpace, const std::vector<int>* nonBasicVarIndex) const {
  // Ensure bound is finite
  if (isInfinity(std::abs(bound), solver->getInfinity())) {
    return solver->getInfinity();
  }

  const double rayVal = dotProductWithHplane(rayOut, solver, true, inNBSpace,
      nonBasicVarIndex);

  return distanceToHplane(varVal, rayVal, solver);
} /* distanceToHplane 5.1 */

double Hplane::distanceToHplane(const double varVal,
    const CoinPackedVector& rayOut, const PHASolverInterface* const solver,
    /*const bool inNBSpace,*/ const SolutionInfo& solnInfo) const {
  // Ensure bound is finite
  if (isInfinity(std::abs(bound), solver->getInfinity())) {
    return solver->getInfinity();
  }

  const double rayVal = dotProductWithHplane(rayOut, solnInfo, true);

  return distanceToHplane(varVal, rayVal, solver);
} /* distanceToHplane 5.2 */

double Hplane::dotProductWithHplane(const CoinPackedVector& vec,
    const PHASolverInterface* solver, const bool vecIsRay,
    const bool inNBSpace, const std::vector<int>* nonBasicVarIndex) const {
  const double VECEPS = vecIsRay ? param.getRAYEPS() : param.getEPS();
  const int numElems = vec.getNumElements();
  const int* vecIndex = vec.getIndices();
  const double* vecVal = vec.getElements();

//  const int solverFactorizationStatus = solver->canDoSimplexInterface();
//  if (solverFactorizationStatus == 0) {
  solver->enableFactorization();
//  }
  std::vector<double> rowOfBInvAStruct(solver->getNumCols()), rowOfBInvASlack(
      solver->getNumRows());
  if (row >= 0) {
    solver->getBInvARow(row, &rowOfBInvAStruct[0], &rowOfBInvASlack[0]);
  }
//  if (solverFactorizationStatus == 0) {
  solver->disableFactorization();
//  }
  if (inNBSpace) {
    if (nonBasicVarIndex == NULL) {
      error_msg(errstr, "In non-basic space and set of non-basic variables to do conversion is null.\n");
      writeErrorToII(errstr, GlobalVariables::log_file);
      exit(1);
    }
  }

  double dot_prod = 0.0;
  for (int el = 0; el < numElems; el++) {
    const int vec_ind = vecIndex[el]; // If in the NB space, this is a non-basic index
    if (isZero(vecVal[el], VECEPS)) {
      continue;
    }

    const int struct_ind = (inNBSpace) ? (*nonBasicVarIndex)[vec_ind] : vec_ind;
    if (!inNBSpace && !isRayVar(solver, struct_ind)) {
      continue;
    }

    const int rayMult = (isNonBasicUBVar(solver, struct_ind)) ? -1 : 1;

    if (row >= 0) {
      const double rowOfBInvAVal =
          (struct_ind < solver->getNumCols()) ?
              rowOfBInvAStruct[struct_ind] :
              rowOfBInvASlack[struct_ind - solver->getNumCols()];
      const double rayVal = rayMult * -1 * rowOfBInvAVal;
      if (!isZero(rayVal, param.getRAYEPS())) {
        dot_prod += vecVal[el] * rayVal;
      }
    } else if (struct_ind == var) {
      dot_prod += rayMult *  vecVal[el];
      break; // var appears at most once in vec
    }
  }
  return dot_prod;
} /* dotProductWithHplane */

double Hplane::dotProductWithHplane(const CoinPackedVector& vec,
    const SolutionInfo& solnInfo, const bool vecIsRay) const {
  const double VECEPS = vecIsRay ? CoinMin(param.getRAYEPS(), solnInfo.EPS) : CoinMin(param.getEPS(), solnInfo.EPS);
  const int numElems = vec.getNumElements();
  const int* vecIndex = vec.getIndices();
  const double* vecVal = vec.getElements();

  double dot_prod = 0.0;
  for (int el = 0; el < numElems; el++) {
    const int vec_ind = vecIndex[el]; // If in the NB space, this is a non-basic index
    if (isZero(vecVal[el], VECEPS)) {
      continue;
    }

    if (row >= 0) {
      const double rayVal = solnInfo.raysOfC1[vec_ind][row];
      if (!isZero(rayVal, CoinMin(param.getRAYEPS(), solnInfo.EPS)))
        dot_prod += vecVal[el] * rayVal;
    } else if (solnInfo.nonBasicVarIndex[vec_ind] == var) {
      dot_prod += vecVal[el]; // The rayVal component is 1.0 in comp nb space
      break; // var appears at most once in vec
    }
  }
  return dot_prod;
} /* dotProductWithHplane */

double Hplane::getRowActivity(const CoinPackedVector& vertex,
    const PHASolverInterface* solver, const bool inNBSpace,
    const std::vector<int>* nonBasicVarIndex) const {
  const int numVertexElems = vertex.getNumElements();
  const int* vertexIndex = vertex.getIndices();
  const double* vertexVal = vertex.getElements();

  const CoinShallowPackedVector currRow = solver->getMatrixByRow()->getVector(row);
  const int numRowElems = currRow.getNumElements();
  const int* rowIndex = currRow.getIndices();
  const double* rowCoeff = currRow.getElements();

  int row_el = 0;
  double dot_prod = 0.0;
  for (int vert_el = 0; vert_el < numVertexElems; vert_el++) {
    const int vert_ind =
        (inNBSpace) ?
            (*nonBasicVarIndex)[vertexIndex[vert_el]] : vertexIndex[vert_el];
    if (vert_ind > solver->getNumCols()) {
      return dot_prod;
    }
    const double vert_val = vertexVal[vert_el];

    double row_val = 0.0;
    if (row_el < numRowElems) {
      int row_ind = rowIndex[row_el];

      while (row_ind < vert_ind) {
        row_el++;
        if (row_el < numRowElems) {
          row_ind = rowIndex[row_el];
        } else {
          break;
        }
      }

      if (row_ind == vert_ind) {
        row_val = rowCoeff[row_ind];
      }
    }

    if (!isZero(vert_val) && !isZero(row_val)) {
      dot_prod += row_val * vert_val;
    }
  }
  return dot_prod;
} /* getRowActivity */
