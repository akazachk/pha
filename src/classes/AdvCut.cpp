//============================================================================
// Name        : AdvCut.cpp
// Author      : Aleksandr M. Kazachkov
// Version     : 0.2014.03.10
// Copyright   : Your copyright notice
// Description : Advanced cut class that extends OsiRowCut to allow storage of
//         other helpful data, such as the cut in NB space.
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#include <cmath> // Remember to put std::abs
#include "AdvCut.hpp"
#include "Utility.hpp"

AdvCut::AdvCut(bool NB, const int num_splits_involved) :
    OsiRowCut() {
  this->postCutObj = -1.0;
  this->NBSpace = NB;
//  this->comp_rhs = -3.1415;
  this->num_coeff = -1;
  this->cgsIndex = -1;
  this->splitVarIndex.resize(num_splits_involved,-1);
  this->cutHeur = CutHeuristics::NUM_CUT_HEUR;
}

/**
 * Copy constructor
 */
AdvCut::AdvCut(const AdvCut& rhs) :
  OsiRowCut(rhs),
  postCutObj(rhs.postCutObj),
  NBSpace(rhs.NBSpace),
  num_coeff(rhs.num_coeff),
  cgsIndex(rhs.cgsIndex),
  cgsName(rhs.cgsName),
  splitVarIndex(rhs.splitVarIndex),
  cutHeur(rhs.cutHeur) 
{
  // Nothing to do here
}

/**
 * Sets the cut based on the input temporary cut
 */
void AdvCut::assign(const AdvCut & tmpCut) {
  this->NBSpace = tmpCut.NBSpace;
  this->num_coeff = tmpCut.num_coeff;
//  this->comp_rhs = tmpCut.comp_rhs;
  setLb(tmpCut.rhs());
  this->cgsIndex = tmpCut.cgsIndex;
  this->cgsName = tmpCut.cgsName;
  this->splitVarIndex = tmpCut.splitVarIndex;
  setRow(tmpCut.row());
  if (!tmpCut.distAlongRay.empty())
    this->distAlongRay = tmpCut.distAlongRay;
}

/**
 * Assignment operator
 */
AdvCut & AdvCut::operator=(const AdvCut & rhs) {
  if (this != &rhs) {
    OsiRowCut::operator=(rhs);
    this->NBSpace = rhs.NBSpace;
    this->num_coeff = rhs.num_coeff;
    setLb(rhs.rhs());
    this->cgsIndex = rhs.cgsIndex;
    this->cgsName = rhs.cgsName;
    this->splitVarIndex = rhs.splitVarIndex;
    setRow(rhs.row());
    if (!rhs.distAlongRay.empty())
      this->distAlongRay = rhs.distAlongRay;
  }
  return *this;
}

double AdvCut::getTwoNormSquared() const {
  const double twoNorm = this->row().twoNorm();
  return twoNorm * twoNorm;;
}

double AdvCut::getTwoNorm() const {
  return this->row().twoNorm();
}

/**
 * The cut is stored as alpha x >= beta.
 */
void AdvCut::setOsiRowCut(const std::vector<double>& coeff, const double rhs,
    const double EPS) {
  setOsiRowCut((int) coeff.size(), coeff.data(), rhs, EPS);
}

void AdvCut::setOsiRowCut(const OsiRowCut& tmpCut) {
  OsiRowCut::operator=(tmpCut);
  // num_coeff is not set
}

/**
 * The cut is stored as alpha x >= beta.
 */
void AdvCut::setOsiRowCut(const int num_coeff, const double* coeff,
    const double rhs, const double EPS) {
//  if (NBSpace) // should only save in basic space?
//    return;
  this->num_coeff = num_coeff;
  std::vector<int> indices;
  std::vector<double> sparse_coeff;
  indices.reserve(num_coeff);
  sparse_coeff.reserve(num_coeff);

  for (int i = 0; i < num_coeff; i++) {
    if (!isZero(coeff[i], EPS)) {
      indices.push_back(i);
      sparse_coeff.push_back(coeff[i]);
    }
  }
  OsiRowCut::setLb(rhs);
  OsiRowCut::setRow((int) indices.size(), indices.data(),
      sparse_coeff.data());
}

/**
 * Depends crucially on vec being sorted identically to row() for the cut
 */
double AdvCut::getActivity(const CoinPackedVector& vec) const {
  const int numVecElem = vec.getNumElements();
  if (numVecElem == 0)
    return 0.0;

  const int* vecIndex = vec.getIndices();
  const double* vecElem = vec.getElements();

  const int numCutElem = row().getNumElements();
  const int* cutIndex = row().getIndices();
  const double* cutElem = row().getElements();

  double activity = 0.0;

  int vec_el = 0;
  for (int cut_el = 0; cut_el < numCutElem; cut_el++) {
    const int cut_ind = cutIndex[cut_el];

    int vec_ind = vecIndex[vec_el];
    while (vec_ind < cut_ind) {
      vec_el++;
      // If we have reached end of vec, rest is zero
      if (vec_el >= numVecElem) {
        return activity;
      }
      vec_ind = vecIndex[vec_el];
    }

    if (vec_ind == cut_ind) {
      if (!isZero(cutElem[cut_el]) && !isZero(vecElem[vec_el])) {
        activity += cutElem[cut_el] * vecElem[vec_el];
      }
    }
  }

  return activity;
}

double AdvCut::getOrthogonality(const CoinPackedVector& vec) const {
  const double scalarprod = sparseDotProduct(vec, this->row());
  const double parallelism = std::abs(scalarprod)
      / (vec.twoNorm() * this->row().twoNorm());
  return 1. - parallelism;
}

/***********************************************************************/
/**
 * @brief Set current cut by converting it from non-basic coeffs given in coeff
 */
void AdvCut::setCutFromJSpaceCoefficients(const double* const coeff,
    const double beta, const PHASolverInterface* const solver,
    const std::vector<int> &nonBasicVarIndex, const bool inCompSpace,
    const double EPS) {
  // Useful info
  int numcols = solver->getNumCols();

  // Allocate space
  this->NBSpace = false;
  this->num_coeff = numcols;
  double StructCutRHS = beta;
  const double* rowRHS = solver->getRightHandSide();
  const char* rowSense = solver->getRowSense();
  const CoinPackedMatrix* mat = solver->getMatrixByRow();

  std::vector<double> struct_coeff(numcols);
  for (int i = 0; i < numcols; i++) {
    const int curr_ind = i;
    const int curr_var = nonBasicVarIndex[curr_ind];
    const double curr_val = coeff[i];

    // If we are in the complemented space,
    // at this point we need to do the complementing
    // ub: y_j = u_j - x_j
    // lb: y_j = x_j - l_j
    if (curr_var < numcols) {
      struct_coeff[curr_var] = curr_val;
      if (inCompSpace) {
        if (solver->getModelPtr()->getColumnStatus(curr_var)
            == ClpSimplex::atUpperBound) {
          StructCutRHS -= curr_val * solver->getColUpper()[curr_var];
          struct_coeff[curr_var] *= -1;
        } else if (solver->getModelPtr()->getColumnStatus(curr_var)
            == ClpSimplex::atLowerBound
            || solver->getModelPtr()->getColumnStatus(curr_var)
                == ClpSimplex::isFixed) {
          StructCutRHS += curr_val * solver->getColLower()[curr_var];
        }
      }
    } else {
      const int curr_row = curr_var - numcols;
      const double mult =
          (inCompSpace && rowSense[curr_row] == 'G') ? 1.0 : -1.0;

      const CoinShallowPackedVector slackRow = mat->getVector(curr_row);
      int slackRowNumElements = slackRow.getNumElements();
      const int* slackRowIndices = slackRow.getIndices();
      const double* slackRowCoeffs = slackRow.getElements();

      StructCutRHS += mult * rowRHS[curr_row] * curr_val;
      for (int j = 0; j < slackRowNumElements; j++) {
        struct_coeff[slackRowIndices[j]] += mult * slackRowCoeffs[j] * curr_val;
      }
    }
  }

  // Normalize
  const double absRHS = std::abs(StructCutRHS);
  if (!isZero(StructCutRHS)) {
    for (int i = 0; i < numcols; i++) {
      struct_coeff[i] = struct_coeff[i] / absRHS;
    }
  }

  if (greaterThanVal(StructCutRHS, 0.0)) {
    StructCutRHS = 1.0;
  } else {
    if (lessThanVal(StructCutRHS, 0.0)) {
      StructCutRHS = -1.0;
    } else {
      StructCutRHS = 0.0;
    }
  }

  this->setOsiRowCut(struct_coeff, StructCutRHS, EPS);
} /* setCutFromJspaceCoefficients */

int AdvCut::cleanCut(const OsiSolverInterface* const solver,
    const SolutionInfo& probData /*, const double beta*/) {
  removeSmallCoefficients(solver, param.getEPS());
  if (badDynamism(probData,
      CoinMax((probData.EPS > 0) ? 1. / probData.EPS : param.getMAXDYN(),
          param.getMAXDYN()), probData.EPS)) {
    return -1 * (CutSolverFails::SCALING_CUTSOLVER_FAIL_IND + 1);
  }
  return 0;
} /* cleanCut */

/**
 * Taken from CglGMI with some minor modifications, including switching to >= cut as we use
 */
void AdvCut::removeSmallCoefficients(const OsiSolverInterface* const solver,
    const double EPS) {
  CoinPackedVector vec = this->mutableRow();
  int cutNz = vec.getNumElements();
  int* cutIndex = vec.getIndices();
  double* cutElem = vec.getElements();
  double cutRhs = this->rhs();
  const double* colLower = solver->getColLower();
  const double* colUpper = solver->getColUpper();

  double coeff, absval;
  int currPos = 0;
  int col;
  for (int i = 0; i < cutNz; ++i) {
    col = cutIndex[i];
    coeff = cutElem[i];
    absval = std::abs(coeff);
    if (isZero(absval, EPS)
        && isZero(absval * solver->getObjCoefficients()[col], EPS)) {
      continue; // throw away this element
    }
    if (absval <= GlobalVariables::param.getEPS_COEFF()) {
      // small coefficient: remove and adjust rhs if possible
      if ((coeff < 0.0) && !isNegInfinity(colLower[col])) {
        cutRhs -= coeff * colLower[col];
      } else if ((coeff > 0.0) && !isInfinity(colUpper[col])) {
        cutRhs -= coeff * colUpper[col];
      }
    } else if (absval > GlobalVariables::param.getEPS_COEFF()) {
      if (currPos < i) {
        cutElem[currPos] = cutElem[i];
        cutIndex[currPos] = cutIndex[i];
      }
      currPos++;
    }
  }
  cutNz = currPos;
  this->setRow(cutNz, cutIndex, cutElem);
  if (isZero(cutRhs, EPS)) {
    cutRhs = 0.;
  }
  this->setLb(cutRhs);
} /* removeSmallCoefficients */

bool AdvCut::badViolation(const OsiSolverInterface* const solver) {
  // Check the cut
  const double currCutNorm = this->getTwoNorm();
  const double violation = this->violated(solver->getColSolution());
  //    const bool cuttingSolutionFlagAbs = (inNBSpace)
  //        ? lessThanVal(cutSolver->getObjValue() - beta, 0.0, param.getMIN_VIOL_ABS())
  //        : greaterThanVal(currCut.violated(origSolver->getColSolution()), 0.0, param.getMIN_VIOL_ABS());
  //    const bool cuttingSolutionFlagRel = (inNBSpace)
  //        ? lessThanVal(cutSolver->getObjValue() - beta, 0.0, currCutNorm * param.getMIN_VIOL_REL())
  //        : greaterThanVal(currCut.violated(origSolver->getColSolution()), 0.0, currCutNorm * param.getMIN_VIOL_REL());
  const bool cuttingSolutionFlagAbs = //(inNBSpace) ||
      greaterThanVal(violation, 0.0, param.getMIN_VIOL_ABS());
  const bool cuttingSolutionFlagRel = //(inNBSpace) ||
      greaterThanVal(violation, 0.0, currCutNorm * param.getMIN_VIOL_REL());
  if (!cuttingSolutionFlagAbs || !cuttingSolutionFlagRel) {
    return true;
  }
  return false;
} /* badViolation */

bool AdvCut::badDynamism(const SolutionInfo& probData, const double SENSIBLE_MAX, const double EPS) {
  const CoinPackedVector vec = this->row();
  const int num_el = vec.getNumElements();
  const double* el = vec.getElements();

  double minAbsElem = 0.0, maxAbsElem = 0.0;
  for (int i = 0; i < num_el; i++) {
    const double absCurr = std::abs(el[i]);
    if (isZero(absCurr, EPS)) { // make sure min and max are not zero
      continue;
    }

    if ((minAbsElem > absCurr) || isZero(minAbsElem, EPS)) {
      minAbsElem = absCurr;
    }
    if (maxAbsElem < absCurr) {
      maxAbsElem = absCurr;
    }
  }
  const double absCurr = std::abs(this->rhs());
  if (absCurr > maxAbsElem) {
    maxAbsElem = absCurr;
  }
  if (isZero(maxAbsElem, EPS)
      // || isZero(minAbsElem / maxAbsElem)
      || greaterThanVal(maxAbsElem / minAbsElem, SENSIBLE_MAX, EPS)
      || (!isZero(probData.minAbsCoeff, EPS)
          && greaterThanVal(std::abs(maxAbsElem / probData.minAbsCoeff),
              std::abs(SENSIBLE_MAX * probData.maxAbsCoeff / probData.minAbsCoeff), EPS))) {
    return true;
  }
  return false;
} /* badDynamism */

/***********************************************************************/
/**
 * @brief Converts a cut from the non-basic space to the structural space.
 *
 * @param structCut         ::  OUTPUT: [split var index, rhs, coeffs]
 * @param nonBasicOrigVarIndex    ::  Indices of NB structural variables
 * @param nonBasicSlackVarIndex   ::  Indices of NB slack variables
 * @param solver          ::  Solver
 * @param JspaceCut         ::  [split var index, rhs, coeffs]
 */
void AdvCut::convertCutFromJSpaceToStructSpace(AdvCut& structCut,
    const PHASolverInterface* const solver,
    const std::vector<int> &nonBasicVarIndex, const AdvCut& JspaceCut,
    const bool inCompSpace, const double EPS) {
  // Useful info
  int numcols = solver->getNumCols();

  // Allocate space
  structCut.NBSpace = false;
  structCut.num_coeff = numcols;
  structCut.cgsIndex = JspaceCut.cgsIndex;
  structCut.cgsName = JspaceCut.cgsName;
  structCut.splitVarIndex = JspaceCut.splitVarIndex;

  double StructCutRHS = JspaceCut.rhs();
  const double* rowRHS = solver->getRightHandSide();
  const char* rowSense = solver->getRowSense();
  const CoinPackedMatrix* mat = solver->getMatrixByRow();

  const CoinPackedVector vec = JspaceCut.row();
  const int numElem = vec.getNumElements();
  const int* vecInd = vec.getIndices();
  const double* vecEl = vec.getElements();

  std::vector<double> struct_coeff(numcols);
  for (int i = 0; i < numElem; i++) {
    const int curr_ind = vecInd[i];
    const int curr_var = nonBasicVarIndex[curr_ind];
    const double curr_val = vecEl[i];

    // If we are in the complemented space,
    // at this point we need to do the complementing
    // ub: y_j = u_j - x_j
    // lb: y_j = x_j - l_j
    if (curr_var < numcols) {
      struct_coeff[curr_var] = curr_val;
      if (inCompSpace) {
        if (solver->getModelPtr()->getColumnStatus(curr_var)
            == ClpSimplex::atUpperBound) {
          StructCutRHS -= curr_val * solver->getColUpper()[curr_var];
          struct_coeff[curr_var] *= -1;
        } else if (solver->getModelPtr()->getColumnStatus(curr_var)
            == ClpSimplex::atLowerBound
            || solver->getModelPtr()->getColumnStatus(curr_var)
                == ClpSimplex::isFixed) {
          StructCutRHS += curr_val * solver->getColLower()[curr_var];
        }
      }
    } else {
      const int curr_row = curr_var - numcols;
      const double mult =
          (inCompSpace && rowSense[curr_row] == 'G') ? 1.0 : -1.0;

      const CoinShallowPackedVector slackRow = mat->getVector(curr_row);
      int slackRowNumElements = slackRow.getNumElements();
      const int* slackRowIndices = slackRow.getIndices();
      const double* slackRowCoeffs = slackRow.getElements();

      StructCutRHS += mult * rowRHS[curr_row] * curr_val;
      for (int j = 0; j < slackRowNumElements; j++) {
        struct_coeff[slackRowIndices[j]] += mult * slackRowCoeffs[j]
            * curr_val;
      }
    }
  }

  // Normalize
  const double absRHS = std::abs(StructCutRHS);
  if (!isZero(StructCutRHS)) {
    for (int i = 0; i < numcols; i++) {
      struct_coeff[i] = struct_coeff[i] / absRHS;
    }
  }

  if (greaterThanVal(StructCutRHS, 0.0)) {
    StructCutRHS = 1.0;
  } else {
    if (lessThanVal(StructCutRHS, 0.0)) {
      StructCutRHS = -1.0;
    } else {
      StructCutRHS = 0.0;
    }
  }

  structCut.setOsiRowCut(struct_coeff, StructCutRHS, EPS);
} /* convertCutFromJSpaceToStructSpace */

/***********************************************************************/
/**
 * @brief Converts a cut from the non-basic space to the structural space.
 *
 * @param structCut         ::  OUTPUT: [split var index, rhs, coeffs]
 * @param nonBasicOrigVarIndex    ::  Indices of NB structural variables
 * @param nonBasicSlackVarIndex   ::  Indices of NB slack variables
 * @param solver          ::  Solver
 * @param JspaceCut         ::  [split var index, rhs, coeffs]
 */
void AdvCut::convertCutCoeffFromJSpaceToStructSpaceNonPacked(
    std::vector<double>& structCut, const PHASolverInterface* const solver,
    const std::vector<int>& nonBasicVarIndex, const double* JspaceCut,
    const double JspaceRHS, const bool inCompSpace) {
  // Useful info
  const int numcols = solver->getNumCols();
  const int numNB = nonBasicVarIndex.size();
  structCut.clear();
  structCut.resize(numcols);
  double StructCutRHS = JspaceRHS;

  const double* rowRHS = solver->getRightHandSide();
  const char* rowSense = solver->getRowSense();
  const CoinPackedMatrix* mat = solver->getMatrixByRow();

//  std::vector<double> struct_coeff(numcols);
  for (int curr_ind = 0; curr_ind < numNB; curr_ind++) {
    const int curr_var = nonBasicVarIndex[curr_ind];
    const double curr_val = JspaceCut[curr_ind];

    // If we are in the complemented space,
    // at this point we need to do the complementing
    // ub: y_j = u_j - x_j
    // lb: y_j = x_j - l_j
    if (curr_var < numcols) {
      structCut[curr_var] = curr_val;
      if (inCompSpace) {
        if (solver->getModelPtr()->getColumnStatus(curr_var)
            == ClpSimplex::atUpperBound) {
          StructCutRHS -= curr_val * solver->getColUpper()[curr_var];
          structCut[curr_var] *= -1;
        } else if (solver->getModelPtr()->getColumnStatus(curr_var)
            == ClpSimplex::atLowerBound
            || solver->getModelPtr()->getColumnStatus(curr_var)
                == ClpSimplex::isFixed) {
          StructCutRHS += curr_val * solver->getColLower()[curr_var];
        }
      }
    } else {
      const int curr_row = curr_var - numcols;
      const double mult =
          (inCompSpace && rowSense[curr_row] == 'G') ? 1.0 : -1.0;

      const CoinShallowPackedVector slackRow = mat->getVector(curr_row);
      int slackRowNumElements = slackRow.getNumElements();
      const int* slackRowIndices = slackRow.getIndices();
      const double* slackRowCoeffs = slackRow.getElements();

      StructCutRHS += mult * rowRHS[curr_row] * curr_val;
      for (int j = 0; j < slackRowNumElements; j++) {
        structCut[slackRowIndices[j]] += mult * slackRowCoeffs[j]
            * curr_val;
      }
    }
  }

  // Normalize
  const double absRHS = std::abs(StructCutRHS);
  if (!isZero(StructCutRHS)) {
    for (int i = 0; i < numcols; i++) {
      structCut[i] = structCut[i] / absRHS;
    }
  }

  if (greaterThanVal(StructCutRHS, 0.0)) {
    StructCutRHS = 1.0;
  } else {
    if (lessThanVal(StructCutRHS, 0.0)) {
      StructCutRHS = -1.0;
    } else {
      StructCutRHS = 0.0;
    }
  }
}

/********************************************************************************/
/**
 * @brief Decide if two cuts are the same.
 * @return 0: seem same in coeff and rhs, +/-1: seem same in coeff, but diff in rhs (-1: cut1 better, +1: cut2 better), 2: seem different in coeff and rhs.
 *
 * The value of eps will be used to determine whether a cut coefficient is zero or not.
 * We are roughly checking that whether cut1 = ratio * cut2 for some ratio.
 * Let ratio = rhs1 / rhs2. (If both are non-zero. Typically both are = 1.)
 * We say the cuts are "truly different" if |coeff^1_j - ratio * coeff^2_j| >= diffeps.
 * However, might be that the cut coefficients are scaled versions, but the rhs values are different.
 * So we also compute ratio_coeff = coeff^1_j / coeff^2_j for the first j where both are non-zero.
 * If |coeff^1_j - ratio * coeff^2_j| >= diffeps for some j,
 * but |coeff^1_j - ratio_coeff * coeff^2_j| < diffeps for all j,
 * then the return code will be 1 (essentially, one cut will dominate the other).
 *
 * By the way, this is all essentially checking orthogonality...
 */
int AdvCut::cutsDifferent(const OsiRowCut &cut1, const OsiRowCut &cut2,
    const double eps) {
  return isRowDifferent(cut1.row(), cut1.rhs(), cut2.row(), cut2.rhs(), eps);
} /* cutsDifferent */

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/* AdvCuts */
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

AdvCuts::AdvCuts(bool NB) :
    OsiCuts() {
  this->NBSpace = NB;
}

/**
 * Copy constructor
 */
AdvCuts::AdvCuts(const AdvCuts& rhs) :
  OsiCuts(rhs),
  cuts(rhs.cuts),
  score(rhs.score),
  orthogonality(rhs.orthogonality),
  NBSpace(rhs.NBSpace)
{
  // Nothing to do here
}

void AdvCuts::setOsiCuts() {
//  if (NBSpace) // should only save in basic space?
//    return;

  for (int cutnum = OsiCuts::sizeCuts(); cutnum < (int) cuts.size(); cutnum++) {
    OsiCuts::insert(cuts[cutnum]);
  }
}

/**
 * Adds cut if it is not a duplicate
 * @return True if cut added, false otherwise
 */
bool AdvCuts::addNonDuplicateCut(const AdvCut & tmpCut,
    bool checkForDuplicate) {
  if (!checkForDuplicate || !isDuplicateCut(tmpCut)) {
    cuts.push_back(tmpCut);
    return true;
  } else {
    return false;
  }
}

/**
 * Checks if cut is not a duplicate
 * @return True if cut is duplicate, false otherwise
 */
bool AdvCuts::isDuplicateCut(const AdvCut & tmpCut) const {
  return (duplicateCutIndex(tmpCut) >= 0);
}

/**
 * Checks if cut is not a duplicate
 * @return Index if cut is duplicate, -1 otherwise
 */
int AdvCuts::duplicateCutIndex(const AdvCut & tmpCut) const {
  for (int i = 0; i < (int) this->cuts.size(); i++) {
    if (tmpCut == this->cuts[i])
      return i;
  }
  // If we have reached the end, then this cut was not found
  return -1;
}

/**
 * Checks not only if cut is duplicate, but also orthogonality
 */
int AdvCuts::howDuplicate(const AdvCut& tmpCut, int& duplicateCutIndex,
    int& minOrthoIndex, double& minOrtho, const double eps) const {
  duplicateCutIndex = -1;
  minOrthoIndex = -1;
  minOrtho = 1.;
  for (int i = 0; i < (int) this->cuts.size(); i++) {
//    if (GlobalVariables::param.getMINORTHOGONALITY() > 0) {
      const double this_ortho = tmpCut.getOrthogonality(this->cuts[i].row());
      if (this_ortho < minOrtho) {
        minOrtho = this_ortho;
        minOrthoIndex = i;
      }
//    }

      const int howDifferent = AdvCut::cutsDifferent(this->cuts[i], tmpCut, eps);
    if (howDifferent != 2) {
      duplicateCutIndex = i;
      return howDifferent;
    }
  }
  return 0;
}

/***********************************************************************/
/**
 * @brief Useful to set cuts if cuts were inserted using OsiCuts interface
 */
void AdvCuts::setCuts() {
  cuts.resize(sizeCuts());
  for (int cut_ind = 0; cut_ind < sizeCuts(); cut_ind++) {
    cuts[cut_ind].setOsiRowCut(rowCut(cut_ind));
  }
}

/***********************************************************************/
/**
 * @brief Assignment operator
 */
AdvCuts& AdvCuts::operator=(const AdvCuts& rhs) {
  if (this != &rhs) {
    OsiCuts::operator=(rhs);
    cuts = rhs.cuts;
    score = rhs.score;
    orthogonality = rhs.orthogonality;
    NBSpace = rhs.NBSpace;
  }
  return *this;
}

/**
 * @brief This assignment operator does not set NBSpace, as it does not have it available
 */
/*
AdvCuts& AdvCuts::operator=(const OsiCuts& rhs) {
  try {
    return AdvCuts::operator=(dynamic_cast<AdvCuts&>(rhs));
  } catch (std::exception& e) {
    if (this != &rhs) {
      OsiCuts::operator=(rhs);
      cuts.resize(rhs.sizeCuts());
      for (int cut_ind = 0; cut_ind < rhs.sizeCuts(); cut_ind++) {
        cuts[cut_ind] = rhs.rowCut(cut_ind);
      }
    }
    return *this;
  }
}
*/


void AdvCuts::printCuts(std::string cut_type,
    const OsiSolverInterface* const solver,
    const SolutionInfo& probData) const {
  char filename[300];
  char headerText[300];

  // Print cuts in NB space
  if (NBSpace) {
    snprintf(filename, sizeof(filename) / sizeof(char), "%s-NB_%s.csv",
        GlobalVariables::out_f_name_stub.c_str(), cut_type.c_str());
    snprintf(headerText, sizeof(headerText) / sizeof(char),
        "\n## Printing %s in non-basic space. ##\n", cut_type.c_str());
#ifdef TRACE
    printf("%s", headerText);
#endif
    writeCutsInNBSpace(headerText, cuts, probData, solver, filename);
  } else {
    // Print cuts in structural space
    snprintf(filename, sizeof(filename) / sizeof(char), "%s-%s.csv",
        GlobalVariables::out_f_name_stub.c_str(), cut_type.c_str());
    snprintf(headerText, sizeof(headerText) / sizeof(char),
        "\n## Printing %s in structural space. ##\n", cut_type.c_str());
#ifdef TRACE
    printf("%s", headerText);
#endif
    writeCutsInStructSpace(headerText, cuts, probData, solver, filename);
  }
}

void AdvCuts::printCuts(std::string cut_type,
    const OsiSolverInterface* const solver, const bool packedForm) const {
  char filename[300];
  char headerText[300];

//  // Print cuts in NB space
//  if (NBSpace) {
//    snprintf(filename, sizeof(filename) / sizeof(char), "%s-NB_%s.csv",
//        GlobalVariables::out_f_name_stub.c_str(), cut_type.c_str());
//    snprintf(headerText, sizeof(headerText) / sizeof(char),
//        "\n## Printing %s in non-basic space. ##\n", cut_type.c_str());
//#ifdef TRACE
//    printf("%s", headerText);
//#endif
//    writeCutsInNBSpace(headerText, cuts, probData, solver, filename);
//  } else {
  // Print cuts in structural space
  snprintf(filename, sizeof(filename) / sizeof(char), "%s-%s.csv",
      GlobalVariables::out_f_name_stub.c_str(), cut_type.c_str());
  snprintf(headerText, sizeof(headerText) / sizeof(char),
      "\n## Printing %s in structural space. ##\n", cut_type.c_str());
#ifdef TRACE
  printf("%s", headerText);
#endif
  writeCutsInStructSpace(headerText, *this, solver, filename, packedForm);
//  }
}

/***********************************************************************/
/**
 * @brief Writes sets of cuts that are stored in NB space.
 */
void AdvCuts::writeCutsInNBSpace(std::string headerText,
    const std::vector<AdvCut> &NBCut, const SolutionInfo &probData,
    const OsiSolverInterface* const solver, char* filename) {
  FILE* NBCutOut = fopen(filename, "w");
  if (NBCutOut == NULL) {
    error_msg(errorstring,
        "Unable to open file to write cuts in non-basic space for filename %s.\n",
        filename);
    writeErrorToII(errorstring, GlobalVariables::log_file);
    exit(1);
  }
  fprintf(NBCutOut, "%s", headerText.c_str());
  fprintf(NBCutOut, "%s,%s,%s,%s,%s,%s,%s\n", "Cut #", "Cut heur", "Facet",
      "Facet name", "Vars", "Rhs", "Coeff");
  fprintf(NBCutOut, ",,,,,,");
  for (int nbcol = 0; nbcol < (int) probData.nonBasicOrigVarIndex.size();
      nbcol++) {
    int currvar = probData.nonBasicOrigVarIndex[nbcol];
    fprintf(NBCutOut, "%d,", currvar);
  }
  for (int nbcol = 0; nbcol < (int) probData.nonBasicSlackVarIndex.size();
      nbcol++) {
    int currvar = probData.nonBasicSlackVarIndex[nbcol];
    fprintf(NBCutOut, "%d,", currvar);
  }
  fprintf(NBCutOut, "\n");
  fprintf(NBCutOut, ",,,,,,");
  for (int nbcol = 0; nbcol < (int) probData.nonBasicOrigVarIndex.size();
      nbcol++) {
    int currvar = probData.nonBasicOrigVarIndex[nbcol];
    fprintf(NBCutOut, "%s,", solver->getColName(currvar).c_str());
  }
  for (int nbcol = 0; nbcol < (int) probData.nonBasicSlackVarIndex.size();
      nbcol++) {
    int currvar = probData.nonBasicSlackVarIndex[nbcol];
    fprintf(NBCutOut, "%s,",
        solver->getRowName(currvar - probData.numCols).c_str());
  }
  fprintf(NBCutOut, "\n");
  for (int currcut = 0; currcut < (int) NBCut.size(); currcut++) {
    fprintf(NBCutOut, "%d,", currcut);
    fprintf(NBCutOut, "%d,", NBCut[currcut].cutHeur);
    fprintf(NBCutOut, "%d,", NBCut[currcut].cgsIndex);
    fprintf(NBCutOut, "%s,", NBCut[currcut].cgsName.c_str());
    fprintf(NBCutOut, "%s,", NBCut[currcut].splitVarIndexString().c_str());
    fprintf(NBCutOut, "%f,", NBCut[currcut].rhs());
    for (int c = 0; c < (int) NBCut[currcut].num_coeff; c++) {
      fprintf(NBCutOut, "%f,", NBCut[currcut][c]);
    }
    fprintf(NBCutOut, "\n");
  }
  fprintf(NBCutOut, "\n");

  fclose(NBCutOut);
}

/***********************************************************************/
/**
 * @brief Writes sets of cuts that are stored in structural space.
 */
void AdvCuts::writeCutsInStructSpace(std::string headerText,
    const std::vector<AdvCut> &structCut, const SolutionInfo &probData,
    const OsiSolverInterface* const solver, char* filename) {
  FILE* StructCutOut = fopen(filename, "w");
  if (StructCutOut == NULL) {
    error_msg(errorstring,
        "Unable to open file to write cuts in structural space for filename %s.\n",
        filename);
    writeErrorToII(errorstring, GlobalVariables::log_file);
    exit(1);
  }
  fprintf(StructCutOut, "%s", headerText.c_str());
  fprintf(StructCutOut, "%s,%s,%s,%s,%s,%s,%s\n", "Cut #", "Cut heur", "Facet",
      "Facet name", "Vars", "Rhs", "Coeff");
  fprintf(StructCutOut, ",,,,,,");
  for (int col = 0; col < probData.numCols; col++) {
    fprintf(StructCutOut, "%d,", col);
  }
  fprintf(StructCutOut, "\n");
  fprintf(StructCutOut, ",,,,,,");
  for (int col = 0; col < probData.numCols; col++) {
    fprintf(StructCutOut, "%s,", solver->getColName(col).c_str());
  }
  fprintf(StructCutOut, "\n");
  for (int currcut = 0; currcut < (int) structCut.size(); currcut++) {
    fprintf(StructCutOut, "%d,", currcut);
    fprintf(StructCutOut, "%d,", structCut[currcut].cutHeur);
    fprintf(StructCutOut, "%d,", structCut[currcut].cgsIndex);
    fprintf(StructCutOut, "%s,", structCut[currcut].cgsName.c_str());
    fprintf(StructCutOut, "%s,", structCut[currcut].splitVarIndexString().c_str());
    fprintf(StructCutOut, "%f,", structCut[currcut].rhs());
    for (int c = 0; c < (int) structCut[currcut].num_coeff; c++) {
      fprintf(StructCutOut, "%f,", structCut[currcut][c]);
//      if (isnan(structCut[currcut].coeff[c])) {
//        error_msg(errstr, "Found NAN! Coeff %d of cut %d is %f.\n", c,
//            currcut, structCut[currcut].coeff[c]);
//        writeErrorToII(errstr, LiftGICsParam::inst_info_out);
//        exit(1);
//      }
    }
    fprintf(StructCutOut, "\n");
  }
  fprintf(StructCutOut, "\n");

  fclose(StructCutOut);
}

/***********************************************************************/
/**
 * @brief Writes sets of cuts that are stored in structural space.
 */
void AdvCuts::writeCutsInStructSpace(std::string headerText,
    const OsiCuts &structCut, const OsiSolverInterface* const solver,
    char* filename, const bool packedForm) {
  FILE* StructCutOut = fopen(filename, "w");
  if (StructCutOut == NULL) {
    error_msg(errorstring,
        "Unable to open file to write cuts in structural space for filename %s.\n",
        filename);
    writeErrorToII(errorstring, GlobalVariables::log_file);
    exit(1);
  }

  fprintf(StructCutOut, "%s", headerText.c_str());
  fprintf(StructCutOut, "%s,%s,%s\n", "Cut #", "Rhs", "Coeff");
  if (!packedForm) {
    fprintf(StructCutOut, ",,");
    for (int col = 0; col < solver->getNumCols(); col++) {
      fprintf(StructCutOut, "%d,", col);
    }
    fprintf(StructCutOut, "\n");
    fprintf(StructCutOut, ",,");
    for (int col = 0; col < solver->getNumCols(); col++) {
      fprintf(StructCutOut, "%s,", solver->getColName(col).c_str());
    }
    fprintf(StructCutOut, "\n");
  }
  for (int curr_cut_ind = 0; curr_cut_ind < structCut.sizeCuts();
      curr_cut_ind++) {
    const OsiRowCut* currRowCut = structCut.rowCutPtr(curr_cut_ind);
    const double mult = currRowCut->sense() == 'G' ? 1.0 : -1.0;
    fprintf(StructCutOut, "%d,", curr_cut_ind);
    fprintf(StructCutOut, "%f,", mult * currRowCut->rhs());

    const int numCutElements = currRowCut->row().getNumElements();
    const int* cutIndex = currRowCut->row().getIndices();
    const double* cutElement = currRowCut->row().getElements();
    if (packedForm) {
      for (int c = 0; c < numCutElements; c++) {
        const int curr_var = cutIndex[c];
        fprintf(StructCutOut, "%f,%s,", mult * cutElement[c],
            solver->getColName(curr_var).c_str());
      }
    } else {
      int curr_ind = 0;
      for (int c = 0; c < numCutElements; c++) {
        const int curr_var = cutIndex[c];
        while (curr_var > curr_ind) {
            fprintf(StructCutOut, "0,");
            curr_ind++;
        }
        fprintf(StructCutOut, "%f,", mult * cutElement[c]);
        curr_ind++;
      }
      while (curr_ind < solver->getNumCols()) {
        fprintf(StructCutOut, "0,");
        curr_ind++;
      }
    }
    fprintf(StructCutOut, "\n");
  }
  fprintf(StructCutOut, "\n");

  fclose(StructCutOut);
}
