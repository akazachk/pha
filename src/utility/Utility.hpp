//============================================================================
// Name        : Utility.hpp
// Author      : akazachk
// Version     : 2018.11
// Description : Header for Utility.hpp, which has miscellaneous useful functions.
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#pragma once
#include <string>
#include <vector>

#include "Point.hpp"
#include "Vertex.hpp"
#include "Ray.hpp"
#include "Hplane.hpp"
#include "IntersectionInfo.hpp"

#include "typedefs.hpp"
#include "SolutionInfo.hpp"
#include "AdvCut.hpp"

#include <cstdlib>  // strtol
#include <limits>

/**
 * Calculate strong branching lower bound
 */
double calculateStrongBranchingLB(double& worstOptVal, double& bestOptVal,
    OsiSolverInterface* solver, const SolutionInfo& probData);

/**
 * Checks whether a solver is optimal
 * Something that can go wrong (e.g., bc1 -64 sb5 tmp_ind = 14):
 * Solver is declared optimal for the scaled problem but there are primal or dual infeasibilities in the unscaled problem
 * In this case, secondary status will be 2 or 3
 * Sometimes cleanup helps
 * Sometimes things break after enabling factorization (e.g., secondary status = 0, but sum infeasibilities > 0)
 */
bool checkSolverOptimality(OsiSolverInterface* const solver,
    const bool exitOnDualInfeas,
    const double timeLimit = std::numeric_limits<double>::max(),
    const int maxNumResolves = 2);

/**
 * Enable factorization, and check whether some cleanup needs to be done
 */
bool enableFactorization(OsiSolverInterface* const solver, const int resolveFlag = 2);

/**
 * @brief Parses int from string using strtol
 */
bool parseInt(const char *str, int &val);

/**
 * @brief Parses double from string using strtod
 */
bool parseDouble(const char *str, double &val);

/**
 * @brief Universal way to check whether we reached the limit for the number of cuts for each split
 * This allows us to change between restricting number of cuts per split and total number of cuts easily
 */
bool reachedCutLimit(const int num_cuts_this_cgs, const int num_cuts_total);

/**
 * @brief Universal way to check whether we reached the limit for the time for some subroutine
 */
bool reachedTimeLimit(const Stats& timer, const std::string& timeName,
    const double max_time);

/**
 * Compare two double arrays for equality
 * No effort is made to check that both have the same number of elements
 */
bool vectorsEqual(const int n, const double* vec1, const double* vec2);
int isRowDifferent(const CoinPackedVector& cut1Vec, const double cut1rhs,
    const CoinPackedVector& cut2Vec, const double cut2rhs,
    const double eps = param.getDIFFEPS());

double getOrthogonality(const CoinPackedVector& vec1,
    const CoinPackedVector& vec2);

void getOrthogonalityStatsForRound(const OsiCuts& cs,
    double& minOrtho, double& maxOrtho, double& avgOrtho,
    const std::vector<int>* const roundWhenCutApplied = NULL,
    const int this_round = -1);

double getMinOrthogonality(const CoinPackedVector& vec, const OsiCuts& cs);

/**
 * @brief Set full vector coordinates from packed vector
 *
 * No effort is made to reset the indices that are not present in packed vector
 * so these should be initialized as desired prior to calling this function
 *
 * No check is performed for indices existing,
 * so this should be ensured before calling this function
 */
void setFullVecFromPackedVec(std::vector<double>& fullVec,
    const CoinPackedVector& vec);

/**
 * @brief Any indices that appear in vec will be put to zero in fullVec
 *
 * No effort is made to reset the indices that are not present in packed vector
 * so these should be initialized as desired prior to calling this function
 *
 * No check is performed for indices existing,
 * so this should be ensured before calling this function
 */
void clearFullVecFromPackedVec(std::vector<double>& fullVec,
    const CoinPackedVector& vec);

/** From CglLandP: return the coefficients of the intersection cut */
inline double unstrengthenedIntersectionCutCoeff(double abar, double f0) {
  if (abar > 0) {
    //return alpha_i * (1 - beta);
    return (abar / f0);
  } else {
    //return -alpha_i * beta;
    return (-abar / (1 - f0));
  }
}

/** From CglLandP compute the modularized row coefficient for an integer variable */
inline double modularizedCoeff(double abar, double f0) {
  double f_i = abar - floor(abar);
  if (f_i <= f0) {
    return f_i;
  } else {
    return f_i - 1;
  }
}

/** Adapted from CglLandP: return the coefficients of the strengthened intersection cut */
inline double strengthenedIntersectionCutCoeff(double abar, double f0,
    const OsiSolverInterface* const solver = NULL, const int currFracVar = -1) {
  if ((currFracVar >= 0)
      && ((currFracVar >= solver->getNumCols())
          || !(solver->isInteger(currFracVar)))) {
    return unstrengthenedIntersectionCutCoeff(abar, f0);
  } else {
//    return modularizedCoeff(abar, f0);
    double f_i = abar - std::floor(abar);
    if (f_i < f0) {
//      return f_i * (1 - f0);
      return f_i / f0;
    } else {
//      return (1 - f_i) * f0;
      return (1 - f_i) / (1 - f0);
    }
  }
}

inline double intersectionCutCoeff(double abar, double f0,
    const OsiSolverInterface* const solver = NULL, const int currFracVar = -1,
    const bool strengthen = false) {
  if (!strengthen || (currFracVar >= solver->getNumCols())
      || !(solver->isInteger(currFracVar))) {
    return unstrengthenedIntersectionCutCoeff(abar, f0);
  } else {
    return strengthenedIntersectionCutCoeff(abar, f0);
  }
}

/**
 * @brief returns lb of col (recall complementing of slacks)
 */
double getVarLB(const OsiSolverInterface* const solver, const int col);
/**
 * @brief returns ub of col (recall complementing of slacks)
 */
double getVarUB(const OsiSolverInterface* const solver, const int col);
/**
 * @brief returns value of col (recall complementing of slacks)
 */
double getVarVal(const PHASolverInterface* const solver, const int col);
double getCNBVarVal(const PHASolverInterface* const solver, const int col);
double getCNBVarVal(const PHASolverInterface* const solver, const int col, const bool complement);

/**
 * For working with non-basic variables
 */
inline bool isBasicCol(const std::vector<int> &cstat, const int col) {
  return (cstat[col] == ClpSimplex::Status::basic);
}

inline bool isNonBasicUBCol(const std::vector<int> &cstat, const int col) {
  return (cstat[col] == ClpSimplex::Status::atUpperBound);
}

inline bool isNonBasicLBCol(const std::vector<int> &cstat, const int col) {
  return (cstat[col] == ClpSimplex::Status::atLowerBound);
}

inline bool isBasicSlack(const std::vector<int> &rstat, const int row) {
  return (rstat[row] == ClpSimplex::Status::basic);
}

inline bool isNonBasicUBSlack(const std::vector<int> &rstat, const int row) {
  return (rstat[row] == ClpSimplex::Status::atUpperBound);
}

inline bool isNonBasicLBSlack(const std::vector<int> &rstat, const int row) {
  return (rstat[row] == ClpSimplex::Status::atLowerBound);
}

/**
 * Status of structural variables
 */
inline bool isNonBasicFreeCol(const OsiClpSolverInterface* const solver,
    const int col) {
  return (solver->getModelPtr()->getColumnStatus(col) ==  ClpSimplex::Status::isFree);
}

inline bool isBasicCol(const OsiClpSolverInterface* const solver, const int col) {
  return (solver->getModelPtr()->getColumnStatus(col) ==  ClpSimplex::Status::basic);
}

inline bool isNonBasicUBCol(const OsiClpSolverInterface* const solver,
    const int col) {
  const ClpSimplex::Status status = solver->getModelPtr()->getColumnStatus(col);
  return ((status == ClpSimplex::Status::atUpperBound)
      || (status == ClpSimplex::Status::isFixed
          && lessThanVal(solver->getReducedCost()[col], 0.)));
}

inline bool isNonBasicLBCol(const OsiClpSolverInterface* const solver,
    const int col) {
  const ClpSimplex::Status status = solver->getModelPtr()->getColumnStatus(col);
  return ((status == ClpSimplex::Status::atLowerBound)
      || (status == ClpSimplex::Status::isFixed
          && !lessThanVal(solver->getReducedCost()[col], 0.)));
}

inline bool isNonBasicFixedCol(const OsiClpSolverInterface* const solver,
    const int col) {
  return (solver->getModelPtr()->getColumnStatus(col) == ClpSimplex::Status::isFixed);
}

/**
 * Status of slack variables
 * The only tricky aspect is that the model pointer row status for slack variables
 * at ub is actually that they are at lb, and vice versa.
 * Basically, when Clp says a row is atLowerBound,
 * then it means it is a >= row, so that rowActivity == rowLower.
 * In Osi, the status of the corresponding slack is at its UPPER bound,
 * because Osi keeps all slacks with a +1 coefficient in the row.
 * When it says a row is atUpperBound,
 * then the row is a <= row and rowActivity == rowUpper.
 * In Osi, the status of the corresponding slack is at its LOWER bound,
 * since this is the regular slack we know and love.
 */
inline bool isNonBasicFreeSlack(const OsiClpSolverInterface* const solver,
    const int row) {
  return (solver->getModelPtr()->getRowStatus(row) == ClpSimplex::Status::isFree);
}

inline bool isBasicSlack(const OsiClpSolverInterface* const solver, const int row) {
  return (solver->getModelPtr()->getRowStatus(row) == ClpSimplex::Status::basic);
}

// NOTE: We do at *lower* bound on purpose because of comment above:
// the underlying model flips how it holds these with respect to Clp
inline bool isNonBasicUBSlack(const OsiClpSolverInterface* const solver,
    const int row) {
  const ClpSimplex::Status status = solver->getModelPtr()->getRowStatus(row);
  return ((status == ClpSimplex::Status::atLowerBound)
      || (status == ClpSimplex::Status::isFixed
          && greaterThanVal(solver->getRowPrice()[row], 0.)));
}

// NOTE: We do at *upper* bound on purpose because of comment above:
// the underlying model flips how it holds these with respect to Clp
inline bool isNonBasicLBSlack(const OsiClpSolverInterface* const solver,
    const int row) {
  const ClpSimplex::Status status = solver->getModelPtr()->getRowStatus(row);
  return ((status == ClpSimplex::Status::atUpperBound)
      || (status == ClpSimplex::Status::isFixed
          && !greaterThanVal(solver->getRowPrice()[row], 0.)));
}

inline bool isNonBasicFixedSlack(const OsiClpSolverInterface* const solver,
    const int row) {
  return (solver->getModelPtr()->getRowStatus(row) == ClpSimplex::Status::isFixed);
}

inline bool isBasicVar(const OsiClpSolverInterface* const solver, const int col) {
  return (solver->getModelPtr()->getStatus(col) ==  ClpSimplex::Status::basic);
}

/**
 * The variables that lead to rays are non-basic structural variables
 * and non-basic slacks not coming from equality constraints.
 */
bool isRayVar(const OsiClpSolverInterface* const solver, const int var);

inline bool isNonBasicLBVar(const OsiClpSolverInterface* const solver,
    const int var) {
  return
      (var < solver->getNumCols()) ?
          isNonBasicLBCol(solver, var) :
          isNonBasicLBSlack(solver, var - solver->getNumCols());
}

inline bool isNonBasicUBVar(const OsiClpSolverInterface* const solver,
    const int var) {
  return
      (var < solver->getNumCols()) ?
          isNonBasicUBCol(solver, var) :
          isNonBasicUBSlack(solver, var - solver->getNumCols());
}

inline bool isNonBasicFreeVar(const OsiClpSolverInterface* const solver,
    const int var) {
  return
      (var < solver->getNumCols()) ?
          isNonBasicFreeCol(solver, var) :
          isNonBasicFreeSlack(solver, var - solver->getNumCols());
}

inline void setMessageHandler(PHASolverInterface* const solver,
    const double max_time = std::numeric_limits<double>::max()) {
#ifndef TRACE
  solver->messageHandler()->setLogLevel(0);
  solver->getModelPtr()->messageHandler()->setLogLevel(0);
#endif
  solver->getModelPtr()->setMaximumSeconds(max_time);
}

void copySolver(PHASolverInterface* const solverCopy,
    const PHASolverInterface* const solverOld);

/**
 * Returns the index in raysCutByHplane of the first ray that is cut by this hyperplane in this activation
 */
int getStartIndex(const Hplane& hplane, //const std::vector<Hplane>& hplaneStore, const int hplane_ind,
    const int split_ind, const int act_ind);

/**
 * Returns the index in raysCutByHplane of the last ray that is cut by this hyperplane in this activation
 */
int getEndIndex(const Hplane& hplane, //const std::vector<Hplane>& hplaneStore, const int hplane_ind,
    const int split_ind, const int act_ind);

bool duplicatePoint(const CoinPackedMatrix& interPtsAndRays,
    const std::vector<double>& rhs, const Point& newPoint);

bool duplicateRay(const CoinPackedMatrix& interPtsAndRays,
    const std::vector<double>& rhs, const Ray& newRay);

int findDuplicateRay(const CoinPackedMatrix& interPtsAndRays,
    const std::vector<double>& rhs, const Ray& newRay);

bool rayCutByHplaneForSplit(const int oldHplaneIndex,
    const std::vector<Hplane>& hplaneStore, const int ray_ind,
    const std::vector<Ray>& rayStore, const int split_ind);

/*
 * Method checks if a hyperplane has been previous activated.
 * Only checks indices specified
 * @return the index if yes and -1 otherwise
 */
int hasHplaneBeenActivated(const std::vector<Hplane>& hplanesPrevAct,
    const std::vector<int> actHplaneIndex, const Hplane& tmpHplane,
    const bool returnIndexWIActivatedVec);

/*
 * Method checks if a hyperplane has been previous activated.
 * @return the index if yes and -1 otherwise
 */
int hasHplaneBeenActivated(const std::vector<Hplane>& hplanesPrevAct,
    const Hplane& tmpHplane, const int start_ind = 0);

/**
 * Method checks if a vertex has been previously created
 * @return the index if yes and -1 otherwise
 */
int hasVertexBeenActivated(const std::vector<Vertex>& verticesPrevAct,
    const Vertex& tmpVertex);

/********************/
/* RAY FUNCTIONS */
/********************/
bool calcIntersectionPointWithSplit(Point& intPt, double& dist,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
    const Vertex& vertex, const Ray& ray, const int split_var, const bool recalc_dist = true);
bool calcIntersectionPointWithSplit(Point& intPt,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
    const Vertex& vertex, const Ray& ray, const int split_var);

/**
 * Starting at origin, we proceed along ray r^j until we intersect H_h.
 * This leads to vertex v^{jh}.
 * Thus, we want to perform the pivot in which j enters the basis and h leaves it.
 * @param rayOutIndex :: the ray index r^j we want to pivot out (one component, in jth index, of 1.0 value)
 *
 * @return Ray \ell from pivot in which j leaves basis and h enters.
 */
void calcNewRayCoordinatesRank1(CoinPackedVector& rayIn, //const PointCutsSolverInterface* const solver,
    const SolutionInfo& solnInfo, const int rayOutIndex,
    const int newRayIndex, //const std::vector<Ray>& rayStore,
    const Hplane& hplane, const double RAYEPS);

/********************/
/* HPLANE FUNCTIONS */
/********************/

/**
 * @brief Reimplemented from AdvCuts class because of precision issues with taking inverses
 *
 * @param split_var :: split
 * @param vec :: must be sorted in increasing order (in terms of indices) and be in NB space
 */
double getSICActivity(const int split_var,
    const std::vector<Vertex>& vertexStore,
    const std::vector<Ray>& rayStore, const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const CoinPackedVector& vec);

/**
 * Finds coordinates of new vertex obtained by proceeding along rayOut from vertex a distance of dist
 */
void calcNewVectorCoordinates(CoinPackedVector& newVector, const double dist,
    const CoinPackedVector& vertex, const CoinPackedVector& rayOut, const double EPS = param.getEPS());

double distanceToSplit(const CoinPackedVector& vertex,
    const CoinPackedVector& rayOut, const int split_var,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo);

double dotProductWithOptTableau(const CoinPackedVector& vec, const int row_ind,
    const SolutionInfo& solnInfo);

/********************/
/* POINT STATISTICS */
/********************/

bool pointInP(const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const IntersectionInfo& interPtsAndRays,
    const int pt_ind);

bool pointInP(const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const int numElmts, const int* elmtIndices,
    const double* elmtVals);

inline bool pointInP(const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const CoinPackedVector& vec) {
  return pointInP(solver, solnInfo, vec.getNumElements(), vec.getIndices(),
      vec.getElements());
}

bool isRayFinal(const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const int ray_ind);
bool isRayFinal(const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const IntersectionInfo& interPtsAndRays,
    const int pt_ind);
bool isRayFinal(const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const int numElmts, const int* elmtIndices,
    const double* elmtVals);
inline bool isRayFinal(const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const CoinPackedVector& vec) {
  return isRayFinal(solver, solnInfo, vec.getNumElements(), vec.getIndices(),
      vec.getElements());
}

/*********************/
/* GENERAL FUNCTIONS */
/*********************/

/***********************************************************************/
/**
 * @brief Get environment variable by key.
 * @param key :: Key to look up.
 */
std::string get_env_var(std::string const & key);

/***********************************************************************/
/**
 * @brief Find value in a std::vector and return where it is in the std::vector.
 * @param key :: Key to look up.
 */
int find_val(const int &key, const std::vector<int> &vec);

/***********************************************************************/
/**
 * @brief Find value in a std::vector and return where it is in the std::vector.
 * @param key :: Key to look up.
 */
int find_val(const double &key, const std::vector<double> &vec);

/***********************************************************************/
/**
 * @brief Find value in a std::vector and return where it is in the std::vector.
 * @param key :: Key to look up.
 */
int find_val(const std::string &key, const std::vector<std::string> &vec);

/****************************************************************************/
/**
 * @brief Check if a given file exists in the system.
 *
 * @returns True if exists, False otherwise
 */
bool fexists(const char *filename);

double computeMinCutEucl(const AdvCuts &cutinfo);
double computeMaxCutEucl(const AdvCuts &cutinfo);
double computeAverageCutEucl(const AdvCuts &cutinfo);
double computeMinCutObj(const AdvCuts &cutinfo);
double computeMaxCutObj(const AdvCuts &cutinfo);
double computeAverageCutObj(const AdvCuts &cutinfo);
//
//inline bool compareAsc(ptSortInfo *tmp1, ptSortInfo *tmp2) {
//  if ((*tmp1).cost <= (*tmp2).cost - param.getEPS())
//    return true;
//  else
//    return false;
//}
//
//inline bool compareDesc(ptSortInfo *tmp1, ptSortInfo *tmp2) {
//  if ((*tmp1).cost >= (*tmp2).cost + param.getEPS())
//    return true;
//  else
//    return false;
//}

template<class T> struct index_cmp {

  const T arr; // array used to be sorted

  // Constructor
  index_cmp(const T _arr) :
      arr(_arr) {
  }

  // Comparator
  bool operator()(const size_t a, const size_t b) const {
    return arr[a] > arr[b] + param.getEPS();
  }
};

template<class T> struct index_cmp_asc {

  const T arr; // array used to be sorted

  // Constructor
  index_cmp_asc(const T _arr) :
      arr(_arr) {
  }

  // Comparator
  bool operator()(const size_t a, const size_t b) const {
    return arr[a] < arr[b] - param.getEPS();
  }
};

template<class T1, class T2> struct index_cmp2 {

  const T1 arr1; // array used to be sorted
  const T2 arr2;
  bool flag;
  // Constructor
  index_cmp2(const T1 _arr1, const T2 _arr2) :
      arr1(_arr1), arr2(_arr2) {
    flag = false;
  }

  // Comparator
  bool operator()(const size_t a, const size_t b) {
    flag = false;
    if (arr1[a] - arr1[b] >= param.getEPS())
      flag = true;
    else {
      if ((arr1[a] - arr1[b]) > -param.getEPS() && (arr1[a] - arr1[b]) < param.getEPS()) {
        if (arr2[a] - arr2[b] <= -param.getEPS())
          flag = true;
      }
    }
    return flag;
  }
};

template<class T1, class T2> struct index_cmp3 {

  const T1 arr1; // array used to be sorted
  const T2 arr2;
  bool flag;
  // Constructor
  index_cmp3(const T1 _arr1, const T2 _arr2) :
      arr1(_arr1), arr2(_arr2) {
    flag = false;
  }

  // Comparator
  bool operator()(const size_t a, const size_t b) {
    flag = false;
    if (arr1[a] - arr1[b] <= -param.getEPS())
      flag = true;
    else {
      if ((arr1[a] - arr1[b]) > -param.getEPS() && (arr1[a] - arr1[b]) < param.getEPS()) {
        if (arr2[a] - arr2[b] >= param.getEPS())
          flag = true;
      }
    }
    return flag;
  }
};

template<class T1, class T2> struct index_cmp4 {

  const T1 arr1; // array used to be sorted
  const T2 arr2;
  bool flag;
  // Constructor
  index_cmp4(const T1 _arr1, const T2 _arr2) :
      arr1(_arr1), arr2(_arr2) {
    flag = false;
  }

  // Comparator
  bool operator()(const size_t a, const size_t b) {
    flag = false;
    if (arr1[a] - arr1[b] >= param.getEPS())
      flag = true;
    else {
      if ((arr1[a] - arr1[b]) > -param.getEPS() && (arr1[a] - arr1[b]) < param.getEPS()) {
        if (arr2[a] - arr2[b] >= param.getEPS())
          flag = true;
      }
    }
    return flag;
  }
};

template<class T1, class T2, class T3, class T4> struct index_cmp5 {

  const T1 arr1; // array used to be sorted
  const T2 arr2;
  const T3 arr3;
  const T4 arr4;
  bool flag;
  // Constructor
  index_cmp5(const T1 _arr1, const T2 _arr2, const T3 _arr3, const T4 _arr4) :
      arr1(_arr1), arr2(_arr2), arr3(_arr3), arr4(_arr4) {
    flag = false;
  }

  // Comparator
  bool operator()(const size_t a, const size_t b) {
    flag = false;
    if (arr1[a] - arr1[b] >= param.getEPS())
      flag = true;
    else {
      if ((arr1[a] - arr1[b]) > -param.getEPS() && (arr1[a] - arr1[b]) < param.getEPS()) {
        if (arr2[a] - arr2[b] >= param.getEPS())
          flag = true;
        else {
          if ((arr2[a] - arr2[b]) > -param.getEPS()
              && (arr2[a] - arr2[b]) < param.getEPS()) {
            if (arr3[a] - arr3[b] >= param.getEPS())
              flag = true;
            else {
              if ((arr3[a] - arr3[b]) > -param.getEPS()
                  && (arr3[a] - arr3[b]) < param.getEPS()) {
                if (arr4[a] - arr4[b] >= param.getEPS())
                  flag = true;
              }
            }
          }
        }
      }

    }
    return flag;
  }
};

double computeAverage(const std::vector<double> &vec);

double computeAverage(const std::vector<int> &vec);

double computeAverage(const std::vector<double> &vec, const int start_index,
    const int end_index);

double computeAverage(const std::vector<int> &vec, const int start_index,
    const int end_index);

double computeAverage(const std::vector<double> &vec,
    const std::vector<int> &index);

double computeAverage(const std::vector<int> &vec,
    const std::vector<int> &index);

double computeAverageDiff(const std::vector<double> &vec1,
    const std::vector<double> &vec2);

double computeAverageDiff(const std::vector<int> &vec1,
    const std::vector<int> &vec2);

double computeAverageDiff(const std::vector<double> &vec1,
    const std::vector<double> &vec2, const int start_index,
    const int end_index);

double computeAverageDiff(const std::vector<int> &vec1,
    const std::vector<int> &vec2, const int start_index,
    const int end_index);

double computeAverageDiff(const std::vector<double> &vec1,
    const std::vector<double> &vec2, const std::vector<int> &index);

double computeAverageDiff(const std::vector<int> &vec1,
    const std::vector<int> &vec2, const std::vector<int> &index);

std::vector<double>::const_iterator computeMin(const std::vector<double> &vec);

std::vector<int>::const_iterator computeMin(const std::vector<int> &vec);

double computeMinAbs(const int num_el, const double* vec);

double computeMaxAbs(const int num_el, const double* vec);

std::vector<double>::const_iterator computeMin(const std::vector<double> &vec,
    const std::vector<int> &index);

std::vector<int>::const_iterator computeMin(const std::vector<int> &vec,
    const std::vector<int> &index);

std::vector<double>::const_iterator computeMax(const std::vector<double> &vec);

std::vector<int>::const_iterator computeMax(const std::vector<int> &vec);

std::vector<double>::const_iterator computeMax(const std::vector<double> &vec,
    const std::vector<int> &index);

std::vector<int>::const_iterator computeMax(const std::vector<int> &vec,
    const std::vector<int> &index);

double computeMinDiff(const std::vector<double> &vec1,
    const std::vector<double> &vec2);

int computeMinDiff(const std::vector<int> &vec1, const std::vector<int> &vec2);

double computeMinDiff(const std::vector<double> &vec1,
    const std::vector<double> &vec2, const std::vector<int> &index);

int computeMinDiff(const std::vector<int> &vec1, const std::vector<int> &vec2,
    const std::vector<int> &index);

double computeMaxDiff(const std::vector<double> &vec1,
    const std::vector<double> &vec2);

int computeMaxDiff(const std::vector<int> &vec1, const std::vector<int> &vec2);

double computeMaxDiff(const std::vector<double> &vec1,
    const std::vector<double> &vec2, const std::vector<int> &index);

int computeMaxDiff(const std::vector<int> &vec1, const std::vector<int> &vec2,
    const std::vector<int> &index);

/********************/
/* CONVERSION FUNCS */
/********************/

/***********************************************************************/
/**
 * The solution to the LP relaxation of the problem has all non-basic
 *  variables are at zero, but these are not all zero when we compute
 *  rays that intersect with bd S. Since this computation is done all
 *  in the non-basic space, but we are using the structural space for
 *  computing the lifting, we need to find the coordinates of the
 *  intersection points in the structural space.
 *
 *  @param solver       ::  Solver we are using
 *  @param rowOfVar     ::  Row in which variable is basic (indexed from 0 to numcols)
 *  @param newpt        ::  [point][value]
 *  @param nbinfo       ::  [point][index of nb ray][value]
 */
void compNBCoor(std::vector<double>& NBSolution, const double* structSolution,
    const PHASolverInterface* const tmpSolver,
    const PHASolverInterface* const origSolver,
    const SolutionInfo& origSolnInfo);

/********************/
/* STRING FUNCTIONS */
/********************/

/***********************************************************************/
/**
 * Parses string line by the set of delimiters given.
 * (Only works well with one for now.)
 * Parses string based on delimiters provided and puts output into vector.
 */
void stringparse(std::string &line, std::string delims,
    std::vector<std::string> &strings);

/****************************************************************************/
double compute_condition_number_norm1(OsiSolverInterface *solver);

/****************************************************************************/
double compute_condition_number_norm2(const OsiSolverInterface* const solver);

/****************************************************************************/
double compute_condition_number(OsiSolverInterface *solver, const int norm);

std::string toLowerString(const std::string str);
std::string toUpperString(const std::string str);
std::string toLowerStringNoUnderscore(const std::string str);
std::string toUpperStringNoUnderscore(const std::string str);
std::vector<std::string> lowerCaseStringVector(const std::vector<std::string>& strVec);
double getObjValueFromFile(const char* filename, const std::string& this_inst_name);
void getOptSolFromFile(const char* filename, std::vector<double>& opt_sol);

/** From Giacomo Nannicini */
double dotProduct(const double* a, const double* b, int dimension);
double dotProduct(int sizea, const int* indexa, const double* a, const double* b);
double dotProduct(int sizea, const int* indexa, const double* a, int sizeb,
    const int* indexb, const double* b);
