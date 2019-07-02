//============================================================================
// Name        : Hplane.hpp
// Author      : Aleksandr M. Kazachkov
// Version     : 0.2014.05.23
// Copyright   : Your copyright notice
// Description : Hyperplane class
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#pragma once

#include <vector>
#include "CoinPackedVector.hpp"
#include "typedefs.hpp"
#include "SolutionInfo.hpp"

enum class HplaneVarFlag {
  none = -1
};

enum class HplaneRowFlag {
  none = -2,
  NB
};

enum class HplaneBoundFlag {
  none = -3,
  toLBNB,
  toLBSlack,
  toLB,
  toUB,
  toUBNB
};

class Hplane {
public:
  // Constructors
  // variable index, row index, UB flag, bound
  Hplane() :
      Hplane(HplaneVarFlag::none, HplaneRowFlag::none, HplaneBoundFlag::none,
          0.) {
  }
  ;
  Hplane(const int numSplits) :
      Hplane(HplaneVarFlag::none, HplaneRowFlag::none, HplaneBoundFlag::none,
          0., numSplits) {
  }
  ;
  Hplane(const HplaneVarFlag tmpVar, const HplaneRowFlag tmpRow,
      const HplaneBoundFlag tmpUBFlag, const int numSplits = 0) :
      Hplane(static_cast<int>(tmpVar), static_cast<int>(tmpRow),
          static_cast<int>(tmpUBFlag), 0., numSplits) {
  }
  Hplane(const int tmpVar, const int tmpRow, const int tmpUBFlag, const int numSplits = 0) :
    Hplane(tmpVar, tmpRow, tmpUBFlag, 0., numSplits) {
  }
  Hplane(const HplaneVarFlag tmpVar, const HplaneRowFlag tmpRow,
      const HplaneBoundFlag tmpUBFlag, const double tmpBound,
      const int numSplits = 0) :
      Hplane(static_cast<int>(tmpVar), static_cast<int>(tmpRow),
          static_cast<int>(tmpUBFlag), tmpBound, numSplits) {
  }
  Hplane(const int tmpVar, const int tmpRow, const int tmpUBFlag,
      const double tmpBound, const int numSplits = 0);

  ~Hplane();

  // var: variable index
  // row: row in which it is basic (-1 if in none; i.e., is nb)
  // ubflag: 2 if nb to ub (ub - lb), 1 if ub, 0 if lb,
  //        -1 if slack (so to lb = 0), -2 if nb to lb (= 0)
  int var, row, ubflag;
  double bound;

  // [split][index] = vert/ray index
  // To get the ray we need the vertex + ray index together,
  // which give the index of the corresponding ray coordinates in rayStore.
  std::vector<std::vector<int> > rayToBeCutByHplane;
  std::vector<std::vector<int> > vertOfRayToBeCutByHplane;

  // [split][index] = vert index
  // Index of vertex created by that ray being cut by this hplane
  std::vector<std::vector<int> > vertCreatedByRayToBeCutByHplane;

  // These are *all* the rays that the hyperplane intersects
  // Useful for quick computation of new rays
  std::vector<std::vector<int> > allRaysIntByHplane;
  std::vector<std::vector<int> > vertOfAllRaysIntByHplane;

  // [split][act] = num hplanes cut in all activations <= act
  // For keeping where next activation starts
  std::vector<std::vector<int> > numRaysCutByHplaneInclPrevAct;

  inline bool operator==(const Hplane& hplane2) const {
    return (this->var == hplane2.var) && (this->ubflag == hplane2.ubflag);
  }

  void resizeVectors(const int numSplits = 0);

  void setHplaneBoundFromUBFlag(const PHASolverInterface* const solver);

  /**
   * distanceToHplane 0
   * Given current value varVal and the value of the ray rayVal for this hplane,
   * calculate the distance to the current hyperplane
   */
  double distanceToHplane(const double varVal, const double rayVal,
      const PHASolverInterface* const solver) const;

  /**
   * distanceToHplane 1
   * Assume starting at optimal vertex
   * The ray is given by the variable (in the structural space)
   */
  double distanceToHplane(const int structRayOutVar,
      const PHASolverInterface* const solver) const;

  /**
   * distanceToHplane 2
   * Given current value varVal of the hyperplane's variable
   * and the ray given by the variable (in the struct space),
   * calculate the distance to the current hyperplane
   */
  double distanceToHplane(const double varVal, const int structRayOutVar,
      const PHASolverInterface* const solver) const;

  /**
   * distanceToHplane 3
   * Assume starting at optimal vertex
   * The ray is given by the CoinPackedVector rayOut
   * This could be in the structural space or the non-basic space
   * If it is in the non-basic space, then the set of non-basic variables needs to be provided
   * In addition, the solver should be optimal with the cobasis being the same as the nb space
   */
  double distanceToHplane(const CoinPackedVector& rayOut,
      const PHASolverInterface* const solver, const bool inNBSpace,
      const std::vector<int>* nonBasicVarIndex = NULL) const;

  /**
   * distanceToHplane 4
   * Starting at given vertex, which could be provided in the struct or non-basic space
   * The ray is given by the CoinPackedVector rayOut
   * This also could be in either the struct or non-basic space
   * If it is in the non-basic space, then the set of non-basic variables needs to be provided
   * In addition, the solver should be optimal with the cobasis being the same as the nb space
   */
  double distanceToHplane(const CoinPackedVector &vertex,
      const CoinPackedVector& rayOut,
      const PHASolverInterface* const solver, const bool inNBSpace,
      const std::vector<int>* nonBasicVarIndex = NULL) const;
  double distanceToHplane(const CoinPackedVector &vertex,
        const CoinPackedVector& rayOut,
        const PHASolverInterface* const solver, const bool inNBSpace,
        const SolutionInfo& solnInfo) const;

  /**
   * distanceToHplane 5
   * Given current value varVal of the hyperplane's variable
   * and the ray given by the CoinPackedVector rayOut,
   * which could be in either the struct or non-basic space,
   * calculate the distance to the current hyperplane
   * If it is in the non-basic space, then the set of non-basic variables needs to be provided
   * In addition, the solver should be optimal with the cobasis being the same as the nb space
   */
  double distanceToHplane(const double varVal, const CoinPackedVector& rayOut,
      const PHASolverInterface* const solver, const bool inNBSpace,
      const std::vector<int>* nonBasicVarIndex = NULL) const;
  double distanceToHplane(const double varVal, const CoinPackedVector& rayOut,
        const PHASolverInterface* const solver, /*const bool inNBSpace,*/
        const SolutionInfo& solnInfo) const;

  /**
   * dotProductWithHplane
   */
  double dotProductWithHplane(const CoinPackedVector& vec,
      const PHASolverInterface* solver, const bool vecIsRay,
      const bool inNBSpace, const std::vector<int>* nonBasicVarIndex = NULL) const;

  double dotProductWithHplane(const CoinPackedVector& vec,
      const SolutionInfo& solnInfo,
      const bool vecIsRay/*, const bool inNBSpace*/) const;

  /**
   * getRowActivity
   */
  double getRowActivity(const CoinPackedVector& vertex,
      const PHASolverInterface* solver, const bool inNBSpace,
      const std::vector<int>* nonBasicVarIndex = NULL) const;

  inline void print() const {
    printf("(var: %d, row: %d, ubflag: %d, bound: %1.6f)\n", this->var,
        this->row, this->ubflag, this->bound);
  }

  inline void printAll() const {
      printf("(var: %d, row: %d, ubflag: %d, bound: %1.6f)\n", this->var,
          this->row, this->ubflag, this->bound);
      printf("\tCut rays: \n");
      for (int split_ind = 0; split_ind < (int) rayToBeCutByHplane.size(); split_ind++) {
        printf("\t\tSplit: %d\tRays:", split_ind);
        for (auto ray_ind : rayToBeCutByHplane[split_ind]) {
          printf(" %d", ray_ind);
        }
        printf("\n");
      }
    }

  inline static void printHplaneStore(const std::vector<Hplane>& store) {
    for (int h = 0; h < (int) store.size(); h++) {
      printf("%d: ", h);
      store[h].print();
    }
  }
}; /* Hplane */
