//============================================================================
// Name        : Ray.hpp
// Author      : Aleksandr M. Kazachkov
// Version     : 2018.08.22
// Description : Class Ray and related methods
//============================================================================

#pragma once

#include "OsiSolverInterface.hpp"
#include "CoinPackedVector.hpp"
#include "SolutionInfo.hpp"

#include "CoinPackedVector.hpp"
#include <string>
#include <vector>

class Ray: public CoinPackedVector {
public:
  // For each split, a quick way to check whether we have cut the ray or not
  std::vector<int> cutFlag;

  Ray() {};
  ~Ray();

  /**
   * Overloading the standard CoinPackedVectorBase method since that one
   * uses std::equal, which is fine for the indices, but for the double vec
   * it is checking exact equality, and we are okay with eps difference.
   */
  inline bool operator==(const CoinPackedVector& ray2) const {
    return !raysDifferent(*this, ray2);
  }

  inline bool operator!=(const Ray& ray2) const {
    return !((*this) == ray2);
  }

  static bool raysDifferent(const Ray& ray1, const CoinPackedVector& ray2);

#ifdef TRACE
  void printSparseVector() const;
  void printSparseVector(FILE* fptr) const;
  const std::string vectorString() const;
#endif
};
