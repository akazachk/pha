//============================================================================
// Name        : Point.hpp
// Author      : Aleksandr M. Kazachkov
// Version     : 2018.08.22
// Description : Class Point and related methods
//============================================================================

#pragma once

#include "OsiSolverInterface.hpp"
#include "CoinPackedVector.hpp"
#include "SolutionInfo.hpp"

#include <string>
#include <vector>

class Point: public CoinPackedVector {
public:
  // These are the indices, in the vector of rays, of all rays emanating from this point
  Point() {};
  ~Point();

  /**
   * Overloading the standard CoinPackedVectorBase method since that one
   * uses std::equal, which is fine for the indices, but for the double vec
   * it is checking exact equality, and we are okay with eps difference.
   */
  inline bool operator==(const CoinPackedVector& pt2) const {
    return !pointsDifferent(*this, pt2);
  }

  inline bool operator!=(const Point& pt2) const {
    return !((*this) == pt2);
  }

  static bool pointsDifferent(const Point& pt1, const CoinPackedVector& pt2);

#ifdef TRACE
  void printSparseVector() const;
  void printSparseVector(FILE* fptr) const;
#endif
  const std::string vectorString() const;
};
