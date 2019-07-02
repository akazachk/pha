//============================================================================
// Name        : Vertex.hpp
// Author      : Aleksandr M. Kazachkov
// Version     : 0.2014.05.23
// Copyright   : Your copyright notice
// Description : Class Point and related methods
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#pragma once

#include "CoinPackedVector.hpp"
#include <string>
#include <vector>

class Vertex: public CoinPackedVector {
public:
// These are the indices, in the vector of rays, of all rays emanating from this point
  std::vector<int> rayIndices;

  // Stores whether each ray should still be considered
  std::vector<bool> rayActiveFlag;

  Vertex() {};
  ~Vertex();

  Vertex(const std::vector<double>& coeff);
  Vertex(int num_coeff, const double* coeff);

  void setPackedVector(int num_coeff, const double* coeff);

  /**
   * Overloading the standard CoinPackedVectorBase method since that one
   * uses std::equal, which is fine for the indices, but for the double vec
   * it is checking exact equality, and we are okay with eps difference.
   */
  inline bool operator==(const Vertex& pt2) const {
    return !verticesDifferent(*this, pt2);
  }

  inline bool operator!=(const Vertex& pt2) const {
    return !((*this) == pt2);
  }

  static bool verticesDifferent(const Vertex& pt1, const Vertex& pt2);

#ifdef TRACE
  void printSparseVector() const;
  void printSparseVector(FILE* fptr) const;
  const std::string vectorString() const;
#endif
};
