//============================================================================
// Name        : Ray.cpp
// Author      : Aleksandr M. Kazachkov
// Version     : 0.2014.03.08
// Copyright   : Your copyright notice
// Description : Class Ray and related methods
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#include "GlobalConstants.hpp"
#include <cmath>
#include "Ray.hpp"
#include "Utility.hpp"

Ray::~Ray() {
}

/**
 * Input has to have the indices sorted in the same order (ideally increasing)
 * Also all zeroes should have been eliminated
 * Returns true if one point is same as the other (up to scaling for the ray)
 */
bool Ray::raysDifferent(const Ray& ray1, const CoinPackedVector& ray2) {
  return isRowDifferent(ray1, 0., ray2, 0., param.getEPS()) != 0;
} /* raysDifferent */

#ifdef TRACE
void Ray::printSparseVector() const {
  printSparseVector(stdout);
}

void Ray::printSparseVector(FILE* fptr) const {
  int numElems = CoinPackedVector::getNumElements();
  const int* index = CoinPackedVector::getIndices();
  const double* element = CoinPackedVector::getElements();
  printf("Num elements is %d.\n", numElems);
  for (int i = 0; i < numElems; i++) {
    fprintf(fptr, "\t%d: %f\n", index[i], element[i]);
  }
}

const std::string Ray::vectorString() const {
  const int numElems = this->getNumElements();
  const int* vecIndex = this->getIndices();
  const double* vecVal = this->getElements();
  char vec_string[200];
  if (numElems == 0) {
    snprintf(vec_string, sizeof(vec_string) / sizeof(char),
        "(zero vector)");
  } else if (numElems == 1) {
    snprintf(vec_string, sizeof(vec_string) / sizeof(char), "(%d: %s)",
        vecIndex[0], stringValue(vecVal[0]).c_str());
  } else if (numElems == 2) {
    snprintf(vec_string, sizeof(vec_string) / sizeof(char),
        "(%d: %s, %d: %s)", vecIndex[0], stringValue(vecVal[0]).c_str(),
        vecIndex[1], stringValue(vecVal[1]).c_str());
  }
  std::string tmp(vec_string);
  return tmp;
}
#endif
