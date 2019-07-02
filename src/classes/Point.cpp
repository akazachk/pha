//============================================================================
// Name        : Point.cpp
// Author      : Aleksandr M. Kazachkov
// Version     : 0.2014.03.08
// Copyright   : Your copyright notice
// Description : Class Point and related methods
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#include "GlobalConstants.hpp"
#include <cmath>
#include "Point.hpp"
#include "Utility.hpp"

Point::~Point() {
}

/**
 * Input has to have the indices sorted in the same order (ideally increasing)
 * Returns true if one point is exactly same as the other (no scaling for the point)
 */
bool Point::pointsDifferent(const Point& pt1, const CoinPackedVector& pt2) {
  return isRowDifferent(pt1, 1., pt2, 1., param.getEPS()) != 0;
} /* pointsDifferent */

const std::string Point::vectorString() const {
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

#ifdef TRACE
void Point::printSparseVector() const {
  printSparseVector(stdout);
}

void Point::printSparseVector(FILE* fptr) const {
  int numElems = CoinPackedVector::getNumElements();
  const int* index = CoinPackedVector::getIndices();
  const double* element = CoinPackedVector::getElements();
  printf("Num elements is %d.\n", numElems);
  for (int i = 0; i < numElems; i++) {
    fprintf(fptr, "\t%d: %f\n", index[i], element[i]);
  }
}
#endif
