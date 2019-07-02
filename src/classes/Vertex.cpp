//============================================================================
// Name        : Vertex.cpp
// Author      : Aleksandr M. Kazachkov
// Version     : 2018.08.22
// Description : Class Point and related methods
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#include "GlobalConstants.hpp"
#include <cmath>
#include "GlobalConstants.hpp"
#include "Vertex.hpp"
#include "Utility.hpp"

Vertex::Vertex(const std::vector<double>& coeff) {
  setPackedVector(coeff.size(), coeff.data());
}

Vertex::Vertex(int num_coeff, const double* coeff) {
  setPackedVector(num_coeff, coeff);
}

Vertex::~Vertex() {
}

void Vertex::setPackedVector(int num_coeff, const double* coeff) {
  std::vector<int> index;
  std::vector<double> val;
  for (int i = 0; i < num_coeff; i++) {
    if (!isZero(coeff[i])) {
      index.push_back(i);
      val.push_back(coeff[i]);
    }
  }

  setVector((int) index.size(), index.data(), val.data());
}

/**
 * Input has to have the indices sorted in the same order (ideally increasing)
 * Returns true if one point is exactly same as the other (no scaling for the point)
 */
bool Vertex::verticesDifferent(const Vertex& pt1, const Vertex& pt2) {
  return isRowDifferent(pt1, 1., pt2, 1., param.getEPS()) != 0;
} /* verticesDifferent */

#ifdef TRACE
void Vertex::printSparseVector() const {
  printSparseVector(stdout);
}

void Vertex::printSparseVector(FILE* fptr) const {
  int numElems = CoinPackedVector::getNumElements();
  const int* index = CoinPackedVector::getIndices();
  const double* element = CoinPackedVector::getElements();
  printf("Num elements is %d.\n", numElems);
  for (int i = 0; i < numElems; i++) {
    fprintf(fptr, "\t%d: %f\n", index[i], element[i]);
  }
}

const std::string Vertex::vectorString() const {
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
