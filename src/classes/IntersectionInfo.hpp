//============================================================================
// Name        : IntersectionInfo.hpp
// Author      : Aleksandr M. Kazachkov
// Version     : 0.2014.06.12
// Copyright   : Your copyright notice
// Description : Class IntesectionInfo and related methods
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#pragma once

#include <vector>
#include <CoinPackedMatrix.hpp>

class IntersectionInfo: public CoinPackedMatrix {
public:
  // As many entries as activations
  std::vector<int> numPointsRemoved, numPointsAdded;
  std::vector<double> avgRemovedDepth, avgAddedDepth;

  std::vector<double> RHS;

  // Hplane information
  // As many entries as there are intersection points / rays (rows)
  // Gives index of hyperplane (if one exists) that cuts ray
  // -1 if the hplane does not exist
  // -2 if from the PCut procedure (unless we start keeping track of the activated hplanes)
  std::vector<int> hplane;
  std::vector<int> hplaneIndexForSplit;
  // Also we keep the corresponding created vertex index
  // -2 if from the PCut procedure 
  std::vector<int> vertex;
  // Also, we want to keep the unique set of indices of hplanes activated for the split
  std::vector<std::vector<int> > hplaneIndexByAct; // [act_ind][index] = hplane
  std::vector<int> allHplaneToAct; // [hplane] = hplane
  std::vector<std::vector<int> > indexOfHplaneInOverallOrderForSplit; // Needed to activate in proper order

  // Ray information
  // As many entries as there are intersection points / rays (rows)
  // cutRay is ray of C1 that is cut
  // newRay is along which coordinate direction in the nb space
  // we move after the cutting of cutRay
  // If no ray was cut (so this is a ray of C1), then cutRay = -1, and newRay is the ray of C1
  // If this is from the PCut procedure, both cutRay and newRay will be -2 for now (unless we start keeping track of these rays)
  std::vector<int> cutRay;
  std::vector<int> newRay;

  // Which facet the intersection point was created on
  // As many entries as there are intersection points / rays (rows)
  // -1 if it is a ray that is parallel to the split and generated from a vertex internal to the split
  // 0 if it is on the floor of the split
  // 1 if it is on the ceiling of the split
  std::vector<int> facetIndex;
  std::vector<int> multiplicity;

  // Information about point depth
  // Will not be deleted if the point is deleted; has to be manually controlled
  // Should only be set after all intersection points are generated so none are deleted
  std::vector<bool> finalFlag; // As many entries as there are intersection points / rays (rows)

  int numPoints, numFinalPoints, numRays, numFinalRays;
//  std::vector<int> pointRowIndex; // Only entries for the intersection points, not rays

  std::vector<double> SICDepth;
  std::vector<double> objDepth;
  double minSICDepth, maxSICDepth, totalSICDepth;
  double minObjDepth, maxObjDepth, totalObjDepth;
  int indexMinSICDepth, indexMaxSICDepth, indexMinObjDepth, indexMaxObjDepth;

  // CONSTRUCTOR
  inline IntersectionInfo() :
      numPoints(0), numFinalPoints(0), numRays(0), numFinalRays(0), minSICDepth(0), maxSICDepth(
          0), totalSICDepth(0), minObjDepth(0), maxObjDepth(0), totalObjDepth(
          0), indexMinSICDepth(-1), indexMaxSICDepth(-1), indexMinObjDepth(-1), indexMaxObjDepth(
          -1) {
  }

  IntersectionInfo & operator=(const IntersectionInfo& rhs) {
    if (this != &rhs) {
      CoinPackedMatrix::operator=(rhs);
      this->numPointsRemoved = rhs.numPointsRemoved;
      this->numPointsAdded = rhs.numPointsAdded;
      this->avgRemovedDepth = rhs.avgRemovedDepth;
      this->avgAddedDepth = rhs.avgAddedDepth;
      this->RHS = rhs.RHS;
      this->hplane = rhs.hplane;
      this->hplaneIndexForSplit = rhs.hplaneIndexForSplit;
      this->vertex = rhs.vertex;
      this->hplaneIndexByAct = rhs.hplaneIndexByAct;
      this->allHplaneToAct = rhs.allHplaneToAct;
      this->indexOfHplaneInOverallOrderForSplit =
          rhs.indexOfHplaneInOverallOrderForSplit;
      this->cutRay = rhs.cutRay;
      this->newRay = rhs.newRay;
      this->facetIndex = rhs.facetIndex;
      this->multiplicity = rhs.multiplicity;
      this->finalFlag = rhs.finalFlag;
      this->numPoints = rhs.numPoints;
      this->numFinalPoints = rhs.numFinalPoints;
      this->numRays = rhs.numRays;
      this->numFinalRays = rhs.numFinalRays;
      this->SICDepth = rhs.SICDepth;
      this->objDepth = rhs.objDepth;
      this->minSICDepth = rhs.minSICDepth;
      this->maxSICDepth = rhs.maxSICDepth;
      this->totalSICDepth = rhs.totalSICDepth;
      this->minObjDepth = rhs.minObjDepth;
      this->maxObjDepth = rhs.maxObjDepth;
      this->totalObjDepth = rhs.totalObjDepth;
      this->indexMinSICDepth = rhs.indexMinSICDepth;
      this->indexMaxSICDepth = rhs.indexMaxSICDepth;
      this->indexMinObjDepth = rhs.indexMinObjDepth;
      this->indexMaxObjDepth = rhs.indexMaxObjDepth;
    }
    return *this;
  }

  inline void reserveSpace(const int n) {
    this->RHS.reserve(n);
    this->hplane.reserve(n);
    this->hplaneIndexForSplit.reserve(n);
    this->vertex.reserve(n);
    this->cutRay.reserve(n);
    this->newRay.reserve(n);
    this->facetIndex.reserve(n);
    this->multiplicity.reserve(n);
//    this->pointRowIndex.reserve(n/2);
    this->finalFlag.reserve(n);
    this->SICDepth.reserve(n);
    this->objDepth.reserve(n);
  }

  inline void resizeInfoByActivationRound(const int num_act) {
    hplaneIndexByAct.resize(num_act);
    indexOfHplaneInOverallOrderForSplit.resize(num_act);
    numPointsAdded.resize(num_act);
    numPointsRemoved.resize(num_act);
    avgAddedDepth.resize(num_act);
    avgRemovedDepth.resize(num_act);
  }

  inline void addPointOrRayToIntersectionInfo(const CoinPackedVectorBase& vec,
      const double rhs, const int ray_to_cut, const int new_ray_ind,
      const int hplane_ind, const int hplane_ind_in_split, const int vert_ind,
      const int facet_ind, const double currSICDepth, const double currObjDepth,
      const bool currFinalFlag) {
    this->appendRow(vec);
    this->RHS.push_back(rhs);
    this->hplane.push_back(hplane_ind); // Store hplane associated with this point
    this->hplaneIndexForSplit.push_back(hplane_ind_in_split); // Index in allHplaneToAct
    this->vertex.push_back(vert_ind);
    this->cutRay.push_back(ray_to_cut);
    this->newRay.push_back(new_ray_ind);
    this->facetIndex.push_back(facet_ind);
    this->multiplicity.push_back(1); // First time it is added
    this->finalFlag.push_back(currFinalFlag);
    if (rhs < -1e-3 || rhs > 1e-3) {
      this->numPoints++;
      this->SICDepth.push_back(currSICDepth);
      this->objDepth.push_back(currObjDepth);
      this->totalSICDepth += currSICDepth;
      if (currSICDepth < this->minSICDepth) {
        this->minSICDepth = currSICDepth;
        this->indexMinSICDepth = this->getNumRows() - 1;
      }
      if (currSICDepth > this->maxSICDepth) {
        this->maxSICDepth = currSICDepth;
        this->indexMaxSICDepth = this->getNumRows() - 1;
      }
      this->totalObjDepth += currObjDepth;
      if (currObjDepth < this->minObjDepth) {
        this->minObjDepth = currObjDepth;
        this->indexMinObjDepth = this->getNumRows() - 1;
      }
      if (currObjDepth > this->maxObjDepth) {
        this->maxObjDepth = currObjDepth;
        this->indexMaxObjDepth = this->getNumRows() - 1;
      }
      this->numFinalPoints += currFinalFlag;
    } else {
      this->numRays++;
      this->SICDepth.push_back(0.0);
      this->objDepth.push_back(0.0);
      this->numFinalRays += currFinalFlag;
    }
  }

  inline void increaseMultiplicity(const int row) {
    this->multiplicity[row]++;
    if (this->RHS[row] < -1e-3 || this->RHS[row] > 1e-3) {
      numPoints++;
      numFinalPoints += finalFlag[row];
      totalSICDepth += SICDepth[row];
      totalObjDepth += objDepth[row];
    } else {
      numRays++;
      numFinalRays += finalFlag[row];
    }
  }

  inline void decreaseMultiplicity(const int row) {
    this->multiplicity[row]--;
    if (this->RHS[row] > -1e-3 || this->RHS[row] < 1e-3) {
      numPoints--;
      numFinalPoints -= finalFlag[row];
      totalSICDepth -= SICDepth[row];
      totalObjDepth -= objDepth[row];
    } else {
      numRays--;
      numFinalRays -= finalFlag[row];
    }
  }

  inline void setMultiplicity(const int row, const int num) {
    const int oldMult = this->multiplicity[row];
    const int delta = num - oldMult;
    this->multiplicity[row] = num;
    if (this->RHS[row] > -1e-3 || this->RHS[row] < 1e-3) {
      numPoints += delta;
      numFinalPoints += finalFlag[row] * delta;
      totalSICDepth += SICDepth[row] * delta;
      totalObjDepth += objDepth[row] * delta;
    } else {
      numRays += delta;
      numFinalRays += finalFlag[row] * delta;
    }
  }

  inline void deleteRowInfo(const int numDel, const int *indDel,
      const bool checkMultiplicity = false) {
    int realNumberToDelete = numDel;
    std::vector<int> indicesToDelete;

    if (checkMultiplicity) {
      indicesToDelete.reserve(numDel);
      for (int i = 0; i < numDel; i++) {
        const int row = indDel[i];
        if (this->multiplicity[row] > 1) {
          this->multiplicity[row]--;
          if (this->RHS[row] < -1e-3 || this->RHS[row] > 1e-3) {
            numPoints--;
            numFinalPoints -= finalFlag[row];
            totalSICDepth -= SICDepth[row];
            totalObjDepth -= objDepth[row];
          } else {
            numRays--;
            numFinalRays -= finalFlag[row];
          }
        } else {
          indicesToDelete.push_back(row);
        }
      }
      indDel = indicesToDelete.data(); // it is int const *, not int const * const!
      realNumberToDelete = indicesToDelete.size();
    }

    this->deleteRows(realNumberToDelete, indDel);

    for (int p = realNumberToDelete - 1; p > -1; p--) {
      const int row = indDel[p];
      if (this->RHS[row] < -1e-3 || this->RHS[row] > 1e-3) {
        numPoints--;
        numFinalPoints -= finalFlag[row];
        totalSICDepth -= SICDepth[row];
        totalObjDepth -= objDepth[row];
        if (indexMinSICDepth == row) {
          indexMinSICDepth = -1;
        }
        if (indexMaxSICDepth == row) {
          indexMaxSICDepth = -1;
        }
        if (indexMinObjDepth == row) {
          indexMinObjDepth = -1;
        }
        if (indexMaxObjDepth == row) {
          indexMaxObjDepth = -1;
        }
      } else {
        numRays--;
        numFinalRays -= finalFlag[row];
      }
      this->RHS.erase(this->RHS.begin() + indDel[p]);
      this->hplane.erase(this->hplane.begin() + indDel[p]);
      this->hplaneIndexForSplit.erase(
          this->hplaneIndexForSplit.begin() + indDel[p]);
      this->cutRay.erase(this->cutRay.begin() + indDel[p]);
      this->newRay.erase(this->newRay.begin() + indDel[p]);
      this->facetIndex.erase(this->facetIndex.begin() + indDel[p]);
      this->multiplicity.erase(this->multiplicity.begin() + indDel[p]);
      this->SICDepth.erase(this->SICDepth.begin() + indDel[p]);
      this->objDepth.erase(this->objDepth.begin() + indDel[p]);
      this->finalFlag.erase(this->finalFlag.begin() + indDel[p]);
    }
  }


  inline const CoinShallowPackedVector operator[](const int index) {
    return this->getVector(index);
  }

  inline double getMinSICDepth() {
    if (indexMinSICDepth < 0) {
      // Need to go through all of the rows
      if (this->getNumRows() == 0) {
        return 0.0;
      }
      minSICDepth = SICDepth[0];
      indexMinSICDepth = 0;
      for (int row = 1; row < this->getNumRows(); row++) {
        if (SICDepth[row] < minSICDepth) {
          minSICDepth = SICDepth[row];
          indexMinSICDepth = row;
        }
      }
    }
    return minSICDepth;
  }

  inline double getMaxSICDepth() {
    if (indexMaxSICDepth < 0) {
      // Need to go through all of the rows
      if (this->getNumRows() == 0) {
        return 0.0;
      }
      maxSICDepth = SICDepth[0];
      indexMaxSICDepth = 0;
      for (int row = 1; row < this->getNumRows(); row++) {
        if (SICDepth[row] > maxSICDepth) {
          maxSICDepth = SICDepth[row];
          indexMaxSICDepth = row;
        }
      }
    }
    return maxSICDepth;
  }

  inline double getMinObjDepth() {
    if (indexMinObjDepth < 0) {
      // Need to go through all of the rows
      if (this->getNumRows() == 0) {
        return 0.0;
      }
      minObjDepth = objDepth[0];
      indexMinObjDepth = 0;
      for (int row = 1; row < this->getNumRows(); row++) {
        if (objDepth[row] < minObjDepth) {
          minObjDepth = objDepth[row];
          indexMinObjDepth = row;
        }
      }
    }
    return minObjDepth;
  }

  inline double getMaxObjDepth() {
    if (indexMaxObjDepth < 0) {
      // Need to go through all of the rows
      if (this->getNumRows() == 0) {
        return 0.0;
      }
      maxObjDepth = objDepth[0];
      indexMaxObjDepth = 0;
      for (int row = 1; row < this->getNumRows(); row++) {
        if (objDepth[row] > maxObjDepth) {
          maxObjDepth = objDepth[row];
          indexMaxObjDepth = row;
        }
      }
    }
    return maxObjDepth;
  }

  inline const double* denseVector(const int index, const int numIndices) const {
    return this->getVector(index).denseVector(numIndices);
  }

  inline void denseMatrix(const int numIndices) const {
    for (int i = 0; i < this->getNumRows(); i++) {
      const double* currRow = this->getVector(i).denseVector(numIndices);
      for (int c = 0; c < numIndices; c++) {
        std::cout << currRow[c] << ",";
      }
      std::cout << std::endl;
      if (currRow) {
        delete currRow;
      }
    }
  }

  inline void clear() {
    CoinPackedMatrix::clear();
    RHS.clear();
    hplane.clear();
    hplaneIndexForSplit.clear();
    cutRay.clear();
    newRay.clear();
    facetIndex.clear();
    multiplicity.clear();
    hplaneIndexByAct.resize(0);
    allHplaneToAct.resize(0);
    indexOfHplaneInOverallOrderForSplit.resize(0);
    numPointsAdded.resize(0);
    numPointsRemoved.resize(0);
    avgAddedDepth.resize(0);
    avgRemovedDepth.resize(0);
    finalFlag.clear();
    numPoints = 0;
    numFinalPoints = 0;
    numRays = 0;
    numFinalRays = 0;
    SICDepth.clear();
    objDepth.clear();
    minSICDepth = 0.0;
    maxSICDepth = 0.0;
    totalSICDepth = 0.0;
    minObjDepth = 0.0;
    maxObjDepth = 0.0;
    totalObjDepth = 0.0;
    indexMinSICDepth = -1;
    indexMaxSICDepth = -1;
    indexMinObjDepth = -1;
    indexMaxObjDepth = -1;
  }
};
