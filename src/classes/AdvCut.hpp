//============================================================================
// Name        : AdvCut.hpp
// Author      : Aleksandr M. Kazachkov
// Version     : 0.2015.03
// Description : Advanced cut class that extends OsiRowCut to allow storage of
//         other helpful data, such as the cut in NB space.
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#pragma once

#include "GlobalConstants.hpp" // For CutHeuristics
#include "SolutionInfo.hpp"
#include "typedefs.hpp"
#include "OsiRowCut.hpp"
#include "OsiCuts.hpp"

class AdvCut: public OsiRowCut {
public:
  double postCutObj;
  bool NBSpace; // True if the cut is stored in NB space

  int num_coeff;
  //  double rhs; // if NB, we complement so that NB space means var at zero
  //  std::vector<double> coeff;
  //  std::vector<double> comp_coeff; // comp_coeff is if we complement so that NB space means var at zero

  // If in NB space
  std::vector<double> distAlongRay;

  // This is the name and index of the split for which this cut was generated
  int cgsIndex;
  std::string cgsName;
  std::vector<int> splitVarIndex;

  // This is the cut heuristic that generated this cut
  enum CutHeuristics cutHeur;

  // Constructors
  AdvCut() :
      AdvCut(false, 1) {
  }
  ;
  AdvCut(bool NB, const int num_splits_involved);
  AdvCut(const AdvCut& rhs);

  inline std::string splitVarIndexString() const {
    std::string tmp;
    for (int i = 0; i < (int) splitVarIndex.size(); i++) {
      if (i > 0) {
        tmp += NAME_SEPARATOR;
      }
      tmp += std::to_string(splitVarIndex[i]);
    }
     return tmp;
  }

  // Operators and such
  inline size_t size() const {
    return num_coeff;
  }

  inline double operator[](int i) const {
    return row()[i];
  }

  /**
   * Overloads OsiColCut operator
   * Does not check simply whether cut components are equal,
   * but also if one cut is a scaled version of the other.
   */
  inline bool operator==(const AdvCut& cut2) const {
    return !cutsDifferent(*this, cut2);
  }

  /**
   * Overloads OsiColCut operator
   * Does not check simply whether cut components are equal,
   * but also if one cut is a scaled version of the other.
   */
  inline bool operator!=(const AdvCut& cut2) const {
    return !((*this) == cut2);
  }

//  /**
//   * Overloads OsiColCut operator
//   * Does not check simply whether cut components are equal,
//   * but also if one cut is a scaled version of the other.
//   */
//  inline bool addCuts(const AdvCut& cut2, const double mutliplier) {
//    this->row_ += cut2.row();
//  }

  // Assignment operator
  AdvCut & operator=(const AdvCut& rhs);

  double getTwoNorm() const;
  double getTwoNormSquared() const;

  void assign(const AdvCut & tmpCut);

  void setOsiRowCut(const std::vector<double>& coeff, const double rhs, const double EPS);

  void setOsiRowCut(const int num_coeff, const double* coeff, const double rhs, const double EPS);

  void setOsiRowCut(const OsiRowCut& tmpCut);

  double getActivity(const CoinPackedVector& vec) const;

  double getOrthogonality(const CoinPackedVector& vec) const;

  void setCutFromJSpaceCoefficients(const double* const coeff,
      const double beta, const PHASolverInterface* const solver,
      const std::vector<int> &nonBasicVarIndex, const bool inCompSpace,
      const double EPS);

  int cleanCut(const OsiSolverInterface* const solver,
          const SolutionInfo& probData /*, const double beta*/);
  void removeSmallCoefficients(const OsiSolverInterface* const solver,
      const double EPS = param.getEPS());
  bool badViolation(const OsiSolverInterface* const solver);
  bool badDynamism(const SolutionInfo& probData, const double SENSIBLE_MAX =
      param.getMAXDYN(), const double EPS = 1e-20); // TODO remove hard-coded eps

  static void convertCutFromJSpaceToStructSpace(AdvCut& structCut,
      const PHASolverInterface* const solver,
      const std::vector<int> &nonBasicVarIndex, const AdvCut& JspaceCut,
      const bool inCompSpace, const double EPS);

  static void convertCutCoeffFromJSpaceToStructSpaceNonPacked(
      std::vector<double>& structCut,
      const PHASolverInterface* const solver,
      const std::vector<int>& nonBasicVarIndex, const double* JspaceCut,
      const double JspaceRHS, const bool inCompSpace = true);

  static int cutsDifferent(const OsiRowCut &cut1, const OsiRowCut &cut2,
      const double eps = param.getDIFFEPS());

  static constexpr const char NAME_SEPARATOR = ';';
};

class AdvCuts: public OsiCuts {
public:
  std::vector<AdvCut> cuts;
  std::vector<double> score;
  std::vector<double> orthogonality;
  bool NBSpace;

  AdvCuts() :
      AdvCuts(false) {
  }
  ;
  AdvCuts(bool NB);
  // Copy constructor
  AdvCuts(const AdvCuts& rhs);

  const std::vector<AdvCut>& getCuts() const {
    return cuts;
  }

  void setCuts();

  inline void setCuts(const AdvCuts& tmpCuts) {
    OsiCuts::operator=(tmpCuts);
    this->NBSpace = tmpCuts.NBSpace;
    this->cuts = tmpCuts.cuts;
    this->score = tmpCuts.score;
    this->orthogonality = tmpCuts.orthogonality;
  }

  bool addNonDuplicateCut(const AdvCut & tmpCut, bool checkForDuplicate = true);
  bool isDuplicateCut(const AdvCut& tmpCut) const;
  int duplicateCutIndex(const AdvCut & tmpCut) const;
  int howDuplicate(const AdvCut& tmpCut, int& duplicateCutIndex,
      int& minOrthoIndex, double& minOrtho, const double eps) const;

  void setOsiCuts();

  AdvCuts & operator=(const AdvCuts& rhs);

  inline size_t size() const {
    return cuts.size();
  }

  inline AdvCut operator[](int i) const {
    return cuts[i];
  }

  inline void resize(const int n, const int num_splits_involved = 1) {
    cuts.resize(n);
    for (int i = 0; i < n; i++) {
      cuts[i].splitVarIndex.resize(num_splits_involved, -1);
    }
    score.resize(n, 0);
    orthogonality.resize(n, 1);
  }

  inline void reserve(int n) {
    cuts.reserve(n);
  }

  inline const AdvCut* getCutPointer(const int i) const {
    return (i < (int) cuts.size()) ? &(cuts[i]) : NULL;
  }

  void printCuts(std::string cut_type, const OsiSolverInterface* const solver,
      const SolutionInfo& probData) const;
  void printCuts(std::string cut_type, const OsiSolverInterface* const solver,
      const bool packedForm) const;

  inline void printCuts() const {
    OsiCuts::printCuts();
  }

  inline void printCuts(const int start, const int numCuts) const {
    for (int cut_ind = start; cut_ind < start + numCuts; cut_ind++) {
      this->rowCut(cut_ind).print();
    }
  }

  inline void printCuts(const int numCuts) const {
    printCuts(0, numCuts);
  }

  static void writeCutsInNBSpace(std::string headerText,
      const std::vector<AdvCut> &NBCut, const SolutionInfo &probData,
      const OsiSolverInterface* const solver, char* filename);
  static void writeCutsInStructSpace(std::string headerText,
      const std::vector<AdvCut> &NBCut, const SolutionInfo &probData,
      const OsiSolverInterface* const solver, char* filename);
  static void writeCutsInStructSpace(std::string headerText,
      const OsiCuts& structCut, const OsiSolverInterface* const solver,
      char* filename, const bool packedForm);
};
