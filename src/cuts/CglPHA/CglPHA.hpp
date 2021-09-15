#pragma once

#include "CglGICParam.hpp"
#include "CglCutGenerator.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinFactorization.hpp"
#include "hplaneActivation.hpp"

// (the AdvCut and AdvCuts classes could probably be removed eventually,
// but they have some useful features wrt OsiCuts)
#include "AdvCut.hpp"
#include "SolutionInfo.hpp"
#include "Hplane.hpp"

class CglPHA : public CglCutGenerator {

  friend void CglPHAUnitTest(const OsiSolverInterface * siP,
				const std::string mpdDir);
public:

  /**@name generateCuts */
  //@{
  /** Generate cuts for the model of the solver
      interface si.

      Insert the generated cuts into OsiCuts cs.
  */
  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
      const SolutionInfo& probData,
      const AdvCuts& structMSICs, const AdvCuts& NBMSICs);

  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo());

  /// For compatibility with CglCutGenerator (const method)
  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo()) const;

  /// Return true if needs optimal basis to do cuts (will return true)
  virtual bool needsOptimalBasis() const;
  //@}
  
  /**@name Public Methods */
  //@{

  // Set the parameters to the values of the given CglGMIParam object.
  void setParam(const CglGICParam &source); 
  // Return the CglGMIParam object of the generator. 
  inline CglGICParam getParam() const {return param;}
  inline int getNumObjTried() const {return num_obj_tried;}
  inline const OsiCuts getStructSICs() const {return structSICs;}

  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglPHA();

  // Constructor with specified parameters
  CglPHA(const CglGICParam&);

  /// Copy constructor 
  CglPHA(const CglPHA &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglPHA & operator=(const CglPHA& rhs);
  
  /// Destructor 
  virtual ~CglPHA();

  /// Copy our stuff
  virtual void copyOurStuff(const CglPHA* const rhs = NULL);

  /// Set SICs
  inline void setSICs(AdvCuts& structSICs, AdvCuts& NBSICs) {
    this->structSICs = structSICs;
    this->NBSICs = NBSICs;
  }

  //@}
    
private:
  
  // Private member methods

/**@name Private member methods */

  //@{
  // Method generating the cuts after all CglGIC members are properly set.
  void startPHAGeneration(AdvCuts &structMGICs, const SolutionInfo& probData,
      const AdvCuts &structSICs, const AdvCuts &NBSICs);
  void algorithm1(int& totalNumPointsRemoved, double& avgRemovedDepth,
      int& totalNumPointsAdded, double& avgAddedDepth,
      IntersectionInfo& interPtsAndRays, const int split_ind,
      const int hplane_ind, const int hplane_ind_in_split,
      const std::vector<Hplane>& hplaneStore, const Hplane& hplane,
      std::vector<Ray>& rayStore, const std::vector<Vertex>& vertexStore,
      const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
      const AdvCut* NBSIC, const double NBSICnorm, const AdvCut& objCut,
      const double objNorm, const int act_ind);
  void algorithm2(AdvCuts &structPHA, const SolutionInfo& solnInfo,
      const AdvCuts &structSICs, const AdvCuts &NBSICs,
      Stats& phaTimeStats, std::vector<int>& num_gics_per_split,
      AdvCut& NBObj, const double objNorm,
      std::vector<std::vector<double> >& distToSplitAlongRay,
      std::vector<Hplane>& hplaneStore, std::vector<Vertex>& vertexStore,
      std::vector<Ray>& rayStore, std::vector<int>& sortedRaysIndex,
//      std::vector<int>& allRayIndices,
      std::vector<int>& allRaysCutAcrossSplits,
      std::vector<std::vector<int> >& raysCutForSplit);
  void algorithm3(AdvCuts &structPHA, const SolutionInfo& solnInfo,
      const AdvCuts &structSICs, const AdvCuts &NBSICs,
      Stats& phaTimeStats, std::vector<int>& num_gics_per_split,
      AdvCut& NBObj, const double objNorm,
      std::vector<std::vector<double> >& distToSplitAlongRay,
      std::vector<Hplane>& hplaneStore, std::vector<Vertex>& vertexStore,
      std::vector<Ray>& rayStore, std::vector<int>& sortedRaysIndex,
//      std::vector<int>& allRayIndices,
      std::vector<int>& allRaysCutAcrossSplits,
      std::vector<std::vector<int> >& raysCutForSplit);
  void runAnalyses(const SolutionInfo* const solnInfoPtr,
      std::vector<Hplane>& hplaneStore, std::vector<int>& num_SIC_points,
      std::vector<int>& num_SIC_final_points, 
      std::vector<int>& num_SIC_final_rays, 
      std::vector<int>& num_SIC_rays,
      const AdvCuts& structSICs, const AdvCuts& NBSICs,
      const AdvCuts& structPHA);

  void PHAHelper(const int indexOfFirstRayToCut, const int numRaysToBeCut,
      const SolutionInfo& solnInfo, const AdvCuts& structSICs,
      const AdvCuts &NBSICs, const AdvCut& NBObj, const double objNorm,
      const std::vector<std::vector<double> >& distToSplitAlongRay,
      std::vector<Hplane>& hplaneStore, std::vector<Vertex>& vertexStore,
      std::vector<Ray>& rayStore, const std::vector<int>& sortIndex,
//      std::vector<int>& allRayIndices,
      std::vector<int> &allRaysCutAcrossSplits,
      std::vector<std::vector<int> > &raysCutForSplit,
      std::vector<IntersectionInfo>& interPtsAndRays, const int act_ind,
      const int total_num_act, const bool isAlg2,
      const HplaneScoringFn option, AdvCuts& structPHA,
      std::vector<int>& num_gics_per_split, int& num_total_pha,
      Stats& phaTimeStats);

  bool chooseHplaneByRay(//const std::vector<int>& allRayIndices,
      const int originatingVertexIndex, const Ray& rayOut, const int ray_ind,
      const PHASolverInterface* const solver,
      const SolutionInfo& solnInfo, std::vector<Hplane>& hplaneStore,
      std::vector<Vertex>& vertexStore, std::vector<Ray>& rayStore,
      const AdvCuts& NBSIC,
      const std::vector<std::vector<double> >& distToSplitAlongRay,
      std::vector<int>& allRaysCutAcrossSplits,
      std::vector<std::vector<int> >& raysCutForSplit,
      std::vector<IntersectionInfo>& interPtsAndRays, const int act_ind,
      const bool selectByDist, const bool selectByDepth,
      const bool selectByFinalPoints, const bool selectByMostRemoved);
  bool chooseHplaneBySplit(//const std::vector<int>& allRayIndices,
      const int originatingVertexIndex, const int numRaysToCut,
      const std::vector<int>& sortIndex,
      const PHASolverInterface* const solver,
      const SolutionInfo& solnInfo, std::vector<Hplane>& hplaneStore,
      std::vector<Vertex>& vertexStore, std::vector<Ray>& rayStore,
      const AdvCuts& NBSIC,
      const std::vector<std::vector<double> >& distToSplitAlongRay,
      std::vector<int>& allRaysCutAcrossSplits,
      std::vector<std::vector<int> >& raysCutForSplit,
      std::vector<IntersectionInfo>& interPtsAndRays, const int act_ind,
      const bool selectByDist, const bool selectByDepth,
      const bool selectByFinalPoints, const bool selectByMostRemoved);
  bool chooseHplaneHelper(const int split_ind, const double distToSplit,
      const bool isVertexOnSplit, const int ray_ind,
      std::vector<Hplane>& hplane, const Hplane& tmpHplane,
      const double tmpDistToHplane, double& minDistOverall,
      double& distToHplaneChosenForSplit,
      double& addedDepthByHplaneChosenForSplit,
      int& numFinalPointsByHplaneChosenForSplit,
      int& numPointsCutByHplaneChosenForSplit,
      const Vertex& vertCreated, double& addedDepthForThisHplane,
      int& numFinalPointsThisHplane, int& numPointsCutThisHplane,
      /*const std::vector<int>& allRayIndices,*/
      const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
      const IntersectionInfo& interPtsAndRays, const std::vector<Vertex>& vertexStore,
      const std::vector<Ray>& rayStore, const AdvCuts& NBSIC,
      const std::vector<std::vector<int> >& raysCutForSplit, const bool selectByDist,
      const bool selectByDepth, const bool selectByFinalPoints,
      const bool calcPointsRemoved, const bool isNBBound,
      const bool chooseByRay);
  void fakeActivateHplane(int& numFinalPointsAdded, const bool calcFinalPoints,
      double& avgAddedDepth, const bool calcDepth, /*const double distToHplane,*/
      const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
      const int split_ind, const std::vector<Vertex>& vertexStore,
      const Vertex& vertCreated, const std::vector<Ray>& rayStore,
      const int ray_ind, const Hplane& hplane, const AdvCut* NBSIC,
      /*const std::vector<int>& allRayIndices,*/ IntersectionInfo& interPtsAndRays,
      const IntersectionInfo& oldInterPtsAndRays,
      const std::vector<int>& raysPrevCut);

  void pointsFromRayViolatingHplane(int& numPointsCutByHplane,
      std::vector<int>& pointsToRemove,
      const int rayInd,
//      const std::vector<Vertex>& vertexStore, const std::vector<Ray>& rayStore,
      const Hplane& hplane, const IntersectionInfo& interPtsAndRays,
//      const int splitsAffByHplane,
      const PHASolverInterface* const solver, const SolutionInfo& solnInfo);

  bool addInterPtFromUncutRayOfC1(IntersectionInfo& interPtsAndRays,
      double& distToSplitAlongRay, const SolutionInfo& solnInfo,
      const int split_ind, const std::vector<Ray>& rayStore, const int ray_ind,
      const std::vector<Vertex>& vertexStore, const bool recalc_dist = true);

  void performActivations(const SolutionInfo& solnInfo,
      const std::vector<Vertex>& vertexStore, std::vector<Ray>& rayStore,
      const std::vector<Hplane>& hplaneStore,
      std::vector<IntersectionInfo>& interPtsAndRays, const AdvCuts& NBSIC,
      const AdvCut& objCut, const double objNorm, const int act_ind);

    void tryObjectives(AdvCuts& structPHA, std::vector<int>& num_cuts_per_cgs,
        int& num_total_gics, int& num_obj_tried, const SolutionInfo& solnInfo,
//        const int dim,
        const std::vector<IntersectionInfo>& interPtsAndRays,
//        const std::vector<Ray>& rayStore,
        const std::vector<Vertex>& vertexStore,
        const std::vector<Hplane>& hplaneStore, const AdvCuts& structSICs,
        Stats& phaTimeStats);
  void tryObjectivesForCgs(PHASolverInterface* cutSolver, AdvCuts& structPHA,
      int& num_cuts_per_cgs, int& num_cuts_generated, int& num_obj_tried,
      const SolutionInfo& probData, //const int dim,
      const int cgs_ind, const std::string& cgsName,
//      const IntersectionInfo& interPtsAndRays,
      const std::vector<Vertex>& vertexStore, const AdvCuts& structSICs);

  bool tieBreakingForChooseHplane(const bool selectByDist,
      const bool selectByDepth, const bool selectByFinalPoints,
      const double distToHplaneChosenForSplit,
      const double addedDepthByHplaneChosenForSplit,
      const int numFinalPointsByHplaneChosenForSplit,
      const double tmpDistToHplane,
      const double& addedDepthForThisHplane,
      const int numFinalPointsThisHplane) const;

	void storeNewRayNBIndices(int nNonBasicCols,
			std::vector<int> &newRayColIndices,
			const std::vector<int> &allRaysCutByHplane,
			std::vector<bool> &raysCutFlag) const;

	void generateSICs(AdvCuts& NBSICs, AdvCuts& structSICs,
			SolutionInfo& solnInfo);

  //@}
  
  // Private member data

/**@name Private member data */

  //@{

  /// Object with CglGICParam members. 
  CglGICParam param;

  /// Pointer on solver. Reset by each call to generateCuts().
  OsiSolverInterface *solver;

  /// Pointer to cast verson of solver
  PHASolverInterface* castsolver; // Currently has to be Clp because of some methods checking Clp status

  // Intersection points and rays
  std::vector<IntersectionInfo> interPtsAndRays;

  /// Pointer on point to separate. Reset by each call to generateCuts().
  const double *xlp;

  /// Pointer on matrix of coefficient ordered by rows. 
  /// Reset by each call to generateCuts().
  const CoinPackedMatrix *byRow;

  /// Pointer on matrix of coefficient ordered by columns. 
  /// Reset by each call to generateCuts().
  const CoinPackedMatrix *byCol;

  // AdvCuts versions of cuts we generate
  AdvCuts nb_advcs, struct_advcs;

  // If SICs are not provided, assume there are none
  AdvCuts structSICs, NBSICs;

  // For calculating beta scaling
  double min_nb_obj_val;

  // Number of objectives tried
  int num_obj_tried = 0;

  // Number of GICs generated
  int num_total_gics = 0;

  //@}
};

//#############################################################################
/** A function that tests the methods in the CglGIC class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void CglPHAUnitTest(const OsiSolverInterface * siP,
			 const std::string mpdDir );
