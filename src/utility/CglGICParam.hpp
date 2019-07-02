// Name:     CglGICParam.hpp
// Author:   Aleksandr M. Kazachkov 
//           based on CglGMIParam.hpp by Giacomo Nannicini
//           and on CglRedSplitParam.hpp by Francois Margot
// Date:     02/02/16
//-----------------------------------------------------------------------------
#pragma once

#include <string>
#include <vector>
#include "CglParam.hpp"
#include <limits>
#include <cmath> // sqrt

namespace CglGICParamNamespace {
  enum ParamIndices {
    SICS_PARAM_IND,
    PHA_PARAM_IND,
    HPLANE_SCORING_FN_PARAM_IND, // -1: skip, 0: first, 1: best, 2: best_final, 3: most_removed
    NUM_ALG2_ROUNDS_PARAM_IND, // # rounds of alg 2 (one hyperplane activated per split per round); if negative, skip alg 3 (tilting)
    NUM_RAYS_CUT_PARAM_IND, // 0 means as many as possible
    MAX_FRAC_VAR_PARAM_IND,
    CUT_LIMIT_PARAM_IND, // if > 0, replaces num_cuts_per_split parameter for determining max num cuts
    USE_SPLIT_SHARE_PARAM_IND, // should the split share be used to generate cuts (may be slow): 0 no, >0 yes first, <0 yes last (value = max # points to try)
    NUM_CUTS_ITER_BILINEAR_PARAM_IND, // 0: do not use the iterative bilinear cut generation, 1+: num cuts to try to generate
    USE_UNIT_VECTORS_HEUR_PARAM_IND, // 0: do not use, 1+: num to try, <0: abs(val) * sqrt(n)
    USE_CUT_VERT_HEUR_PARAM_IND, // 0: do not use the cut vertices cut selection heuristic
    USE_TIGHT_POINTS_HEUR_PARAM_IND, // 0: do not use the tight on points/rays cut heuristic, 1+: num points/rays to try
    CUT_PRESOLVE_PARAM_IND, // use presolve in cutSolver?
    ROUNDS_PARAM_IND, // number rounds
    NUM_PARAMS
  };

  const std::vector<std::string> ParamName {
      "SICS", "PHA",
      "HPLANE_SCORING_FN",
      "NUM_ALG2_ROUNDS",
      "NUM_RAYS_CUT", "MAX_FRAC_VAR", "CUT_LIMIT",
      "USE_SPLIT_SHARE", "NUM_CUTS_ITER_BILINEAR",
      "USE_UNIT_VECTORS_HEUR", "USE_CUT_VERT_HEUR", "USE_TIGHT_POINTS_HEUR",
      "CUT_PRESOLVE", "ROUNDS" };

  const std::vector<int> ParamDefaults {
      1,     //    SICS_DEFAULT
      0,     //    PHA_DEFAULT
      2,     //    HPLANE_SCORING_FN_DEFAULT
      0,     //    NUM_ALG2_ROUNDS_DEFAULT
      0,     //    NUM_RAYS_CUT_DEFAULT
      0,     //    MAX_FRAC_VAR_DEFAULT
      100,   //    CUT_LIMIT_DEFAULT
      1000,  //    USE_SPLIT_SHARE_DEFAULT
      1,     //    NUM_CUTS_ITER_BILINEAR_DEFAULT
      0,     //    USE_UNIT_VECTORS_HEUR_DEFAULT
      1,     //    USE_CUT_VERT_HEUR_DEFAULT
      0,     //    USE_TIGHT_POINTS_HEUR_DEFAULT
      1,     //    CUT_PRESOLVE_DEFAULT
      1,     //    ROUNDS_DEFAULT
  };

  enum class HplaneScoringFn {
    skip = -1, neighbor, depth, final, mostremoved,
  };

  const std::vector<std::string> HplaneScoringFnName { "neighbor", "depth",
      "final", "mostremoved" };
  const std::vector<std::string> HplaneScoringFnDescription {
      "activate the first hyperplane intersected for each ray",
      "activate the hyperplane leading to best average depth of added points",
      "activate the hyperplane leading to most final points",
      "activate the hyperplane removing the most points" };
} /* namespace CglGICParamNamespace */

using namespace CglGICParamNamespace;

class CglGICParam : public CglParam {
  public:
    std::vector<int> paramVal; // TODO Should eventually be protected

    /**@name Constructors and destructors */
    //@{
    /// Default constructor
    CglGICParam(
        const double inf = COIN_DBL_MAX,
        const double eps = DEFAULT_EPS,
        const double eps_coeff = DEFAULT_EPS_COEFF,
        const int max_supp_abs = DEFAULT_MAX_SUPPORT_ABS);

    /// Copy constructor
    CglGICParam(const CglGICParam& source);

    /// Clone
    virtual CglGICParam* clone() const;

    /// Assignment operator
    virtual CglGICParam& operator=(const CglGICParam &rhs);

    /// Destructor 
    virtual ~CglGICParam();
    //@}

    /**@name Static methods */
    //@{
    /** Check whether the parameter has a legal value */
    bool checkParam(const int paramIndex, const bool given_value = false, int value = -1) const;
    //@}

    /**@name Set/get methods */
    //@{
    /** To set parameter values */
    void setParams(const CglGICParam* const rhs = NULL);
    void setParams(const char* filename);
    bool finalizeParams(const int numNB, const int numSplits);

    /** Command line parameter values */
    virtual bool setParamVal(const int index, const int value);
    // name is assumed to have been standardized already
    virtual bool setParamVal(const std::string name, const double value);
    // name is assumed to have been standardized already
    int getParamIndex(const std::string name) const;
    inline std::string getParamName(const int index) const {return ParamName[index];}
    inline int getParamVal(const int index) const {return paramVal[index];}
    
    /** Aliases for parameter get/set method in the base class CglParam */
    /** Value for Infinity. Default: DBL_MAX */
    inline void setInf(const double value) {setINFINIT(value);}
    inline double getInf() const {return INFINIT;}
    inline void setInfinity(const double value) {setINFINIT(value);}
    inline double getInfinity() const {return INFINIT;}
    inline void setINFINITY(const double value) {setINFINIT(value);}
    inline double getINFINITY() const {return INFINIT;}

    /** Epsilon for comparing numbers. Default: 1.0e-6  */
    inline void setEps(const double value) {setEPS(value);}
    inline double getEps() const {return EPS;}

    /** Epsilon for zeroing out coefficients. Default: 1.0e-5 */
    inline void setEpsCoeff(const double value) {setEPS_COEFF(value);}
    inline double getEpsCoeff() const {return EPS_COEFF;}

    /** Maximum support of the cutting planes. Default: INT_MAX */
    inline void setMaxSupport(const int value) {setMAX_SUPPORT(value);}
    inline int getMaxSupport() const {return MAX_SUPPORT;}
    /** Alias for consistency with our naming scheme */
    inline void setMaxSupportAbs(const int value) {setMAX_SUPPORT(value);}
    inline int getMaxSupportAbs() const {return MAX_SUPPORT;}
    inline int getMAX_SUPPORT_ABS() const {return MAX_SUPPORT;}

    /** Set the value of RAYEPS. Default: 1e-7 */
    virtual bool setRAYEPS(const double value);
    inline double getRAYEPS() const {return RAYEPS;}
    /// Aliases
    inline bool set_rayeps(const bool value) {return setRAYEPS(value);}
    inline double get_rayeps() const {return RAYEPS;}
    inline bool setRayEps(const bool value) {return setRAYEPS(value);}
    inline double getRayEps() const {return RAYEPS;}

    /** Set the value of DIFFEPS. Default: 1e-3 */
    virtual bool setDIFFEPS(const double value);
    inline double getDIFFEPS() const {return DIFFEPS;}
    /// Aliases
    inline bool set_diffeps(const bool value) {return setDIFFEPS(value);}
    inline double get_diffeps() const {return DIFFEPS;}
    inline bool setDiffEps(const bool value) {return setDIFFEPS(value);}
    inline double getDiffEps() const {return DIFFEPS;}

    /** Set the value of AWAY. Default: 1e-3 */
    virtual bool setAWAY(const double value);
    inline double getAWAY() const {return AWAY;}
    /// Aliases
    inline bool set_away(const bool value) {return setAWAY(value);}
    inline double get_away() const {return AWAY;}
    inline bool setAway(const bool value) {return setAWAY(value);}
    inline double getAway() const {return AWAY;}

    /** Set the value of MINORTHOGONALITY. Default: 0.0 */
    virtual bool setMINORTHOGONALITY(const double value);
    inline double getMINORTHOGONALITY() const {return MIN_ORTHOGONALITY;}
    /// Aliases
    inline bool set_minorthogonality(const bool value) {return setMINORTHOGONALITY(value);}
    inline double get_minorthogonality() const {return MIN_ORTHOGONALITY;}
    inline bool setMinOrthogonality(const bool value) {return setMINORTHOGONALITY(value);}
    inline double getMinOrthogonality() const {return MIN_ORTHOGONALITY;}
    inline bool setMin_Orthogonality(const bool value) {return setMINORTHOGONALITY(value);}
    inline double getMin_Orthogonality() const {return MIN_ORTHOGONALITY;}
    inline bool setMIN_ORTHOGONALITY(const bool value) {return setMINORTHOGONALITY(value);}
    inline double getMIN_ORTHOGONALITY() const {return MIN_ORTHOGONALITY;}

    /** Set the value of TIMELIMT. Default: 3600 */
    virtual bool setTIMELIMIT(const double value);
    inline double getTIMELIMIT() const {return TIMELIMIT;}
    /// Aliases
    inline bool set_timelimit(const bool value) {return setTIMELIMIT(value);}
    inline double get_timelimit() const {return TIMELIMIT;}
    inline bool setTimeLimit(const bool value) {return setTIMELIMIT(value);}
    inline double getTimeLimit() const {return TIMELIMIT;}

    /** Set the value of CUTSOLVER_TIMELIMT. Default: 30 */
    virtual bool setCUTSOLVERTIMELIMIT(const double value);
    inline double getCUTSOLVERTIMELIMIT() const {return CUTSOLVER_TIMELIMIT;}
    /// Aliases
    inline bool set_cutsolvertimelimit(const bool value) {return setCUTSOLVERTIMELIMIT(value);}
    inline double get_cutsolvertimelimit() const {return CUTSOLVER_TIMELIMIT;}
    inline bool setCutSolverTimeLimit(const bool value) {return setCUTSOLVERTIMELIMIT(value);}
    inline double getCutSolverTimeLimit() const {return CUTSOLVER_TIMELIMIT;}
    inline bool setCutSolver_TimeLimit(const bool value) {return setCUTSOLVERTIMELIMIT(value);}
    inline double getCutSolver_TimeLimit() const {return CUTSOLVER_TIMELIMIT;}
    inline bool setCUTSOLVER_TIMELIMIT(const bool value) {return setCUTSOLVERTIMELIMIT(value);}
    inline double getCUTSOLVER_TIMELIMIT() const {return CUTSOLVER_TIMELIMIT;}

    /**
     * Set the maximum ratio between largest and smallest non zero 
     * coefficients in a cut. Default: 1e6.
     */
    virtual bool setMAXDYN(const double value);
    /** Get the value of MAXDYN */
    inline double getMAXDYN() const {return MAXDYN;}
    /// Aliases
    inline bool setMaxDyn(const double value) {return setMAXDYN(value);}
    inline double getMaxDyn() const {return MAXDYN;}

    /**
     * Set maximum number of cuts we should create 
     * (overrides num_act_rounds * num_cut_gen_sets)
     */
    virtual bool setMAXCUTS(const int value);
    inline double getMAXCUTS() const {return MAXCUTS;}
    inline bool setMaxCuts(const int value) {return setMAXCUTS(value);}
    inline double getMaxCuts() const {return getMAXCUTS();}

    /**
     * Set maximum number of iterations for the cut solver
     */
    virtual bool setCUTSOLVER_MAX_ITER(const int value);
    inline double getCUTSOLVER_MAX_ITER() const {return CUTSOLVER_MAX_ITER;}
    
    /**
     * Set maximum number of iterations of the bilinear program to try
     */
    virtual bool setITER_PER_CUT_BILINEAR(const int value);
    inline int getITER_PER_CUT_BILINEAR() const {return ITER_PER_CUT_BILINEAR;}

    /**
     * Set maximum number of pivots from optimum (for lifting)
     */
    virtual bool setMAX_PIVOTS(const int value);
    inline double getMAX_PIVOTS() const {return MAX_PIVOTS;}
    inline bool setMaxPivots(const int value) {return setMAX_PIVOTS(value);}
    inline double getMaxPivots() const {return getMAX_PIVOTS();}

    /**
     * Set the value of MIN_VIOL_ABS, the minimum violation for the current
     * basic solution in a generated cut. Default: 1e-7 
     **/
    virtual bool setMIN_VIOL_ABS(const double value);
    /** Get the value of MIN_VIOL_ABS */
    inline double getMIN_VIOL_ABS() const {return MIN_VIOL_ABS;}
    /// Aliases
    inline bool setMinViolAbs(const double value) {return setMIN_VIOL_ABS(value);}
    inline double getMinViolAbs() const {return MIN_VIOL_ABS;}

    /**
     * Set the value of MIN_VIOL_REL, the minimum violation for the current
     * basic solution in a generated cut. Default: 1e-7 
     **/
    virtual bool setMIN_VIOL_REL(const double value);
    /** Get the value of MIN_VIOL_REL */
    inline double getMIN_VIOL_REL() const {return MIN_VIOL_REL;}
    /// Aliases
    inline bool setMinViolRel(const double value) {return setMIN_VIOL_REL(value);}
    inline double getMinViolRel() const {return MIN_VIOL_REL;}

    /** Set the value of CHECK_DUPLICATES. Default: 0 */
    virtual void setCHECK_DUPLICATES(const bool value);
    /** Get the value of CHECK_DUPLICATES */
    inline bool getCHECK_DUPLICATES() const {return CHECK_DUPLICATES;}
    /// Aliases
    inline void setCheckDuplicates(bool value) {setCHECK_DUPLICATES(value);}
    inline bool getCheckDuplicates() const {return CHECK_DUPLICATES;}

    /** Set the value of EXIT_ON_INFEAS. Default: 0 */
    virtual void setEXIT_ON_INFEAS(const bool value);
    /** Get the value of EXIT_ON_INFEAS */
    inline bool getEXIT_ON_INFEAS() const {return EXIT_ON_INFEAS;}
    /// Aliases
    inline void setExitOnInfeas(bool value) {setEXIT_ON_INFEAS(value);}
    inline bool getExitOnInfeas() const {return EXIT_ON_INFEAS;}

    /** Set the value of TEMP. Default: 0 */
    virtual void setTEMP(const double value);
    /** Get the value of EXIT_ON_INFEAS */
    inline double getTEMP() const {return TEMP;}
    /// Aliases
    inline void setTemp(double value) {setTEMP(value);}
    inline double getTemp() const {return TEMP;}
    //@}

    void printParams(FILE* outfile = stdout, const char SEP = ' ') const;

    inline int getMaxNumCgs(const int numSplits) const {
      int max_num_cgs = numSplits;
      return max_num_cgs;
    }

  protected:
    /**@name Parameters */
    //@{
    // Instance parameters

    /** 
     * Value to determine whether a ray is parallel to a hyperplane
     * Default: 1e-7 
     */
    double RAYEPS;

    /** 
     * A more conservative value of epsilon, meaning larger epsilon
     * Default: 1e-3
     */
    double DIFFEPS;
    //double PIVOTEPS;
    double AWAY;
    double MIN_ORTHOGONALITY;

    /**
     * Maximum total time
     */
    double TIMELIMIT;
    double CUTSOLVER_TIMELIMIT;

    /** 
     * Maximum ratio between largest and smallest non zero 
     * coefficients in a cut. Default: 1e6.
     */
    double MAXDYN;

    /**
     * Maximum total number of cuts to create
     */
    int MAXCUTS;

    /**
     * Maximum number of iterations for the cut solver
     */
    int CUTSOLVER_MAX_ITER;

    /**
     * Maximum number of iterations of the bilinear cut generation program
     */
    int ITER_PER_CUT_BILINEAR;

    /**
     * Maximum number of pivots from the optimal solution
     */
    int MAX_PIVOTS;

    /** 
     * Minimum absolute violation for the current basic solution in a generated cut.
     * Default: 1e-7. 
     */
    double MIN_VIOL_ABS;

    /** 
     * Minimum relative violation for the current basic solution in a generated cut.
     * Default: 1e-7. 
     */
    double MIN_VIOL_REL;


    /** Check for duplicates when adding the cut to the collection? */
    bool CHECK_DUPLICATES;

    /** Exit if adding GICs make problem infeasible */
    bool EXIT_ON_INFEAS;

    /**
     * A temporary option to be used for debugging or test
     */
    double TEMP;
    //@}

    // Default values
    // These are defined in CglParam
    static constexpr double  DEFAULT_EPS = 1E-7;
    static constexpr double  DEFAULT_EPS_COEFF = 1E-7;
    static constexpr double  DEFAULT_INFINITY = std::numeric_limits<double>::max(); //1E20;
    static constexpr int     DEFAULT_MAX_SUPPORT_ABS = std::numeric_limits<int>::max();

    // The following are specific for GICs
    static constexpr double  DEFAULT_RAYEPS = 1E-7;
    static constexpr double  DEFAULT_DIFFEPS = 1E-3;
    static constexpr double  DEFAULT_AWAY = 1E-3;
    static constexpr double  DEFAULT_MIN_ORTHOGONALITY = 0;

    static constexpr double  DEFAULT_TIMELIMIT = 3600; //10800; // 3 hours in seconds
    static constexpr double  DEFAULT_CUTSOLVER_TIMELIMIT = 5;

    //static constexpr double  DEFAULT_EPS_PIVOT = 1E-1;
    static constexpr double  DEFAULT_MAXDYN = 1E6; // |alpha_max| / |alpha_min| upper bound
    static constexpr int     DEFAULT_MAXCUTS = 500;
    static constexpr int     DEFAULT_CUTSOLVER_MAX_ITER = 1E6;
    static constexpr int     DEFAULT_ITER_PER_CUT_BILINEAR = 1;
    static constexpr int     DEFAULT_MAX_PIVOTS = 0;

    static constexpr double  DEFAULT_MIN_VIOL_ABS = 1E-7;
    static constexpr double  DEFAULT_MIN_VIOL_REL = 1E-7;
    static constexpr bool    DEFAULT_CHECK_DUPLICATES = true;
    static constexpr bool    DEFAULT_EXIT_ON_INFEAS = false;

    static constexpr double  DEFAULT_TEMP = 0.;
};
