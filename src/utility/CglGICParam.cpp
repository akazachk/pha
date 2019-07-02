// Name:     CglGICParam.cpp
// Author:   Aleksandr M. Kazachkov 
//           based on CglGMIParam.hpp by Giacomo Nannicini
//           and on CglRedSplitParam.hpp by Francois Margot
// Date:     02/02/16
//-----------------------------------------------------------------------------

#include "CglGICParam.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include "Utility.hpp"
#include "GlobalConstants.hpp"

/***********************************************************************/
CglGICParam::CglGICParam(const double inf, const double eps,
    const double eps_coeff, const int max_supp_abs) :
  CglParam(inf, eps, eps_coeff, max_supp_abs) {
  setParams();
} /* constructor */

/***********************************************************************/
CglGICParam::CglGICParam(const CglGICParam& source) :
  CglParam(source) {
  setParams(&source);
} /* copy constructor from CglGICParam */

/***********************************************************************/
CglGICParam* CglGICParam::clone() const {
  return new CglGICParam(*this);
} /* clone */

/***********************************************************************/
CglGICParam& CglGICParam::operator=(const CglGICParam &rhs) {
  if(this != &rhs) {
    CglParam::operator=(rhs);
  }
  setParams(&rhs);
  return *this;
} /* assignment operator */

/***********************************************************************/
CglGICParam::~CglGICParam() {
  paramVal.resize(0);
} /* destructor */

/***********************************************************************/
void CglGICParam::setParams(const CglGICParam* const rhs) {
  // Initialize array
  if ((int) paramVal.size() < ParamIndices::NUM_PARAMS) {
    paramVal.resize(ParamIndices::NUM_PARAMS, 0);
  }
  if (rhs) {
    for (int i = 0; i < (int) CglGICParamNamespace::ParamName.size(); i++) {
      paramVal[i] = rhs->getParamVal(i);
    }
    RAYEPS = rhs->getRAYEPS();
    DIFFEPS = rhs->getDIFFEPS();
    AWAY = rhs->getAWAY();
    MIN_ORTHOGONALITY = rhs->getMINORTHOGONALITY();
    TIMELIMIT = rhs->getTIMELIMIT();
    CUTSOLVER_TIMELIMIT = rhs->getCUTSOLVER_TIMELIMIT();
    MAXDYN = rhs->getMAXDYN();
    MAXCUTS = rhs->getMAXCUTS();
    CUTSOLVER_MAX_ITER = rhs->getCUTSOLVER_MAX_ITER();
    ITER_PER_CUT_BILINEAR = rhs->getITER_PER_CUT_BILINEAR();
    MAX_PIVOTS = rhs->getMAX_PIVOTS();
    MIN_VIOL_ABS = rhs->getMIN_VIOL_ABS();
    MIN_VIOL_REL = rhs->getMIN_VIOL_REL();
    CHECK_DUPLICATES = rhs->getCHECK_DUPLICATES();
    EXIT_ON_INFEAS = rhs->getEXIT_ON_INFEAS();
    TEMP = rhs->getTEMP();
  } else {
    for (int i = 0; i < (int) CglGICParamNamespace::ParamName.size(); i++) {
      paramVal[i] = ParamDefaults[i];
    }
    RAYEPS = DEFAULT_RAYEPS;
    DIFFEPS = DEFAULT_DIFFEPS;
    AWAY = DEFAULT_AWAY;
    MIN_ORTHOGONALITY = DEFAULT_MIN_ORTHOGONALITY;
    TIMELIMIT = DEFAULT_TIMELIMIT;
    CUTSOLVER_TIMELIMIT = DEFAULT_CUTSOLVER_TIMELIMIT;
    MAXDYN = DEFAULT_MAXDYN;
    MAXCUTS = DEFAULT_MAXCUTS;
    CUTSOLVER_MAX_ITER = DEFAULT_CUTSOLVER_MAX_ITER;
    ITER_PER_CUT_BILINEAR = DEFAULT_ITER_PER_CUT_BILINEAR;
    MAX_PIVOTS = DEFAULT_MAX_PIVOTS;
    MIN_VIOL_ABS = DEFAULT_MIN_VIOL_ABS;
    MIN_VIOL_REL = DEFAULT_MIN_VIOL_REL;
    CHECK_DUPLICATES = DEFAULT_CHECK_DUPLICATES;
    EXIT_ON_INFEAS = DEFAULT_EXIT_ON_INFEAS;
    TEMP = DEFAULT_TEMP;
  }
} /* setParams */

/**********************************************************/
void CglGICParam::setParams(const char* filename) {
  if (!filename){
    setParams();
    return;
  }

  std::ifstream infile(filename);
  if (infile.is_open()) {
    std::string line;
    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      if (line.empty()) {
        continue;
      }
      std::string param_name;
      if (!(iss >> param_name)) {
        fprintf(stderr, "*** ERROR: Could not read parameter name. String is %s.\n", line.c_str());
        continue;
      }
      double val;
      if (!(iss >> val)) {
        fprintf(stderr, "*** ERROR: Could not read parameter value. String is %s.\n", line.c_str());
        continue;
      }

      setParamVal(toUpperStringNoUnderscore(param_name), val); // standardize name
    }
    infile.close();
  } else {
    // So far we are exiting disgracefully
    fprintf(stderr, "*** ERROR: Could not open parameter file %s.\n", filename);
    exit(1);
  }
} /* setParams (from file) */

/***********************************************************************/
bool CglGICParam::checkParam(const int paramIndex, const bool given_value, int value) const {
  // First we check if value is set to a dummy value
  if (!given_value) {
    value = paramVal[paramIndex];
  }
  if (paramIndex >= 0 && paramIndex < ParamIndices::NUM_PARAMS) {
    switch (paramIndex) { 
      case SICS_PARAM_IND:
        if (value < 0) {
          printf(
            "*** ERROR: "
            "SICs param needs to be a nonnegative integer. Provided is: %d.\n",
            value);
          return false;
        }
        break;
      case PHA_PARAM_IND:
        if (value < 0 || value > 1) {
          printf(
            "*** ERROR: "
            "Use PHA param needs to be false (0) or true (1). Provided is: %d.\n",
            value);
          return false;
        }
        break;
//      case NB_SPACE_PARAM_IND: // 0: Use NB space when generating cuts, 1: do not use it
//        if (value < 0 || value > 1) {
//          printf(
//            "*** ERROR: "
//            "Use NB space param needs to be false (0) or true (1). Provided is: %d.\n",
//            value);
//          return false;
//        }
//        break;
//      case PIVOT_PARAM_IND: // 0: Pivot to better solutions when generating cuts, 1: do not pivot
//        if (value < 0 || value > 1) {
//          printf(
//            "*** ERROR: "
//            "Pivot param needs to be false (0) or true (1). Provided is: %d.\n",
//            value);
//          return false;
//        }
//        break;
      case HPLANE_SCORING_FN_PARAM_IND: // -1: skip, 0: first, 1: best, 2: best_final
        if (value < -1 || value > 2) {
          printf(
            "*** ERROR: "
            "Hplane heuristic needs to be between -1 and 2. Provided is: %d.\n",
            value);
          return false;
        }
        break;
//      case NEXT_HPLANE_FINAL_PARAM_IND: // If generating more than 1 hplane cutting each ray, how to choose
//        if (value < 0 || value > 1) {
//          printf(
//            "*** ERROR: "
//            "Next hplane final flag needs to be false (0) or true (1). Provided is: %d.\n",
//            value);
//          return false;
//        }
//        break;
      case NUM_ALG2_ROUNDS_PARAM_IND: // <= 1 if next_hplane_final = false, else can be more; negative if no tilting round
//        if (value < 0) {
//          printf(
//            "*** ERROR: "
//            "Number of hyperplanes activated per ray needs to be at least 0. Provided is: %d.\n",
//            value);
//          return false;
//        }
        break;
      case NUM_RAYS_CUT_PARAM_IND: // -2 means sqrt(n)/2 rays, -1 means sqrt(n), 0 means as many as possible
        if (value < -2) {
          printf(
            "*** ERROR: "
            "Number of rays cut needs to be -2 or -1 (if we want sqrt(n)/2 or sqrt(n) rays cut, respectively) or non-negative. Provided is %d.\n",
            value);
          return false;
        }
        break;
      case MAX_FRAC_VAR_PARAM_IND:
        if (value < -2) {
          printf(
              "*** ERROR: "
                  "Number of variables cut needs to be -2 or -1 (if we want sqrt(n)/2 or sqrt(n) considered, respectively) or non-negative. Provided is %d.\n",
              value);
          return false;
        }
        break;
      case CUT_LIMIT_PARAM_IND: // limit on number of cuts; see reachedCutLimit in Utility
        /*if (value < 0) {
          printf(
            "*** ERROR: "
            "Maximum number of cuts needs to be at least 0. Provided is %d.\n",
            value);
          return false;
        }*/
        break;
      case USE_SPLIT_SHARE_PARAM_IND: // should the split share be used to generate cuts (may be slow): 0 no, > 0 yes first, < 0 yes last
        break;
      case NUM_CUTS_ITER_BILINEAR_PARAM_IND: // 0: do not use the iterative bilinear cut generation, 1+: num cuts to try to generate
        if (value < 0) {
          printf(
            "*** ERROR: "
            "Number of cuts for the iterative bilinear program needs to be non-negative. Provided is: %d.\n",
            value);
          return false;
        }
        break;
      case USE_UNIT_VECTORS_HEUR_PARAM_IND: // 0: do not use, 1+: num to try, <0: abs(val) * sqrt(n)
        break;
      case USE_CUT_VERT_HEUR_PARAM_IND: // 0: do not use the cut vertices cut selection heuristic
        if (value < 0 || value > 1) {
          printf(
            "*** ERROR: "
            "Cut vertices heuristic is either 1 (yes) or 0 (no). Provided is: %d.\n",
            value);
          return false;
        }
        break;
      case USE_TIGHT_POINTS_HEUR_PARAM_IND: // 0: do not use the tight on points/rays cut heuristic, 1+: num points/rays to try
        if (value < 0) {
          printf(
            "*** ERROR: "
            "Number of points/rays to try for tight on points/rays heuristic must be at least 0. Provided is: %d.\n",
            value);
          return false;
        }
        break;
      case CUT_PRESOLVE_PARAM_IND: // use presolve in cutSolver?
        if (value < 0 || value > 1) {
          printf(
            "*** ERROR: "
            "Whether to use cut presolve is a boolean. Provided is: %d.\n",
            value);
          return false;
        }
        break;
      case ROUNDS_PARAM_IND: // number rounds
        if (value < 0) {
          printf(
              "*** ERROR: "
                  "The number of rounds needs to be non-negative. Provided is: %d.\n",
              value);
          return false;
        }
        break;
    }
  } else {
    printf("*** WARNING: CglGICParam::checkParam(): index %d out of bounds (needs to be at most %d)\n",
        paramIndex, ParamIndices::NUM_PARAMS);
    return false;
  }
  return true;
}

/***********************************************************************/
bool CglGICParam::setParamVal(const int index, const int value) {
  if (index >= 0 && index < ParamIndices::NUM_PARAMS) {
    paramVal[index] = value;  // We do not check for this being in bounds here
    return true;
  } else {
    printf("*** WARNING: CglGICParam::setParam(): index %d out of bounds (needs to be at most %d)\n",
        index, ParamIndices::NUM_PARAMS);
    return false;
  }
}

/***********************************************************************/
/*
 * name is assumed to have been standardized already
 */
bool CglGICParam::setParamVal(const std::string name, const double value) {
  const int index = getParamIndex(name);
  if (index >= 0) {
    return setParamVal(index, int(value));
  } else if (name.compare("EPS") == 0) {
    setEPS(value);
    return true;
  } else if (name.compare("RAYEPS") == 0) {
    setRAYEPS(value);
    return true;
  }  else if (name.compare("AWAY") == 0) {
    setAWAY(value);
    return true;
  } else if (name.compare("MINORTHOGONALITY") == 0) {
    setMINORTHOGONALITY(value);
    return true;
  } else if (name.compare("MINVIOLABS") == 0) {
    setMIN_VIOL_ABS(value);
    return true;
  } else if (name.compare("MINVIOLREL") == 0) {
    setMIN_VIOL_REL(value);
    return true;
  } else if (name.compare("TIMELIMIT") == 0) {
    setTIMELIMIT(value);
    return true;
  } else {
    printf("*** WARNING: Could not find parameter %s.\n", name.c_str());
    return false;
  }
}

/***********************************************************************/
/*
 * name is assumed to have been standardized already
 */
int CglGICParam::getParamIndex(const std::string name) const {
  for (unsigned i = 0; i < CglGICParamNamespace::ParamName.size(); i++) {
    std::string str1 = toUpperStringNoUnderscore(
        CglGICParamNamespace::ParamName[i]);
    if (str1.compare(name) == 0) {
      return i;
    }
  }
  return -1;
}

/***********************************************************************/
bool CglGICParam::setRAYEPS(const double value) {
  if (value >= 0.0) {
    RAYEPS = value;
    return true;
  } 
  else {
    printf("*** WARNING: CglGICParam::setRAYEPS(): value: %f ignored (needs to be at least 0)\n",
        value);
    return false;
  }
} /* setRAYEPS */

/***********************************************************************/
bool CglGICParam::setDIFFEPS(const double value) {
  if (value >= 0.0) {
    DIFFEPS = value;
    return true;
  } 
  else {
    printf("*** WARNING: CglGICParam::setDIFFEPS(): value: %f ignored (needs to be at least 0)\n",
        value);
    return false;
  }
} /* setDIFFEPS */

/***********************************************************************/
bool CglGICParam::setAWAY(const double value) {
  if (value >= 0.0) {
    AWAY = value;
    return true;
  } else {
    printf(
        "*** WARNING: CglGICParam::setAWAY(): value: %f ignored (needs to be at least 0)\n",
        value);
    return false;
  }
} /* setAWAY */

/***********************************************************************/
bool CglGICParam::setMINORTHOGONALITY(const double value) {
  if (value >= 0.0 && value <= 1.0) {
    MIN_ORTHOGONALITY= value;
    return true;
  } else {
    printf(
        "*** WARNING: CglGICParam::setMINORTHOGONALITY(): value: %f ignored (needs to be at least 0 and at most 1)\n",
        value);
    return false;
  }
} /* setMINORTHOGONALITY */

/***********************************************************************/
bool CglGICParam::setTIMELIMIT(const double value) {
  if (value >= 0.0) {
    TIMELIMIT = value;
    return true;
  }
  else {
    printf("*** WARNING: CglGICParam::setTIMELIMIT(): value: %f ignored (needs to be at least 0)\n",
        value);
    return false;
  }
} /* setTIMELIMIT */

/***********************************************************************/
bool CglGICParam::setCUTSOLVERTIMELIMIT(const double value) {
  if (value >= 0.0) {
    CUTSOLVER_TIMELIMIT = value;
    return true;
  }
  else {
    printf("*** WARNING: CglGICParam::setCUTSOLVERTIMELIMIT(): value: %f ignored (needs to be at least 0)\n",
        value);
    return false;
  }
} /* setCUTSOLVERTIMELIMIT */

/***********************************************************************/
bool CglGICParam::setMAXDYN(const double value) {
  if (value >= 1.0) {
    MAXDYN = value;
    return true;
  } else {
    printf("*** WARNING: CglGICParam::setMAXDYN(): value: %f ignored (needs to be at least 1)\n",
        value);
    return false;
  }
} /* setMAXDYN */

/***********************************************************************/
bool CglGICParam::setMAXCUTS(const int value) {
  if (value >= 0.0) {
    MAXCUTS = value;
    return true;
  } else {
    printf("*** WARNING: CglGICParam::setMAXCUTS(): value: %d ignored (needs to be at least 0)\n",
        value);
    return false;
  }
} /* setMAXCUTS */

/***********************************************************************/
bool CglGICParam::setCUTSOLVER_MAX_ITER(const int value) {
  if (value >= 0.0) {
    CUTSOLVER_MAX_ITER = value;
    return true;
  } else {
    printf("*** WARNING: CglGICParam::setCUTSOLVER_MAX_ITER(): value: %d ignored (needs to be at least 0)\n",
        value);
    return false;
  }
} /* setCUTSOLVER_MAX_ITER */

/***********************************************************************/
bool CglGICParam::setITER_PER_CUT_BILINEAR(const int value) {
  if (value >= 0.0) {
    ITER_PER_CUT_BILINEAR = value;
    return true;
  } else {
    printf("*** WARNING: CglGICParam::setITER_PER_CUT_BILINEAR (): value: %d ignored (needs to be at least 0)\n",
        value);
    return false;
  }
} /* setCUT_SOLVER_MAX_ITER */

/***********************************************************************/
bool CglGICParam::setMAX_PIVOTS(const int value) {
  if (value >= 0.0) {
    MAXCUTS = value;
    return true;
  } else {
    printf("*** WARNING: CglGICParam::setMAXCUTS(): value: %d ignored (needs to be at least 0)\n",
        value);
    return false;
  }
} /* setMAXDYN */

/***********************************************************************/
bool CglGICParam::setMIN_VIOL_ABS(const double value) {
  if (value >= 0.0) {
    MIN_VIOL_ABS = value;
    return true;
  } 
  else {
    printf("*** WARNING: CglGICParam::setMIN_VIOL_ABS(): value: %f ignored (needs to be at least 0)\n",
        value);
    return false;
  }
} /* setMIN_VIOL_ABS */

/***********************************************************************/
bool CglGICParam::setMIN_VIOL_REL(const double value) {
  if (value >= 0.0) {
    MIN_VIOL_REL = value;
    return true;
  } 
  else {
    printf("*** WARNING: CglGICParam::setMIN_VIOL_REL(): value: %f ignored (needs to be at least 0)\n",
        value);
    return false;
  }
} /* setMIN_VIOL_REL */

/***********************************************************************/
void CglGICParam::setCHECK_DUPLICATES(const bool value) {
  CHECK_DUPLICATES = value;
} /* setCHECK_DUPLICATES */

/***********************************************************************/
void CglGICParam::setEXIT_ON_INFEAS(const bool value) {
  EXIT_ON_INFEAS = value;
} /* setEXIT_ON_INFEAS */

/***********************************************************************/
void CglGICParam::setTEMP(const double value) {
  TEMP = value;
} /* setTEMP */

/***********************************************************************/
// Check that the parameters are okay given the number of rays and splits
// If returns false, we will be exiting out
bool CglGICParam::finalizeParams(const int numNB, const int numSplits) {
  // Ensure parameters given to the problem are feasible
  if (paramVal[ParamIndices::NUM_RAYS_CUT_PARAM_IND] == 0) {
    paramVal[ParamIndices::NUM_RAYS_CUT_PARAM_IND] = numNB;
  }

  if (paramVal[ParamIndices::NUM_RAYS_CUT_PARAM_IND] < 0) {
    paramVal[ParamIndices::NUM_RAYS_CUT_PARAM_IND] = std::ceil(std::sqrt(numNB))
        / (-1 * paramVal[ParamIndices::NUM_RAYS_CUT_PARAM_IND]);
  }

  if (paramVal[ParamIndices::NUM_RAYS_CUT_PARAM_IND] > numNB) {
    error_msg(errorstring,
        "The maximum number of rays that can be cut is %d. The value specified of %d is too high. Changing to the maximum number of rays.\n",
        numNB,
        paramVal[ParamIndices::NUM_RAYS_CUT_PARAM_IND]);
    writeErrorToII(errorstring, GlobalVariables::log_file);
    paramVal[ParamIndices::NUM_RAYS_CUT_PARAM_IND] = numNB;
  }

  // If number of rays to be cut is STILL 0, there's a problem.
  // Note that there might be some strange things happening if
  // a variable is restricted to be integer, but it has a
  // bound that is fractional. We should clean the problem maybe
  // to ensure that we don't have non-basic integer-restricted
  // variables that are fractional.
  if (paramVal[ParamIndices::NUM_RAYS_CUT_PARAM_IND] == 0) {
    error_msg(errorstring,
        "The number of non-basic columns is %d. There are no rays to cut.\n",
        numNB);
    writeErrorToII(errorstring, GlobalVariables::log_file);
    return false;
  }

  // Check that the proper number of splits is available and the number of cut-generating sets is set correctly
  if (numSplits < 1) {
    error_msg(errorstring,
        "The number of splits is %d. No cut-generating sets are possible.\n",
        numSplits);
    writeErrorToII(errorstring, GlobalVariables::log_file);
    return false;
  }

  if (paramVal[ParamIndices::MAX_FRAC_VAR_PARAM_IND] == 0) {
    paramVal[ParamIndices::MAX_FRAC_VAR_PARAM_IND] = numSplits;
  }

  if (paramVal[ParamIndices::MAX_FRAC_VAR_PARAM_IND] < 0) {
    paramVal[ParamIndices::MAX_FRAC_VAR_PARAM_IND] = std::ceil(
        std::sqrt(numSplits))
        / (-1 * paramVal[ParamIndices::MAX_FRAC_VAR_PARAM_IND]);
    }

  if (paramVal[ParamIndices::MAX_FRAC_VAR_PARAM_IND] > numSplits) {
    warning_msg(warnstring,
        "Number of requested fractional vars is %d, but the maximum number we can choose is %d. Reducing the requested number of frac vars.\n",
        paramVal[ParamIndices::MAX_FRAC_VAR_PARAM_IND], numSplits);
    paramVal[ParamIndices::MAX_FRAC_VAR_PARAM_IND] = numSplits;
  }

  if (!paramVal[ParamIndices::SICS_PARAM_IND]
      && paramVal[ParamIndices::PHA_PARAM_IND]) {
    if (paramVal[ParamIndices::HPLANE_SCORING_FN_PARAM_IND] == static_cast<int>(HplaneScoringFn::depth)) {
      error_msg(errorstring, "To use average depth to select hyperplanes, need to generate SICs.\n");
      writeErrorToII(errorstring, GlobalVariables::log_file);
      return false;
    }
  }

  return true;
} /* finalizeParams */

/***********************************************************************/
void CglGICParam::printParams(FILE* outfile, const char SEP) const {
  for (int i = 0; i < ParamIndices::NUM_PARAMS; ++i) {
    fprintf(outfile, "%s%c%d\n", toLowerString(ParamName[i].c_str()).c_str(), SEP, paramVal[i]);
  }
  fprintf(outfile, "eps%c%e\n", SEP, EPS);
  fprintf(outfile, "rayeps%c%e\n", SEP, RAYEPS);
  fprintf(outfile, "away%c%e\n", SEP, AWAY);
  fprintf(outfile, "min_orthogonality%c%e\n", SEP, MIN_ORTHOGONALITY);
  fprintf(outfile, "min_viol_abs%c%e\n", SEP, MIN_VIOL_ABS);
  fprintf(outfile, "min_viol_rel%c%e\n", SEP, MIN_VIOL_REL);
  fprintf(outfile, "timelimit%c%e\n", SEP, TIMELIMIT);
  fflush(outfile);
}
