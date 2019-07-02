#pragma once

#include "CglGICParam.hpp"
#include "Utility.hpp"
#include "optionparser.h"

enum optionIndex { 
  UNKNOWN, HELP, 
  INST_NAME, OUT_DIR_NAME, LOG_FILE, PARAM_FILE, OPT_FILE,
  SICS, PHA,
  HPLANE_SCORING_FN, NUM_ALG2_ROUNDS,
  NUM_RAYS_CUT,
  MAX_FRAC_VAR, CUT_LIMIT,
  USE_SPLIT_SHARE, NUM_CUTS_ITER_BILINEAR,
  USE_UNIT_VECTORS_HEUR,
  USE_CUT_VERT_HEUR, USE_TIGHT_POINTS_HEUR, CUT_PRESOLVE,
  ROUNDS,
  EPS, RAYEPS, TIMELIMIT,
  MIN_ORTHOGONALITY, MIN_VIOL_ABS, MIN_VIOL_REL,
  TEMP,
  NUM_OPTIONS
};

const std::vector<std::string> optionName {
  "UNKNOWN", "HELP", 
    "INST_NAME", "OUT_DIR_NAME", "LOG_FILE", "PARAM_FILE", "OPT_FILE",
    "SICS", "PHA",
    "HPLANE_SCORING_FN", "NUM_ALG2_ROUNDS",
    "NUM_RAYS_CUT",
    "MAX_FRAC_VAR", "CUT_LIMIT",
    "USE_SPLIT_SHARE", "NUM_CUTS_ITER_BILINEAR",
    "USE_UNIT_VECTORS_HEUR",
    "USE_CUT_VERT_HEUR", "USE_TIGHT_POINTS_HEUR", "CUT_PRESOLVE",
    "ROUNDS",
    "EPS", "RAYEPS", "TIMELIMIT",
    "MIN_ORTHOGONALITY", "MIN_VIOL_ABS", "MIN_VIOL_REL",
    "TEMP",
    "NUM_OPTIONS"
};

inline const char* helpText(const std::string shortOpt, const int option_ind, const std::string description) {
  std::string str;
  if (!shortOpt.empty())
    str += "  -" + shortOpt;  
  str += "\t--";
  str += toLowerString(optionName[option_ind].c_str());
  if (!description.empty())
    str += "\t" + description;
  return str.c_str();
}

struct Arg: public option::Arg {
  static void printError(const char* msg1, const option::Option& opt, const char* msg2)
  {
    fprintf(stderr, "%s", msg1);
    fwrite(opt.name, opt.namelen, 1, stderr);
    fprintf(stderr, "%s", msg2);
  }

  static option::ArgStatus Unknown(const option::Option& option, bool msg)
  {
    if (msg) printError("Unknown option '", option, "'\n");
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus Required(const option::Option& option, bool msg)
  {
    if (option.arg != 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires an argument\n");
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus NonEmpty(const option::Option& option, bool msg)
  {
    if (option.arg != 0 && option.arg[0] != 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires a non-empty argument\n");
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus Numeric(const option::Option& option, bool msg)
  {
    double val;
    if (option.arg != 0 && parseDouble(option.arg, val))
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires a numeric argument\n");
    return option::ARG_ILLEGAL;
  }
  
  static option::ArgStatus NumericNonNeg(const option::Option& option, bool msg)
  {
    double val = -1.0;
    if (option.arg != 0 && parseDouble(option.arg, val) && val >= 0.0)
      return option::ARG_OK;
    
    if (msg) printError("Option '", option, "' requires a non-negative numeric argument\n");
    return option::ARG_ILLEGAL;
  }
  
  static option::ArgStatus Binary(const option::Option& option, bool msg)
  {
    int val;
    if (option.arg != 0 && parseInt(option.arg, val)) {
      if (val == 0 || val == 1) {
        return option::ARG_OK;
      }
    }
    if (msg) printError("Option '", option, "' requires a binary argument (0 or 1)\n");
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus BinaryOptional(const option::Option& option, bool msg)
  {
    if (option.arg != 0) {
      int val;
      if (parseInt(option.arg, val) 
          && (val == 0 || val == 1)) {
        return option::ARG_OK;
      }
      if (msg) printError("Option '", option, "' takes an optional binary argument (0 or 1)\n");
      return option::ARG_ILLEGAL;
    } else {
      return option::ARG_IGNORE;
    }
  }
}; /* struct Arg */

void processArgs(int argc, char** argv, CglGICParam& param,
    std::string& tmppath, std::string& tmpname, std::string& tmpoutdir,
    std::string& opt_file, /*std::string& sol_file,*/ std::string& instdir);
