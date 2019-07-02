#include "optionhelper.hpp"

/**********************************************************
 * Implement arg parsing using getopt and use CglGIC class
 *********************************************************/
// To prevent memory leaks from toLower
const std::vector<std::string> lowerCaseNames = lowerCaseStringVector(optionName);

//const std::string tmp[] = {toLowerString(optionName[optionIndex::UNKNOWN]), toLowerString(optionName[optionIndex::HELP])};
const option::Descriptor usage[] = {
    { optionIndex::UNKNOWN, 0, "", "", Arg::Unknown,
        "USAGE: example [options]\n\n"
            "Options:" },
    { optionIndex::HELP, 0, "h", 
        lowerCaseNames[optionIndex::HELP].c_str(),
        Arg::None,
        "  -h, --help \tPrint usage and exit." },
    // REQUIRED OPTIONS
    { optionIndex::UNKNOWN, 0, "", "", Arg::None, "\n  REQUIRED OPTIONS:"},
    { optionIndex::INST_NAME, 0, "i",
        lowerCaseNames[optionIndex::INST_NAME].c_str(),
        Arg::Required,
        "  -i <arg>, --inst_name=<arg>  \tInstance (/dir/name), including extension." },
    { optionIndex::OUT_DIR_NAME, 0, "o", 
        lowerCaseNames[optionIndex::OUT_DIR_NAME].c_str(),
        Arg::Required,
        "  -o <arg>, --out_dir_name=<arg> \tDirectory to which output is written." },
    // PARAMETERS AND OUTPUT FILE
    { optionIndex::UNKNOWN, 0, "", "", Arg::None, "\n  OTHER GENERAL OPTIONS:"},
    { optionIndex::LOG_FILE, 0, "",
        lowerCaseNames[optionIndex::LOG_FILE].c_str(),
         Arg::Optional,
         "  --log_file=<arg> \tWhere to compile instance information (usually across instances)." },
    { optionIndex::PARAM_FILE, 0, "p", 
        lowerCaseNames[optionIndex::PARAM_FILE].c_str(),
        Arg::Required,
        "  -p <arg>, --param_file=<arg> \tWhere to find the parameter file." },
		{ optionIndex::OPT_FILE, 0, "",
				lowerCaseNames[optionIndex::OPT_FILE].c_str(),
				Arg::Optional,
				"  --opt_file=<arg> \tWhere to find the file with optimal values." },
    // WHICH CUTTING PLANE TECHNIQUES TO USE
    { optionIndex::UNKNOWN, 0, "", "", Arg::None, "\n  CUTTING PLANE OPTIONS:"},
    { optionIndex::SICS, 0, "", 
        lowerCaseNames[optionIndex::SICS].c_str(),
        Arg::Numeric,
        "  --sics=<arg>"
            "\tNumber of rounds of standard intersection cuts to generate." },
    { optionIndex::PHA, 0, "",
        lowerCaseNames[optionIndex::PHA].c_str(),
        Arg::Binary,
        "  --pha=<arg>"
            "\tGenerate GICs by partial hyperplane activation?" },
    // COMMON CUTTING PLANE OPTIONS
    { optionIndex::UNKNOWN, 0, "", "", Arg::None, "\n  COMMON CUTTING PLANE OPTIONS:"},
    { optionIndex::NUM_RAYS_CUT, 0, "", 
        lowerCaseNames[optionIndex::NUM_RAYS_CUT].c_str(),
        Arg::Numeric,
        "  --num_rays_cut=<arg>"
            "\tNumber of rays to be cut. A value of "
            "0: cut as many as possible; -1: cut sqrt(n); -2: cut sqrt(n)/2." },
    { optionIndex::MAX_FRAC_VAR, 0, "", 
        lowerCaseNames[optionIndex::MAX_FRAC_VAR].c_str(),
        Arg::Numeric,
        "  --max-frac-core=<arg>"
            "\tNumber of fractional variables to consider. "
            "A value of "
            "0: cut as many as possible; "
            "-1: cut sqrt(n); "
            "-2: cut sqrt(n)/2. " },
    { optionIndex::CUT_LIMIT, 0, "", 
        lowerCaseNames[optionIndex::CUT_LIMIT].c_str(),
        Arg::Numeric,
        "  --cut_limit=<arg>"
            "\tLimit on number of PHA cuts generated. 0: unlimited; > 0: overall limit; < 0: limit set per split." },
    { optionIndex::ROUNDS, 0, "",
        lowerCaseNames[optionIndex::ROUNDS].c_str(),
        Arg::Numeric,
        "  --rounds=<arg>"
            "\tNumber of cutting plane rounds? "
            "(0: turn off all cuts, "
            "> 0: do this many rounds)" },
    { optionIndex::EPS, 0, "", 
        lowerCaseNames[optionIndex::EPS].c_str(),
        Arg::Numeric,
        "  --eps=<arg>"
            "\tValue of EPS." },
    { optionIndex::RAYEPS, 0, "", 
        lowerCaseNames[optionIndex::RAYEPS].c_str(),
        Arg::Numeric,
        "  --rayeps=<arg>"
            "\tValue of RAYEPS (zero threshold for ray components)." },
    { optionIndex::MIN_ORTHOGONALITY, 0, "",
        lowerCaseNames[optionIndex::MIN_ORTHOGONALITY].c_str(),
        Arg::Numeric,
        "  --min_orthogonality=<arg>"
            "\tValue of MIN_ORTHOGONALITY, min orthogonality to accept a cut." },
    { optionIndex::MIN_VIOL_ABS, 0, "", 
        lowerCaseNames[optionIndex::MIN_VIOL_ABS].c_str(),
        Arg::Numeric,
        "  --minviolabs=<arg>"
            "\tValue of MIN_VIOL_ABS, absolute minimum by which cut should remove the interior vertex." },
    { optionIndex::MIN_VIOL_REL, 0, "", 
        lowerCaseNames[optionIndex::MIN_VIOL_REL].c_str(),
        Arg::Numeric,
        "  --minviolrel=<arg>"
            "\tValue of MIN_VIOL_REL, relative (Euclidean distance) "
            "minimum by which cut should remove the interior vertex." },
    { optionIndex::TIMELIMIT, 0, "", 
        lowerCaseNames[optionIndex::TIMELIMIT].c_str(),
        Arg::Numeric,
        "  --timelimit=<arg>"
            "\tSuggestion of maximum of time that can be spent (but more complicated timing stuff happens)." },
    { optionIndex::TEMP, 0, "", 
        lowerCaseNames[optionIndex::TEMP].c_str(),
        Arg::Numeric,
        "  --temp=<arg>"
            "\tTemporary option that can be set for debugging or testing purposes." },
    // PHA OPTIONS
    { optionIndex::UNKNOWN, 0, "", "", Arg::None, "\n  PHA OPTIONS:"},
    { optionIndex::HPLANE_SCORING_FN, 0, "", 
        lowerCaseNames[optionIndex::HPLANE_SCORING_FN].c_str(),
        Arg::Required,
        "  --pha_act_option=<arg>"
            "\tHyperplane heuristic (-1: skip, 0: first, 1: best_avg_depth, 2: best_final_pts)." },
    { optionIndex::NUM_ALG2_ROUNDS, 0, "", 
        lowerCaseNames[optionIndex::NUM_ALG2_ROUNDS].c_str(),
        Arg::Numeric,
        "  --num_alg2_rounds=<arg>"
            "\tNumber hyperplanes to activate per split; negative means no algorithm 3 (tilting) done. " },
    // OBJECTIVE FUNCTION OPTIONS
    { optionIndex::UNKNOWN, 0, "", "", Arg::None, "\n  OBJECTIVE FUNCTION OPTIONS:"},
//    { optionIndex::USE_ALL_ONES_HEUR, 0, "",
//        lowerCaseNames[optionIndex::USE_ALL_ONES_HEUR].c_str(),
//        Arg::Binary,
//        "  --use_all_ones_heur=<arg>"
//            "\t(Objective function) Use the all ones objective?" },
    { optionIndex::USE_UNIT_VECTORS_HEUR, 0, "",
        lowerCaseNames[optionIndex::USE_UNIT_VECTORS_HEUR].c_str(),
        Arg::Numeric,
        "  --use_unit_vectors_heur=<arg>"
            "\t(Objective function) Should we try the axis directions?" },
    { optionIndex::USE_SPLIT_SHARE, 0, "", 
        lowerCaseNames[optionIndex::USE_SPLIT_SHARE].c_str(),
        Arg::Numeric,
        "  --use_split_share=<arg>"
            "\t(Objective function) Should the split share be used to generate cuts (may be slow): "
            "0 no, 1 yes first, 2 yes last." },
    { optionIndex::NUM_CUTS_ITER_BILINEAR, 0, "", 
        lowerCaseNames[optionIndex::NUM_CUTS_ITER_BILINEAR].c_str(),
        Arg::Numeric,
        "  --num_cuts_iter_bilinear=<arg>"
            "\t(Objective function) Number of cuts to try for the iterative bilinear cut generation." },
    { optionIndex::USE_CUT_VERT_HEUR, 0, "", 
        lowerCaseNames[optionIndex::USE_CUT_VERT_HEUR].c_str(),
        Arg::Binary,
        "  --use_cut_vert_heur=<arg>"
            "\t(Objective function) Use cut vertices cut heuristic?" },
    { optionIndex::USE_TIGHT_POINTS_HEUR, 0, "", 
        lowerCaseNames[optionIndex::USE_TIGHT_POINTS_HEUR].c_str(),
        Arg::Numeric,
        "  --use_tight_points_heur=<arg>"
            "\t(Objective function) How many rows to use in the tight on points and rays cut heuristic." },
    { optionIndex::CUT_PRESOLVE, 0, "", 
        lowerCaseNames[optionIndex::CUT_PRESOLVE].c_str(),
        Arg::Binary,
        "  --cut_presolve=<arg>"
            "\tUse presolve in cut generation (cutSolver)?" },
    { optionIndex::UNKNOWN, 0, "", "", Arg::None, "\nExamples:\n"
            "  example --unknown -- --this_is_no_option\n"
            "  example -unk --plus -ppp file1 file2\n" }, { 0, 0, 0, 0, 0, 0 }
}; /* usage */

void processArgs(int argc, char** argv, CglGICParam& param,
    std::string& tmppath, std::string& tmpname, std::string& tmpoutdir,
    std::string& opt_file, /*std::string& sol_file,*/ std::string& instdir) {
  // Process parameters
  argc -= (argc > 0);
  argv += (argc > 0); // skip program name argv[0] if present
  option::Stats stats(usage, argc, argv);
  option::Option* options = new option::Option[stats.options_max];
  option::Option* buffer = new option::Option[stats.buffer_max];
  option::Parser parse(usage, argc, argv, options, buffer);

  if (parse.error()) {
    fprintf(stderr, "*** ERROR: Error parsing parameters. Exiting.\n");
    delete[] options;
    delete[] buffer;
    exit(1);
  }

  if (options[optionIndex::HELP] || argc < 0) {
    int columns = getenv("COLUMNS") ? atoi(getenv("COLUMNS")) : 82;
    option::printUsage(fwrite, stdout, usage, columns);
    delete[] options;
    delete[] buffer;
    exit(0);
  }

  if (options[optionIndex::PARAM_FILE]) {
    std::string param_file = options[optionIndex::PARAM_FILE].arg;
    printf("## Setting parameters from file %s. ##\n", param_file.c_str());
    param.setParams(param_file.c_str());
  }

  if (argc < 2) {
    fprintf(stderr,
        "*** ERROR: Not enough arguments provided. Need at least 2 (inst name and out dir). Exiting.\n");
    delete[] options;
    delete[] buffer;
    exit(1);
  }

  for (int i = 0; i < parse.optionsCount(); ++i) {
    option::Option& opt = buffer[i];
    switch (opt.index()) {
      case optionIndex::HELP: // not possible, because handled further above and exits the program
      case optionIndex::PARAM_FILE: { // Handled above
        break;
      }
			case optionIndex::OPT_FILE: { // Decide where to get optimal solution value info
				opt_file = options[optionIndex::OPT_FILE].arg;
				break;
			}
      case optionIndex::INST_NAME: {
        tmppath = opt.arg;
        // Extract the overall directory as well as the instance directory and filename stub
        size_t found_slash = tmppath.find_last_of("/\\");
        if ((found_slash != std::string::npos)
            && (found_slash < tmppath.length())) {
          instdir = tmppath.substr(0, found_slash + 1);
          tmpname = tmppath.substr(found_slash + 1);
        } else { // Otherwise, no directory was specified, so the parent directory is where to save instance info
          tmpname = tmppath;
          instdir = "./";
        }
        break;
      }
      case optionIndex::OUT_DIR_NAME: {
        tmpoutdir = opt.arg;
        break;
      }
      case optionIndex::LOG_FILE: {
        GlobalVariables::LOG_FILE = opt.arg;
        break;
      }
      case optionIndex::SICS: {
        const int paramIndex = ParamIndices::SICS_PARAM_IND;
        if (opt.arg) {
          parseInt(opt.arg, param.paramVal[paramIndex]);
        } else {
          param.paramVal[paramIndex] = ParamDefaults[paramIndex];
        }
        if (!param.checkParam(paramIndex)) {
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        break;
      }
      case optionIndex::PHA: {
        const int paramIndex = ParamIndices::PHA_PARAM_IND;
        if (opt.arg) {
          parseInt(opt.arg, param.paramVal[paramIndex]);
        } else {
          param.paramVal[paramIndex] = ParamDefaults[paramIndex];
        }
        if (!param.checkParam(paramIndex)) {
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        break;
      }
      case optionIndex::HPLANE_SCORING_FN: {
        const int paramIndex = ParamIndices::HPLANE_SCORING_FN_PARAM_IND;
        parseInt(opt.arg, param.paramVal[paramIndex]);
        if (!param.checkParam(paramIndex)) {
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        break;
      }
      case optionIndex::NUM_ALG2_ROUNDS: {
        const int paramIndex = ParamIndices::NUM_ALG2_ROUNDS_PARAM_IND;
        parseInt(opt.arg, param.paramVal[paramIndex]);
        if (!param.checkParam(paramIndex)) {
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        break;
      }
      case optionIndex::NUM_RAYS_CUT: {
        const int paramIndex = ParamIndices::NUM_RAYS_CUT_PARAM_IND;
        parseInt(opt.arg, param.paramVal[paramIndex]);
        if (!param.checkParam(paramIndex)) {
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        break;
      }
      case optionIndex::MAX_FRAC_VAR: {
        const int paramIndex = ParamIndices::MAX_FRAC_VAR_PARAM_IND;
        parseInt(opt.arg, param.paramVal[paramIndex]);
        if (!param.checkParam(paramIndex)) {
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        break;
      }
      case optionIndex::CUT_LIMIT: {
        const int paramIndex = ParamIndices::CUT_LIMIT_PARAM_IND;
        parseInt(opt.arg, param.paramVal[paramIndex]);
        if (!param.checkParam(paramIndex)) {
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        break;
      }
      case optionIndex::USE_UNIT_VECTORS_HEUR: {
        const int paramIndex = ParamIndices::USE_UNIT_VECTORS_HEUR_PARAM_IND;
        parseInt(opt.arg, param.paramVal[paramIndex]);
        if (!param.checkParam(paramIndex)) {
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        break;
      }
      case optionIndex::USE_SPLIT_SHARE: {
        const int paramIndex = ParamIndices::USE_SPLIT_SHARE_PARAM_IND;
        parseInt(opt.arg, param.paramVal[paramIndex]);
        if (!param.checkParam(paramIndex)) {
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        break;
      }
      case optionIndex::NUM_CUTS_ITER_BILINEAR: {
        const int paramIndex = ParamIndices::NUM_CUTS_ITER_BILINEAR_PARAM_IND;
        parseInt(opt.arg, param.paramVal[paramIndex]);
        if (!param.checkParam(paramIndex)) {
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        break;
      }
      case optionIndex::USE_CUT_VERT_HEUR: {
        const int paramIndex = ParamIndices::USE_CUT_VERT_HEUR_PARAM_IND;
        parseInt(opt.arg, param.paramVal[paramIndex]);
        if (!param.checkParam(paramIndex)) {
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        break;
      }
      case optionIndex::USE_TIGHT_POINTS_HEUR: {
        const int paramIndex = ParamIndices::USE_TIGHT_POINTS_HEUR_PARAM_IND;
        parseInt(opt.arg, param.paramVal[paramIndex]);
        if (!param.checkParam(paramIndex)) {
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        break;
      }
      case optionIndex::CUT_PRESOLVE: {
        const int paramIndex = ParamIndices::CUT_PRESOLVE_PARAM_IND;
        parseInt(opt.arg, param.paramVal[paramIndex]);
        if (!param.checkParam(paramIndex)) {
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        break;
      }
      case optionIndex::ROUNDS: {
        const int paramIndex = ParamIndices::ROUNDS_PARAM_IND;
        parseInt(opt.arg, param.paramVal[paramIndex]);
        if (!param.checkParam(paramIndex)) {
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        break;
      }
      case optionIndex::EPS: {
        double tmpEps;
        parseDouble(opt.arg, tmpEps);
        if (tmpEps <= 0 || tmpEps > 1) {
          error_msg(errstr,
              "EPS needs to be a positive (small) value. Provided is %e.\n",
              tmpEps);
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        param.setEPS(tmpEps);
        break;
      }
      case optionIndex::RAYEPS: {
        double tmpEps;
        parseDouble(opt.arg, tmpEps);
        if (tmpEps <= 0 || tmpEps > 1) {
          error_msg(errstr,
              "RAYEPS needs to be a positive (small) value. Provided is %e.\n",
              tmpEps);
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        param.setRAYEPS(tmpEps);
        break;
      }
      case optionIndex::TIMELIMIT: {
        double tmpEps;
        parseDouble(opt.arg, tmpEps);
        if (tmpEps <= 0) {
          error_msg(errstr,
              "TOTAL_TIME needs to be a positive value. Provided is %e.\n",
              tmpEps);
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        param.setTIMELIMIT(tmpEps);
        break;
      }
      case optionIndex::MIN_ORTHOGONALITY: {
        double tmpVal;
        parseDouble(opt.arg, tmpVal);
        if (tmpVal < 0 || tmpVal > 1) {
          error_msg(errstr,
              "MIN_ORTHOGONALITY needs to be a value in [0,1]. Provided is %e.\n",
              tmpVal);
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        param.setMINORTHOGONALITY(tmpVal);
        break;
      }
      case optionIndex::MIN_VIOL_ABS: {
        double tmpVal;
        parseDouble(opt.arg, tmpVal);
        if (tmpVal <= 0 || tmpVal > 1) {
          error_msg(errstr,
              "MIN_VIOL_ABS needs to be a positive (small) value. Provided is %e.\n",
              tmpVal);
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        param.setMIN_VIOL_ABS(tmpVal);
        break;
      }
      case optionIndex::MIN_VIOL_REL: {
        double tmpVal;
        parseDouble(opt.arg, tmpVal);
        if (tmpVal <= 0 || tmpVal > 1) {
          error_msg(errstr,
              "MIN_VIOL_REL needs to be a positive (small) value. Provided is %e.\n",
              tmpVal);
          delete[] options;
          delete[] buffer;
          exit(1);
        }
        param.setMIN_VIOL_REL(tmpVal);
        break;
      }
      case optionIndex::TEMP: {
        double tmpVal;
        parseDouble(opt.arg, tmpVal);
        param.setTEMP(tmpVal);
        break;
      }
      case optionIndex::UNKNOWN: // not possible because Arg::Unknown returns ARG_ILLEGAL// which aborts the parse with an error
        break;
    }
  }

  delete[] options;
  delete[] buffer;
} /* processArgs */
