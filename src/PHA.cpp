//============================================================================
// Name        : PHA.cpp
// Author      : Aleksandr M. Kazachkov (based on code by Selvaprabu Nadarajah)
// Version     : 2018.11
// Description : Partial hyperplane activation
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#include <string>
#include <ctime>
#include <chrono>

#include "optionhelper.hpp"

// COIN-OR
#include "CoinPackedMatrix.hpp"
#include "typedefs.hpp"
#include "OsiCuts.hpp"

// My files
#include "GlobalConstants.hpp"
#include "Utility.hpp"
#include "Output.hpp"

#include "SolutionInfo.hpp"
#include "Vertex.hpp"
#include "Point.hpp"
#include "Ray.hpp"
#include "Hplane.hpp"
#include "CglPHA.hpp"
#include "CglSIC.hpp"
#include "CutHelper.hpp"
#ifdef SHOULD_USE_GUROBI
#include "BBHelper.hpp"
#endif

#ifdef TRACE
#include "Debug.hpp"
#endif

// Catch abort signal if it ever gets sent
#include <csignal>

void signal_handler_with_error_msg(int signal_number) {
  /*Your code goes here. You can output debugging info.
   If you return from this function, and it was called
   because abort() was called, your program will exit or crash anyway
   (with a dialog box on Windows).
   */
  error_msg(errorstring, "Abort or seg fault message received. Signal number: %d.\n", signal_number);
  writeErrorToII(errorstring, GlobalVariables::log_file);
  exit(1);
} /* signal_handler_with_error_msg */

/*************/
/* VARIABLES */
/*************/
std::string instdir;
std::string f_name;
std::string in_file_ext;
std::string opt_file;
std::string FINAL_LOG_NAME;
std::string compare_cuts_output;
int num_obj_tried = 0;
CglGICParam saved_param; // this is to save parameters as they were inputted (because we change them throughout the code for various reasons)

//std::vector<double> opt_sol;

// For timing related names
time_t start_time_t, end_time_t;
struct tm* timeinfo;

/**********************************************************/
void introMessage() {
  printf("===\n");
  printf("%s\n", "Partial hyperplane activation, version 2018.11");
  printf("Start time: %s", asctime(timeinfo));
  printf("===\n\n");
  return;
} /* introMessage */

int init(int argc, char** argv, CglGICParam &param);
int wrapUp(int returncode);

/********************************************************************************/
/********************************************************************************/
/*******************************  MAIN  *****************************************/
/********************************************************************************/
/********************************************************************************/
/**
 * NOTE:   This is meant to be run on *NIX systems. For any other OS, the code
 *         may need to be modified before it will compile.
 */
int main(int argc, char **argv) {
  std::signal(SIGABRT, signal_handler_with_error_msg);
	std::signal(SIGSEGV, signal_handler_with_error_msg);

  /***********************************************************************************
   * Read inputs (if any) and set directory / file information
   ***********************************************************************************/
  GlobalVariables::start_time = std::chrono::system_clock::now();
  time(&start_time_t);
  timeinfo = localtime(&start_time_t);

  introMessage();
  init(argc, argv, param);
  saved_param = param;
	const int savedMaxFracVarParamVal = saved_param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND);

  /***********************************************************************************
   * Initialize solver and initially solve
   ***********************************************************************************/
  GlobalVariables::timeStats.start_timer(GlobalConstants::TOTAL_TIME);

  GlobalVariables::timeStats.start_timer(INIT_SOLVE_TIME);

  PHASolverInterface* solver = new PHASolverInterface;
  setMessageHandler(solver);

  if (in_file_ext.compare("lp") == 0) {
    // Read LP file
#ifdef TRACE
    printf("\n## Reading LP file. ##\n");
#endif
    solver->readLp(f_name.c_str());
  } else {
    if (in_file_ext.compare("mps") == 0) {
      // Read MPS file
#ifdef TRACE
      printf("\n## Reading MPS file. ##\n");
#endif
      solver->readMps(f_name.c_str());
    } else {
      error_msg(errorstring, "Unrecognized extension: %s.\n",
          in_file_ext.c_str());
      writeErrorToII(errorstring, GlobalVariables::log_file);
      exit(1);
    }
  }

  // Calculate objective value and compare against stored one (if it exists)
//  if (opt_sol.size() > 0) {
//    if ((int) opt_sol.size() != solver->getNumCols()) {
//      error_msg(errorstring, "Incorrect number of variables read from solution file: cols: %d, vars read: %d\n",
//          solver->getNumCols(), (int) opt_sol.size());
//      exit(1);
//    }
//
//    const double opt_val = dotProduct(opt_sol.data(), solver->getObjCoefficients(), solver->getNumCols());
//    if (!isInfinity(std::abs(GlobalVariables::bestObjValue))) {
//      if (!isVal(opt_val, GlobalVariables::bestObjValue)) {
//        error_msg(errorstring, "Optimal IP values do not match. From sol: %e. From opt file: %e.\n",
//          opt_val, GlobalVariables::bestObjValue);
//        exit(1);
//      }
//    }
//    GlobalVariables::bestObjValue = opt_val; // if wishing to test feasible instead of optimal solutions, comment out this line and error above
//    printf("Optimal solution read in and value set to %e.\n", opt_val);
//  } /* calculate IP obj value from solution when it exists */

  // Make sure we are doing a minimization problem; this is just to make later
  // comparisons simpler (i.e., a higher LP obj after adding the cut is better).
  if (solver->getObjSense() < param.getEPS()) {
    printf(
        "\n## Detected maximization problem. Negating objective function to make it minimization. ##\n");
    solver->setObjSense(1.0);
    const double* obj = solver->getObjCoefficients();
    for (int col = 0; col < solver->getNumCols(); col++) {
      solver->setObjCoeff(col, -1. * obj[col]);
    }
  }

#ifdef TRACE
  printf("\n## Performing initial solve. ##\n");
#endif

  solver->initialSolve();

#ifdef TRACE
  printf("Saving LP relaxation solution information.\n");
#endif
  // We'll be using the ray info, so we should not use the brief SolutionInfo
  SolutionInfo solnInfoOrig(solver,
      param.paramVal[ParamIndices::MAX_FRAC_VAR_PARAM_IND], false);
#ifdef SHOULD_WRITE_BRIEF_SOLN
  solnInfoOrig.printBriefSolutionInfo(solver, GlobalVariables::out_f_name_stub);
#endif

#ifdef SHOULD_WRITE_SOLN
  writeSoln(solver, GlobalVariables::out_f_name_stub);
#endif

  GlobalVariables::timeStats.end_timer(INIT_SOLVE_TIME);

  /***********************************************************************************
   * Finalize parameters
   ***********************************************************************************/
#ifdef TRACE
  printf("\n## Finalizing parameters. ##\n");
#endif

  if (solnInfoOrig.numSplits == 0) {
    error_msg(errorstring, "solnInfoOrig.numSplits: %d. Should not be zero.\n",
        solnInfoOrig.numSplits);
    writeErrorToII(errorstring, GlobalVariables::log_file);
    exit(1);
  }

  if (!param.finalizeParams(solnInfoOrig.numNB, solnInfoOrig.numSplits)) {
    error_msg(errstr,
        "Issue finalizing parameters (no rays to cut or issue with parameter combinations).\n");
    writeErrorToII(errstr, GlobalVariables::log_file);
    exit(1);
  }

  AdvCuts structSICs(false), structPHA(false);
  AdvCuts NBSICs(true);

  if (param.getParamVal(ParamIndices::SICS_PARAM_IND) >= 1) {
    /***********************************************************************************
     * Generate SICs
     ***********************************************************************************/
#ifdef TRACE
    printf("\n## Generating SICs for the main problem. ##\n");
#endif

    GlobalVariables::timeStats.start_timer(GEN_MSICS_TIME);
    if (param.getParamVal(ParamIndices::SICS_PARAM_IND) > 1) {
      printf("\n############### Starting round 1 of SIC generation. ###############\n");
    }
    // Only feasible ones are returned
    CglSIC SICGen(param);
    SICGen.generateCuts(*solver, NBSICs, structSICs, solnInfoOrig);

		// Check if we should do more rounds
    PHASolverInterface* tmpSolver = NULL;
    if (param.getParamVal(ParamIndices::SICS_PARAM_IND) > 1) {
      tmpSolver = dynamic_cast<PHASolverInterface*>(solver->clone());
			int oldNumSICs = 0;
			for (int round_ind = 1;
					round_ind < param.getParamVal(ParamIndices::SICS_PARAM_IND);
					round_ind++) {
				printf("\n############### Starting round %d of SIC generation. ###############\n", round_ind + 1);
				applyCuts(tmpSolver, structSICs, oldNumSICs, -1);
				checkSolverOptimality(tmpSolver, false);

				oldNumSICs = structSICs.sizeCuts();

				OsiCuts currStructSICs;
				CglGICParam tmpParam(saved_param);
				tmpParam.setParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND, 0); // will cause number cuts added per round to not be the same
				CglSIC CurrSICGen(tmpParam);
				CurrSICGen.generateCuts(*tmpSolver, currStructSICs);

				// Check if we should stop early
				if (currStructSICs.sizeCuts() == 0) {
					break;
				}

				structSICs.insert(currStructSICs);
			} /* rounds of cuts */
			if (tmpSolver) {
				delete tmpSolver;
			}
		} /* check if we should do rounds of SICs */

		GlobalVariables::numCutsFromHeur[CutHeuristics::SIC_CUT_GEN] = structSICs.sizeCuts();
		GlobalVariables::timeStats.end_timer(GEN_MSICS_TIME);

#ifdef TRACE
    printf("***** SICs generated: %d.\n", structSICs.sizeCuts());
#endif

#ifdef SHOULD_WRITE_SICS
    NBSICs.printCuts("SICs", solver, solnInfoOrig);
    structSICs.printCuts("SICs", solver, solnInfoOrig);
#endif

    if (structSICs.sizeCuts() == 0) {
      error_msg(errorstring, "Number of SICs is 0.\n");
      writeErrorToII(errorstring, GlobalVariables::log_file);
      exit(1);
    }
  } /* generate SICs? */

  /***********************************************************************************
   * Check parameters once more after SIC generation
   ***********************************************************************************/
  if (param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND)
      != solnInfoOrig.numFeasSplits) {
    error_msg(errorstring,
        "Number of *feasible* splits %d is not equal to the number of requested splits %d, "
        "and this is not in subspace, so should not happen unless problem is infeasible.\n",
        solnInfoOrig.numFeasSplits,
        param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND));
    writeErrorToII(errorstring, GlobalVariables::log_file);
    exit(1);
  }

#ifdef TRACE
  // Count intersecting rays for the chosen splits
  solnInfoOrig.countIntersectingRays(solnInfoOrig.feasSplitVar);
#endif

#if defined TRACE && defined SHOULD_WRITE_SPLIT_FEAS_INFO
  printSplitInfo(solver, solnInfoOrig);
#endif

  /***********************************************************************************
   * Get range of objective values attainable from any one split
   ***********************************************************************************/
  {
    // The below is so that this check does not mess anything up
    PHASolverInterface* copySolver = dynamic_cast<PHASolverInterface*>(solver->clone());
    calculateStrongBranchingLB(GlobalVariables::worstSplitObjValue,
        GlobalVariables::bestSplitObjValue, copySolver, solnInfoOrig);
    if (copySolver) {
      delete copySolver;
    }
#ifdef TRACE
    printf(
        "\n## Best bound attainable from any one split ranges from %.2f to %.2f. ##\n",
        GlobalVariables::worstSplitObjValue,
        GlobalVariables::bestSplitObjValue);
#endif
  }

  // Having finalized everything, write parameter information to file
  writeParamInfoToII(saved_param, GlobalVariables::log_file);
  writeSolnInfoToII(solnInfoOrig, GlobalVariables::log_file);

  /***********************************************************************************
   * Generate GICs
   ***********************************************************************************/
  if (param.getParamVal(ParamIndices::PHA_PARAM_IND)) {
    GlobalVariables::timeStats.start_timer(GEN_PHA_TIME);
    if (param.getParamVal(ParamIndices::ROUNDS_PARAM_IND) > 1) {
      printf("\n############### Starting round 1 of PHA cut generation. ###############\n");
    }
    CglPHA GICGen(param);
    GICGen.generateCuts(*solver, structPHA, solnInfoOrig, structSICs, NBSICs);
    num_obj_tried += GICGen.getNumObjTried();

		// Check if we should do more rounds
    PHASolverInterface* tmpSolver = NULL;
    if (param.getParamVal(ParamIndices::ROUNDS_PARAM_IND) > 1) {
      tmpSolver = dynamic_cast<PHASolverInterface*>(solver->clone());
    }
		int oldNumGICs = 0;
    for (int round_ind = 1;
        round_ind < param.getParamVal(ParamIndices::ROUNDS_PARAM_IND);
        round_ind++) {
			printf("\n############### Starting round %d of PHA cut generation. ###############\n", round_ind + 1);
			applyCuts(tmpSolver, structPHA, oldNumGICs, -1);
			checkSolverOptimality(tmpSolver, false);

      oldNumGICs = structPHA.sizeCuts();
			
      // Generate new cuts but not with new SICs
      //const int savedSICParamVal = saved_param.getParamVal(ParamIndices::SICS_PARAM_IND);
      //saved_param.setParamVal(ParamIndices::SICS_PARAM_IND, 0);

      CglGICParam tmpParam(saved_param);
      CglPHA CurrPHAGen(tmpParam);
      CurrPHAGen.setSICs(structSICs, NBSICs);
      CurrPHAGen.generateCuts(*tmpSolver, structPHA);
      num_obj_tried += CurrPHAGen.getNumObjTried();
      saved_param.setParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND, savedMaxFracVarParamVal);

      // Check if we should stop early
      if ((int) structPHA.size() == oldNumGICs) {
        break;
      }

      // Add new cuts to the total collection
      //for (int i = 0; i < (int) currCuts.size(); i++) {
      //  structPHA.addNonDuplicateCut(currCuts.cuts[i], false);
      //}
      //structPHA.setOsiCuts();
    } /* rounds of cuts */
    if (tmpSolver) {
      delete tmpSolver;
    }
    GlobalVariables::timeStats.end_timer(GEN_PHA_TIME);
  } /* check if pha flag is set to true */

  /***********************************************************************************
   * Compare cuts
   ***********************************************************************************/
  printf("\n## Finished cut generation.");
  if (structSICs.sizeCuts() > 0) {
    printf("  SIC: %d.", structSICs.sizeCuts());
  }
  if (structPHA.sizeCuts() > 0) {
    printf("  PHA: %d.", structPHA.sizeCuts());
  }
  printf(" ##\n");
  const int numCutsToAddPerRound =
//      (structSICs.sizeCuts() > 0) ? structSICs.sizeCuts() :
          param.getParamVal(ParamIndices::MAX_FRAC_VAR_PARAM_IND);
  compareCuts(num_obj_tried, solver, structSICs, structPHA,
      numCutsToAddPerRound, compare_cuts_output);

#ifdef SHOULD_USE_GUROBI
  /***********************************************************************************
   * Check that cuts do not remove IP opt
   ***********************************************************************************/
  if (param.getTEMP() == 1) {
    int num_errors = 0;
    double max_diff = 0.;
    printf("\n## Checking validity of cuts wrt integer-optimal solution. ##\n");
    // Get the original solution
    // First try to find it in the solution folder
    const std::string solution_folder;
    BBInfo tmp_bb_info;
    std::vector<double> solution;
    doBranchAndBoundWithGurobi(f_name.c_str(), tmp_bb_info, &solution);

    // Check cuts
    for (int cut_ind = 0; cut_ind < structPHA.sizeCuts(); cut_ind++) {
      const double rhs = structPHA.cuts[cut_ind].rhs();
      const int num_el = structPHA.cuts[cut_ind].row().getNumElements();
      const int* ind = structPHA.cuts[cut_ind].row().getIndices();
      const double* el = structPHA.cuts[cut_ind].row().getElements();
      const double activity = dotProduct(num_el, ind, el, solution.data());

      if (lessThanVal(activity, rhs)) {
        warning_msg(warnstring, "Cut %d removes optimal solution. Activity: %.10f. Rhs: %.10f.\n", cut_ind, activity, rhs);
        num_errors++;
        max_diff = CoinMax(max_diff, rhs - activity);
      }
    }
    printf("Number errors: %d.\tMax diff: %e.\n", num_errors, max_diff);
    fprintf(GlobalVariables::log_file, "%d,", num_errors);
    fprintf(GlobalVariables::log_file, "%.10e,", max_diff);
  } // checking cuts for violating the IP opt 
#endif

  /***********************************************************************************
   * Wrap up
   ***********************************************************************************/
  if (solver) {
    delete solver;
  }
  return wrapUp(0);
} /* main */

/**********************************************************/
/**
 * Initialization helper
 */
int init(int argc, char** argv, CglGICParam &param) {
  // These will prove useful
  std::string tmppath, tmpname, tmpoutdir;

  // Process the input arguments
  processArgs(argc, argv, param, tmppath, tmpname, tmpoutdir, opt_file, /*sol_file,*/ instdir);

  // FINISHED READING INPUT
  // NOW WE PROCESS IT
  GlobalVariables::initialize();

  // If nothing entered, assume the out directory is the same place as the instance location
  if (tmpoutdir.empty()) {
    tmpoutdir = instdir;
  }

  // Check to make sure that instdir has a '/' on the end
  if (!instdir.empty() && instdir.at(instdir.size() - 1) != '/')
    instdir += '/';

  // Check to make sure that tmpoutdir has a "/" on the end
  if (!tmpoutdir.empty() && tmpoutdir.at(tmpoutdir.size() - 1) != '/')
    tmpoutdir += '/';

  const std::string TEST_NAME = tmpname;
  tmpname.clear();

  // Open instances info output file
  const size_t found_slash_ii = GlobalVariables::LOG_FILE.find_last_of(
      "/\\");
  const size_t found_dot_ii = GlobalVariables::LOG_FILE.find_last_of('.');
  const std::string FINAL_II_OUT_DIR =
      (found_slash_ii != std::string::npos)
          && (found_slash_ii < GlobalVariables::LOG_FILE.length()) ?
          "" : tmpoutdir;
  const std::string FINAL_II_OUT_EXT =
      (found_dot_ii != std::string::npos)
          && (found_dot_ii < GlobalVariables::LOG_FILE.length())
          && ((found_slash_ii == std::string::npos)
              || (found_slash_ii < found_dot_ii)) ?
          "" : II_OUT_EXT; // The last condition is to ensure we have not inputted something like "../stub"
  FINAL_LOG_NAME = (FINAL_II_OUT_DIR
      + GlobalVariables::LOG_FILE + FINAL_II_OUT_EXT);
  bool ii_exists = fexists(FINAL_LOG_NAME.c_str());
  GlobalVariables::log_file = fopen(FINAL_LOG_NAME.c_str(), "a");
  if (GlobalVariables::log_file == NULL) {
    error_msg(errstr, "Cannot open file %s to write instance info.\n",
        FINAL_LOG_NAME.c_str());
    exit(1);
  }

  // Get file name stub
  size_t found_dot = TEST_NAME.find_last_of('.');

  // Put string after last '.' into string in_file_ext
  if (found_dot >= TEST_NAME.length()) { // Previously checked found_dot >= 0
    std::cerr
        << "*** ERROR: Cannot find the file extension (no '.' in input file name).\n"
        << std::endl;
    exit(1);
  }

  GlobalVariables::prob_name = TEST_NAME.substr(0, found_dot);

  std::string tmp_file_ext(TEST_NAME.substr(found_dot + 1));
  unsigned found_dot_tmp = found_dot;

  if (tmp_file_ext.compare("gz") == 0) {
    found_dot_tmp = GlobalVariables::prob_name.find_last_of('.');

    // Put string after last '.' into string in_file_ext
    if (found_dot_tmp >= GlobalVariables::prob_name.length()) {
      std::cerr
          << "*** ERROR: Other than gz, cannot find the file extension (no '.' in input file name).\n"
          << std::endl;
      exit(1);
    }

    tmp_file_ext = GlobalVariables::prob_name.substr(found_dot_tmp + 1);
    GlobalVariables::prob_name = GlobalVariables::prob_name.substr(0,
        found_dot_tmp);
  }

  in_file_ext.assign(tmp_file_ext);

  std::cout << "Directory for instance has been set to \"" + instdir + "\"."
      << std::endl;

  std::cout << "The filename is \"" + TEST_NAME + "\"." << std::endl;

  std::cout << "Output directory has been set to \"" + tmpoutdir + "\"."
      << std::endl;

  std::cout
      << "Summary instance information will be sent to \"" + FINAL_LOG_NAME
          + "\"." << std::endl;

  if (!opt_file.empty()) { 
    std::cout << "Reading objective information from \"" + opt_file + "\"."
        << std::endl;
    GlobalVariables::bestObjValue = getObjValueFromFile(opt_file.c_str(), prob_name);
    if (isInfinity(std::abs(GlobalVariables::bestObjValue))) {
      error_msg(errorstring, "Could not find instance %s in obj value file.\n", prob_name.c_str());
      writeErrorToII(errorstring, GlobalVariables::log_file);
      exit(1);
    }
    std::cout << "Best known IP objective value (read from file): "
        << GlobalVariables::bestObjValue << "." << std::endl;
  }

//  if (!sol_file.empty()) {
//    printf("Reading solution from %s.\n", sol_file.c_str());
//    getOptSolFromFile(sol_file.c_str(), opt_sol);
//  }

  // If the instance info file did not exist, then we want to add a header so that we know what each column is
  if (!ii_exists) {
    writeHeaderToII(GlobalVariables::log_file);
  }

  // Set f_name to directory+TEST_NAME. Same for the stub.
  f_name = instdir + TEST_NAME;
//  f_name_stub = dir + GlobalVariables::prob_name;
  // Do the same for the out stub
  GlobalVariables::out_f_name_stub = tmpoutdir + GlobalVariables::prob_name;

  // Check that the file exists
  bool f_exists = fexists(f_name.c_str());
  if (!f_exists) {
    // Check if maybe it is a compressed file and the ".gz" extension was forgotten
    f_exists = fexists((f_name + ".gz").c_str());
    if (f_exists) {
#ifdef TRACE
      std::cout << "Inputted filename was missing .gz extension. "
          << "As no such file exists, but an archived version does exist, working with the archived file."
          << std::endl;
#endif
      f_name = f_name + ".gz";
    } else {
      error_msg(errorstring, "File %s does not exist!\n", f_name.c_str());
      writeErrorToII(errorstring, GlobalVariables::log_file);
      exit(1);
    }
  }

  // Write the problem name into the instances info file
  writeProbNameToII(GlobalVariables::prob_name, GlobalVariables::log_file);

  return 0;
} /* init */

/**********************************************************/
/**
 * Wrap everything up and return summary statistics before exiting
 */
int wrapUp(int returncode) {
  char start_time_string[25];
  snprintf(start_time_string, sizeof(start_time_string) / sizeof(char), "%s",
      asctime(timeinfo));

  time(&end_time_t);
  timeinfo = localtime(&end_time_t);
  char end_time_string[25];
  snprintf(end_time_string, sizeof(end_time_string) / sizeof(char), "%s",
      asctime(timeinfo));

  GlobalVariables::timeStats.end_timer(GlobalConstants::TOTAL_TIME);
  char filename[300];
  snprintf(filename, sizeof(filename) / sizeof(char), "%s-Times.csv",
      GlobalVariables::out_f_name_stub.c_str());
#ifdef SHOULD_WRITE_TIME_INFO
#  ifdef TRACE
  printf("\n## Saving time information to file. ##\n");
#  endif
  writeTimeInfo(GlobalVariables::timeStats, GlobalVariables::time_name,
      filename, end_time_string);
#endif
  writeTimeInfoToII(GlobalVariables::timeStats, GlobalVariables::time_name,
      end_time_string, GlobalVariables::log_file);

#if defined TRACE && defined SHOULD_WRITE_CUTSOLVER_ERRORS
  if (!NO_GEN_CUTS) {
    printf("\n## Saving errors in cut solver to file. ##\n");
    snprintf(filename, sizeof(filename) / sizeof(char),
        "%s-CutSolverErrors.csv", GlobalVariables::out_f_name_stub.c_str());
    writeCutSolverErrors(num_obj_tried, filename);
  }
#endif

  writeEntryToII(prob_name, GlobalVariables::log_file);
  writeEntryToII("DONE!", GlobalVariables::log_file, '\n');
//  fprintf(GlobalVariables::inst_info_out, "%s,", end_time_string);
//  fprintf(GlobalVariables::inst_info_out, "DONE!\n");

  if (GlobalVariables::log_file != NULL)
    fclose(GlobalVariables::log_file);

  // For reference print parameters that were used
  printf("\n## Parameters that were used. ##\n");
#ifdef SHOULD_WRITE_PARAMS
  snprintf(filename, sizeof(filename) / sizeof(char), "%s-params.txt",
      GlobalVariables::out_f_name_stub.c_str());
  FILE* paramOut = fopen(filename, "w");
  saved_param.printParams(paramOut);
  fclose(paramOut);
#endif
  printf("PROB_NAME: %s\n", GlobalVariables::prob_name.c_str());
  printf("OUT_F_NAME_STUB: %s\n", GlobalVariables::out_f_name_stub.c_str());
  printf("INSTANCES_INFO: %s\n", FINAL_LOG_NAME.c_str());
  for (int param_ind = 0; param_ind < NUM_PARAMS; param_ind++) {
    printf("%s: %d\n", ParamName[param_ind].c_str(), param.paramVal[param_ind]);
  }
  printf("EPS: %e\n", param.getEPS());
  printf("RAYEPS: %e\n", param.getRAYEPS());
  printf("MIN_ORTHOGONALITY: %e\n", param.getMINORTHOGONALITY());
  printf("MIN_VIOL_ABS: %e\n", param.getMIN_VIOL_ABS());
  printf("MIN_VIOL_REL: %e\n", param.getMIN_VIOL_REL());
  printf("TIMELIMIT: %e\n", param.getTIMELIMIT());

  printf("%s", compare_cuts_output.c_str());

  printf("\n");
  printf("Start time: %s\n", start_time_string);
  printf("Finish time: %s\n", end_time_string);
  printf("Elapsed time: %.f seconds\n", difftime(end_time_t, start_time_t));
  printf("\n## DONE! ##\n");
  return returncode;
} /* wrapUp */
