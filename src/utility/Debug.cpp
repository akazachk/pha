#include <fstream>
#include <cmath> // isnan
#include <algorithm> // For max_element, min_element
#include "Utility.hpp"
#include "Debug.hpp"

/*
 *
 */
void writeCondNumToII(double LGICBasisCondNumInit, double LGICBasisCondNumAll,
    double LSICBasisCondNumInit, double LSICBasisCondNumAll,
    double OSICBasisCondNum, double LGICAndLSICBasisCondNumInit,
    double LGICAndLSICBasisCondNumAll, double LGICAndOSICBasisCondNumInit,
    double LGICAndOSICBasisCondNumAll, FILE* myfile) {
  if (myfile == NULL)
    return;
  //////////
  // Basis cond num
  //////////
  fprintf(myfile, "%2.3f%c", LGICBasisCondNumInit, SEP); // 1
  fprintf(myfile, "%2.3f%c", LGICBasisCondNumAll, SEP); // 2
  fprintf(myfile, "%2.3f%c", LSICBasisCondNumInit, SEP); // 3
  fprintf(myfile, "%2.3f%c", LSICBasisCondNumAll, SEP); // 4
  fprintf(myfile, "%2.3f%c", OSICBasisCondNum, SEP); // 5
  fprintf(myfile, "%2.3f%c", LGICAndLSICBasisCondNumInit, SEP); // 6
  fprintf(myfile, "%2.3f%c", LGICAndLSICBasisCondNumAll, SEP); // 7
  fprintf(myfile, "%2.3f%c", LGICAndOSICBasisCondNumInit, SEP); // 7
  fprintf(myfile, "%2.3f%c", LGICAndOSICBasisCondNumAll, SEP); // 7
} /* writeCondNumToII */

/***********************************************************************/
void printSplitInfo(const OsiSolverInterface* const solverMain,
    SolutionInfo& solnInfoMain) {
//#ifdef TRACE
//  printf(
//      "\n## Store split feasibility info as well as num cuts generated per each feasible split. ##\n");
//#endif
  char filename[300];
  int tmpindex1 = 0;
  snprintf(filename, sizeof(filename) / sizeof(char),
      "%sMain-splitFeasInfo.csv", GlobalVariables::out_f_name_stub.c_str());
  FILE* splitInfoOuttmp = fopen(filename, "w");
  fprintf(splitInfoOuttmp,
      "Split Feasibility Info (Main, num splits: %d, num feasible: %d)\n",
      solnInfoMain.numSplits, solnInfoMain.numFeasSplits);
  fprintf(splitInfoOuttmp,
      "Split,Split var index,Split name,Feas sides > 0 && post-SIC feas && chosen?,"
          "Num feas sides,Floor feas?,Ceil feas?,Num intersecting rays,Num rays int floor,Num rays int ceil,"
          "Percent intersecting rays,Num GICs generated\n");
  for (int split = 0; split < solnInfoMain.numSplits; split++) {
    int currsplit = solnInfoMain.fractionalCore[split];
    fprintf(splitInfoOuttmp, "%d,", split);
    fprintf(splitInfoOuttmp, "%d,", currsplit);
    fprintf(splitInfoOuttmp, "%s,",
        solverMain->getColName(currsplit).c_str());
    fprintf(splitInfoOuttmp, "%d,",
        (int) solnInfoMain.splitFeasFlag[split]);
    fprintf(splitInfoOuttmp, "%d,",
        (int) (solnInfoMain.fracCoreFloorFeasFlag[split]
            + solnInfoMain.fracCoreCeilFeasFlag[split]));
    fprintf(splitInfoOuttmp, "%d,",
        (int) solnInfoMain.fracCoreFloorFeasFlag[split]);
    fprintf(splitInfoOuttmp, "%d,",
        (int) solnInfoMain.fracCoreCeilFeasFlag[split]);
    if (solnInfoMain.splitFeasFlag[split]) {
      fprintf(splitInfoOuttmp, "%d,",
          solnInfoMain.count_intersecting_rays[tmpindex1]);
      fprintf(splitInfoOuttmp, "%d,",
          solnInfoMain.count_intersecting_rays_floor[tmpindex1]);
      fprintf(splitInfoOuttmp, "%d,",
          solnInfoMain.count_intersecting_rays_ceil[tmpindex1]);
      fprintf(splitInfoOuttmp, "%2.3f,",
          (double) solnInfoMain.count_intersecting_rays[tmpindex1]
              / solnInfoMain.numNB);
      tmpindex1++;
    }
    fprintf(splitInfoOuttmp, "\n");
  }
  fclose(splitInfoOuttmp);
} /* printSplitInfo */

/***********************************************************************/
/**
 * @brief Writes euclidean distance by which initial vertex is cut,
 * and the objective function after adding the single cut
 */
void writeCutSummary(std::string headerText, AdvCuts &structCut, char* filename,
    FILE* inst_info_out) {
  FILE* CutOut = fopen(filename, "w");
  if (CutOut == NULL) {
    error_msg(errorstring,
        "Unable to open file to write cut summary info for filename %s.\n",
        filename);
    writeErrorToII(errorstring, inst_info_out);
    exit(1);
  }
  fprintf(CutOut, "%s", headerText.c_str());
  fprintf(CutOut, "Cut#,Facet,Facetname,Euclideandist,LPobj\n");
  for (int currcut = 0; currcut < (int) structCut.size(); currcut++) {
//    fprintf(CutOut, "%d,%d,%s,%2.3f,%2.3f\n", currcut, currsplit,
//        solver->getColName(currsplit).c_str(), structCut[currcut].eucl,
//        structCut[currcut].obj);
    fprintf(CutOut, "%d,%d,%s,%2.3f,%2.3f\n", currcut,
        structCut[currcut].cgsIndex, structCut[currcut].cgsName.c_str(),
        -1.0, structCut[currcut].postCutObj);
  }
  fclose(CutOut);
} /* writeCutSummary */

/***********************************************************************/
/**
 * @brief Writes euclidean distance by which initial vertex is cut,
 * and the objective function after adding the single cut
 */
void writeCutSummary(std::string headerText, std::vector<AdvCuts> &structCut, char* filename,
    FILE* inst_info_out) {
  FILE* CutOut = fopen(filename, "w");
  if (CutOut == NULL) {
    error_msg(errorstring,
        "Unable to open file to write cut summary info for filename %s.\n",
        filename);
    writeErrorToII(errorstring, inst_info_out);
    exit(1);
  }
  fprintf(CutOut, "%s", headerText.c_str());
  fprintf(CutOut, "Cut#,Cgs,Facet,Facetname,Euclideandist,LPobj\n");
  for (int cgs_ind = 0; cgs_ind < (int) structCut.size(); cgs_ind++) {
    for (int currcut = 0; currcut < (int) structCut[cgs_ind].size();
        currcut++) {
//      fprintf(CutOut, "%d,%d,%s,%2.3f,%2.3f\n", currcut, currsplit,
//          solver->getColName(currsplit).c_str(),
//          structCut[splitind][currcut].eucl,
//          structCut[splitind][currcut].obj);
      fprintf(CutOut, "%d,%d,%d,%s,%2.3f,%2.3f\n", currcut, cgs_ind,
          structCut[cgs_ind][currcut].cgsIndex, structCut[cgs_ind][currcut].cgsName.c_str(),
          -1.0, structCut[cgs_ind][currcut].postCutObj);
    }
  }
  fclose(CutOut);
} /* writeCutSummary */

/***********************************************************************/
/**
 * @brief Writes sets of cuts that are stored in structural space.
 */
void writeCutsInStructSpace(std::string headerText, AdvCuts &structCut,
    SolutionInfo &solnInfo, const PHASolverInterface* const solver, char* filename,
    FILE* inst_info_out) {
  FILE* StructCutOut = fopen(filename, "w");
  if (StructCutOut == NULL) {
    error_msg(errorstring,
        "Unable to open file to write cuts in structural space for filename %s.\n",
        filename);
    writeErrorToII(errorstring, inst_info_out);
    exit(1);
  }
  fprintf(StructCutOut, "%s", headerText.c_str());
  fprintf(StructCutOut, "%s,%s,%s,%s,%s\n", "Cut #", "Facet", "Facet name", "Vars"
      "Rhs", "Coeff");
  fprintf(StructCutOut, ",,,,,");
  for (int col = 0; col < solnInfo.numCols; col++) {
    fprintf(StructCutOut, "%d,", col);
  }
  fprintf(StructCutOut, "\n");
  fprintf(StructCutOut, ",,,,,");
  for (int col = 0; col < solnInfo.numCols; col++) {
    fprintf(StructCutOut, "%s,", solver->getColName(col).c_str());
  }
  fprintf(StructCutOut, "\n");
  for (int currcut = 0; currcut < (int) structCut.size(); currcut++) {
    fprintf(StructCutOut, "%d,", currcut);
    fprintf(StructCutOut, "%d,", structCut[currcut].cgsIndex);
    fprintf(StructCutOut, "%s,", structCut[currcut].cgsName.c_str());
    fprintf(StructCutOut, "%f,", structCut[currcut].rhs());
    for (int coeff = 0; coeff < (int) structCut[currcut].num_coeff;
        coeff++) {
//      if (isZero(structCut[currcut][coeff])) {
//        structCut[currcut][coeff] = 0.0;
//      }
      fprintf(StructCutOut, "%f,", structCut[currcut][coeff]);
      if (std::isnan(structCut[currcut][coeff])) {
        error_msg(errorstring, "Found NAN! Coeff %d of cut %d is %f.\n",
            coeff, currcut, structCut[currcut][coeff]);
        writeErrorToII(errorstring, inst_info_out);
        exit(1);
      }
    }
    fprintf(StructCutOut, "\n");
  }
  fprintf(StructCutOut, "\n");

  fclose(StructCutOut);
}

/***********************************************************************/
/**
 * @brief Writes sets of cuts that are stored in structural space.
 */
void writeCutsInStructSpace(std::string headerText,
    std::vector<AdvCuts> &structCut, SolutionInfo &solnInfo,
    const PHASolverInterface* const solver, char* filename, FILE* inst_info_out) {
  FILE* StructCutOut = fopen(filename, "w");
  if (StructCutOut == NULL) {
    error_msg(errorstring,
        "Unable to open file to write cuts in structural space for filename %s.\n",
        filename);
    writeErrorToII(errorstring, inst_info_out);
    exit(1);
  }
  fprintf(StructCutOut, "%s", headerText.c_str());
  fprintf(StructCutOut, "%s,%s,%s,%s,%s\n", "Cut #", "Split", "Split name",
      "Rhs", "Coeff");
  fprintf(StructCutOut, ",,,,");
  for (int col = 0; col < solnInfo.numCols; col++) {
    fprintf(StructCutOut, "%d,", col);
  }
  fprintf(StructCutOut, "\n");
  fprintf(StructCutOut, ",,,,");
  for (int col = 0; col < solnInfo.numCols; col++) {
    fprintf(StructCutOut, "%s,", solver->getColName(col).c_str());
  }
  fprintf(StructCutOut, "\n");
  for (int cgs_ind = 0; cgs_ind < (int) structCut.size(); cgs_ind++) {
    for (int currcut = 0; currcut < (int) structCut[cgs_ind].size();
        currcut++) {
      fprintf(StructCutOut, "%d-%d,", cgs_ind, currcut);
      fprintf(StructCutOut, "%d,", structCut[cgs_ind][currcut].cgsIndex);
      fprintf(StructCutOut, "%s,", structCut[cgs_ind][currcut].cgsName.c_str());
      fprintf(StructCutOut, "%f,", structCut[cgs_ind][currcut].rhs());
      for (int coeff = 0;
          coeff < (int) structCut[cgs_ind][currcut].num_coeff;
          coeff++) {
//        if (std::abs(structCut[cgs_ind][currcut][coeff]) < param.getEPS()) {
//          structCut[cgs_ind][currcut][coeff] = 0.0;
//        }
        fprintf(StructCutOut, "%f,",
            structCut[cgs_ind][currcut][coeff]);
        if (std::isnan(structCut[cgs_ind][currcut][coeff])) {
          error_msg(errorstring,
              "Found NAN! Coeff %d of cut %d is %f.\n", coeff,
              currcut, structCut[cgs_ind][currcut][coeff]);
          writeErrorToII(errorstring, inst_info_out);
          exit(1);
        }
      }
      fprintf(StructCutOut, "\n");
    }
  }
  fprintf(StructCutOut, "\n");

  fclose(StructCutOut);
}

/***********************************************************************/
/**
 * @brief Writes sets of cuts that are stored in NB space.
 */
void writeCutsInNBSpace(std::string headerText, AdvCuts &NBCut,
    SolutionInfo &solnInfo, const PHASolverInterface* const solver, char* filename,
    FILE* inst_info_out) {
  FILE* NBCutOut = fopen(filename, "w");
  if (NBCutOut == NULL) {
    error_msg(errorstring,
        "Unable to open file to write cuts in non-basic space for filename %s.\n",
        filename);
    writeErrorToII(errorstring, inst_info_out);
    exit(1);
  }
  fprintf(NBCutOut, "%s", headerText.c_str());
  fprintf(NBCutOut, "%s,%s,%s,%s,%s\n", "Cut #", "Split", "Split name", "Rhs",
      "Coeff");
  fprintf(NBCutOut, ",,,,");
  for (int nbcol = 0; nbcol < (int) solnInfo.nonBasicOrigVarIndex.size();
      nbcol++) {
    int currvar = solnInfo.nonBasicOrigVarIndex[nbcol];
    fprintf(NBCutOut, "%d,", currvar);
  }
  for (int nbcol = 0; nbcol < (int) solnInfo.nonBasicSlackVarIndex.size();
      nbcol++) {
    int currvar = solnInfo.nonBasicSlackVarIndex[nbcol];
    fprintf(NBCutOut, "%d,", currvar);
  }
  fprintf(NBCutOut, "\n");
  fprintf(NBCutOut, ",,,,");
  for (int nbcol = 0; nbcol < (int) solnInfo.nonBasicOrigVarIndex.size();
      nbcol++) {
    int currvar = solnInfo.nonBasicOrigVarIndex[nbcol];
    fprintf(NBCutOut, "%s,", solver->getColName(currvar).c_str());
  }
  for (int nbcol = 0; nbcol < (int) solnInfo.nonBasicSlackVarIndex.size();
      nbcol++) {
    int currvar = solnInfo.nonBasicSlackVarIndex[nbcol];
    fprintf(NBCutOut, "%s,",
        solver->getRowName(currvar - solnInfo.numCols).c_str());
  }
  fprintf(NBCutOut, "\n");
  for (int currcut = 0; currcut < (int) NBCut.size(); currcut++) {
    fprintf(NBCutOut, "%d,", currcut);
    fprintf(NBCutOut, "%d,", NBCut[currcut].cgsIndex);
    fprintf(NBCutOut, "%s,", NBCut[currcut].cgsName.c_str());
    fprintf(NBCutOut, "%f,", NBCut[currcut].rhs());
    for (int coeff = 0; coeff < (int) NBCut[currcut].num_coeff; coeff++) {
//      if (std::abs(NBCut[currcut].coeff[coeff]) < param.getEPS()) {
//        NBCut[currcut].coeff[coeff] = 0.0;
//      }
      fprintf(NBCutOut, "%f,", NBCut[currcut][coeff]);
      if (std::isnan(NBCut[currcut][coeff])) {
        error_msg(errorstring, "Found NAN! Coeff %d of cut %d is %f.\n",
            coeff, currcut, NBCut[currcut][coeff]);
        writeErrorToII(errorstring, inst_info_out);
        exit(1);
      }
    }
    fprintf(NBCutOut, "\n");
  }
  fprintf(NBCutOut, "\n");

  fclose(NBCutOut);
}

/***********************************************************************/
/**
 * @brief Writes sets of cuts that are stored in NB space.
 */
void writeCutsInNBSpace(std::string headerText, std::vector<AdvCuts> &NBCut,
    SolutionInfo &solnInfo, const PHASolverInterface* const solver, char* filename,
    FILE* inst_info_out) {
  FILE* NBCutOut = fopen(filename, "w");
  if (NBCutOut == NULL) {
    error_msg(errorstring,
        "Unable to open file to write cuts in non-basic space for filename %s.\n",
        filename);
    writeErrorToII(errorstring, inst_info_out);
    exit(1);
  }
  fprintf(NBCutOut, "%s", headerText.c_str());
  fprintf(NBCutOut, "%s,%s,%s,%s,%s\n", "Cut #", "Split", "Split name", "Rhs",
      "Coeff");
  fprintf(NBCutOut, ",,,,");
  for (int nbcol = 0; nbcol < (int) solnInfo.nonBasicOrigVarIndex.size();
      nbcol++) {
    int currvar = solnInfo.nonBasicOrigVarIndex[nbcol];
    fprintf(NBCutOut, "%d,", currvar);
  }
  for (int nbcol = 0; nbcol < (int) solnInfo.nonBasicSlackVarIndex.size();
      nbcol++) {
    int currvar = solnInfo.nonBasicSlackVarIndex[nbcol];
    fprintf(NBCutOut, "%d,", currvar);
  }
  fprintf(NBCutOut, "\n");
  fprintf(NBCutOut, ",,,,");
  for (int nbcol = 0; nbcol < (int) solnInfo.nonBasicOrigVarIndex.size();
      nbcol++) {
    int currvar = solnInfo.nonBasicOrigVarIndex[nbcol];
    fprintf(NBCutOut, "%s,", solver->getColName(currvar).c_str());
  }
  for (int nbcol = 0; nbcol < (int) solnInfo.nonBasicSlackVarIndex.size();
      nbcol++) {
    int currvar = solnInfo.nonBasicSlackVarIndex[nbcol];
    fprintf(NBCutOut, "%s,",
        solver->getRowName(currvar - solnInfo.numCols).c_str());
  }
  fprintf(NBCutOut, "\n");
  for (int cgs_ind = 0; cgs_ind < (int) NBCut.size(); cgs_ind++) {
    for (int currcut = 0; currcut < (int) NBCut[cgs_ind].size();
        currcut++) {
      fprintf(NBCutOut, "%d,", currcut);
      fprintf(NBCutOut, "%d,", NBCut[cgs_ind][currcut].cgsIndex);
      fprintf(NBCutOut, "%s,", NBCut[cgs_ind][currcut].cgsName.c_str());
      fprintf(NBCutOut, "%f,", NBCut[cgs_ind][currcut].rhs());
      for (int coeff = 0;
          coeff < (int) NBCut[cgs_ind][currcut].num_coeff; coeff++) {
//        if (std::abs(NBCut[splitind][currcut].coeff[coeff]) < param.getEPS()) {
//          NBCut[splitind][currcut].coeff[coeff] = 0.0;
//        }
        fprintf(NBCutOut, "%f,", NBCut[cgs_ind][currcut][coeff]);
        if (std::isnan(NBCut[cgs_ind][currcut][coeff])) {
          error_msg(errorstring,
              "Found NAN! Coeff %d of cut %d is %f.\n", coeff,
              currcut, NBCut[cgs_ind][currcut][coeff]);
          writeErrorToII(errorstring, inst_info_out);
          exit(1);
        }
      }
      fprintf(NBCutOut, "\n");
    }
  }
  fprintf(NBCutOut, "\n");

  fclose(NBCutOut);
}

/***********************************************************************/
/**
 * @brief Writes how far along each NB ray we go for each cut.
 */
void writeRayDistInfo(std::string headerText, const PHASolverInterface* const solver,
    SolutionInfo &solnInfo, AdvCuts &structCut,
    std::vector<std::vector<double> > &distAlongRay, char* filename,
    FILE* inst_info_out) {
  FILE* NBDist = fopen(filename, "w");
  if (NBDist == NULL) {
    error_msg(errorstring,
        "Unable to open file to write ray distance info for filename %s.\n",
        filename);
    writeErrorToII(errorstring, inst_info_out);
    exit(1);
  }
  fprintf(NBDist, "%s", headerText.c_str());
  fprintf(NBDist, "Ray,Ray name,Cut\n");
  fprintf(NBDist, ",,");
  for (int currcut = 0; currcut < (int) structCut.size(); currcut++) {
    fprintf(NBDist, "%d,", structCut[currcut].cgsIndex);
  }
  fprintf(NBDist, "\n");
  fprintf(NBDist, ",,");
  for (int currcut = 0; currcut < (int) structCut.size(); currcut++) {
    fprintf(NBDist, "%s,", structCut[currcut].cgsName.c_str());
  }
  fprintf(NBDist, "\n");
  for (int ray = 0;
      ray
          < (int) (solnInfo.nonBasicOrigVarIndex.size()
              + solnInfo.nonBasicSlackVarIndex.size()); ray++) {
    std::string rayname;
    if (ray < (int) solnInfo.nonBasicOrigVarIndex.size()) {
      int currvar = solnInfo.nonBasicOrigVarIndex[ray];
      rayname = solver->getColName(currvar);
      if (solnInfo.cstat[currvar] == 2) {
        // Then we are going in the negative direction
        rayname = "-" + rayname;
      }
    } else {
      rayname = solver->getRowName(
          solnInfo.nonBasicSlackVarIndex[ray
              - (int) solnInfo.nonBasicOrigVarIndex.size()]
              - solnInfo.numCols);
    }

    fprintf(NBDist, "%d,%s,", ray, rayname.c_str());
    for (int currcut = 0; currcut < (int) structCut.size(); currcut++) {
      double currdist = distAlongRay[currcut][ray];
      if (isInfinity(currdist, solver->getInfinity(), param.getEPS())) {
        fprintf(NBDist, "inf,");
      } else {
        fprintf(NBDist, "%2.3f,", currdist);
      }
    }
    fprintf(NBDist, "\n");
  }

  fprintf(NBDist, ",SICEuclDist,");
  for (int currcut = 0; currcut < (int) structCut.size(); currcut++) {
//    fprintf(NBDist, "%2.3f,", structCut[currcut].eucl);
    fprintf(NBDist, "not_impl,");
  }
  fprintf(NBDist, "\n");

  fprintf(NBDist, ",SICObj,");
  for (int currcut = 0; currcut < (int) structCut.size(); currcut++) {
//    int indexOfSplitInFeasSplits = structCut[currcut].splitIndex;
////    printf("Split index is %d.\n", indexOfSplitInFeasSplits);
//    int indexOfSplitInFracCore = solnInfo.feasSplit[indexOfSplitInFeasSplits];
////    printf("Split index in frac core is %d.\n", indexOfSplitInFracCore);
//    if (solnInfo.splitFeas[indexOfSplitInFracCore]) {
//    if (structCut[currcut].feas) {
    double currobj = structCut[currcut].postCutObj;
    if (!isInfinity(currobj, solver->getInfinity(), param.getEPS()))
      fprintf(NBDist, "%2.3f,", currobj);
    else
      fprintf(NBDist, "error,");
//    } else {
//      fprintf(NBDist, "infeas,");
//    }
  }
  fprintf(NBDist, "\n");

  fclose(NBDist);
}

/***********************************************************************/
/**
 * @brief Writes how far along each NB ray we go for each cut.
 */
void writeRayDistInfo(std::string headerText, const PHASolverInterface* const solver,
    SolutionInfo &solnInfo, std::vector<AdvCuts> &structCut,
    std::vector<std::vector<std::vector<double> > > &distAlongRay,
    char* filename, FILE* inst_info_out) {
  FILE* NBDist = fopen(filename, "w");
  if (NBDist == NULL) {
    error_msg(errorstring,
        "Unable to open file to write ray distance info for filename %s.\n",
        filename);
    writeErrorToII(errorstring, inst_info_out);
    exit(1);
  }
  fprintf(NBDist, "%s", headerText.c_str());
  fprintf(NBDist, "Ray,Ray name,Cut\n");
  fprintf(NBDist, ",,");
  for (int splitind = 0; splitind < (int) structCut.size(); splitind++) {
    for (int currcut = 0; currcut < (int) structCut[splitind].size();
        currcut++) {
      fprintf(NBDist, "%d,", structCut[splitind][currcut].cgsIndex);
    }
  }
  fprintf(NBDist, "\n");
  fprintf(NBDist, ",,");
  for (int splitind = 0; splitind < (int) structCut.size(); splitind++) {
    for (int currcut = 0; currcut < (int) structCut[splitind].size();
        currcut++) {
      fprintf(NBDist, "%s,", structCut[splitind][currcut].cgsName.c_str());
    }
  }
  fprintf(NBDist, "\n");
  for (int ray = 0;
      ray
          < (int) (solnInfo.nonBasicOrigVarIndex.size()
              + solnInfo.nonBasicSlackVarIndex.size()); ray++) {
    std::string rayname;
    if (ray < (int) solnInfo.nonBasicOrigVarIndex.size()) {
      int currvar = solnInfo.nonBasicOrigVarIndex[ray];
      rayname = solver->getColName(currvar);
      if (solnInfo.cstat[currvar] == 2) {
        // Then we are going in the negative direction
        rayname = "-" + rayname;
      }
    } else {
      rayname = solver->getRowName(
          solnInfo.nonBasicSlackVarIndex[ray
              - (int) solnInfo.nonBasicOrigVarIndex.size()]
              - solnInfo.numCols);
    }

    fprintf(NBDist, "%d,%s,", ray, rayname.c_str());
    for (int splitind = 0; splitind < (int) structCut.size(); splitind++) {
      for (int currcut = 0; currcut < (int) structCut[splitind].size();
          currcut++) {
        double currdist = distAlongRay[splitind][currcut][ray];
        if (isInfinity(currdist, solver->getInfinity(), param.getEPS())) {
          fprintf(NBDist, "inf,");
        } else {
          fprintf(NBDist, "%2.3f,", currdist);
        }
      }
    }
    fprintf(NBDist, "\n");
  }

  fprintf(NBDist, ",SICEuclDist,");
  for (int splitind = 0; splitind < (int) structCut.size(); splitind++) {
    for (int currcut = 0; currcut < (int) structCut[splitind].size();
        currcut++) {
//      fprintf(NBDist, "%2.3f,", structCut[splitind][currcut].eucl);
      fprintf(NBDist, "not_impl,");
    }
  }
  fprintf(NBDist, "\n");

  fprintf(NBDist, ",SICObj,");
  for (int splitind = 0; splitind < (int) structCut.size(); splitind++) {
    for (int currcut = 0; currcut < (int) structCut[splitind].size();
        currcut++) {
//    int indexOfSplitInFeasSplits = structCut[currcut].splitIndex;
////    printf("Split index is %d.\n", indexOfSplitInFeasSplits);
//    int indexOfSplitInFracCore = solnInfo.feasSplit[indexOfSplitInFeasSplits];
////    printf("Split index in frac core is %d.\n", indexOfSplitInFracCore);
//    if (solnInfo.splitFeas[indexOfSplitInFracCore]) {
//      if (structCut[splitind][currcut].feas) {
      double currobj = structCut[splitind][currcut].postCutObj;
      if (!isInfinity(currobj, solver->getInfinity(), param.getEPS()))
        fprintf(NBDist, "%2.3f,", currobj);
      else
        fprintf(NBDist, "error,");
//      } else {
//        fprintf(NBDist, "infeas,");
//      }
    }
  }
  fprintf(NBDist, "\n");

  fclose(NBDist);
}

void writeInterPtsAndRaysInJSpace(FILE* fptr, const SolutionInfo &solnInfo,
    const std::vector<IntersectionInfo>& interPtsAndRays,
    const std::vector<Hplane>& hplaneStore) {
  int numSplits = interPtsAndRays.size();

  fprintf(fptr,
      "Contains information about the intersection points and rays created during PCuts\n\n");

  fprintf(fptr,
      "splitIndx,numPoints,numFinalPoints,numRays,Total\n");
  for (int s = 0; s < numSplits; s++) {
    fprintf(fptr, "%d,%d,%d,%d,%d,\n",
        solnInfo.feasSplitVar[s], interPtsAndRays[s].numPoints,
        interPtsAndRays[s].numFinalPoints, interPtsAndRays[s].numRays,
        interPtsAndRays[s].numPoints + interPtsAndRays[s].numRays);
  }
  fprintf(fptr, "\n\n");

//  int numElmts;
  //, count;
  fprintf(fptr,
      "splitIndx,No,cutRayNBColIndx,cutRayVarIndx,hplaneRowIndx,newRayNBColIndx,newRayVarIndx,PtFlag,indx1,val1,indx2,val2,finalFlag\n");
  //  fprintf(fptr, ",,,");
  //  for (int j = 0; j < numCols; j++)
  //    fprintf(fptr, "%d,", j);
  //  fprintf(fptr, "\n");
  std::vector<bool> flag(solnInfo.numNB);
  std::vector<int> index(solnInfo.numNB);

  for (int s = 0; s < numSplits; s++) {
    const int numPointsOrRays = interPtsAndRays[s].getNumRows();
    for (int r = 0; r < numPointsOrRays; r++) {
      bool pointOrRayFlag = !isZero(interPtsAndRays[s].RHS[r]);
      const int hplane_ind = interPtsAndRays[s].hplane[r];
      const int cut_ray_ind = interPtsAndRays[s].cutRay[r];
      const int new_ray_ind = interPtsAndRays[s].newRay[r];

      const CoinShallowPackedVector tmpVector =
          interPtsAndRays[s].getVector(r);
//      numElmts = tmpVector.getNumElements();
//      const int* elmtIndices = tmpVector.getIndices();
//      const double* elmtVals = tmpVector.getElements();

      if (cut_ray_ind >= 0) { // Then necessarily new_ray_ind >= 0
        fprintf(fptr, "%d,%d,%d,%d,%d,%d,%d,%d,",
            solnInfo.feasSplitVar[s], r,
            cut_ray_ind,
            solnInfo.nonBasicVarIndex[cut_ray_ind],
            hplaneStore[hplane_ind].row,
            new_ray_ind,
            solnInfo.nonBasicVarIndex[new_ray_ind],
            pointOrRayFlag);
      } else { 
        if (new_ray_ind >= 0) {
          fprintf(fptr, "%d,%d,-1,-1,-1,%d,%d,%d,",
              solnInfo.feasSplitVar[s], r,
              interPtsAndRays[s].newRay[r],
              solnInfo.nonBasicVarIndex[interPtsAndRays[s].newRay[r]],
              pointOrRayFlag);
        } else {
          fprintf(fptr, "%d,%d,-1,-1,-1,-2,-2,%d,",
              solnInfo.feasSplitVar[s], r,
            pointOrRayFlag);
        }
      }

      if (cut_ray_ind >= 0) {
        fprintf(fptr, "%d,%f,", cut_ray_ind, tmpVector[cut_ray_ind]);
      } else {
        fprintf(fptr, "%d,0.0,", cut_ray_ind);
      }
      if (new_ray_ind >= 0) {
        fprintf(fptr, "%d,%f,", new_ray_ind, tmpVector[new_ray_ind]);
      } else {
        fprintf(fptr, "%d,0.0,", new_ray_ind);
      }

//      for (int k = 0; k < numElmts; k++) {
//        fprintf(fptr, "%d,%lf,", elmtIndices[k], elmtVals[k]);
//      }

      const int final_flag = interPtsAndRays[s].finalFlag[r];
      fprintf(fptr, "%d,", final_flag);
      fprintf(fptr, "\n");
    }
  }
}

void writeInterPtsAndRaysInJSpace(FILE* fptr, const SolutionInfo &solnInfo,
    const std::vector<std::vector<IntersectionInfo> >& interPtsAndRays,
    const std::vector<std::string>& cgsName,
    const std::vector<Hplane>& hplaneStore) {
  const int num_cgs = interPtsAndRays.size();

  fprintf(fptr,
      "Contains information about the intersection points and rays created\n\n");

  fprintf(fptr, "cgsName,numPoints,numFinalPoints,numRays,Total\n");
  for (int cgs_ind = 0; cgs_ind < num_cgs; cgs_ind++) {
    int numPoints = 0, numFinalPoints = 0, numRays = 0;
    for (int f = 0; f < (int) interPtsAndRays[cgs_ind].size(); f++) {
      numPoints += interPtsAndRays[cgs_ind][f].numPoints;
      numFinalPoints += interPtsAndRays[cgs_ind][f].numFinalPoints;
      numRays += interPtsAndRays[cgs_ind][f].numRays;
    }
    fprintf(fptr, "%s,%d,%d,%d,%d,\n", cgsName[cgs_ind].c_str(), numPoints,
        numFinalPoints, numRays, numPoints + numRays);
  }
  fprintf(fptr, "\n\n");

//  int numElmts;
  //, count;
  fprintf(fptr,
      "cgsName,No,cutRayNBColIndx,cutRayVarIndx,hplaneRowIndx,newRayNBColIndx,newRayVarIndx,PtFlag,indx1,val1,indx2,val2,finalFlag\n");
  //  fprintf(fptr, ",,,");
  //  for (int j = 0; j < numCols; j++)
  //    fprintf(fptr, "%d,", j);
  //  fprintf(fptr, "\n");
  std::vector<bool> flag(solnInfo.numNB);
  std::vector<int> index(solnInfo.numNB);

  for (int cgs_ind = 0; cgs_ind < num_cgs; cgs_ind++) {
    for (int f = 0; f < (int) interPtsAndRays[cgs_ind].size(); f++) {
      for (int r = 0; r < interPtsAndRays[cgs_ind][f].getNumRows(); r++) {
        bool pointOrRayFlag = !isZero(interPtsAndRays[cgs_ind][f].RHS[r]);
        const int hplane_ind = interPtsAndRays[cgs_ind][f].hplane[r];
        const int cut_ray_ind = interPtsAndRays[cgs_ind][f].cutRay[r];
        const int new_ray_ind = interPtsAndRays[cgs_ind][f].newRay[r];

        const CoinShallowPackedVector tmpVector =
            interPtsAndRays[cgs_ind][f].getVector(r);

        if (cut_ray_ind >= 0) { // Then necessarily new_ray_ind >= 0
          fprintf(fptr, "%s,%d,%d,%d,%d,%d,%d,%d,", cgsName[cgs_ind].c_str(), r,
              cut_ray_ind, solnInfo.nonBasicVarIndex[cut_ray_ind],
              hplaneStore[hplane_ind].row, new_ray_ind,
              solnInfo.nonBasicVarIndex[new_ray_ind], pointOrRayFlag);
        } else {
          if (new_ray_ind >= 0) {
            fprintf(fptr, "%s,%d,-1,-1,-1,%d,%d,%d,", cgsName[cgs_ind].c_str(),
                r, interPtsAndRays[cgs_ind][f].newRay[r],
                solnInfo.nonBasicVarIndex[interPtsAndRays[cgs_ind][f].newRay[r]],
                pointOrRayFlag);
          } else {
            fprintf(fptr, "%s,%d,-1,-1,-1,-2,-2,%d,", cgsName[cgs_ind].c_str(),
                r, pointOrRayFlag);
          }
        }

        if (cut_ray_ind >= 0) {
          fprintf(fptr, "%d,%f,", cut_ray_ind, tmpVector[cut_ray_ind]);
        } else {
          fprintf(fptr, "%d,0.0,", cut_ray_ind);
        }
        if (new_ray_ind >= 0) {
          fprintf(fptr, "%d,%f,", new_ray_ind, tmpVector[new_ray_ind]);
        } else {
          fprintf(fptr, "%d,0.0,", new_ray_ind);
        }

        //      for (int k = 0; k < numElmts; k++) {
        //        fprintf(fptr, "%d,%lf,", elmtIndices[k], elmtVals[k]);
        //      }

        const int final_flag = interPtsAndRays[cgs_ind][f].finalFlag[r];
        fprintf(fptr, "%d,", final_flag);
        fprintf(fptr, "\n");
      }
    }
  }
}

void writeInterPtsAndRaysInJSpace(FILE* fptr, const SolutionInfo &solnInfo,
    const std::vector<IntersectionInfo>& interPtsAndRays,
    const std::vector<Hplane>& hplaneStore,
    const std::vector<int>& num_SIC_points,
    const std::vector<int>& num_SIC_final_points,
    const std::vector<int>& num_SIC_rays) {
  const int numSplits = solnInfo.feasSplitVar.size();

  fprintf(fptr,
      "Contains information about the intersection points and rays created from partial hyperplane activation\n\n");

  fprintf(fptr,
      "splitIndx,numPoints,numFinalPoints,numRays,Total,,SICPoints,SICFinalPoints,SICRays\n");
//  int numPoints, numRays;
  for (int s = 0; s < numSplits; s++) {
//    numPoints = 0;
//    numRays = 0;
//    const int numPointsOrRays = interPtsAndRays[s].getNumRows();
//    for (int r = 0; r < numPointsOrRays; r++) {
//      if (std::abs(interPtsAndRays[s].RHS[r] - 1.0) < param.getEPS())
//        numPoints++;
//      else
//        numRays++;
//    }
    fprintf(fptr, "%d,%d,%d,%d,%d,,%d,%d,%d\n",
        solnInfo.feasSplitVar[s], interPtsAndRays[s].numPoints,
        interPtsAndRays[s].numFinalPoints, interPtsAndRays[s].numRays,
        interPtsAndRays[s].numPoints + interPtsAndRays[s].numRays,
        num_SIC_points[s], num_SIC_final_points[s], num_SIC_rays[s]);
  }
  fprintf(fptr, "\n\n");

//  int numElmts;
  //, count;
  fprintf(fptr,
      "splitIndx,No,cutRayNBColIndx,cutRayVarIndx,hplaneRowIndx,newRayNBColIndx,newRayVarIndx,PtFlag,indx1,val1,indx2,val2,finalFlag\n");
  //  fprintf(fptr, ",,,");
  //  for (int j = 0; j < numCols; j++)
  //    fprintf(fptr, "%d,", j);
  //  fprintf(fptr, "\n");
  std::vector<bool> flag(solnInfo.numNB);
  std::vector<int> index(solnInfo.numNB);

  for (int s = 0; s < numSplits; s++) {
    const int numPointsOrRays = interPtsAndRays[s].getNumRows();
    for (int r = 0; r < numPointsOrRays; r++) {
      bool pointOrRayFlag = !isZero(interPtsAndRays[s].RHS[r]);
      const int hplane_ind = interPtsAndRays[s].hplane[r];

      const CoinShallowPackedVector tmpVector =
          interPtsAndRays[s].getVector(r);
//      numElmts = tmpVector.getNumElements();
//      const int* elmtIndices = tmpVector.getIndices();
//      const double* elmtVals = tmpVector.getElements();

      if (interPtsAndRays[s].cutRay[r] >= 0) {
        fprintf(fptr, "%d,%d,%d,%d,%d,%d,%d,%d,",
            solnInfo.feasSplitVar[s], r,
            interPtsAndRays[s].cutRay[r],
            solnInfo.nonBasicVarIndex[interPtsAndRays[s].cutRay[r]],
            hplaneStore[hplane_ind].row,
            interPtsAndRays[s].newRay[r],
            solnInfo.nonBasicVarIndex[interPtsAndRays[s].newRay[r]],
            pointOrRayFlag);
      } else {
        fprintf(fptr, "%d,%d,-1,-1,-1,%d,%d,%d,",
            solnInfo.feasSplitVar[s], r,
            interPtsAndRays[s].newRay[r],
            solnInfo.nonBasicVarIndex[interPtsAndRays[s].newRay[r]],
            pointOrRayFlag);
      }

      fprintf(fptr, "%d,%f,", interPtsAndRays[s].cutRay[r],
          tmpVector[interPtsAndRays[s].cutRay[r]]);
      fprintf(fptr, "%d,%f,", interPtsAndRays[s].newRay[r],
          tmpVector[interPtsAndRays[s].newRay[r]]);

//      for (int k = 0; k < numElmts; k++) {
//        fprintf(fptr, "%d,%lf,", elmtIndices[k], elmtVals[k]);
//      }

      const int final_flag = interPtsAndRays[s].finalFlag[r];
      fprintf(fptr, "%d,", final_flag);
      fprintf(fptr, "\n");
    }
  }
}

/***********************************************************************/
/**
 * @brief Writes cut solver errors to separate file from instances info.
 */
void writeCutSolverErrors(int num_obj_tried, char* filename) {
  const bool fileExists = fexists(filename);
  FILE* CutErrorsOut = fopen(filename, "a");
  if (CutErrorsOut == NULL) {
    error_msg(errorstring,
        "Unable to open file to write time info for filename %s.\n",
        filename);
    writeErrorToII(errorstring, GlobalVariables::log_file);
    exit(1);
  }
  if (!fileExists) {
    fprintf(CutErrorsOut, "\n## Writing cut solver errors. ##\n");
    fprintf(CutErrorsOut, "%s%c", "OBJ", SEP);
    for (int fail_ind = 0; fail_ind < CutSolverFails::NUM_CUTSOLVER_FAILS;
        fail_ind++) {
      fprintf(CutErrorsOut, "%s%c", CutSolverFailName[fail_ind].c_str(), SEP);
    }
    fprintf(CutErrorsOut, "\n");
  }

  fprintf(CutErrorsOut, "%d%c", num_obj_tried, SEP);
  printf("OBJ TRIED: %d\n", num_obj_tried);
  for (int fail_ind = 0; fail_ind < CutSolverFails::NUM_CUTSOLVER_FAILS;
      fail_ind++) {
    fprintf(CutErrorsOut, "%d%c",
        GlobalVariables::numCutSolverFails[fail_ind], SEP);
    printf("%s: %d\n", CutSolverFailName[fail_ind].c_str(),
        GlobalVariables::numCutSolverFails[fail_ind]);
  }
  fclose(CutErrorsOut);
}

/***********************************************************************/
/**
 * Prints our the basis corresponding the COIN instance solver
 */
void printSimplexTableauWithNames(FILE* fptr, PHASolverInterface* const solver) {
  solver->enableFactorization();
  fprintf(fptr, "Optimal,simplex,tableau,of,P:\n");
  int numCols = solver->getNumCols();
  int numRows = solver->getNumRows();
  std::vector<double> a0(numRows);
  std::vector<int> tableauRowToVarIndex(numRows);
  solver->getBasics(&tableauRowToVarIndex[0]);
  const double* rowActivity = solver->getRowActivity();
  const double* rowRHS = solver->getRightHandSide();
  const double* colSoln = solver->getColSolution();

  std::vector<int> rstat(numRows), cstat(numCols);

  for (int i = 0; i < numRows; i++) {
    if (tableauRowToVarIndex[i] < numCols) {
      a0[i] = colSoln[tableauRowToVarIndex[i]];
    } else {
      a0[i] = -rowActivity[tableauRowToVarIndex[i] - numCols]
          + rowRHS[tableauRowToVarIndex[i] - numCols];
    }
  }

  int i, j;
  double* basisBasicRow = new double[numCols];
  double* basisSlackRow = new double[numRows];

  const char* rowSense = solver->getRowSense();

  fprintf(fptr, ",,");
  for (j = 0; j < numCols; j++)
    fprintf(fptr, "%d,", j);
  for (j = 0; j < numRows; j++)
    fprintf(fptr, "%d,", numCols + j);
  fprintf(fptr, "\n");
  fprintf(fptr, ",,");
  for (j = 0; j < numCols; j++)
    fprintf(fptr, "%s,", solver->getColName(j).c_str());
  for (j = 0; j < numRows; j++)
    fprintf(fptr, "%s,", solver->getRowName(j).c_str());
  fprintf(fptr, "LP solution (not RHS)\n");
  for (i = 0; i < numRows; i++) {
    fprintf(fptr, "%d,", i);
    if (tableauRowToVarIndex[i] < numCols)
      fprintf(fptr, "%s,",
          solver->getColName(tableauRowToVarIndex[i]).c_str());
    else
      fprintf(fptr, "%s,",
          solver->getRowName(tableauRowToVarIndex[i] - numCols).c_str());
    solver->getBInvARow(i, basisBasicRow, basisSlackRow);
    for (j = 0; j < numCols; j++)
      fprintf(fptr, "%2.4lf,", basisBasicRow[j]); //cout << basisBasicRow[j] << " ";
    for (j = 0; j < numRows; j++) {
      if (rowSense[j] == 'G')
        fprintf(fptr, "%2.4lf,", -basisSlackRow[j]); //cout << basisSlackRow[j] << " ";
      else {
        fprintf(fptr, "%2.4lf,", basisSlackRow[j]); //cout << basisSlackRow[j] << " ";
      }
    }
    //fprintf(fptr,"%1.2f ",primalSoln[i]);
    fprintf(fptr, "%2.4lf,", a0[i]);
    fprintf(fptr, "\n");

  }
  fprintf(fptr, ",Row0,");
  const double* redCostStruct = solver->getReducedCost();
  const double* redCostSlacks = solver->getRowPrice();
  solver->getBasisStatus(&cstat[0], &rstat[0]);
  for (j = 0; j < numCols; j++) {
    fprintf(fptr, "%2.4lf,", redCostStruct[j]);
    //    printf("%2.4lf\n", redCostStruct[j]);
  } //cout << basisBasicRow[j] << " ";
  for (j = 0; j < numRows; j++) {

    if (rowSense[j] == 'G')
      fprintf(fptr, "%2.4lf,", redCostSlacks[j]);
    else
      fprintf(fptr, "%2.4lf,", -redCostSlacks[j]);
  }
  fprintf(fptr, "%2.4lf\n", solver->getObjValue());
  // cout << "\n";
//  if (solverFactorizationStatus == 0) {
  solver->disableFactorization();
//  }
  delete[] basisBasicRow;
  delete[] basisSlackRow;
} /* printSimplexTableauWithNames */

/***********************************************************************/
/**
 * Prints our the basis corresponding the COIN instance solver; only non-basic non-artificial columns printed.
 */
void printNBSimplexTableauWithNames(FILE* fptr, PHASolverInterface* const solver) {
  //Print optimal basis
//  const int solverFactorizationStatus = solver->canDoSimplexInterface();
//  if (solverFactorizationStatus == 0) {
  solver->enableFactorization();
//  }
  fprintf(fptr, "NB, optimal,simplex,tableau,of,P:\n");
  int numCols = solver->getNumCols();
  int numRows = solver->getNumRows();
  std::vector<double> a0(numRows);
  std::vector<int> tableauRowToVarIndex(numRows);
  solver->getBasics(&tableauRowToVarIndex[0]);
  const double* rowActivity = solver->getRowActivity();
  const double* rowRHS = solver->getRightHandSide();
  const double* colSoln = solver->getColSolution();

  std::vector<int> rstat(numRows), cstat(numCols);
  solver->getBasisStatus(&cstat[0], &rstat[0]);

  for (int i = 0; i < numRows; i++) {
    if (tableauRowToVarIndex[i] < numCols) {
      a0[i] = colSoln[tableauRowToVarIndex[i]];
    } else {
      a0[i] = -rowActivity[tableauRowToVarIndex[i] - numCols]
          + rowRHS[tableauRowToVarIndex[i] - numCols];
    }
  }

  int i, j;
  double* basisBasicRow = new double[numCols];
  double* basisSlackRow = new double[numRows];
  fprintf(fptr, ",,");
  for (j = 0; j < numCols; j++) {
    if (cstat[j] == 2 || cstat[j] == 3)
      fprintf(fptr, "%d,", j);
  }
  for (j = 0; j < numRows; j++) {
    if ((rstat[j] == 2 || rstat[j] == 3) && solver->getRowSense()[j] != 'E')
      fprintf(fptr, "%d,", numCols + j);
  }
  fprintf(fptr, "\n");
  fprintf(fptr, ",,");
  for (j = 0; j < numCols; j++)
    if (cstat[j] == 2 || cstat[j] == 3)
      fprintf(fptr, "%s,", solver->getColName(j).c_str());
  for (j = 0; j < numRows; j++)
    if ((rstat[j] == 2 || rstat[j] == 3) && solver->getRowSense()[j] != 'E')
      fprintf(fptr, "%s,", solver->getRowName(j).c_str());
  fprintf(fptr, "RHS\n");
  for (i = 0; i < numRows; i++) {
    fprintf(fptr, "%d,", i);
    if (tableauRowToVarIndex[i] < numCols)
      fprintf(fptr, "%s,",
          solver->getColName(tableauRowToVarIndex[i]).c_str());
    else
      fprintf(fptr, "%s,",
          solver->getRowName(tableauRowToVarIndex[i] - numCols).c_str());
    solver->getBInvARow(i, basisBasicRow, basisSlackRow);
    for (j = 0; j < numCols; j++)
      if (cstat[j] == 2 || cstat[j] == 3)
        fprintf(fptr, "%2.7lf,", basisBasicRow[j]); //cout << basisBasicRow[j] << " ";
    for (j = 0; j < numRows; j++)
      if ((rstat[j] == 2 || rstat[j] == 3)
          && solver->getRowSense()[j] != 'E')
        fprintf(fptr, "%2.7lf,", basisSlackRow[j]); //cout << basisSlackRow[j] << " ";
    //fprintf(fptr,"%1.2f ",primalSoln[i]);
    fprintf(fptr, "%2.7lf,", a0[i]);
    fprintf(fptr, "\n");

  }
  fprintf(fptr, ",Row0,");
  const double* redCostStruct = solver->getReducedCost();
  const double* redCostSlacks = solver->getRowPrice();
  const char* rowSense = solver->getRowSense();
  for (j = 0; j < numCols; j++) {
    if (cstat[j] == 2 || cstat[j] == 3)
      fprintf(fptr, "%2.7lf,", redCostStruct[j]);
    //    printf("%2.4lf\n", redCostStruct[j]);
  } //cout << basisBasicRow[j] << " ";
  for (j = 0; j < numRows; j++) {
    if ((rstat[j] == 2 || rstat[j] == 3)
        && solver->getRowSense()[j] != 'E') {
      if (rowSense[j] == 'G')
        fprintf(fptr, "%2.7lf,", redCostSlacks[j]);
      else
        fprintf(fptr, "%2.7lf,", -redCostSlacks[j]);
    }
  }
  fprintf(fptr, "%2.7lf\n", solver->getObjValue());
  // cout << "\n";
//  if (solverFactorizationStatus == 0) {
  solver->disableFactorization();
//  }
  delete[] basisBasicRow;
  delete[] basisSlackRow;
} /* printNBSimplexTableauWithNames */

/***********************************************************************/
/**
 * @brief Writes solution (both primal and dual) to file.
 * @param solver          ::  The solver whose solution we are writing.
 * @param soln_name       ::  The name we are outputting to.
 * @param out_stream      ::  Out stream.
 * @param width_varname   ::  Width of variable name field.
 * @param width_val       ::  Width of variable value field.
 * @param width_activity  ::  Width of row activity field.
 * @param after_dec       ::  Number of digits after the decimal.
 * @param EXT             ::  Extension added to end of file output.
 */
void writeSoln(const PHASolverInterface* const solver, std::string soln_name,
    int width_varname, int width_val, int width_stat, int width_rc,
    int width_activity, int width_rhs, int after_dec,
    const std::string EXT) {
#ifdef TRACE
  printf("\n## Writing solution for %s. ##\n", soln_name.c_str());
#endif
  FILE *myfile;
  myfile = fopen((soln_name + EXT).c_str(), "w");
  if (myfile != NULL) {
    std::vector<int> cstat(solver->getNumCols()), rstat(
        solver->getNumRows());
    solver->getBasisStatus(&cstat[0], &rstat[0]);

    fprintf(myfile, "Obj: %s with value %f.\n",
        solver->getObjName().c_str(), solver->getObjValue());
    fprintf(myfile, "\n");

    fprintf(myfile, "%-*s\t%-*s\t%-*s\t%-*s\t\n", width_varname,
        "Variable name", width_val, "Value", width_stat, "cstat",
        width_rc, "Reduced Cost");
    for (int col = 0; col < solver->getNumCols(); col++) {
      fprintf(myfile, "%-*.*s\t% -*.*f\t%-*d\t%-*f\t\n", width_varname,
          width_varname, solver->getColName(col).c_str(), width_val,
          after_dec, solver->getColSolution()[col], width_stat,
          cstat[col], width_rc, solver->getReducedCost()[col]);
    }
    fprintf(myfile, "\n");

    fprintf(myfile, "%-*s\t%-*s\t%-*s\t%-*s\t%-*s\t%-*s\t\n", width_varname,
        "Row name", width_val, "Dual Value", width_activity, "Activity",
        width_rhs, ">=", width_rhs, "<=", width_stat, "rstat");
    for (int row = 0; row < solver->getNumRows(); row++) {
      char rowsense = solver->getRowSense()[row];
      if (rowsense == 'E' || rowsense == 'R')
        fprintf(myfile,
            "%-*.*s\t% -*.*f\t% -*.*f\t% -*.*f\t% -*.*f\t%-*d\t\n",
            width_varname, width_varname,
            solver->getRowName(row).c_str(), width_val, after_dec,
            solver->getRowPrice()[row], width_activity, after_dec,
            solver->getRowActivity()[row], width_rhs, after_dec,
            solver->getRowLower()[row], width_rhs, after_dec,
            solver->getRowUpper()[row], width_stat, rstat[row]);
      else if (rowsense == 'G')
        fprintf(myfile,
            "%-*.*s\t% -*.*f\t% -*.*f\t% -*.*f\t%-*s\t%-*d\t\n",
            width_varname, width_varname,
            solver->getRowName(row).c_str(), width_val, after_dec,
            solver->getRowPrice()[row], width_activity, after_dec,
            solver->getRowActivity()[row], width_rhs, after_dec,
            solver->getRowLower()[row], width_rhs, "", width_stat,
            rstat[row]);
      else if (rowsense == 'L')
        fprintf(myfile,
            "%-*.*s\t% -*.*f\t% -*.*f\t%-*s\t% -*.*f\t%-*d\t\n",
            width_varname, width_varname,
            solver->getRowName(row).c_str(), width_val, after_dec,
            solver->getRowPrice()[row], width_activity, after_dec,
            solver->getRowActivity()[row], width_rhs, "", width_rhs,
            after_dec, solver->getRowUpper()[row], width_stat,
            rstat[row]);
      else
        fprintf(myfile, "%-*.*s\t% -*.*f\t% -*.*f\t%-*s\t%-*s\t%-*d\t\n",
            width_varname, width_varname,
            solver->getRowName(row).c_str(), width_val, after_dec,
            solver->getRowPrice()[row], width_activity, after_dec,
            solver->getRowActivity()[row], width_rhs, "?",
            width_rhs, "?", width_stat, rstat[row]);
    }
    fclose(myfile);
  }
} /* writeSoln */

double rayPostPivot(const int component, const int nb_var_in1,
    const int struct_var_out1, const int nb_var_in2,
   // const PHASolverInterface* const solver, 
    const SolutionInfo& solnInfo) {
  const int hplane_row = solnInfo.rowOfVar[struct_var_out1];

  if (hplane_row < 0) { // This has to be the same var as nb_var_in1
    if (nb_var_in2 == nb_var_in1) {
      if (component == 1)
        return -1.0;
      else
        return 0.0;
    } else {
      if (component == 1)
        return 0.0;
      else
        return 1.0;
    }
  }

  if (component == 2)
    return 1.0;
  else
    return -1 * solnInfo.raysOfC1[nb_var_in2][hplane_row]
        / solnInfo.raysOfC1[nb_var_in1][hplane_row];
}

double vertexPostSinglePivot(const int component, const double rayDirnToHplane1,
    const double init_val_out1, const double final_val_out1) {
  if (component == 2)
    return 0.0;

  const double numer = final_val_out1 - init_val_out1;
  const double denom = rayDirnToHplane1;

  return numer / denom;
}

double pointPostDoublePivot(const int component, const double rayDirnToHplane1,
    const double init_val_out1, const double final_val_out1,
    const double rayDirnToSplit1, const double rayDirnToSplit2,
    const double rayDirnToHplane2, const double init_val_out2,
    const double final_val_out2) {
  double numer;
  if (component == 1) {
    numer = init_val_out1 * rayDirnToSplit2
        - init_val_out2 * rayDirnToHplane2
        - final_val_out1 * rayDirnToSplit2
        + final_val_out2 * rayDirnToHplane2;
  } else {
    numer = -init_val_out1 * rayDirnToSplit1
        + init_val_out2 * rayDirnToHplane1
        + final_val_out1 * rayDirnToSplit1
        - final_val_out2 * rayDirnToHplane1;
  }
  const double denom = rayDirnToSplit1 * rayDirnToHplane2
      - rayDirnToSplit2 * rayDirnToHplane1;

  return numer / denom;
}

double pointPostPivot(const int component, const int nb_var_in1,
    const int struct_var_out1, const int nb_var_in2,
    const int struct_var_out2, const PHASolverInterface* const solver,
    const SolutionInfo& solnInfo, const bool computeVertexFlag) {
  // If we are only pivoting in one variable, and not computing a vertex,
  // then we need to compute the intersection point of this ray with
  // the boundary of the split (which will be the struct_var_out1)
  if (nb_var_in2 < 0) {
    if (component == 2) {
      return 0.0; // No second nb component
    }

    if (!computeVertexFlag) {
      const int split_row = solnInfo.rowOfVar[struct_var_out1];
      const double split_val = solnInfo.a0[split_row];
      const double rayDirnToSplit =
          solnInfo.raysOfC1[nb_var_in1][split_row];

      double bound = 0.0;
      if (isZero(rayDirnToSplit, param.getRAYEPS())) {
        return solver->getInfinity();
      } else if (greaterThanVal(rayDirnToSplit, 0.0, param.getRAYEPS())) {
        bound = std::ceil(split_val);
      } else {
        bound = std::floor(split_val);
      }

      return (bound - split_val) / rayDirnToSplit;
    }
  }

  const int hplane_row = solnInfo.rowOfVar[struct_var_out1];

  // If the ray direction is positive, we intersect the hyperplane ub,
  // if it is negative, we intersect the hyperplane lb,
  // and otherwise, we do not intersect it at all
  double init_val_out1 = +0.0;
  double final_val_out1 = +0.0;

  const double rayDirnToHplane1 =
      (hplane_row >= 0) ? solnInfo.raysOfC1[nb_var_in1][hplane_row] : 1.0;

  // First we check if a non-basic var is being pivoted out
  // (in which case ray will intersect the
  // "other" bound, since we are assuming all non-basic vars are at lb = 0)
  if (hplane_row < 0) {
    const double UB = getVarUB(solver, struct_var_out1);
    const double LB = getVarLB(solver, struct_var_out1);
    if (!isInfinity(UB, solver->getInfinity())
        && !isNegInfinity(LB, solver->getInfinity())) {
      final_val_out1 = (UB - LB);

      if (component == 1) // Same if nb_var_in2 < 0 or not
        return final_val_out1;
    } else {
      return solver->getInfinity();
    }
  } else {
    // Otherwise it is a basic variable being pivoted out
    init_val_out1 = solnInfo.a0[hplane_row];

    // If it is a slack variable, the ray better have a negative direction
    // (There is only one bound on slacks, which is a lower bound)
    if (struct_var_out1 >= solver->getNumCols()) {
      if (!lessThanVal(rayDirnToHplane1, 0.0, param.getRAYEPS())) {
        return solver->getInfinity();
      }
      final_val_out1 = +0.0;
    } else { // Otherwise it is a structural variable being pivoted out
      if (isZero(rayDirnToHplane1, param.getRAYEPS())) {
        return solver->getInfinity();
      } else if (greaterThanVal(rayDirnToHplane1, 0.0, param.getRAYEPS())) { // Ray intersects ub
        final_val_out1 = getVarUB(solver, struct_var_out1);
      } else { // Ray intersects lb
        final_val_out1 = getVarLB(solver, struct_var_out1);
      }
    }
  }

  // Check if this is after a single pivot and we want the vertex value
  if (nb_var_in2 < 0 || computeVertexFlag) {
    return vertexPostSinglePivot(component, rayDirnToHplane1, init_val_out1,
        final_val_out1);
  }

  // Get values for the second pivot
  const int split_row = solnInfo.rowOfVar[struct_var_out2]; // Note: this HAS to be a non-negative number

  // If the ray direction is positive, we intersect the split ceiling,
  // if it is negative, we intersect the split floor,
  // and otherwise, we do not intersect it at all
  double init_val_out2 = solnInfo.a0[split_row];

  const double rayDirnToSplit1 = solnInfo.raysOfC1[nb_var_in1][split_row];
  const double rayDirnToSplit2 = solnInfo.raysOfC1[nb_var_in2][split_row];

  // If the first ray intersected its own bound, then the next ray,
  // (which will necessarily be a different ray) will not change the first
  // component at all (handled above), and the second
  // component will be as in the "one variable pivoted out" case
  if (hplane_row < 0) {
    const double split_val = solnInfo.a0[split_row]
        + final_val_out1 * rayDirnToSplit1; // Note: must not go past an integer value

    double bound = 0.0;
    if (isZero(rayDirnToSplit2, param.getRAYEPS())) {
      return solver->getInfinity();
    } else if (greaterThanVal(rayDirnToSplit2, 0.0, param.getRAYEPS())) {
      bound = std::ceil(split_val);
    } else {
      bound = std::floor(split_val);
    }

    return (bound - split_val) / rayDirnToSplit2;
  }

  const double rayDirnToHplane2 = solnInfo.raysOfC1[nb_var_in2][hplane_row]; // Note: we know hplane_row >= 0

  const double checkNewRayDirnToSplit = rayDirnToSplit2
      - rayDirnToHplane2 * rayDirnToSplit1 / rayDirnToHplane1;

  double final_val_out2;
  if (isZero(checkNewRayDirnToSplit, param.getRAYEPS())) {
    return solver->getInfinity();
  } else if (greaterThanVal(checkNewRayDirnToSplit, 0.0, param.getRAYEPS())) {
    final_val_out2 = std::ceil(init_val_out2);
  } else {
    final_val_out2 = std::floor(init_val_out2);
  }

  return pointPostDoublePivot(component, rayDirnToHplane1, init_val_out1,
      final_val_out1, rayDirnToSplit1, rayDirnToSplit2, rayDirnToHplane2,
      init_val_out2, final_val_out2);
}


/** Assume sorted vector **/
void packedShallowVectorSum(CoinPackedVector& sum, const double mult1,
    const CoinShallowPackedVector& vec1, const double mult2,
    const CoinShallowPackedVector& vec2) {
  //CoinShallowPackedVector sum;
  if (isZero(mult1)) {
    if (isZero(mult2)) {
      return;
    }
    else {
      sum = vec2;
      return;
    }
  } else if (isZero(mult2)) {
    sum = vec1;
    return;
  }

  const int numVec1Elems = vec1.getNumElements();
  const int numVec2Elems = vec2.getNumElements();
  const int* vec1Index = vec1.getIndices();
  const int* vec2Index = vec2.getIndices();
  const double* vec1Val = vec1.getElements();
  const double* vec2Val = vec2.getElements();

  std::vector<int> sumIndex;
  std::vector<double> sumVal;

  int vec2_el = 0;
  for (int el = 0; el < numVec1Elems; el++) {
    const int vec1_ind = vec1Index[el];
    const double vec1_val = vec1Val[el];

    double vec2_val = 0.0;
    while (vec2_el < numVec2Elems) {
      const int vec2_ind = vec2Index[vec2_el];
      if (vec2_ind < vec1_ind) {
        if (!isZero(vec2Val[vec2_el])) {
          sumIndex.push_back(vec2_ind);
          sumVal.push_back(mult2 * vec2Val[vec2_el]);
        }
        vec2_el++;
      } else if (vec2_ind == vec1_ind) {
        vec2_val = vec2Val[vec2_el];
        vec2_el++;
        break;
      } else {
        vec2_el++;
        break;
      }
    }

    const double new_val = mult1 * vec1_val + mult2 * vec2_val;
    if (!isZero(new_val)) {
      sumIndex.push_back(vec1_ind);
      sumVal.push_back(new_val);
    }
  }

  while (vec2_el < numVec2Elems) {
    if (!isZero(vec2Val[vec2_el])) {
      sumIndex.push_back(vec2Index[vec2_el]);
      sumVal.push_back(mult2 * vec2Val[vec2_el]);
    }
    vec2_el++;
  }

  sum.setVector(sumIndex.size(), sumIndex.data(), sumVal.data(), false);
} /* packedSumVector */

/** Assume sorted vector **/
double packedSum(const CoinPackedVector& vec1,
    const CoinPackedVector& vec2) {
  const int numVec1Elems = vec1.getNumElements();
  const int numVec2Elems = vec2.getNumElements();
  const double* vec1Val = vec1.getElements();
  const double* vec2Val = vec2.getElements();

  double sum = 0.0;

  for (int el = 0; el < numVec1Elems; el++) {
    sum += vec1Val[el];
  }
  for (int el = 0; el < numVec2Elems; el++) {
    sum += vec2Val[el];
  }

  return sum;
}

double packedDotProduct(const CoinPackedVector& vec1,
    const CoinPackedVector& vec2) {
  const int numVec1Elems = vec1.getNumElements();
  const int numVec2Elems = vec2.getNumElements();
  const int* vec1Index = vec1.getIndices();
  const int* vec2Index = vec2.getIndices();
  const double* vec1Val = vec1.getElements();
  const double* vec2Val = vec2.getElements();

  double dot_prod = 0.0;

  int vec2_el = 0;
  for (int el = 0; el < numVec1Elems; el++) {
    const int vec1_ind = vec1Index[el];
    const double vec1_val = vec1Val[el];

    double vec2_val = 0.0;
    if (vec2_el < numVec2Elems) {
      int vec2_ind = vec2Index[vec2_el];

      while (vec2_ind < vec1_ind) {
        vec2_el++;
        if (vec2_el < numVec2Elems) {
          vec2_ind = vec2Index[vec2_el];
        } else {
          break;
        }
      }

      if (vec2_ind == vec1_ind) {
        vec2_val = vec2Val[vec2_el];
      }
    }

    if (!isZero(vec1_val) && !isZero(vec2_val))
      dot_prod += vec1_val * vec2_val;
  }

  return dot_prod;
}

/********************************************************************************/
/**
 * Check if two solvers have the same everything:
 * 1. colSolution
 * 2. rowPrice
 * 3. slackVals
 * 4. cstat and rstat
 * 5. reducedCosts
 */
bool solversEqual(const PHASolverInterface* const clp1, const PHASolverInterface* const clp2) {
  int numrows1 = clp1->getNumRows();
  int numcols1 = clp1->getNumCols();
  int numrows2 = clp2->getNumRows();
  int numcols2 = clp2->getNumCols();

  // First check that numrows and numcols are the same
  if (numrows1 != numrows2 || numcols1 != numcols2) {
#ifdef TRACE
    printf("*** ERROR: numrows (%d vs %d) or numcols (%d vs %d) not same.",
        numrows1, numrows2, numcols1, numcols2);
#endif
    return false;
  }

  // Set up cstat, rstat
  std::vector<int> cstat1(numcols1), rstat1(numrows1);
  std::vector<int> cstat2(numcols2), rstat2(numrows2);

  clp1->getBasisStatus(&cstat1[0], &rstat1[0]);
  clp2->getBasisStatus(&cstat2[0], &rstat2[0]);

  for (int col = 0; col < numcols1; col++) {
    double val1 = clp1->getColSolution()[col];
    double val2 = clp2->getColSolution()[col];
    double rc1 = clp1->getReducedCost()[col];
    double rc2 = clp2->getReducedCost()[col];
    if ((std::abs(val1 - val2) > param.getEPS()) || (std::abs(rc1 - rc2) > param.getEPS())
        || (cstat1[col] != cstat2[col])) {
      char errorstring[300];
      snprintf(errorstring, sizeof(errorstring) / sizeof(char),
          "*** ERROR: New col values not the same as old values. Variable %d:\n\t"
              "val1: %f \t val2: %f \t Diff: %e\n\t"
              "rc1: %f \t rc2: %f \t Diff: %e\n\t"
              "cstat1: %d \t cstat2: %d\n", col, val1, val2,
          val1 - val2, rc1, rc2, rc1 - rc2, cstat1[col], cstat2[col]);
      std::cerr << errorstring << std::endl;
//      writeErrorToII(errorstring, inst_info_out);
      return false;
    }
  }
  for (int row = 0; row < numrows1; row++) {
    double sv1 = std::abs(
        clp1->getRightHandSide()[row] - clp1->getRowActivity()[row]);
    double sv2 = std::abs(
        clp2->getRightHandSide()[row] - clp2->getRowActivity()[row]);
    double rp1 = clp1->getRowPrice()[row];
    double rp2 = clp2->getRowPrice()[row];
    if ((std::abs(sv1 - sv2) > param.getEPS()) || (std::abs(rp1 - rp2) > param.getEPS())
        || (rstat1[row] != rstat2[row])) {
      char errorstring[300];
      snprintf(errorstring, sizeof(errorstring) / sizeof(char),
          "*** ERROR: New row values not the same as old values. Variable %d:\n\t"
              "sv1: %f \t sv2: %f \t Diff: %e\n\t"
              "rp1: %f \t rp2: %f \t Diff: %e\n\t"
              "rstat1: %d \t rstat2: %d\n", row, sv1, sv2,
          sv1 - sv2, rp1, rp2, rp1 - rp2, rstat1[row], rstat2[row]);
      std::cerr << errorstring << std::endl;
//      writeErrorToII(errorstring, inst_info_out);
      return false;
    }

  }
  return true;
}
