/*
 * hplaneActivation.cpp
 *
 *  Created: May 25, 2014
 *  Edited: November 30, 2015
 *  Author: akazachk (based on snadarajah work)
 */

#include <hplaneActivation.hpp>
#include "GlobalConstants.hpp"
#include "Utility.hpp"

/**
 * @return rayIndicesIntByHplane [vertex][ray ind] std::vector keeping all rays cut by this hplane
 */
void addRaysIntByHplaneBeforeBdS(std::vector<int>& vertOfRayIndicesIntByHplane,
    std::vector<int>& rayIndicesIntByHplane,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
    const Hplane& hplane,
    const std::vector<Vertex>& vertexStore, const int init_vert_ind,
    const std::vector<Ray>& rayStore, const int split_ind) {
  // For the vertex, check the rays emanating from it (that still remain to be cut)
  const int numRaysForVertex = vertexStore[init_vert_ind].rayIndices.size();
  for (int r = 0; r < numRaysForVertex; r++) {
    if (!vertexStore[init_vert_ind].rayActiveFlag[r]) {
      continue;
    }

    const int ray_ind = vertexStore[init_vert_ind].rayIndices[r];
    const int split_var = solnInfo.feasSplitVar[split_ind];

    const double distToSplit = distanceToSplit(vertexStore[init_vert_ind],
        rayStore[ray_ind], split_var, solver, solnInfo);

//    const double distToHplane = hplane.distanceToHplane(
//        vertexStore[init_vert_ind], rayStore[ray_ind], solver,
//        true,
//        &(solnInfo.nonBasicVarIndex));
    const double distToHplane = hplane.distanceToHplane(
        vertexStore[init_vert_ind], rayStore[ray_ind], solver,
        true, solnInfo);

//      if (distToHplane >= solver->getInfinity() - param.getEPS()) {
    if (isInfinity(distToHplane)) {
      continue; // hplane does not intersect this ray
//      } else if (distToHplane <= distToSplit - param.getEPS()) {
    } else if (lessThanVal(distToHplane, distToSplit, CoinMin(solnInfo.EPS, param.getEPS()))) {
      vertOfRayIndicesIntByHplane.push_back(init_vert_ind);
      rayIndicesIntByHplane.push_back(ray_ind);
    }
  }
} /* addRaysIntByHplaneBeforeBdS */

int removePtsAndRaysCutByHplane(const PHASolverInterface* const solver,
    const SolutionInfo &solnInfo, const std::vector<bool> &raysCutFlag,
    const int split_var, const Hplane& hplane, const int hplane_ind_in_split,
//    const int hplane_ind, const std::vector<Hplane>& hplaneStore,
    const std::vector<Vertex>& vertexStore, const std::vector<Ray>& rayStore,
    IntersectionInfo& interPtsAndRays, const AdvCut* NBSIC, const double norm,
    double& sumRemovedDepth, const bool fake_act) {
  std::vector<int> delIndices;

  const double* UB = solver->getColUpper();
  const double* LB = solver->getColLower();

  // Only remove points that correspond to rays that are cut by this hplane
  // Also the points should not be those created by hplanes activated "later"
  for (int row_ind = 0; row_ind < interPtsAndRays.getNumRows(); row_ind++) {
    const int curr_hplane_ind_in_split =
        interPtsAndRays.hplaneIndexForSplit[row_ind];
    if (curr_hplane_ind_in_split >= hplane_ind_in_split) {
      continue;
    }

    const int ray_ind =
        (interPtsAndRays.cutRay[row_ind] >= 0) ?
            interPtsAndRays.cutRay[row_ind] :
            interPtsAndRays.newRay[row_ind];
    if (!raysCutFlag[ray_ind]) {
      continue;
    }
    const bool currRowIsRay = isZero(interPtsAndRays.RHS[row_ind]);

    double SIC_activity = +0.0;
    //This if condition ensures that intersection points are removed only
    //if they emanate from a vertex of a ray that is cut
    // In PHA 1.1, this should never be true, since we never cut a ray twice

    //      if ((*intPtOrRayInfo[r]).cutRayColIndx >= 0) { // Uncomment this if we want to remove all points violating hplane
    const CoinShallowPackedVector tmpVector = interPtsAndRays.getVector(
        row_ind);
    const int numElmts = tmpVector.getNumElements();
    const int* elmtIndices = tmpVector.getIndices();
    const double* elmtVals = tmpVector.getElements();
    double activity = 0.0;
    for (int e = 0; e < numElmts; e++) {
      // Either hplaneRowIndx >= 0, in which case we activated a basic var
      // or it is -1, and we activated an upper or lower bound of a nb var
      if (hplane.row >= 0) {
        activity += elmtVals[e]
            * solnInfo.raysOfC1[elmtIndices[e]][hplane.row];
      } else {
        const int nb_ind = (-1 - solnInfo.rowOfVar[hplane.var]);
        if (elmtIndices[e] == nb_ind)
          activity = elmtVals[e];
      }
//      SIC_activity += NBSIC[elmtIndices[e]] * elmtVals[e];
      // NBSIC dist recalculated here to avoid some numerical issues from taking inverses
      SIC_activity += elmtVals[e]
          / distanceToSplit(vertexStore[0], rayStore[elmtIndices[e]],
              split_var, solver, solnInfo);
    }

    // Note that if hplaneRowIndx == -1, then we activated a bound of a nb var,
    // and it can't be a bound of a slack variable, since they are at zero in the cobasis,
    // and that is their only bound.
    if (hplane.var >= solnInfo.numCols) {
      // Recall that all slack variables are non-negative in this
      // setting, since they were complemented in SolutionInfo
      if (!currRowIsRay
          && lessThanVal(solnInfo.a0[hplane.row] + activity, 0.0)) {
        delIndices.push_back(row_ind);
      } else if (currRowIsRay && lessThanVal(activity, 0.0)) {
        delIndices.push_back(row_ind);
      }
    } else {
      // If this was a basic variable, no bounds were changed
      // Otherwise things may have shifted
      double lb, ub, newValue;
      if (hplane.row >= 0) {
        lb = LB[hplane.var];
        ub = UB[hplane.var];
        newValue = solnInfo.a0[hplane.row] + activity;
      } else {
        lb = 0.0;
        ub = UB[hplane.var] - LB[hplane.var];
        newValue = activity;
      }
      if (hplane.ubflag == static_cast<int>(HplaneBoundFlag::toUB)) {
        if (!currRowIsRay && greaterThanVal(newValue, ub)) {
          delIndices.push_back(row_ind);
        } else if (currRowIsRay && greaterThanVal(activity, 0.0)) {
          delIndices.push_back(row_ind);
        }
      } else {
        if (!currRowIsRay && lessThanVal(newValue, lb)) {
          delIndices.push_back(row_ind);
        } else if (currRowIsRay && lessThanVal(activity, 0.0)) {
          delIndices.push_back(row_ind);
        }
      }
    }

    // If we are deleting this point, then count it in deleted depth calculations
    if (!currRowIsRay && (delIndices.size() > 0)
        && (delIndices[delIndices.size() - 1] == row_ind)) {
      double curr_normed_viol = (NBSIC) ? (SIC_activity - NBSIC->rhs()) / norm : 0.0;
      if (isZero(curr_normed_viol)) {
        curr_normed_viol = +0.0;
      }
      sumRemovedDepth += curr_normed_viol;
      //        cout << "sumRemoved" << ": " << sumRemovedDepth << endl;

      if (lessThanVal(curr_normed_viol, 0.0)) {
        error_msg(errstr,
            "Depth %f is negative! Point cannot be above the SIC.\n",
            curr_normed_viol);
        writeErrorToII(errstr, GlobalVariables::log_file);
        exit(1);
      }

    }
  }

  if (!fake_act) {
    // Sort deleted indices
    sort(delIndices.begin(), delIndices.end());
    interPtsAndRays.deleteRowInfo(delIndices.size(), &delIndices[0]);
  }
  return delIndices.size();
} /* removePtsAndRaysCutByHplane */

/**
 * New version, 06/07/2016
 * The ray ray_ind is being cut, creating the new vertex tmpVertex
 */
void updateHplaneStore(std::vector<IntersectionInfo>& interPtsAndRays,
    std::vector<Vertex>& vertexStore, const int init_vert_ind,
    std::vector<Ray>& rayStore,
//    std::vector<std::vector<int> >& allRaysCut,
    std::vector<int>& allRaysCutAcrossSplits,
//    std::vector<std::vector<std::vector<int> > >& raysToBeCut,
    std::vector<std::vector<int> >& raysCutForSplit,
    std::vector<Hplane>& hplaneStore,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
    const Vertex& newVertex, const Hplane& tmpHplane,
//    const std::vector<int>& splitsAffByHplane,
    const int split_ind,
    const int ray_ind,
    const int act_ind, const bool checkIfRayPrevCut) {

//  // If the hyperplane activated affects at least one split, then this ray is cut for some split.
//  // (I.e., if ||splitsEffByHplane|| = 0, then for all splits,
//  // the ray intersects the boundary of the split before it intersects
//  // the hplane chosen.)
//  if (splitsAffByHplane.size() == 0) {
//    return; // If hplane affects no splits, then we do not activate this hplane
//  }

  // The current ray will be cut by a hplane prior to bd S (for at least one split)
  bool rayPrevCutAcrossSplits = false;
  if (checkIfRayPrevCut) {
    rayPrevCutAcrossSplits = (find_val(ray_ind, allRaysCutAcrossSplits) >= 0);
  }
  if (!rayPrevCutAcrossSplits) {
    allRaysCutAcrossSplits.push_back(ray_ind);
  }

  // Add the new hplane
  const int oldHplaneIndex = hasHplaneBeenActivated(hplaneStore, tmpHplane);
  const int newHplaneIndex =
      (oldHplaneIndex == -1) ? hplaneStore.size() : oldHplaneIndex;
  if (oldHplaneIndex == -1) {
    hplaneStore.push_back(tmpHplane);
    if (hplaneStore[newHplaneIndex].rayToBeCutByHplane.size() > 0) {
      for (int s = 0; s < (int) hplaneStore[newHplaneIndex].rayToBeCutByHplane.size(); s++) {
        hplaneStore[newHplaneIndex].rayToBeCutByHplane[s].clear();
      }
    }
    hplaneStore[newHplaneIndex].resizeVectors(solnInfo.numFeasSplits);
  }

  // Add the new vertex
  const int newVertexIndex = addVertex(vertexStore, newVertex);
//      newHplaneIndex);

//  // Store that this ray is "first" cut by the proper hplane / vertex index
//  rayStore[ray_ind].firstHplaneOverall = newHplaneIndex;
//  rayStore[ray_ind].firstVertexOverall = newVertexIndex;

  // Only iterate over splits where the hyperplane intersects the ray before the bdS
//  for (unsigned s1 = 0; s1 < splitsAffByHplane.size(); s1++) {
//    const int split_ind = splitsAffByHplane[s1];

    // This ray is cut for this split
    bool rayPrevCutForSplit = false;
    if (checkIfRayPrevCut) {
      rayPrevCutForSplit = (find_val(ray_ind, raysCutForSplit[split_ind]) >= 0);
    }
    if (!rayPrevCutForSplit) {
      raysCutForSplit[split_ind].push_back(ray_ind);
    }

    // Stores only those rays of C1 that have to be cut by the new hplane for validity of the procedure
    const int prevTotalNumRaysCutByHplane =
        hplaneStore[newHplaneIndex].rayToBeCutByHplane[split_ind].size();
    hplaneStore[newHplaneIndex].rayToBeCutByHplane[split_ind].push_back(
        ray_ind); // Clearly this ray, r^j, is cut
    hplaneStore[newHplaneIndex].vertOfRayToBeCutByHplane[split_ind].push_back(
        init_vert_ind);
    hplaneStore[newHplaneIndex].vertCreatedByRayToBeCutByHplane[split_ind].push_back(
        newVertexIndex);

    // Check if activated hyperplane is new or has been previously activated
    // The vector actHplaneIndexForSplit is needed for performing activations later
    const int oldHplaneForActAndSplit = hasHplaneBeenActivated(hplaneStore,
        interPtsAndRays[split_ind].hplaneIndexByAct[act_ind], tmpHplane, true);
    const int oldHplaneForSplitOverall = hasHplaneBeenActivated(hplaneStore,
        interPtsAndRays[split_ind].allHplaneToAct, tmpHplane, true);

    // If the hplane is new for this split, we need to add it to activated hplanes
    // And we need to get the other rays affected by this hplane for this split
    if (oldHplaneForSplitOverall == -1) {
      // Store the index of this hplane in hplanesToAct
      interPtsAndRays[split_ind].hplaneIndexByAct[act_ind].push_back(
          newHplaneIndex);
      interPtsAndRays[split_ind].indexOfHplaneInOverallOrderForSplit[act_ind].push_back(
          interPtsAndRays[split_ind].allHplaneToAct.size());
      interPtsAndRays[split_ind].allHplaneToAct.push_back(newHplaneIndex);
      // This stores *all* the rays of C1 intersected by the new hyperplane before bd S
      addRaysIntByHplaneBeforeBdS(
          hplaneStore[newHplaneIndex].vertOfAllRaysIntByHplane[split_ind],
          hplaneStore[newHplaneIndex].allRaysIntByHplane[split_ind], solver,
          solnInfo, hplaneStore[newHplaneIndex], vertexStore, init_vert_ind,
          rayStore, split_ind);

      // We need to include in R* (rays cut by the hplane) all rays previously cut
      // This is sufficient for the validity of the procedure
      // TODO might be better to sort here instead of creating this rayIntersectedFlag
      std::vector<bool> rayIntersectedFlag(solnInfo.numNB, false);
      const int numRaysIntByHplane =
          hplaneStore[newHplaneIndex].allRaysIntByHplane[split_ind].size();
      for (int r = 0; r < numRaysIntByHplane; r++) {
        rayIntersectedFlag[hplaneStore[newHplaneIndex].allRaysIntByHplane[split_ind][r]] =
            true;
      }

      // The hiccup is that this hyperplane may have been included for the split on a previous activation
      // Some of the rays cut by other hyperplanes may already be cut by this hyperlane previously
//      // We remedy this by ignoring all rays cut by this hyperplane in a previous activation
//      const int numRaysCutByHplane =
//          hplaneStore[newHplaneIndex].rayToBeCutByHplane[split_ind].size();
//      for (int r = 0; r < numRaysCutByHplane; r++) {
//        rayIntersectedFlag[hplaneStore[newHplaneIndex].rayToBeCutByHplane[split_ind][r]] =
//            false;
//      }

      // For all rays previously cut for this split, check whether the hplane intersects it
      for (int r = 0; r < (int) raysCutForSplit[split_ind].size() - 1; r++) {
        const int curr_ray = raysCutForSplit[split_ind][r];
        if (rayIntersectedFlag[curr_ray]) {
          newVertexHelper(vertexStore, init_vert_ind, rayStore,
              hplaneStore, solver, solnInfo,
              curr_ray, split_ind, newHplaneIndex, true);
        }
      }
    } else { /* oldHplaneForSplitOverall != -1 */
      // Insert this hplane in the proper order in hplaneIndexByAct
      // Also need to update indexOfHplaneInOverallOrder
      if (oldHplaneForActAndSplit == -1) {
        std::vector<int>::iterator it =
            std::lower_bound(
                interPtsAndRays[split_ind].indexOfHplaneInOverallOrderForSplit[act_ind].begin(),
                interPtsAndRays[split_ind].indexOfHplaneInOverallOrderForSplit[act_ind].end(),
                oldHplaneForSplitOverall);
        const int tmp_h =
            (it
                - interPtsAndRays[split_ind].indexOfHplaneInOverallOrderForSplit[act_ind].begin());

        interPtsAndRays[split_ind].hplaneIndexByAct[act_ind].insert(
            interPtsAndRays[split_ind].hplaneIndexByAct[act_ind].begin()
                + tmp_h, newHplaneIndex);
        interPtsAndRays[split_ind].indexOfHplaneInOverallOrderForSplit[act_ind].insert(
            it, oldHplaneForSplitOverall);
      }
      /* finished inserting hplane into proper place in this activation */

      // Simply need to add r^j to all the hyperplanes activated after
      const int numHplanesActForSplit =
          interPtsAndRays[split_ind].allHplaneToAct.size();
      for (int h = oldHplaneForSplitOverall + 1; h < numHplanesActForSplit;
          h++) {
        const int curr_hplane_ind = interPtsAndRays[split_ind].allHplaneToAct[h];
        const bool rayIntByCurrHplane = (find_val(ray_ind,
            hplaneStore[curr_hplane_ind].allRaysIntByHplane[split_ind]) >= 0);
        const bool rayCutByCurrHplane = (find_val(ray_ind,
            hplaneStore[curr_hplane_ind].rayToBeCutByHplane[split_ind]) >= 0);

        // If the hplane intersects the ray prior to bd S, but has not cut it yet
        if (rayIntByCurrHplane && !rayCutByCurrHplane) {
          newVertexHelper(vertexStore, init_vert_ind, rayStore, hplaneStore,
              solver, solnInfo, ray_ind, split_ind, curr_hplane_ind);
        }

        // Need to update numRaysCutByHplaneInclPrevAct
        // Ensure the vectors are properly sized
        const int old_num_act_rounds_for_hplane =
            (int) hplaneStore[curr_hplane_ind].numRaysCutByHplaneInclPrevAct[split_ind].size();
        const int newTotalNumRaysCutByHplane =
            hplaneStore[curr_hplane_ind].rayToBeCutByHplane[split_ind].size();
        const int prevTotal =
            (old_num_act_rounds_for_hplane > 0) ?
                hplaneStore[curr_hplane_ind].numRaysCutByHplaneInclPrevAct[split_ind][old_num_act_rounds_for_hplane
                    - 1] :
                0;
        if (old_num_act_rounds_for_hplane < act_ind + 1) {
          hplaneStore[curr_hplane_ind].numRaysCutByHplaneInclPrevAct[split_ind].resize(
              act_ind + 1, prevTotal);
        }
        hplaneStore[curr_hplane_ind].numRaysCutByHplaneInclPrevAct[split_ind][act_ind] =
            newTotalNumRaysCutByHplane;

        // The hplane has been activated before, but perhaps not for this activation round
        const int oldCurrHplaneForActAndSplit = hasHplaneBeenActivated(
            hplaneStore, interPtsAndRays[split_ind].hplaneIndexByAct[act_ind],
            hplaneStore[curr_hplane_ind], true);
        if (oldCurrHplaneForActAndSplit == -1) {
          std::vector<int>::iterator it =
              std::lower_bound(
                  interPtsAndRays[split_ind].indexOfHplaneInOverallOrderForSplit[act_ind].begin(),
                  interPtsAndRays[split_ind].indexOfHplaneInOverallOrderForSplit[act_ind].end(),
                  h);
          const int tmp_h =
              (it
                  - interPtsAndRays[split_ind].indexOfHplaneInOverallOrderForSplit[act_ind].begin());

          interPtsAndRays[split_ind].hplaneIndexByAct[act_ind].insert(
              interPtsAndRays[split_ind].hplaneIndexByAct[act_ind].begin()
                  + tmp_h, curr_hplane_ind);
          interPtsAndRays[split_ind].indexOfHplaneInOverallOrderForSplit[act_ind].insert(
              it, h);
        }
      }
    } /* else oldHplaneForSplit != -1 */

    // Update numRaysCutByHplaneInclPrevAct, including all activations "missed" by this hyperplane
    // Ensure the vectors are properly sized
    const int old_num_act_rounds_for_hplane =
        (int) hplaneStore[newHplaneIndex].numRaysCutByHplaneInclPrevAct[split_ind].size();
    const int newTotalNumRaysCutByHplane =
            hplaneStore[newHplaneIndex].rayToBeCutByHplane[split_ind].size();
    if (old_num_act_rounds_for_hplane < act_ind + 1) {
      hplaneStore[newHplaneIndex].numRaysCutByHplaneInclPrevAct[split_ind].resize(
          act_ind + 1, prevTotalNumRaysCutByHplane);
    }
    hplaneStore[newHplaneIndex].numRaysCutByHplaneInclPrevAct[split_ind][act_ind] =
        newTotalNumRaysCutByHplane;
//  } /* iterating over splitsAffByHplane */
} /* updateHplaneStore (new) */

/**
 * @brief Adds a new vertex that has not been calculated yet
 *
 * @param check_split :: If true, ensures distToSplit > distToHplane
 */
void newVertexHelper(std::vector<Vertex>& vertexStore, const int init_vert_ind,
    const std::vector<Ray>& rayStore, std::vector<Hplane>& hplaneStore,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
    const int ray_ind,
    const int split_ind, const int hplane_ind, const bool check_split) {
  // We have to create the vertex from this new intersection
  Vertex newVertex;
  const double distToHplane = hplaneStore[hplane_ind].distanceToHplane(
      vertexStore[init_vert_ind], rayStore[ray_ind], solver,
      true, solnInfo);

  if (check_split) {
    const double distToSplit = distanceToSplit(vertexStore[init_vert_ind],
        rayStore[ray_ind], solnInfo.feasSplitVar[split_ind],
        solver, solnInfo);

    // If the ray intersects the split before hitting the hyperplane,
    // then it does not lead to a vertex for this split
    if (!lessThanVal(distToHplane, distToSplit)) {
      return;
    }
  }

  calcNewVectorCoordinates(newVertex, distToHplane,
      vertexStore[init_vert_ind], rayStore[ray_ind]);

  // Add this new vertex
  newVertexHelper(vertexStore, init_vert_ind, hplaneStore, newVertex, ray_ind,
      split_ind, hplane_ind);
} /* newVertexHelper */

/**
 * @brief Adds a new vertex that has already been calculated
 */
void newVertexHelper(std::vector<Vertex>& vertexStore, const int init_vert_ind,
    std::vector<Hplane>& hplaneStore, const Vertex& newVertex,
    const int ray_ind, const int split_ind, const int hplane_ind) {
  // Add this new vertex
  const int vertex_ind = addVertex(vertexStore, newVertex);//, hplane_ind);

  // Add ray to set of rays cut by the hyperplane
  hplaneStore[hplane_ind].rayToBeCutByHplane[split_ind].push_back(ray_ind);
  hplaneStore[hplane_ind].vertOfRayToBeCutByHplane[split_ind].push_back(
      init_vert_ind);
  hplaneStore[hplane_ind].vertCreatedByRayToBeCutByHplane[split_ind].push_back(
      vertex_ind);
} /* newVertexHelper */

int addVertex(std::vector<Vertex>& vertexStore, const Vertex& tmpVertex) {
//    const int newHplaneIndex
  const int oldVertexIndex = hasVertexBeenActivated(vertexStore, tmpVertex);
  const int newVertexIndex =
      (oldVertexIndex == -1) ? vertexStore.size() : oldVertexIndex;
  if (oldVertexIndex == -1) {
    vertexStore.push_back(tmpVertex);
  }
//  vertexStore[newVertexIndex].hplaneIndices.push_back(newHplaneIndex);
  return newVertexIndex;
} /* addVertex */

bool pointIsCutByHplane(const PHASolverInterface* const solver,
    const SolutionInfo &solnInfo, const CoinPackedVector& row,
    const double RHS, const Hplane& hplane) {
  const double* UB = solver->getColUpper();
  const double* LB = solver->getColLower();

  const bool currRowIsRay = isZero(RHS);

  const int numElmts = row.getNumElements();
  const int* elmtIndices = row.getIndices();
  const double* elmtVals = row.getElements();
  double activity = 0.0;
  for (int e = 0; e < numElmts; e++) {
    // Either hplaneRowIndx >= 0, in which case we activated a basic var
    // or it is -1, and we activated an upper or lower bound of a nb var
    if (hplane.row >= 0) {
      activity += elmtVals[e] * solnInfo.raysOfC1[elmtIndices[e]][hplane.row];
    } else {
      const int nb_ind = (-1 - solnInfo.rowOfVar[hplane.var]);
      if (elmtIndices[e] == nb_ind)
        activity = elmtVals[e];
    }
  }

  // Note that if hplaneRowIndx == -1, then we activated a bound of a nb var,
  // and it can't be a bound of a slack variable, since they are at zero in the cobasis,
  // and that is their only bound.
  if (hplane.var >= solnInfo.numCols) {
    // Recall that all slack variables are non-negative in this
    // setting, since they were complemented in SolutionInfo
    if (!currRowIsRay && lessThanVal(solnInfo.a0[hplane.row] + activity, 0.0)) {
      return true;
    } else if (currRowIsRay && lessThanVal(activity, 0.0)) {
      return true;
    }
  } else {
    // If this was a basic variable, no bounds were changed
    // Otherwise things may have shifted
    double lb, ub, newValue;
    if (hplane.row >= 0) {
      lb = LB[hplane.var];
      ub = UB[hplane.var];
      newValue = solnInfo.a0[hplane.row] + activity;
    } else {
      lb = 0.0;
      ub = UB[hplane.var] - LB[hplane.var];
      newValue = activity;
    }
    if (hplane.ubflag == static_cast<int>(HplaneBoundFlag::toUB)) {
      if (!currRowIsRay && greaterThanVal(newValue, ub)) {
        return true;
      } else if (currRowIsRay && greaterThanVal(activity, 0.0)) {
        return true;
      }
    } else {
      if (!currRowIsRay && lessThanVal(newValue, lb)) {
        return true;
      } else if (currRowIsRay && lessThanVal(activity, 0.0)) {
        return true;
      }
    }
  }
  return false;
} /* pointIsCutByHplane */

#ifdef SHOULD_WRITE_HPLANES
void printActivatedHplaneDetails(FILE* fptr,
    const SolutionInfo& solnInfo,
    const std::vector<Hplane>& hplaneStore,
    const std::vector<IntersectionInfo>& interPtsAndRays,
    const int act_ind) {
  fprintf(fptr,
      "This file contains details on the partially activated hyperplanes. UB flag: -1 = hplane not a bound. 0 = lower-bound. 1 = upper-bound.\n\n");
  fprintf(fptr,
      "(UBlag = -1 implies the hyperplane is structural and not a variable bound)\n\n");
  fprintf(fptr,
      "SplitVarIndex,No,hplaneVarIndex,hplaneRowIndex,UBFlag,RaysToBeCut...\n");
  for (int s = 0; s < solnInfo.numFeasSplits; s++) {
    for (unsigned h = 0;
        h < interPtsAndRays[s].hplaneIndexByAct[act_ind].size();
        h++) {
      const int hplane_ind =
          interPtsAndRays[s].hplaneIndexByAct[act_ind][h];
      fprintf(fptr, "%d,%d,%d,%d,%d,", solnInfo.feasSplitVar[s], h,
          hplaneStore[hplane_ind].var, hplaneStore[hplane_ind].row,
          hplaneStore[hplane_ind].ubflag);

      const int start_ind = getStartIndex(hplaneStore[hplane_ind], s,
          act_ind);
      const int end_ind = getEndIndex(hplaneStore[hplane_ind], s, act_ind);
      for (int r = start_ind; r < end_ind; r++) {
        fprintf(fptr, "%d,",
            hplaneStore[hplane_ind].rayToBeCutByHplane[s][r]);
      }
      fprintf(fptr, "\n");
    }
  }

  fprintf(fptr, "\nAll rays to be cut\n");
  for (int s = 0; s < solnInfo.numFeasSplits; s++) {
    for (unsigned h = 0;
        h < interPtsAndRays[s].hplaneIndexByAct[act_ind].size();
        h++) {
      const int hplane_ind =
          interPtsAndRays[s].hplaneIndexByAct[act_ind][h];
      fprintf(fptr, "%d,%d,%d,%d,%d,", solnInfo.feasSplitVar[s], h,
          hplaneStore[hplane_ind].var, hplaneStore[hplane_ind].row,
          hplaneStore[hplane_ind].ubflag);

      const int numAllRaysIntByHplane =
          hplaneStore[hplane_ind].allRaysIntByHplane[s].size();
      for (int r = 0; r < numAllRaysIntByHplane; r++) {
        fprintf(fptr, "%d,",
            hplaneStore[hplane_ind].allRaysIntByHplane[s][r]);
      }
      fprintf(fptr, "\n");
    }
  }
} /* printActivatedHplaneDetails */
#endif
