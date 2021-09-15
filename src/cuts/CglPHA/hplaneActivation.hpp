//============================================================================
// Name        : hplaneActivation.hpp
// Author      : Aleksandr M. Kazachkov
// Version     : 0.2014.05.23
// Copyright   : Your copyright notice
// Description : Methods to activate hyperplanes
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#pragma once

#include "Hplane.hpp"
#include "Vertex.hpp"
#include "Ray.hpp"
#include "SolutionInfo.hpp"
#include "IntersectionInfo.hpp"
#include "AdvCut.hpp"
#include "typedefs.hpp"

void addRaysIntByHplaneBeforeBdS(std::vector<int>& vertOfRayIndicesIntByHplane,
    std::vector<int>& rayIndicesIntByHplane,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
    const Hplane& hplane,
    const std::vector<Vertex>& vertexStore, const int init_vert_ind,
    const std::vector<Ray>& rayStore, const int split_ind);

int removePtsAndRaysCutByHplane(const PHASolverInterface* const solver,
    const SolutionInfo &solnInfo, const std::vector<bool> &raysCutFlag,
    const int split_var, const Hplane& hplane, const int hplane_ind_in_split,
//    const int hplane_ind, const std::vector<Hplane>& hplaneStore,
    const std::vector<Vertex>& vertexStore, const std::vector<Ray>& rayStore,
    IntersectionInfo& interPtsAndRays, const AdvCut* NBSIC, const double norm,
    double& sumRemovedDepth, const bool fake_act = false);

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
    const int split_ind,
    const int ray_ind,
    const int act_ind, const bool checkIfRayPrevCut = false);

/**
 * @brief Adds a new vertex that has not been calculated yet
 *
 * @param check_split :: If true, ensures distToSplit > distToHplane
 */
void newVertexHelper(std::vector<Vertex>& vertexStore, const int init_vert_ind,
    const std::vector<Ray>& rayStore, std::vector<Hplane>& hplaneStore,
    const PHASolverInterface* const solver, const SolutionInfo& solnInfo,
    const int ray_ind,
    const int split_ind, const int hplane_ind, const bool check_split =
        false);

/**
 * @brief Adds a new vertex
 */
void newVertexHelper(std::vector<Vertex>& vertexStore, const int init_vert_ind,
    std::vector<Hplane>& hplaneStore, const Vertex& newVertex,
    const int ray_ind, const int split_ind, const int hplane_ind);

int addVertex(std::vector<Vertex>& vertexStore, const Vertex& tmpVertex);
//    const int newHplaneIndex);

bool pointIsCutByHplane(const PHASolverInterface* const solver,
    const SolutionInfo &solnInfo, const CoinPackedVector& row,
    const double RHS, const Hplane& hplane);

#ifdef SHOULD_WRITE_HPLANES
void printActivatedHplaneDetails(FILE* fptr,
    const SolutionInfo& solnInfo,
    const std::vector<Hplane>& hplaneStore,
    const std::vector<IntersectionInfo>& interPtsAndRays,
    const int act_ind = 0);
#endif
