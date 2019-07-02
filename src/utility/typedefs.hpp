//============================================================================
// Name        : typedefs.hpp
// Author      : akazachk
// Version     : 0.2015.01.15
// Copyright   : Your copyright notice
// Description : Header for typedefs.hpp, which has useful typedefs
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//============================================================================

#pragma once
#include "OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"

// Redefine as whichever type of solver you want
using PHASolverInterface = OsiClpSolverInterface;
