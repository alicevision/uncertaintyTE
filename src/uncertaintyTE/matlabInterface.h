// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/*
* File:   matlabInterface.h
* Author: Michal Polic
*/

#pragma once
#include "compute.h"
#include "bal_problem.h"

#ifdef USE_MATLAB
    #include <mex.h>

    ceres::examples::BALProblem loadSceneMatlab(const mxArray *prhs[]);
    void buildJacobian(ceres::examples::BALProblem *bal_problem, ceres::CRSMatrix &jacobian, cov::Options &options);
#endif