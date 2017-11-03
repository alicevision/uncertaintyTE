#pragma once
#include "compute.h"
#include "bal_problem.h"

#ifdef USE_MATLAB
    #include <mex.h>

    ceres::examples::BALProblem loadSceneMatlab(const mxArray *prhs[]);
    void buildJacobian(ceres::examples::BALProblem *bal_problem, ceres::CRSMatrix &jacobian, cov::Options &options);
#endif