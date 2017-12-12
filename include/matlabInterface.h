#pragma once

#include "bal_problem.h"

ceres::examples::BALProblem loadSceneMatlab(const mxArray *prhs[]);

void buildJacobian(ceres::examples::BALProblem *bal_problem, ceres::CRSMatrix &jacobian, cov::Options &options);