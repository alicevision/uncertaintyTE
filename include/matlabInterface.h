#pragma once

#define MATLAB_INTERFACE_HEADER_INCLUDED

#ifndef CERES_EXAMPLES_BAL_PROBLEM_H_
	#include "bal_problem.h"
#endif // !CERES_EXAMPLES_BAL_PROBLEM_H_


ceres::examples::BALProblem loadSceneMatlab(const mxArray *prhs[]);

void buildJacobian(ceres::examples::BALProblem *bal_problem, ceres::CRSMatrix &jacobian, cov::Options &options);