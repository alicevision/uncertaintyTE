// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/*
* File:   main_uncertainty.cpp
* Author: Michal Polic
*/

#include <uncertaintyTE/compute.h>
#include <uncertaintyTE/auxCmd.h>
#include <uncertaintyTE/FactoryIO.h>
#include <uncertaintyTE/IO.h>
#include <uncertaintyTE/ColmapIO.h>
#include <uncertaintyTE/JacobianIO.h>
#include <uncertaintyTE/JacobianComposer.h>

#include <gflags/gflags.h>

#include <mex.h>
#include <uncertaintyTE/matlabInterface.h>

#ifdef _WIN32
	#define EXECUTABLE_FILE "uncertainty.exe"
#elif __linux__ 
	#define EXECUTABLE_FILE "uncertainty"
#endif

#ifdef _WIN32
    #define MEX_FUNCTION_NAME mexFunction
#elif __linux__
    #define MEX_FUNCTION_NAME __mexFunction
#endif


/**
 * @brief Matlab interface, call in Matlab: unc( .... )
 */
void MEX_FUNCTION_NAME(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	cov::Statistic statistic = cov::Statistic();
	cov::Options options = cov::Options();
	ceres::CRSMatrix jacobian = ceres::CRSMatrix();
	
	// Load and create Jacobian using camera model [angle axis, camera center, focal length, radial distortion]
	ceres::examples::BALProblem bal_problem = loadSceneMatlab( prhs );
	buildJacobian(&bal_problem, jacobian, options);

	// alocate output arrays
	int NcamArr = 0.5 * bal_problem.camera_block_size() * (bal_problem.camera_block_size() + 1) * bal_problem.num_cameras();
	int NptsArr = 6 * bal_problem.num_points();
	double *ptsUnc = (double*)malloc(NptsArr * sizeof(double));
	double *camUnc = (double*)malloc(NcamArr * sizeof(double));
	plhs[0] = mxCreateDoubleMatrix(NptsArr + NcamArr, 1, mxREAL);
	double *outC = mxGetPr(plhs[0]);

	// COMPUTE COVARIANCES
	computeCovariances(options, statistic, jacobian, camUnc, ptsUnc);

	// write results to the output
	for (int i = 0; i < NcamArr; i++)
		outC[i] = camUnc[i];
	for (int i = 0; i < NptsArr; i++)
		outC[NcamArr + i] = ptsUnc[i];

	free(ptsUnc);
	free(camUnc);
	mexPrintf("Computation of covariances ... [done]\n\n");
}
