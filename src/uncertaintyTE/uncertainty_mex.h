// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/*
* File:   uncertainty_mex.h
* Author: Michal Polic
*
* Created on November 5, 2017, 11:30 AM
*/

#ifdef USE_MATLAB
    #include <mex.h>
#endif

#ifdef _WIN32
    #define MEX_FUNCTION_NAME mexFunction
#elif __linux__ 
    #define MEX_FUNCTION_NAME __mexFunction
#endif   

#ifdef USE_MATLAB    
    void MEX_FUNCTION_NAME(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
#endif