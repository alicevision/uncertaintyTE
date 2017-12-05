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