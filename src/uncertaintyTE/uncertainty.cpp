#include "compute.h"
#include "auxCmd.h"
#include "FactoryIO.h"
#include "IO.h"
#include "ColmapIO.h"
#include "JacobianIO.h"
#include "JacobianComposer.h"

/*
Library function for C++ call: getCovariances( ... )
In:
  - options: informations about the reconstruction (numCams, camParams, numPoints, numObs)
  - statistic: blank object for the output statistics
  - jacobian: sparse matrix (form Ceres-solver) with jacobian ( you can use the output of the computation of jacobian in Ceres as the input )
  - h_camUnc: array which contatins covariances for cameras
  - h_ptUnc: array which contatins covariances for points
*/
#ifdef _WIN32
    extern "C" __declspec(dllexport)
#endif
inline void getCovariances(
    cov::Options &options,
    cov::Statistic &statistic,
    ceres::CRSMatrix &jacobian,
    double* points3D,
    double* h_camUnc,
    double* h_ptUnc)
{
    JacobianComposer::findPts2Fix(options, options._numPoints, points3D);
    computeCovariances(options, statistic, jacobian, h_camUnc, h_ptUnc);
}

