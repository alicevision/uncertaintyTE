#pragma once
#include <vector>
#include <iostream>
#include <memory>

#define _USE_MATH_DEFINES
#include <cmath>

#include "ceres/ceres.h"

#define SVD_QR_ITERATION 0
#define SVD_DEVIDE_AND_CONQUER 1
#define TAYLOR_EXPANSION 2

namespace cov{
	struct Options {
	public:
		double _epsilon, _lambda;
		int _algorithm, _numCams, _camParams, _numPoints, _numObs, _svdRemoveN, _maxIterTE;
		int *_pts2fix = NULL;
                bool _debug = false;

		Options() : _lambda(-1), _svdRemoveN(-1), _maxIterTE(-1) {}

		Options(int numCams, int camParams, int numPoints, int numObs) :
			_algorithm(TAYLOR_EXPANSION), _epsilon(1e-10), _lambda(-1), _numCams( numCams ), _camParams(camParams), _numPoints(numPoints), _numObs(numObs), _svdRemoveN(-1), _maxIterTE(-1) {}

		Options(int algorithm, double eps_or_lamb, int numCams, int camParams, int numPoints, int numObs) :
			_algorithm(algorithm), _epsilon(eps_or_lamb), _lambda(eps_or_lamb), _numCams( numCams ), _camParams(camParams), _numPoints(numPoints), _numObs(numObs), _svdRemoveN(-1), _maxIterTE(-1) {}
	
		~Options() {
			free(_pts2fix);
		}
	};

	struct Statistic {
		double timeCreateJ, timeFixJ, timeNormJ, timeMultiplyJJ, timeSplitJJ, timeInvV, timeComposeZ, timeInvZ, timeTE, timePtsUnc, timeAll;
		double lambda;
		int *fixedPts;
		std::vector<double> cycle_change;

		~Statistic(){
			free(fixedPts);
		}
	};
}

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
