#pragma once
#include <vector>
#include <iostream>
#include <memory>

#define _USE_MATH_DEFINES
#include <cmath>

#include "ceres/ceres.h"


namespace cov {

enum EAlgorithm {
    eAlgorithmSvdQrIteration = 0,
    eAlgorithmSvdDeviceAndconquer = 1,
    eAlgorithmSvdTaylorExpansion = 2
};


inline std::string EAlgorithm_enumToString(EAlgorithm alg) {
    switch (alg) {
        case eAlgorithmSvdQrIteration: return "SVD_QR_ITERATION";
        case eAlgorithmSvdDeviceAndconquer: return "SVD_DEVIDE_AND_CONQUER";
        case eAlgorithmSvdTaylorExpansion: return "TAYLOR_EXPANSION";
        default: return "not defined";
    }
}

inline EAlgorithm EAlgorithm_stringToEnum(const std::string& algorithm) {
    if (algorithm == "SVD_QR_ITERATION")
        return eAlgorithmSvdQrIteration;
    if (algorithm == "SVD_DEVIDE_AND_CONQUER")
        return eAlgorithmSvdDeviceAndconquer;
    if (algorithm == "TAYLOR_EXPANSION")
        return eAlgorithmSvdTaylorExpansion;
    throw std::runtime_error(std::string("Unrecognized algorithm: ") + algorithm);
}

struct Options {
public:
    double _epsilon, _lambda;
    EAlgorithm _algorithm;
    int _numCams, _camParams, _numPoints, _numObs, _svdRemoveN, _maxIterTE;
    int *_pts2fix = NULL;
    bool _debug = false;

    Options() : _lambda(-1), _svdRemoveN(-1), _maxIterTE(-1) {}

    Options(int numCams, int camParams, int numPoints, int numObs) :
        _algorithm(eAlgorithmSvdTaylorExpansion), _epsilon(1e-10), _lambda(-1), _numCams( numCams ), _camParams(camParams), _numPoints(numPoints), _numObs(numObs), _svdRemoveN(-1), _maxIterTE(-1) {}

    Options(EAlgorithm algorithm, double eps_or_lamb, int numCams, int camParams, int numPoints, int numObs) :
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

struct Uncertainty {
    std::size_t _nbCovarianceValuePerCam = 0;
    std::vector<double> _camerasUnc;
    std::vector<double> _pointsUnc;

    void init(const Options& options)
    {
        _nbCovarianceValuePerCam = (0.5 * options._camParams * (options._camParams + 1));
        _camerasUnc.resize(options._numCams * _nbCovarianceValuePerCam);
        _pointsUnc.resize(options._numPoints * 6);
    }

    /**
     * @return 21 values per camera
     */
    const std::vector<double>& getCamerasUncRaw() const { return _camerasUnc; }
    /**
     * @return matrix of 6x6 values per point
     */
    const std::vector<double>& getCamerasUncMatrix() const
    {
        // TODO
        return _camerasUnc;
    }
    /**
     * @return 6 values per camera
     */
    const std::vector<double>& getCamerasUncEigenValues() const
    {
        // TODO
        return _camerasUnc;
    }

    /**
     * @return 6 values per point
     */
    const std::vector<double>& getPointsUncRaw() const { return _pointsUnc; }
    /**
     * @return matrix of 3x3 values per point
     */
    const std::vector<double>& getPointsUncMatrix() const
    {
        // TODO
        return _pointsUnc;
    }
    /**
     * @return 3 values per point
     */
    const std::vector<double>& getPointsUncEigenValues() const
    {
        // TODO
        return _pointsUnc;
    }
};
}

/**
 * @param[in] options: informations about the reconstruction (numCams, camParams, numPoints, numObs)
 * @param[out] statistic: blank object for the output statistics
 * @param[in] jacobian: sparse matrix (form Ceres-solver) with jacobian ( you can use the output of the computation of jacobian in Ceres as the input )
 * @param[in] points3D: all 3D points (to select the 3 static points)
 * @param[out] uncertainties: output covariances for cameras and points
 */
#ifdef _WIN32
    extern "C" __declspec(dllexport)
#endif
void getCovariances(
    cov::Options &options,
    cov::Statistic &statistic,
    ceres::CRSMatrix &jacobian,
    double* points3D,
    cov::Uncertainty& uncertainties);
