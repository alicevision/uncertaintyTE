// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/*
* File:   uncertainty.cpp
* Author: Michal Polic
*/

#include "compute.h"
#include "auxCmd.h"
#include "FactoryIO.h"
#include "IO.h"
#include "ColmapIO.h"
#include "JacobianIO.h"
#include "JacobianComposer.h"


void getCovariances(
    cov::Options &options,
    cov::Statistic &statistic,
    ceres::CRSMatrix &jacobian,
    double* points3D,
    cov::Uncertainty& uncertainties)
{
    JacobianComposer::findPts2Fix(options, options._numPoints, points3D);
    uncertainties.init(options);
    computeCovariances(options, statistic, jacobian, uncertainties._camerasUnc.data(), uncertainties._pointsUnc.data());
}



namespace cov {

	void Uncertainty::init(const Options& options)
	{
		_camParams = options._camParams;
		_numCams = options._numCams;
		_numPoints = options._numPoints;
		_nbCovarianceValuePerCam = (0.5 * options._camParams * (options._camParams + 1));
		_camerasUnc.resize(options._numCams * _nbCovarianceValuePerCam);
		_pointsUnc.resize(options._numPoints * 6);
	}

	const std::vector<double> Uncertainty::getCamerasUncRaw() const { return _camerasUnc; }

	const std::vector<double> Uncertainty::getCameraUncMatrix(int id) const {
		std::vector<double> camCov;
		camCov.resize(_camParams * _camParams);
		double *camCovArray = camCov.data();
		int t = 0;
		for (int i = 0; i < _camParams; i++) {
			for (int j = i; j < _camParams; j++) {
				camCovArray[i * _camParams + j] = _camerasUnc.data()[t + id * _nbCovarianceValuePerCam];
				camCovArray[i + j * _camParams] = _camerasUnc.data()[t + id * _nbCovarianceValuePerCam];
				t++;
			}
		}
		return camCov;
	}

	const std::vector<double> Uncertainty::getCamerasUncMatrices() const {
		int numParams = _camParams * _camParams;
		std::vector<double> camsCov;
		camsCov.resize(_numCams * numParams);
		for (int k = 0; k < _numCams; k++){
			std::vector<double> camCov = getCameraUncMatrix(k);
			std::copy(camCov.begin(), camCov.end(), &camsCov.data()[k*numParams]);
		}
		return camsCov;
	}

	const std::vector<double> Uncertainty::getCamerasUncEigenValues() const
	{
		std::vector<double> eigenValues;
		for (int i = 0; i < _numCams; i++) {
			std::vector<double> covMat = getCameraUncMatrix(i);
			double *covArray = (double*)malloc(_camParams*_camParams * sizeof(double));
			std::copy(covMat.begin(), covMat.end(), covArray);

			Eigen::Map<Eigen::MatrixXd> A( covArray, _camParams, _camParams);
			Eigen::VectorXd eigValCamera = A.eigenvalues().real();
			std::sort(eigValCamera.data(), eigValCamera.data() + _camParams, std::greater<double>());
			for (int j = 0; j < eigValCamera.size(); j++)
				eigenValues.push_back(eigValCamera(j));
		}
	
		return eigenValues;
	}




	const std::vector<double> Uncertainty::getPointsUncRaw() const { return _pointsUnc; }

	const std::vector<double> Uncertainty::getPointUncMatrix(int id) const {
		std::vector<double> out;
		out.resize(3 * 3);
		double *tmp = out.data();
		tmp[0] = _pointsUnc[0 + 6 * id];
		tmp[1] = _pointsUnc[1 + 6 * id];
		tmp[2] = _pointsUnc[2 + 6 * id];
		tmp[3] = _pointsUnc[1 + 6 * id];
		tmp[4] = _pointsUnc[3 + 6 * id];
		tmp[5] = _pointsUnc[4 + 6 * id];
		tmp[6] = _pointsUnc[2 + 6 * id];
		tmp[7] = _pointsUnc[4 + 6 * id];
		tmp[8] = _pointsUnc[5 + 6 * id];
		return out;
	}

	const std::vector<double> Uncertainty::getPointsUncMatrices() const {
		std::vector<double> out;
		out.reserve(_numPoints * 3 * 3);
		for (int k = 0; k < _numPoints; k++) {
			std::vector<double> ptUnc = getPointUncMatrix(k);
			out.insert(out.end(), ptUnc.begin(), ptUnc.end());
		}
		return out;
	}

	const std::vector<double> Uncertainty::getPointsUncEigenValues() const	{
		std::vector<double> eigenValues;
		for (int i = 0; i < _numPoints; i++) {
			std::vector<double> ptUnc = getPointUncMatrix(i);
			double *ptCovArray = (double*)malloc(9 * sizeof(double));
			std::copy(ptUnc.begin(), ptUnc.end(), ptCovArray);
			
			Eigen::Map<Eigen::MatrixXd> A(ptCovArray, 3, 3);
			Eigen::VectorXd eigValPt = A.eigenvalues().real();
			std::sort(eigValPt.data(), eigValPt.data() + 3, std::greater<double>());

			for (int j = 0; j < eigValPt.size(); j++)
				eigenValues.push_back(eigValPt(j));
		}
		return eigenValues;
	}

}