#pragma once	

#define _USE_MATH_DEFINES
#include <cmath>
#include "compute.h"


void replaceAll(std::string& str, const std::string& from, const std::string& to);

void loadJacobian(std::ifstream& file, int algorithm, ceres::CRSMatrix& jacobian, cov::Options& options);

void printJacobian(ceres::CRSMatrix &J);
void printJacobianMEX(ceres::CRSMatrix &J);
void setPts2Fix(cov::Options &opt, int N, double *pts);

#ifdef USE_OPENMVG
	int loadSceneOpenMVG(std::string sSfM_Data_Filename_In, openMVG::sfm::SfM_Data &sfm_data);
	void openmvgSfM2Jacobian(openMVG::sfm::SfM_Data &sfm_data, ceres::CRSMatrix &jacobian, cov::Options &opt);
#endif
       
