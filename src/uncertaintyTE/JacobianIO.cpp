// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   JacobianIO.cpp
 * Author: Michal Polic
 * 
 * Created on November 2, 2017, 5:22 PM
 */

#include "JacobianIO.h"

JacobianIO::JacobianIO() {}

JacobianIO::~JacobianIO() {}

int JacobianIO::data_type(){
    return _data_type;
}

bool JacobianIO::read(const std::string& input_dir, Scene& scene){
	std::ifstream in_file(input_dir, std::ios_base::in);   // TODO: move it to the JacobianIO
    return loadJacobian(in_file, cov::eAlgorithmSvdTaylorExpansion, scene._jacobian, scene._options);
}

bool JacobianIO::write(const std::string& output_dir, Scene& scene){
    exit(1);
}

bool JacobianIO::loadJacobian(std::ifstream& file, cov::EAlgorithm algorithm, ceres::CRSMatrix& jacobian, cov::Options& options) {
	if (!file.good()) {
		std::cerr << "\nThe input file doesn't exist.\n";
		return false;
	}

	int numJ;
	int fixPts[3];
	options._algorithm = algorithm;
	options._epsilon = 1e-10;
	options._svdRemoveN = 7;
	options._maxIterTE = -1;

	file >> options._lambda >> options._numCams >> options._camParams >> options._numPoints >> options._numObs;
	file >> fixPts[0] >> fixPts[1] >> fixPts[2];
	if (fixPts[0] >= 0 && fixPts[1] >= 0 && fixPts[2] >= 0)
		options._pts2fix = new int[3]{ fixPts[0], fixPts[1], fixPts[2] };

	file >> jacobian.num_rows >> jacobian.num_cols >> numJ;
	jacobian.rows = std::vector<int>(jacobian.num_rows + 1);
	jacobian.cols = std::vector<int>(numJ);
	jacobian.values = std::vector<double>(numJ);

	for (int i = 0; i <= jacobian.num_rows; ++i)
		file >> jacobian.rows[i];
	for (int i = 0; i < numJ; ++i)
		file >> jacobian.cols[i];
	for (int i = 0; i < numJ; ++i)
		file >> jacobian.values[i];

	return true;
}
