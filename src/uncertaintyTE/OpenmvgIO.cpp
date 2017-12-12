// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   OpenmvgIO.cpp
 * Author: Michal Polic
 * 
 * Created on November 5, 2017, 11:30 AM
 */

#include "OpenmvgIO.h"
#include "auxCmd.h"

OpenmvgIO::OpenmvgIO() {}

OpenmvgIO::~OpenmvgIO() {}

bool OpenmvgIO::read(const std::string& input_file, Scene& scene) {
#ifdef USE_OPENMVG
	std::cout << "Loading a OpenMVG scene: " << input_file << '\n';
	openMVG::sfm::SfM_Data sfm_data;
	loadSceneOpenMVG(input_file, sfm_data);
	scene._jacobian = ceres::CRSMatrix();
	scene._options = cov::Options();
	openmvgSfM2Jacobian(sfm_data, scene._jacobian, scene._options);    // work for just 9 params for camera representation
	return true;
#else
    std::cerr << "Load OpenMVG library to enable its input.";
	exit(1);
#endif
}

bool OpenmvgIO::write(const std::string& output_dir, Scene& scene) {
	std::cerr << "OpenMVG output is not implemented yet.";
	exit(1);
}

int OpenmvgIO::data_type() {
	return -1;
};
