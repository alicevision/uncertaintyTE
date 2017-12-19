// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   AliceVisionIO.cpp
 * Author: Michal Polic
 * 
 * Created on November 5, 2017, 11:30 AM
 */

#include "AliceVisionIO.h"
#include "auxCmd.h"

#ifdef USE_ALICEVISION
#include <aliceVision/sfm/sfm.hpp>
#include <aliceVision/sfm/BundleAdjustmentCeres.hpp>
#endif


AliceVisionIO::AliceVisionIO() {}

AliceVisionIO::~AliceVisionIO() {}

bool AliceVisionIO::read(const std::string& input_file, Scene& scene) {
#ifdef USE_ALICEVISION
	std::cout << "Loading an AliceVision scene: " << input_file << '\n';
    aliceVision::sfm::SfM_Data sfm_data;
    // Load input SfM_Data scene
    if (!aliceVision::sfm::Load(sfm_data, sSfM_Data_Filename_In, aliceVision::sfm::ESfM_Data(aliceVision::sfm::ALL)))
    {
        throw std::runtime_error("The input SfM_Data file \"" + input_file + "\" cannot be read.");
    }

    scene._jacobian = ceres::CRSMatrix();
    scene._options = cov::Options();
    {
      BundleAdjustmentCeres bundleAdjustmentObj;
      BA_Refine refineOptions = BA_REFINE_ROTATION | BA_REFINE_TRANSLATION | BA_REFINE_STRUCTURE;
      bundleAdjustmentObj.createJacobian(sfmData, refineOptions, scene._jacobian);
    }
	return true;
#else
    throw std::runtime_error("Need to be built with AliceVision library to load the scene.");
#endif
}

bool AliceVisionIO::write(const std::string& output_dir, Scene& scene) {
    throw std::runtime_error("AliceVision exporter not implemented.");
}

int AliceVisionIO::data_type() {
	return -1;
}
