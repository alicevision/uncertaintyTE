// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   IO.h
 * Author: Michal Polic
 *
 * Created on October 25, 2017, 10:08 AM
 */

#ifndef IO_H
#define IO_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "Scene.h"

#define SCENE_DATA 0
#define JACOBIAN_DATA 1


bool saveResults(const std::string& filepath, cov::Options& options, cov::Statistic& statistic,
    int num_camera_covar_values, const double* camUnc, const double* ptsUnc);

inline bool saveResults(const std::string& filepath, cov::Options& options, cov::Statistic& statistic, const cov::Uncertainty& uncertainty)
{
    return saveResults(filepath, options, statistic, uncertainty._nbCovarianceValuePerCam, &uncertainty._camerasUnc[0], &uncertainty._pointsUnc[0]);
}

class IO {
public:
    virtual int data_type() = 0;
    virtual bool read(const std::string& input_dir, Scene& scene) = 0;
    virtual bool write(const std::string& output_dir, Scene& scene) = 0;

    bool writeCov2File(const std::string& filepath, Scene& scene, cov::Statistic& statistic);
};

#endif /* IO_H */

