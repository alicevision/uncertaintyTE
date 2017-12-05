/*
 * Copyright (C) 2017 policmic
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * File:   IO.h
 * Author: policmic
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

