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


using namespace std;

class IO {
public:
    virtual int data_type() = 0;
    virtual bool read(const string input_dir, Scene& scene) = 0;
	virtual bool write(const string output_dir, Scene& scene) = 0;

	bool writeCov2File(const string output_dir, Scene& scene, cov::Statistic& statistic);
private:

	bool saveResults(const std::string& out_dir, cov::Options& options, cov::Statistic& statistic,
		int num_camera_covar_values, double* camUnc, double *ptsUnc);

};

#endif /* IO_H */

