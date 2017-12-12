// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   JacobianComposer.hpp
 * Author: Michal Polic
 *
 * Created on November 2, 2017, 4:41 PM
 */

#ifndef JACOBIANCOMPOSER_HPP
#define JACOBIANCOMPOSER_HPP

#include <random>

#include "uncertainty.h"
#include "Scene.h"
#include "snavely_reprojection_error.h"


class JacobianComposer {
public:
    // functions 
    static void scene2Jacobian(std::string cam_model, std::string algorithm, Scene &scene);
    static void  findPts2Fix(cov::Options &opt, int N, double* points3D);
	static void findPts2Fix(cov::Options &opt, int n, std::map<int, Point3D> &pts3D);
	static void findPts2Fix(cov::Options &opt, int N, std::vector<Point3D> &pts3D);

private:
	// constructors 
	JacobianComposer();
	~JacobianComposer();

    // functions 
	static void Scene2Problem(Scene &s, ceres::Problem* problem, std::string cam_model);

};

#endif /* JACOBIANCOMPOSER_HPP */

