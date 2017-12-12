// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   Scene.h
 * Author: Michal Polic
 *
 * Created on October 25, 2017, 10:08 AM
 */

#ifndef SCENE_H
#define SCENE_H

#include <map>
#include "Camera.h"
#include "Image.h"
#include "Point3D.h"
#include "Point2D.h"
#include "uncertainty.h"
#include "ceres/rotation.h"

class Scene {
public:
	// input parameters
    map<int, Camera> _cameras;
    map<int, Image> _images;
    map<int, Point3D> _points3D;

	// internal / input parameters and settings
	ceres::CRSMatrix _jacobian;
	cov::Options _options;

    cov::Uncertainty _uncertainty;

    Scene();
    ~Scene();

};

ostream& operator<< (ostream& out, const Scene& s);

#endif /* SCENE_H */

