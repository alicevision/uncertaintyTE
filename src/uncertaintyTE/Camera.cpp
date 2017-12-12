// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   ColmapCamera.cpp
 * Author: Michal Polic
 * 
 * Created on October 24, 2017, 4:27 PM
 */

#include "Camera.h"

Camera::Camera() {}

Camera::~Camera() {}

ostream& operator<< (ostream& out, const Camera& c){
    out << "> Camera " << c._id << " " << c._model
            << " [focal:" << c._f << ", width:" << c._img_width << ", height:" << c._img_height
            << ", uv:" << c._uv[0] << "," << c._uv[2] 
            << ", rad_dist:" << c._r[0] << "," << c._r[1] << "," << c._r[2] << "]\n";
    return out;
};
