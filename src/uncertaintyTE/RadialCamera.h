// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   RadialCamera.h
 * Author: Michal Polic
 *
 * Created on November 3, 2017, 10:43 AM
 */

#ifndef RADIALCAMERA_H
#define RADIALCAMERA_H

#include <vector>
#include "Camera.h"

class RadialCamera : Camera {
public:
    RadialCamera();
    ~RadialCamera();
    
    // functions 
    int nParameters();
    //double* getParameters();
    
private:
    int _num_parameters = 9;
    std::vector<double> _parameters;
};

#endif /* RADIALCAMERA_H */

