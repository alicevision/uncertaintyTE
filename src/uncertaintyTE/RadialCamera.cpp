// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   RadialCamera.cpp
 * Author: Michal Polic
 * 
 * Created on November 3, 2017, 10:43 AM
 */

#include "RadialCamera.h"

RadialCamera::RadialCamera() {}

RadialCamera::~RadialCamera() {}

int RadialCamera::nParameters(){
    return _num_parameters;
}

//double* RadialCamera::getParameters(){
//
//    // TODO: ... change the camera and image class
//}