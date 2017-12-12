// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   ColmapCamera.h
 * Author: Michal Polic
 *
 * Created on October 24, 2017, 4:27 PM
 */

#ifndef CAMERA_H
#define CAMERA_H

#include <string>
#include <iostream>
using namespace std;

class Camera {
public:
    int _id;
    string _model;
    int _img_width;      // width in px
    int _img_height;     // height in px
    double _uv[2];       // principal point
    double _f;           // focal length
    double _r[2] = {0,0};    // radial distortion 1-2 params
    
    Camera();
    ~Camera();
    
    //virtual int nParameters();
    //virtual double* getParameters();
private:

};

ostream& operator<< (ostream& out, const Camera& c);

#endif /* COLMAPCAMERA_H */

