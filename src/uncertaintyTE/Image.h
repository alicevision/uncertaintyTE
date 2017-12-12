// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   ColmapImage.h
 * Author: Michal Polic
 *
 * Created on October 24, 2017, 4:22 PM
 */

#ifndef IMAGE_H
#define IMAGE_H

#include <vector>
#include <iostream>
#include "Point2D.h"

using namespace std;

class Image {
public:
    int _id;
    int _cam_id;
    double _q[4];    // quaternion (Eigen format)
    double _R[9];    // rotation matrix (Eigen format)
    double _aa[3];   // Euler vector - angle * axis 
    
    double _t[3];    // image translation 
    double _C[3];    // image center
    vector<Point2D> _point2D;
    
    // camera parameters --> replace with proper models
    double _f;      // focal length
    double _r[2];   // radial distortion
    
    Image();
    virtual ~Image();
    
    void qt2aaRC();  // recompute the Colmap input (quatrnion and translation) to the model parameters 
private:

};

ostream& operator<< (ostream& out, const Image& i);

#endif /* COLMAPIMAGE_H */

