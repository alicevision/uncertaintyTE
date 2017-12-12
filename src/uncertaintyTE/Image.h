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
 * File:   ColmapImage.h
 * Author: policmic
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

