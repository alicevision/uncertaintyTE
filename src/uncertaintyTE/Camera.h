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
 * File:   ColmapCamera.h
 * Author: policmic
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

