/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   RadialCamera.h
 * Author: root
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

