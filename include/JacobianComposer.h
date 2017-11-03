/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   JacobianComposer.hpp
 * Author: policmic
 *
 * Created on November 2, 2017, 4:41 PM
 */

#ifndef JACOBIANCOMPOSER_HPP
#define JACOBIANCOMPOSER_HPP

#include "uncertainty.h"
#include "Scene.h"


class JacobianComposer {
public:
    JacobianComposer();
    ~JacobianComposer();
    
    // functions 
    void scene2Jacobian(std::string cam_model, Scene &scene, ceres::CRSMatrix &jacobian, cov::Options &options);
    
private:
    // functions 
    void findPts2Fix(cov::Options &opt, int n, map<int, Point3D> &pts3D);
    void Scene2Problem(Scene &s, ceres::Problem* problem);

};

#endif /* JACOBIANCOMPOSER_HPP */

