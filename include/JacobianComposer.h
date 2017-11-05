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

#include <random>

#include "uncertainty.h"
#include "Scene.h"
#include "snavely_reprojection_error.h"


class JacobianComposer {
public:
    // functions 
    static void scene2Jacobian(std::string cam_model, std::string algorithm, Scene &scene);
    
private:
	// constructors 
	JacobianComposer();
	~JacobianComposer();

    // functions 
	static void findPts2Fix(cov::Options &opt, int n, map<int, Point3D> &pts3D);
	static void Scene2Problem(Scene &s, ceres::Problem* problem, std::string cam_model);

};

#endif /* JACOBIANCOMPOSER_HPP */

