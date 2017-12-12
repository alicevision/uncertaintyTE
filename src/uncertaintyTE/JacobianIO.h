/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   JacobianIO.h
 * Author: policmic
 *
 * Created on November 2, 2017, 5:22 PM
 */

#ifndef JACOBIANIO_H
#define JACOBIANIO_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "IO.h"
#include "Scene.h"


class JacobianIO : public IO {
public:
    JacobianIO();
    ~JacobianIO();
    
    bool read(const std::string& input_dir, Scene& scene);
    bool write(const std::string& output_dir, Scene& scene);
    
    int data_type();
	
private:
    string _jacobian_path;
    int _data_type = JACOBIAN_DATA;

    bool loadJacobian(std::ifstream& file, cov::EAlgorithm algorithm, ceres::CRSMatrix& jacobian, cov::Options& options);
};

#endif /* JACOBIANIO_H */

