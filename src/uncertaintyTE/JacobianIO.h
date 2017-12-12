// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   JacobianIO.h
 * Author: Michal Polic
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

