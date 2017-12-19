// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   FactoryIO.cpp
 * Author: Michal Polic
 * 
 * Created on October 25, 2017, 10:40 AM
 */

#include "FactoryIO.h"


IO* FactoryIO::createIO(const std::string& type){
    std::string COLMAP = string("COLMAP");
    //std::string YASFM = string("YASFM");
    std::string ALICEVISION = string("ALICEVISION");
    std::string JACOBIAN = string("JACOBIAN");
    
    IO* io = NULL;
    if (COLMAP.compare(type) == 0){
        io = new ColmapIO();
    
    }else if(JACOBIAN.compare(type) == 0){
        io = new JacobianIO();
	}
    else if (ALICEVISION.compare(type) == 0) {
		io = new AliceVisionIO();
    }
    else {
        throw std::runtime_error("Failed to create the IO for: " + type);
    }
    
    return io;
}
