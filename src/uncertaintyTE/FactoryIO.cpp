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


IO* FactoryIO::createIO(const string type){
    string COLMAP = string("COLMAP");
    //string YASFM = string("YASFM");
    string OPENMVG = string("OPENMVG");
    string JACOBIAN = string("JACOBIAN");
    
    IO* io = NULL;
    if (COLMAP.compare(type) == 0){
        io = new ColmapIO();
    
    }else if(JACOBIAN.compare(type) == 0){
        io = new JacobianIO();
        
	}
	else if (OPENMVG.compare(type) == 0) {
		io = new OpenmvgIO();

	}
    // ...
    
    return io;
}