// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   OpenmvgIO.h
 * Author: Michal Polic
 *
 * Created on November 5, 2017, 11:30 AM
 */

#ifndef OPENMVGIO_H
#define OPENMVGIO_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "IO.h"
#include "Scene.h"
#include "compute.h"


class OpenmvgIO : public IO {
public:
    OpenmvgIO();
    ~OpenmvgIO();
    
	// Functions and methods
    bool read(const std::string& input_dir, Scene& scene);
    bool write(const std::string& output_dir, Scene& scene);

	int data_type();
private:

};


#endif /* OPENMVGIO_H */

