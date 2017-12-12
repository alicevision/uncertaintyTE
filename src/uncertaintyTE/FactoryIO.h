// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/*
* File:   auxCmd.h
* Author: Michal Polic
*
* Created on October 25, 2017, 10:40 AM
*/

#ifndef FACTORYIO_H
#define FACTORYIO_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>

#include "IO.h"
#include "ColmapIO.h"
#include "JacobianIO.h"
#include "OpenmvgIO.h"

using namespace std;

class FactoryIO {
public:
    static IO* createIO(const string type);
};

#endif /* FACTORYIO_H */

