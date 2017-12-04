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
 * File:   FactoryIO.h
 * Author: policmic
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

