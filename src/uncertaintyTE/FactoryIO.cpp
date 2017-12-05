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
 * File:   FactoryIO.cpp
 * Author: policmic
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