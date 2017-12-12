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
 * File:   ColmapImage.cpp
 * Author: policmic
 * 
 * Created on October 24, 2017, 4:22 PM
 */

#include "Image.h"
#include <math.h>

Image::Image() {}

Image::~Image() {}


void Image::qt2aaRC(){
    // convert quatrnion to rotation matrix 
	_R[0] = 1 - 2 * _q[2] * _q[2] - 2 * _q[3] * _q[3];
	_R[1] = 2 * _q[1] * _q[2] - 2 * _q[0] * _q[3];
	_R[2] = 2 * _q[3] * _q[1] + 2 * _q[0] * _q[2];
	_R[3] = 2 * _q[1] * _q[2] + 2 * _q[0] * _q[3];
	_R[4] = 1 - 2 * _q[1] * _q[1] - 2 * _q[3] * _q[3];
	_R[5] = 2 * _q[2] * _q[3] - 2 * _q[0] * _q[1];
	_R[6] = 2 * _q[3] * _q[1] - 2 * _q[0] * _q[2];
	_R[7] = 2 * _q[2] * _q[3] + 2 * _q[0] * _q[1];
	_R[8] = 1 - 2 * _q[1] * _q[1] - 2 * _q[2] * _q[2];
	
	// convert quatrnion to angle axis
	double sinth = sqrt(_q[1]*_q[1] + _q[2]*_q[2] + _q[3]*_q[3]);
	double angle = 2 * atan2(sinth, _q[0]);
	_aa[0] = angle * _q[1] / sinth;
	_aa[1] = angle * _q[2] / sinth;
	_aa[2] = angle * _q[3] / sinth;

    // convert translation to camera center
	_C[0] = -(_R[0] * _t[0] + _R[3] * _t[1] + _R[6] * _t[2]);
	_C[1] = -(_R[1] * _t[0] + _R[4] * _t[1] + _R[7] * _t[2]);
	_C[2] = -(_R[2] * _t[0] + _R[5] * _t[1] + _R[8] * _t[2]);
}



ostream& operator<< (ostream& out, const Image& i){
    out << "> Image [id:" << i._id << ", nobs:" << i._point2D.size() 
            << ", q:" << i._q[0] << "," << i._q[1] << "," << i._q[2] << "," << i._q[3] 
            << ", t:" << i._t[0] << "," << i._t[1] << "," << i._t[2] << "]\n";
    
    out << ">> Observations: \n" ;
    int k = 0;   
    for ( auto const& p2D : i._point2D ){
        if ( k < 10 ){
            out << p2D;
            k++;
        }else break;
    }
    return out ;
};