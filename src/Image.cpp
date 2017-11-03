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

Image::Image() {}

Image::~Image() {}


void Image::qt2aaRC(){
    
    
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