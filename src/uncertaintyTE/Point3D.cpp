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
 * File:   Point3D.cpp
 * Author: policmic
 * 
 * Created on October 25, 2017, 2:34 PM
 */

#include "Point3D.h"

Point3D::Point3D() {}

Point3D::Point3D(int id, double* X): _id(id){
    _X[0] = X[0];
    _X[1] = X[1];
    _X[2] = X[2];
}

Point3D::~Point3D() {}

ostream& operator<< (ostream& out, const Point3D& p3D){
    out << "> Point3D [" << p3D._id << ", " << 
            p3D._X[0] << ", " << p3D._X[1] << ", " << p3D._X[2] << "]\n";
    return out;
};