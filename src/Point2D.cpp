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
 * File:   Point2D.cpp
 * Author: policmic
 * 
 * Created on October 25, 2017, 3:49 PM
 */

#include "Point2D.h"

Point2D::Point2D() {}

Point2D::~Point2D() {}

bool Point2D::operator < (const Point2D& p) const {
	return (_point3D_id < p._point3D_id);
}

ostream& operator<< (ostream& out, const Point2D& p2D){
    out << ">>> xy: " << p2D._xy[0] << ", " << p2D._xy[1] << "-> p3D: " << p2D._point3D_id << "\n";
    return out;
};
