// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   Point2D.cpp
 * Author: Michal Polic
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
