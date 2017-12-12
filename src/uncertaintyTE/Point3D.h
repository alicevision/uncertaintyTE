// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   Point3D.h
 * Author: Michal Polic
 *
 * Created on October 25, 2017, 2:34 PM
 */

#ifndef POINT3D_H
#define POINT3D_H

#include <iostream>

using namespace std;

class Point3D {
public:
    int _id;
    double _X[3];
    
    Point3D();
    Point3D(int id, double* X);
    ~Point3D();
private:

};

ostream& operator<< (ostream& out, const Point3D& p);

#endif /* POINT3D_H */

