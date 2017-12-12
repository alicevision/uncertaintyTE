// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   Point2D.h
 * Author: Michal Polic
 *
 * Created on October 25, 2017, 3:49 PM
 */

#ifndef POINT2D_H
#define POINT2D_H

#include <iostream>

using namespace std;

class Point2D {
public:
    int _point3D_id;
    double _xy[2];
    
    Point2D();
    ~Point2D();

	bool operator < (const Point2D& p) const;
private:

};

ostream& operator<< (ostream& out, const Point2D& p);

#endif /* POINT2D_H */

