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
 * File:   Point2D.h
 * Author: policmic
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
private:

};

ostream& operator<< (ostream& out, const Point2D& p);

#endif /* POINT2D_H */

