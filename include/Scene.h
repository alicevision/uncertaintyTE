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
 * File:   Scene.h
 * Author: policmic
 *
 * Created on October 25, 2017, 10:08 AM
 */

#ifndef SCENE_H
#define SCENE_H

#include <map>
#include "Camera.h"
#include "Image.h"
#include "Point3D.h"
#include "Point2D.h"

class Scene {
public:
    map<int, Camera> _cameras;
    map<int, Image> _images;
    map<int, Point3D> _points3D;
    
    Scene();
    ~Scene();
private:

};

ostream& operator<< (ostream& out, const Scene& s);

#endif /* SCENE_H */

