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
 * File:   ColmapCamera.cpp
 * Author: policmic
 * 
 * Created on October 24, 2017, 4:27 PM
 */

#include "Camera.h"

Camera::Camera() {}

Camera::~Camera() {}

ostream& operator<< (ostream& out, const Camera& c){
    out << "> Camera " << c._id << " " << c._model
            << " [focal:" << c._f << ", width:" << c._img_width << ", height:" << c._img_height
            << ", uv:" << c._uv[0] << "," << c._uv[2] 
            << ", rad_dist:" << c._r[0] << "," << c._r[1] << "," << c._r[2] << "]\n";
    return out;
};
