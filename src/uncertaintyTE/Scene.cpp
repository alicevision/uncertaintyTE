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
 * File:   Scene.cpp
 * Author: policmic
 * 
 * Created on October 25, 2017, 10:08 AM
 */

#include "Scene.h"

Scene::Scene() {}

Scene::~Scene() {}

ostream& operator<< (ostream& out, const Scene& s){
    out << "Scene [cams:" << " - " << s._cameras.size() 
            << ", images:" << s._images.size() 
            << ", points3D:" << s._points3D.size() << "\n";

    // print first cameras
    int i = 0;
    for(auto const &cam : s._cameras) {
        if ( i < 10 ){
            out << cam.second;
            i++;
        }else break;
    }
    
    // print first images
    i = 0;
    for(auto const &img : s._images) {
        if ( i < 10 ){
            out << img.second;
            i++;
        }else break;
    }
    
    // print first images
    i = 0;
    for(auto const &pt3D : s._points3D) {
        if ( i < 10 ){
            out << pt3D.second;
            i++;
        }else break;
    }
    return out ;
};
