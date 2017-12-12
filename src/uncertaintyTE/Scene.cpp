// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   Scene.cpp
 * Author: Michal Polic
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
