// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   ColmapIO.cpp
 * Author: Michal Polic
 * 
 * Created on October 24, 2017, 2:33 PM
 */

#include <vector>

#include "ColmapIO.h"

ColmapIO::ColmapIO() {}
ColmapIO::~ColmapIO() {}

int ColmapIO::data_type(){
    return _data_type;
}

bool ColmapIO::read(const std::string& input_dir, Scene& scene){
    cout << "Read COLMAP reconstruction from: " << input_dir << "\n";

    if (readCameras(string(input_dir + "/cameras.txt"), scene)  && 
        readImages(string(input_dir + "/images.txt"), scene)    &&
        readPoints3D(string(input_dir + "/points3D.txt"), scene))
    {
        cout << "Reading ... [done]\n";
        return true;
    }
    else
        return false;
}

bool ColmapIO::readCameras(const std::string& file_path, Scene& scene){
    ifstream file(file_path, ios_base::in);
    if (!file.good()) {
        cerr << "\nERROR: The file '" << file_path << "' does not exist!\n";
        return false;
    }
    
    // read lines
    string line = string();
    while ( !file.eof() ) {
        line.clear();
        getline(file,line);
        if (line.size() > 1 && line.at(0) != '#'){
            Camera c = Camera();
            istringstream line_stream(line);
            line_stream >> c._id >> c._model >> c._img_width >> c._img_height >> c._f >> c._uv[0] >> c._uv[1]; 
            
            // additional parameters
            if (c._model.compare("SIMPLE_RADIAL")) 
                line_stream >> c._r[0];
            if (c._model.compare("RADIAL")) 
                line_stream >> c._r[0] >> c._r[1];
            // ...
            
            // save the camera to scene
            scene._cameras[c._id] = c;
        }
    }

    file.close();
    return true;
}

// Images must be readed after the cameras -> camera f,r is included inside for now
bool ColmapIO::readImages(const std::string& file_path, Scene& scene){
    ifstream file(file_path, ios_base::in);
    if (!file.good()) {
        cerr << "\nERROR: The file '" << file_path << "' does not exist!\n";
        return false;
    }
    
    // read lines
    string line = string();
    while ( !file.eof() ) {
        line.clear();
        getline(file,line);
        if (line.size() > 1 && line.at(0) != '#'){
            Image i = Image();
            istringstream line_stream(line);
            
            // read an image - IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, NAME
            line_stream >> i._id >> i._q[0] >> i._q[1] >> i._q[2] >> i._q[3] 
                    >> i._t[0] >> i._t[1] >> i._t[2] >> i._cam_id; 
            
            // image points 2D - POINTS2D[] as (X, Y, POINT3D_ID)
            line.clear();
            getline(file,line);
            istringstream line_stream2(line);
            while (!line_stream2.eof()){
                Point2D p2d = Point2D();
                line_stream2 >> p2d._xy[0] >> p2d._xy[1] >> p2d._point3D_id;
                if (p2d._point3D_id != -1)
                    i._point2D.push_back(p2d);
            }
            
			// TODO: change this model of camera -> e.g. one f, r for more cameras 
			i._f = scene._cameras[i._cam_id]._f;
			i._r[0] = scene._cameras[i._cam_id]._r[0];
			i._r[1] = scene._cameras[i._cam_id]._r[1];

			// Fill values not mentioned in images.txt
			i.qt2aaRC();

            // save image to scene
            scene._images[i._id] = i;
        }
    }

    file.close();
    return true;
}

bool ColmapIO::readPoints3D(const std::string& file_path, Scene& scene){
    ifstream file(file_path, ios_base::in);
    if (!file.good()) {
        cerr << "\nERROR: The file '" << file_path << "' does not exist!\n";
        return false;
    }
    
    // read lines
    string line = string();
    while ( !file.eof() ) {
        line.clear();
        getline(file,line);
        if (line.size() > 1 && line.at(0) != '#'){
            Point3D p = Point3D();
            istringstream line_stream(line);
            
            // load one point3D - POINT3D_ID, X, Y, Z, R, G, B, ERROR, TRACK[] as (IMAGE_ID, POINT2D_IDX)
            line_stream >> p._id >> p._X[0] >> p._X[1] >> p._X[2]; 
            
            // save the point to scene
            scene._points3D[p._id] = p;
        }
    }

    file.close();
    return true;
}


bool ColmapIO::write(const std::string& output_dir, Scene& scene){
    cout << "Write COLMAP reconstruction to: " << output_dir << "\n";
    
    
    cout << "Writing ... [done]\n";
    return false;
}
