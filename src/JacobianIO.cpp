/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   JacobianIO.cpp
 * Author: policmic
 * 
 * Created on November 2, 2017, 5:22 PM
 */

#include "JacobianIO.h"

JacobianIO::JacobianIO() {}

JacobianIO::~JacobianIO() {}

int JacobianIO::data_type(){
    return _data_type;
}

bool JacobianIO::read(const string input_dir, Scene& scene){
    exit(1);
}

bool JacobianIO::write(const string output_dir, Scene& scene){
    exit(1);
}