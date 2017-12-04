/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   OpenmvgIO.h
 * Author: root
 *
 * Created on November 5, 2017, 11:30 AM
 */

#ifndef OPENMVGIO_H
#define OPENMVGIO_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "IO.h"
#include "Scene.h"
#include "compute.h"


class OpenmvgIO : public IO {
public:
    OpenmvgIO();
    ~OpenmvgIO();
    
	// Functions and methods
	bool read(const string input_dir, Scene& scene);
	bool write(const string output_dir, Scene& scene);

	int data_type();
private:

};


#endif /* OPENMVGIO_H */

