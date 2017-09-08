#include "uncertainty.h"

#include <cuda.h>

#include "thrust/sort.h" 
#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_internal.h"
#include "magma_operators.h"
#include "magmasparse.h"
#include "testings.h"

#include <iostream>
#include <fstream>
#include <string>

#include <Eigen/Eigen>
#include <ctime>

#include <mex.h>

#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"


#ifndef COMPUTE_HEADER_INCLUDED
	#define COMPUTE_HEADER_INCLUDED
#endif

#define LEFT 0
#define RIGHT 1
#define checkCudaErrors(val) check( (val), #val, __FILE__, __LINE__)

template<typename T>
void check(T err, const char* const func, const char* const file, const int line) {
	if (err != cudaSuccess) {
		std::cerr << "CUDA error at: " << file << ":" << line << std::endl;
		std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
		exit(1);
	}
}

void computeCovariances(
	cov::Options &options,
	cov::Statistic &statistics,
	ceres::CRSMatrix &jacobian,
	double *camUnc, 
	double *ptUnc);
