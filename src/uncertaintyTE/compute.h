#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include "uncertainty.h"

#include <cuda.h>

#include <iostream>
#include <fstream>
#include <string>

#include <Eigen/Eigen>
#include <ctime>


#ifdef USE_OPENMVG
  #include "openMVG/sfm/sfm_data.hpp"
  #include "openMVG/sfm/sfm_data_io.hpp"
#endif


#include "thrust/sort.h" 
#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_operators.h"
//#include "magmasparse.h"

#ifdef _WIN32
	#include "magma_internal.h"
	#include "testings.h"
	#define MINIMUM min
	#define MAXIMUM max
#else
/*
	#define min
	#define max
	#include "magma_internal.h"
	#include "testings.h"
	#undef min
	#undef max
*/
	#define MINIMUM std::min
	#define MAXIMUM std::max
#define TESTING_CHECK( err )                                                 \
    do {                                                                     \
        magma_int_t err_ = (err);                                            \
        if ( err_ != 0 ) {                                                   \
            fprintf( stderr, "Error: %s\nfailed at %s:%d: error %lld: %s\n", \
                     #err, __FILE__, __LINE__,                               \
                     (long long) err_, magma_strerror(err_) );               \
            exit(1);                                                         \
        }                                                                    \
    } while( 0 )

#endif



class ScaledSparseMatrix;
typedef ScaledSparseMatrix SSM;

//////////////////////////////////////////////////////////////////////////////////////
// TIME MEASURING
#ifdef _WIN32		// ORIGINAL TIME MEASURING
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::time_point<std::chrono::system_clock> s_clock;
typedef std::chrono::steady_clock::time_point tp;

#else      // DISABLE <chrono> for measuring the time - time is alwas 0
typedef double tp;
struct Clock {
	static double now() {
		return 0;
	}
};
#endif
double timeDuration(tp from, tp to);
//////////////////////////////////////////////////////////////////////////////////////



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

void fixPts(tp *s, int *pts, cov::Options &opt, cov::Statistic &statistic, SSM *J);
