#ifndef COMPUTE_HEADER_INCLUDED
	#include "compute.h"
#endif
#ifndef AUXCMD_HEADER_INCLUDED
	#define AUXCMD_HEADER_INCLUDED
#endif

std::string algorihm2str(int alg);

void replaceAll(std::string& str, const std::string& from, const std::string& to);

void loadJacobian(std::ifstream& file, int algorithm, ceres::CRSMatrix& jacobian, cov::Options& options);

void saveResults(std::string& process_file_name, const std::string& current_dir, cov::Options& options, cov::Statistic& statistic,
	int num_camera_covar_values, double* camUnc);