#ifndef COMPUTE_HEADER_INCLUDED
	#include "compute.h"
#endif
#ifndef AUXCMD_HEADER_INCLUDED
	#include "auxCmd.h"
#endif

#ifdef _WIN32
	#define EXECUTABLE_FILE "uncertainty.exe"
#elif __linux__ 
	#define EXECUTABLE_FILE "uncertainty"
#endif

/*
Main function called from command line: unc.exe
  In: 
	- algorithm: 0 = SVD_QR_ITERATION, 1 = SVD_DEVIDE_AND_CONQUER, 2 = TAYLOR_EXPANSION
    - jacobian: path to the file without spaces
  Example: unc.exe 2 input/myFile.jacob
*/
int main(int argc, char* argv[]) {
	ceres::CRSMatrix jacobian = ceres::CRSMatrix();
	cov::Options options = cov::Options();
	cov::Statistic statistic = cov::Statistic();
	options._svdRemoveN = 7;

	// Read reconstruction values from file (jacobian, number of cameras, .... )
	std::string current_exec_name = argv[0];
	std::string process_file_name = argv[2];
	const std::string current_dir = current_exec_name.substr(0, current_exec_name.find(EXECUTABLE_FILE));
	std::ifstream file(current_dir + process_file_name, std::ios_base::in);
	loadJacobian(file, std::stod(argv[1]), jacobian, options);

	// alocate output arrays
	int num_camera_covar_values = 0.5 * options._camParams * (options._camParams + 1);   // save only half of each symmetric matrix
	int camUnc_size = num_camera_covar_values * options._numCams;
	double* camUnc = (double*)malloc(camUnc_size * sizeof(double));
	assert(camUnc != NULL);
	double* ptsUnc = (double*)malloc(6 * options._numPoints * sizeof(double));
	assert(ptsUnc != NULL);

	// COMPUTE COVARIANCES
	computeCovariances(options, statistic, jacobian, camUnc, ptsUnc);

	// write results to the outut file 
	saveResults(process_file_name, current_dir, options, statistic, num_camera_covar_values, camUnc);
	std::cout << "Main function... [done]\n";
	return 0;
}


/*
Library function foc C++ call: getCovariances( ... )
In:
  - options: informations about the reconstruction (numCams, camParams, numPoints, numObs)
  - statistic: blank object for the output statistics
  - jacobian: sparse matrix (form Ceres-solver) with jacobian ( you can use the output of the computation of jacobian in Ceres as the input )
  - h_camUnc: array which contatins covariances for cameras
  - h_ptUnc: array which contatins covariances for points
*/
#ifdef _WIN32
	extern "C" __declspec(dllexport) void getCovariances(
		cov::Options &options,
		cov::Statistic &statistic,
		ceres::CRSMatrix &jacobian,
		double* h_camUnc,
		double* h_ptUnc)
	{
		computeCovariances(options, statistic, jacobian, h_camUnc, h_ptUnc);
	}
#elif __linux__ 
	void getCovariances(
		cov::Options &options,
		cov::Statistic &statistic,
		ceres::CRSMatrix &jacobian,
		double* h_camUnc,
		double* h_ptUnc)
	{
		computeCovariances(options, statistic, jacobian, h_camUnc, h_ptUnc);
	}
#endif


/*
Matlab interface, call in Matlab: unc( .... ) 
TODO: comming soon ...
*/