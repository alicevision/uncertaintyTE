#include "compute.h"
#include "auxCmd.h"
#include "FactoryIO.h"
#include "IO.h"
#include "ColmapIO.h"
#include "JacobianIO.h"
#include "JacobianComposer.h"

#ifdef USE_MATLAB
    #include <mex.h>
    #include "matlabInterface.h"
    #include "uncertainty_mex.h"
#endif
#ifdef _WIN32
	#define EXECUTABLE_FILE "uncertainty.exe"
#elif __linux__ 
	#define EXECUTABLE_FILE "uncertainty"
#endif


#ifndef LOAD_CMD_IO_FLAGS
    #define LOAD_CMD_IO_FLAGS
    DEFINE_string(alg, "/media/policmic/DATA/MichalPolic/data/RedundantSFM/colmap_sparse_reconstruction",
            "algorithm for inversion of Schur complement matrix (SVD_QR_ITERATION, SVD_DEVIDE_AND_CONQUER, TAYLOR_EXPANSION)\n");
    DEFINE_string(in, "/media/policmic/DATA/MichalPolic/data/RedundantSFM/colmap_sparse_reconstruction",
            "path to input scene files (e.g. directory which contains cameras.txt, images.txt, points3D.tx for COLMAP;\n"
            " selected file <path_to_file>.jacob for Jacobian, etc. )\n");
    DEFINE_string(in_form, "COLMAP",
            "the format of input data (e.g. COLMAP, JACOBIAN, OPENMVG ...)\n");
    DEFINE_string(out, "/media/policmic/DATA/MichalPolic/data/RedundantSFM/colmap_sparse_reconstruction",
            "path to output covariance files\n");
    DEFINE_string(cam, "RADIAL",
            "the model of camera ( now just RADIAL ... in future PINHOLE, SIMPLE_RADIAL, OPENCV, FULL_OPENCV, ...)\n");
    
#endif




/*
Main function called from command line: unc.exe
  In: 
	- algorithm: 0 = SVD_QR_ITERATION, 1 = SVD_DEVIDE_AND_CONQUER, 2 = TAYLOR_EXPANSION
    - jacobian/openMVG_scene: path to the file without spaces
  Example Jacobain: unc.exe -in=input/myFile.jacob
  Example OpenMVG:  unc.exe -in=input/myFile.[bin,json,xml]  // TODO ...
*/
int main(int argc, char* argv[]) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    
    ceres::CRSMatrix jacobian = ceres::CRSMatrix();
    cov::Options options = cov::Options();
    cov::Statistic statistic = cov::Statistic();
    
    // read input data
    IO* io = FactoryIO::createIO(FLAGS_in_form);
    if (io->data_type() == SCENE_DATA){
        Scene scene = Scene();
        if ( !io->read( FLAGS_in, scene ) ) 
            exit(1);
        JacobianComposer jc = JacobianComposer();
        jc.scene2Jacobian( FLAGS_cam, scene, jacobian, options );
        
    }else{
        std::ifstream in_file(FLAGS_in, std::ios_base::in);
        loadJacobian(in_file, 2, jacobian, options);       // FLAGS_algorithm = 2
    }
    
    // TODO: add to the loader
    /*#ifdef USE_OPENMVG
    std::cout << "Loading a OpenMVG scene: " << input_file << '\n';
    openMVG::sfm::SfM_Data sfm_data;
    loadSceneOpenMVG(input_file, sfm_data);
    openmvgSfM2Jacobian(sfm_data, jacobian, options);
    #endif // USE_OPENMVG
     */
    
    
    options._svdRemoveN = 7;

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
    saveResults(FLAGS_out, options, statistic, num_camera_covar_values, camUnc, ptsUnc);
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
*/  
#ifdef USE_MATLAB
void MEX_FUNCTION_NAME(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	cov::Statistic statistic = cov::Statistic();
	cov::Options options = cov::Options();
	ceres::CRSMatrix jacobian = ceres::CRSMatrix();
	
	// Load and create Jacobian using camera model [angle axis, camera center, focal length, radial distortion]
	ceres::examples::BALProblem bal_problem = loadSceneMatlab( prhs );
	buildJacobian(&bal_problem, jacobian, options);

	// alocate output arrays
	int NcamArr = 0.5 * bal_problem.camera_block_size() * (bal_problem.camera_block_size() + 1) * bal_problem.num_cameras();
	int NptsArr = 6 * bal_problem.num_points();
	double *ptsUnc = (double*)malloc(NptsArr * sizeof(double));
	double *camUnc = (double*)malloc(NcamArr * sizeof(double));
	plhs[0] = mxCreateDoubleMatrix(NptsArr + NcamArr, 1, mxREAL);
	double *outC = mxGetPr(plhs[0]);

	// COMPUTE COVARIANCES
	computeCovariances(options, statistic, jacobian, camUnc, ptsUnc);

	// write results to the output
	for (int i = 0; i < NcamArr; i++)
		outC[i] = camUnc[i];
	for (int i = 0; i < NptsArr; i++)
		outC[NcamArr + i] = ptsUnc[i];

	free(ptsUnc);
	free(camUnc);
	mexPrintf("Computation of covariances ... [done]\n\n");
}
#endif
