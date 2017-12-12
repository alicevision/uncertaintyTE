#include <uncertaintyTE/compute.h>
#include <uncertaintyTE/auxCmd.h>
#include <uncertaintyTE/FactoryIO.h>
#include <uncertaintyTE/IO.h>
#include <uncertaintyTE/ColmapIO.h>
#include <uncertaintyTE/JacobianIO.h>
#include <uncertaintyTE/JacobianComposer.h>

#include <gflags/gflags.h>

#ifdef USE_MATLAB
    #include <mex.h>
    #include <uncertaintyTE/matlabInterface.h>
    #include <uncertaintyTE/uncertainty_mex.h>
#endif
#ifdef _WIN32
	#define EXECUTABLE_FILE "uncertainty.exe"
#elif __linux__ 
	#define EXECUTABLE_FILE "uncertainty"
#endif


#ifndef LOAD_CMD_IO_FLAGS
    #define LOAD_CMD_IO_FLAGS
    DEFINE_string(alg, "TAYLOR_EXPANSION",
            "algorithm for inversion of Schur complement matrix [SVD_QR_ITERATION, SVD_DEVIDE_AND_CONQUER, TAYLOR_EXPANSION]");
    
    DEFINE_string(in, ".",
                    "path to input scene files (e.g. directory which contains cameras.txt, images.txt, points3D.tx for COLMAP;\n"
                    " or jacobian file <path_to_file>.jacob )");
   
    // TODO: enable the openMVG
    DEFINE_string(in_form, "COLMAP",
        "the format of input data [COLMAP, JACOBIAN, OPENMVG]");

    DEFINE_string(out, ".",
            "path to output covariance file");

    DEFINE_string(cam, "SIMPLE_RADIAL",
            "camera model ( SIMPLE_PINHOLE, PINHOLE, SIMPLE_RADIAL, RADIAL )");
    
    DEFINE_bool(debug, false,
            "Print all the used matrices to txt files.\n");
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
    cov::Statistic statistic = cov::Statistic();
    
    // read input data
    IO* io = FactoryIO::createIO(FLAGS_in_form);
    Scene scene = Scene();
    if ( !io->read( FLAGS_in, scene ) ) 
        exit(1);
        
    // compose Jacobian if not exist
    if( scene._jacobian.values.size() == 0 )
        JacobianComposer::scene2Jacobian( FLAGS_cam, FLAGS_alg, scene );

    // debug -> print matrices
    scene._options._debug = FLAGS_debug;

    scene._uncertainty.init(scene._options);
    
    // COMPUTE COVARIANCES
    computeCovariances(scene._options, statistic, scene._jacobian, &scene._uncertainty._camerasUnc[0], &scene._uncertainty._pointsUnc[0]);

    // write results to the outut file 
    io->writeCov2File(FLAGS_out, scene, statistic );

    std::cout << "Main function... [done]\n";
    return 0;
}



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
