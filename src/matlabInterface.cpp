#ifdef USE_MATLAB
#include <mex.h>
#include "compute.h"
#include "matlabInterface.h"
#include "snavely_reprojection_error.h"
#include "auxCmd.h"


void BuildProblem(ceres::examples::BALProblem *bal_problem, ceres::Problem* problem) {
	const int point_block_size = bal_problem->point_block_size();
	const int camera_block_size = bal_problem->camera_block_size();
	double* points = bal_problem->mutable_points();
	double* points_w = bal_problem->mutable_points_w();
	double* cameras = bal_problem->mutable_cameras();
	double* cameras_w = bal_problem->mutable_cameras_w();

	// Observations is 2*num_observations long array observations =
	// [u_1, u_2, ... , u_n], where each u_i is two dimensional, the x
	// and y positions of the observation.
	const double* observations = bal_problem->observations();

	for (int i = 0; i < bal_problem->num_observations(); ++i) {
		// Each observation correponds to a pair of a camera and a point
		// which are identified by camera_index()[i] and point_index()[i]
		// respectively.
		int cid = bal_problem->camera_index()[i];
		double* camera = cameras + camera_block_size * cid;
		double* camera_w = cameras_w + camera_block_size * cid;
		double* point = points + point_block_size * bal_problem->point_index()[i];
		double* point_w = points_w + point_block_size * bal_problem->point_index()[i];

		ceres::CostFunction* cost_function;
		// Each Residual block takes a point and a camera as input and
		// outputs a 2 dimensional residual.
		switch (bal_problem->bundle_type()[cid]) {
		case 1:
			cost_function = ceres::examples::SnavelyReprojectionErrorWithQuaternions::Create(observations[2 * i + 0], observations[2 * i + 1]);
			break;
		case 2:
			cost_function = ceres::examples::FixFocalRadialReprojectionError::Create(
				observations[2 * i], observations[2 * i + 1], camera[6], camera[7], camera[8]);
			break;
		case 3:
			cost_function = ceres::examples::FixCamsReprojectionError::Create(observations[2 * i], observations[2 * i + 1], camera);
			break;
		case 4:
			cost_function = ceres::examples::FixPtsReprojectionError::Create(observations[2 * i], observations[2 * i + 1], point);
			break;
		case 5:
			cost_function = ceres::examples::FixPtsFocalRadialReprojectionError::Create(observations[2 * i], observations[2 * i + 1], camera, point);
			break;
		case 6:
			cost_function = ceres::examples::FixOneCamCenterCoordinateFocalRadialReprojectionError::Create(observations[2 * i], observations[2 * i + 1], camera);
			break;
		case 7:
			cost_function = ceres::examples::FixFocalReprojectionError::Create(observations[2 * i], observations[2 * i + 1], camera);
			break;
		default:   // = case 0 
			cost_function = ceres::examples::SnavelyReprojectionError::Create(observations[2 * i + 0], observations[2 * i + 1], camera_w[0], camera_w[1],
				camera_w[2], camera_w[3], camera_w[4], camera_w[5], camera_w[6], camera_w[7], camera_w[8], point_w[0], point_w[1], point_w[2]);
			break;
		}

		// The loss function
		ceres::LossFunction* loss_function = NULL; // new ceres::HuberLoss(1.0);

		switch (bal_problem->bundle_type()[cid]) {
		case 1:
			problem->AddResidualBlock(cost_function, loss_function, camera, camera + 4, point);
			break;
		case 2:
			problem->AddResidualBlock(cost_function, loss_function, camera, point);
			break;
		case 3:
			problem->AddResidualBlock(cost_function, loss_function, point);
			break;
		case 4:
			problem->AddResidualBlock(cost_function, loss_function, camera);
			break;
		case 5:
			problem->AddResidualBlock(cost_function, loss_function, camera);
			break;
		case 6:
			problem->AddResidualBlock(cost_function, loss_function, camera, point);
			break;
		case 7:
			problem->AddResidualBlock(cost_function, loss_function, camera, point);
			break;
		default:   // = case 0 
			problem->AddResidualBlock(cost_function, loss_function, camera, point);
			break;
		}
	}
}


ceres::examples::BALProblem loadSceneMatlab(const mxArray *prhs[]) {
	int			n = mxGetNumberOfFields(prhs[0]);
	mxArray		*tmp = NULL;
	int			ncams;
	int			npts;
	double		*parameters = NULL;
	double		*parameters_w = NULL;
	bool		use_quaternions = false;
	long		nobs;
	int			*point_index = NULL;
	int			*camera_index = NULL;
	double		*observations = NULL;
	int			*bundle_type = NULL;
	char		*c_saveJfile = NULL;
	double		eps_or_lamb = -1;
	int			cov_alg = TAYLOR_EXPANSION;
	std::string saveJfile;

	/* read input structure */
	for (int i = 0; i < n; i++) {
		tmp = mxGetFieldByNumber(prhs[0], 0, i);
		std::string name = mxGetFieldNameByNumber(prhs[0], i);
		if (tmp == NULL) {
			mexPrintf("Field %s is empty.\n", name.c_str());
			continue;
		}
		else {
			//mexPrintf("Field %s is filled.\n", name);
		}

		if (name.compare("ncams") == 0) {
			ncams = (int)mxGetScalar(tmp);
		}
		if (name.compare("npts") == 0) {
			npts = (int)mxGetScalar(tmp);
		}
		if (name.compare("parameters") == 0) {
			parameters = mxGetPr(tmp);
		}

		if (name.compare("nobs") == 0) {
			nobs = (long)mxGetScalar(tmp);
		}
		if (name.compare("camera_index") == 0) {
			camera_index = (int *)mxGetData(tmp);
		}
		if (name.compare("point_index") == 0) {
			point_index = (int *)mxGetData(tmp);
		}
		if (name.compare("observations") == 0) {
			observations = mxGetPr(tmp);
		}
		if (name.compare("bundle_type") == 0) {
			bundle_type = (int *)mxGetData(tmp);
		}


		if (name.compare("eps_or_lamb") == 0) {
			eps_or_lamb = *mxGetPr(tmp);
		}
		if (name.compare("parameters_w") == 0) {
			parameters_w = mxGetPr(tmp);
		}
		if (name.compare("save_jacobian2file") == 0) {
			c_saveJfile = mxArrayToString(tmp);
			saveJfile = c_saveJfile;
		}
	}

	// build BA problem object
	return ceres::examples::BALProblem(ncams, npts, parameters, use_quaternions,
		nobs, point_index, camera_index, observations, bundle_type, parameters_w);
}


void buildJacobian(ceres::examples::BALProblem *bal_problem, ceres::CRSMatrix &jacobian, cov::Options &opt) {
	ceres::Problem problem;
	BuildProblem(bal_problem, &problem);

	std::vector<double*> parameter_blocks;
	for(int i = 0; i < bal_problem->num_cameras(); i++)
		parameter_blocks.push_back(bal_problem->mutable_cameras() + bal_problem->camera_block_size() * i);
	for (int i = 0; i < bal_problem->num_points(); i++)
		parameter_blocks.push_back(bal_problem->mutable_points() + 3 * i);

	// Configure Jacobian engine 
	double cost = 0.0;
	ceres::Problem::EvaluateOptions evalOpt;
	evalOpt.parameter_blocks = parameter_blocks;
	evalOpt.num_threads = 8;
	evalOpt.apply_loss_function = true;

	// create Jacobain 
	problem.Evaluate(evalOpt, &cost, NULL, NULL, &jacobian);

	// Configure covarivnce engine ( find the indexes of the most distatnt points etc. )
	setPts2Fix(opt, bal_problem->num_points(), bal_problem->mutable_points());

	opt._numCams = bal_problem->num_cameras();
	opt._camParams = bal_problem->camera_block_size();
	opt._numPoints = bal_problem->num_points();
	opt._numObs = bal_problem->num_observations();
	opt._algorithm = 2;
	opt._epsilon = 1e-10;
	opt._lambda = -1;
	opt._svdRemoveN = -1;
	opt._maxIterTE = -1;
}

#endif