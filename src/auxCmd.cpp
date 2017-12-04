#include "auxCmd.h"
#include "ceres/rotation.h"

#ifdef USE_OPENMVG
  #include "openMVG/sfm/sfm_data_BA_ceres_camera_functor.hpp"
#endif

#ifdef USE_MATLAB
  #include <mex.h>
#endif

#include <random>

#ifndef DBL_MIN
#define DBL_MIN -1e999
#endif

std::string algorihm2str(int alg) {
	switch (alg) {
	case 0: return "SVD_QR_ITERATION";
	case 1: return "SVD_DEVIDE_AND_CONQUER";
	case 2: return "TAYLOR_EXPANSION";
	default: return "not defined";
	}
}

void replaceAll(std::string& str, const std::string& from, const std::string& to) {
	if (from.empty())
		return;
	size_t start_pos = 0;
	while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
}



inline double dist(double *p1, double *p2) {
	return sqrt((p1[0] - p2[0])*(p1[0] - p2[0]) + (p1[1] - p2[1])*(p1[1] - p2[1]) + (p1[2] - p2[2])*(p1[2] - p2[2]));
}

double dist(double *p1, double *p2, double *p3) {
	return dist(p1, p2) + dist(p1, p3) + dist(p2, p3);
}

// find the most distant points
void setPts2Fix(cov::Options &opt, int N, double *pts) {
	//int N = bal_problem->num_points();
	//double *pts = bal_problem->mutable_points();

	// the simplest variant of RANSAC - find triple of most distant points
	double max_dist = DBL_MIN;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, N);
	int i = 0;
	while (i < 100000) {
		int p1 = floor(dis(gen));
		int p2 = floor(dis(gen));
		int p3 = floor(dis(gen));
		if (p1 == p2 || p2 == p3 || p1 == p3) continue;
		double d = dist(pts + 3 * p1, pts + 3 * p2, pts + 3 * p3);
		if (d > max_dist) {
			max_dist = d;
			opt._pts2fix = new int[3]{ p1, p2, p3 };
		}
		i++;
	}
	std::sort(opt._pts2fix, opt._pts2fix + 2);
}

void printJacobian(ceres::CRSMatrix &J) {
	std::cout << "\n\nJ = zeros(" << J.num_rows << ", " << J.num_cols << ");\n";
	for (int i = 0; i < J.num_rows; i++) {
		for (int j = J.rows[i]; j < J.rows[i + 1]; j++) {
			std::cout << "J(" << (i + 1) << "," << (J.cols[j] + 1) << ") = " << J.values[j] << ";";
		}
	}
	std::cout << "\n\n\n";
}

#ifdef USE_MATLAB
void printJacobianMEX(ceres::CRSMatrix &J) {
	mexPrintf("\n\nJ = zeros(%d, %d);\n", J.num_rows, J.num_cols);
	for (int i = 0; i < J.num_rows; i++) {
		for (int j = J.rows[i]; j < J.rows[i + 1]; j++) {
			mexPrintf("J(%d,%d) = %f;", (i + 1), (J.cols[j] + 1), J.values[j]);
		}
	}
	mexPrintf("\n\n\n");
}
#endif


#ifdef USE_OPENMVG
int loadSceneOpenMVG(std::string sSfM_Data_Filename_In, openMVG::sfm::SfM_Data &sfm_data) {
	// Load input SfM_Data scene
	if (!openMVG::sfm::Load(sfm_data, sSfM_Data_Filename_In, openMVG::sfm::ESfM_Data(openMVG::sfm::ALL)))
	{
		std::cerr << std::endl
			<< "The input SfM_Data file \"" << sSfM_Data_Filename_In << "\" cannot be read." << std::endl;
		return EXIT_FAILURE;
	}
	return 1;
}

/// Create the appropriate cost functor according the provided input camera intrinsic model
ceres::CostFunction * IntrinsicsToCostFunction(openMVG::cameras::IntrinsicBase * intrinsic, const openMVG::Vec2 & observation)
{
	switch (intrinsic->getType())
	{
	case openMVG::cameras::PINHOLE_CAMERA:
		return new ceres::AutoDiffCostFunction<openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic, 2, 3, 6, 3>(
			new openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic(observation.data()));
		break;
	case openMVG::cameras::PINHOLE_CAMERA_RADIAL1:
		return new ceres::AutoDiffCostFunction<openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1, 2, 4, 6, 3>(
			new openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1(observation.data()));
		break;
	case openMVG::cameras::PINHOLE_CAMERA_RADIAL3:
		return new ceres::AutoDiffCostFunction<openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3, 2, 6, 6, 3>(
			new openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3(observation.data()));
		break;
	case openMVG::cameras::PINHOLE_CAMERA_BROWN:
		return new ceres::AutoDiffCostFunction<openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2, 2, 8, 6, 3>(
			new openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2(observation.data()));
		break;
	case openMVG::cameras::PINHOLE_CAMERA_FISHEYE:
		return new ceres::AutoDiffCostFunction<openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye, 2, 7, 6, 3>(
			new openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye(observation.data()));
	case openMVG::cameras::PINHOLE_CAMERA_FISHEYE1:
		return new ceres::AutoDiffCostFunction<openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye1, 2, 4, 6, 3>(
			new openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye1(observation.data()));
	default:
		return nullptr;
	}
}

void openmvgSfM2Jacobian(openMVG::sfm::SfM_Data &sfm_data, ceres::CRSMatrix &jacobian, cov::Options &opt){
	ceres::Problem problem;

	// blocks filled in Jacobian ( cameras R,t and 3D points )
	std::vector<double*> parameter_blocks;
	std::vector<double> mutable_points;

	// Data wrapper for refinement:
	openMVG::Hash_Map<openMVG::IndexT, std::vector<double> > map_poses;

	// Setup Poses data & subparametrization  -- UNCOMMENT FOR RUN OPENMVG
	for (openMVG::sfm::Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose){
		const openMVG::IndexT indexPose = itPose->first;

		const openMVG::geometry::Pose3 & pose = itPose->second;
		const openMVG::Mat3 R = pose.rotation();
		const openMVG::Vec3 t = pose.translation();

		double angleAxis[3];
		ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
		map_poses[indexPose].reserve(6); //angleAxis + translation
		map_poses[indexPose].push_back(angleAxis[0]);
		map_poses[indexPose].push_back(angleAxis[1]);
		map_poses[indexPose].push_back(angleAxis[2]);
		map_poses[indexPose].push_back(t(0));
		map_poses[indexPose].push_back(t(1));
		map_poses[indexPose].push_back(t(2));

		double * parameter_block = &map_poses[indexPose][0];
		problem.AddParameterBlock(parameter_block, 6);
		parameter_blocks.push_back(parameter_block);
	}

	// ----------------------------------------------------------------------------------------------------
	// TODO: covarince propagation for diferent camera types ( usually not used in practise )
	// => if the intrinsic parameters differ from camera to camera the covariance computation fail

	openMVG::Hash_Map<openMVG::IndexT, std::vector<double> > map_intrinsics;
	// Setup Intrinsics data & subparametrization
	for (const auto& itIntrinsic : sfm_data.GetIntrinsics()) {
		const openMVG::IndexT idIntrinsics = itIntrinsic.first;

		assert(isValid(itIntrinsic.second->getType()));
		map_intrinsics[idIntrinsics] = itIntrinsic.second->getParams();

		double * parameter_block = &map_intrinsics[idIntrinsics][0];
		problem.AddParameterBlock(parameter_block, map_intrinsics[idIntrinsics].size());

		// Focal length
		if (itIntrinsic.second->initialFocalLengthPix() > 0) {
			// If we have an initial guess, we only authorize a margin around this value.
			assert(map_intrinsics[idIntrinsics].size() >= 1);
			const unsigned int maxFocalErr = 0.2 * (MAXIMUM)(itIntrinsic.second->w(), itIntrinsic.second->h());
			problem.SetParameterLowerBound(parameter_block, 0, (double)itIntrinsic.second->initialFocalLengthPix() - maxFocalErr);
			problem.SetParameterUpperBound(parameter_block, 0, (double)itIntrinsic.second->initialFocalLengthPix() + maxFocalErr);
		}
		else {
			// We don't have an initial guess, but we assume that we use a converging lens, so the focal length should be positive.
			problem.SetParameterLowerBound(parameter_block, 0, 0.0);
		}

		// Optical center - don't refine the optical center
		std::vector<int> vec_constant_params;
		vec_constant_params.push_back(1);
		vec_constant_params.push_back(2);
		ceres::SubsetParameterization *subset_parameterization =
			new ceres::SubsetParameterization(map_intrinsics[idIntrinsics].size(), vec_constant_params);
		problem.SetParameterization(parameter_block, subset_parameterization);
	}
	// ----------------------------------------------------------------------------------------------------

	// Set a LossFunction to be less penalized by false measurements
	//  - set it to NULL if you don't want use a lossFunction.
	ceres::LossFunction * p_LossFunction = NULL; // new ceres::HuberLoss(openMVG::Square(4.0));
	// TODO: make the LOSS function and the parameter an option

	// For all visibility add reprojections errors:
	for (auto& landmarkIt : sfm_data.structure)
	{
		const openMVG::sfm::Observations & observations = landmarkIt.second.observations;
		// Iterate over 2D observation associated to the 3D landmark
		for (const auto& observationIt : observations)
		{
			// Build the residual block corresponding to the track observation:
			const openMVG::sfm::View * view = sfm_data.views.at(observationIt.first).get();

			// Each Residual block takes a point and a camera as input and outputs a 2
			// dimensional residual. Internally, the cost function stores the observed
			// image location and compares the reprojection against the observation.
			ceres::CostFunction* cost_function = IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), observationIt.second.x);

			if (cost_function)
				problem.AddResidualBlock(cost_function,
					p_LossFunction,
					&map_intrinsics[view->id_intrinsic][0],
					&map_poses[view->id_pose][0],
					landmarkIt.second.X.data());   //Do we need to copy 3D point to avoid false motion, if failure ?
		}
		parameter_blocks.push_back(landmarkIt.second.X.data());
		mutable_points.push_back(landmarkIt.second.X(0));
		mutable_points.push_back(landmarkIt.second.X(1));
		mutable_points.push_back(landmarkIt.second.X(2));
	}

	// Configure Jacobian engine 
	double cost = 0.0;
	ceres::Problem::EvaluateOptions evalOpt;
	evalOpt.parameter_blocks = parameter_blocks;
	evalOpt.num_threads = 8;
	evalOpt.apply_loss_function = true;

	// create Jacobain 
	problem.Evaluate(evalOpt, &cost, NULL, NULL, &jacobian);

	// Configure covarivnce engine ( find the indexes of the most distatnt points etc. )
	setPts2Fix(opt, mutable_points.size() / 3, mutable_points.data());
	opt._numCams = sfm_data.views.size();
	opt._camParams = 6;
	opt._numPoints = sfm_data.structure.size();
	opt._numObs = jacobian.num_rows / 2;
	opt._algorithm = 2;
	opt._epsilon = 1e-10;
	opt._lambda = -1;
	opt._svdRemoveN = -1;
	opt._maxIterTE = -1;
}
#endif

