/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   OpenmvgIO.cpp
 * Author: root
 * 
 * Created on November 5, 2017, 11:30 AM
 */

#include "OpenmvgIO.h"
#include "auxCmd.h"

OpenmvgIO::OpenmvgIO() {}

OpenmvgIO::~OpenmvgIO() {}

//int OpenmvgIO::loadSceneOpenMVG(std::string sSfM_Data_Filename_In, openMVG::sfm::SfM_Data &sfm_data) {
//	// Load input SfM_Data scene
//	if (!openMVG::sfm::Load(sfm_data, sSfM_Data_Filename_In, openMVG::sfm::ESfM_Data(openMVG::sfm::ALL)))
//	{
//		std::cerr << std::endl
//			<< "The input SfM_Data file \"" << sSfM_Data_Filename_In << "\" cannot be read." << std::endl;
//		return EXIT_FAILURE;
//	}
//	return 1;
//}
//
///// Create the appropriate cost functor according the provided input camera intrinsic model
//ceres::CostFunction * IntrinsicsToCostFunction(openMVG::cameras::IntrinsicBase * intrinsic, const openMVG::Vec2 & observation)
//{
//	switch (intrinsic->getType())
//	{
//	case openMVG::cameras::PINHOLE_CAMERA:
//		return new ceres::AutoDiffCostFunction<openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic, 2, 3, 6, 3>(
//			new openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic(observation.data()));
//		break;
//	case openMVG::cameras::PINHOLE_CAMERA_RADIAL1:
//		return new ceres::AutoDiffCostFunction<openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1, 2, 4, 6, 3>(
//			new openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1(observation.data()));
//		break;
//	case openMVG::cameras::PINHOLE_CAMERA_RADIAL3:
//		return new ceres::AutoDiffCostFunction<openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3, 2, 6, 6, 3>(
//			new openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3(observation.data()));
//		break;
//	case openMVG::cameras::PINHOLE_CAMERA_BROWN:
//		return new ceres::AutoDiffCostFunction<openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2, 2, 8, 6, 3>(
//			new openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2(observation.data()));
//		break;
//	case openMVG::cameras::PINHOLE_CAMERA_FISHEYE:
//		return new ceres::AutoDiffCostFunction<openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye, 2, 7, 6, 3>(
//			new openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye(observation.data()));
//	case openMVG::cameras::PINHOLE_CAMERA_FISHEYE1:
//		return new ceres::AutoDiffCostFunction<openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye1, 2, 4, 6, 3>(
//			new openMVG::sfm::ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye1(observation.data()));
//	default:
//		return nullptr;
//	}
//}
//
//void openmvgSfM2Jacobian(openMVG::sfm::SfM_Data &sfm_data, ceres::CRSMatrix &jacobian, cov::Options &opt) {
//	ceres::Problem problem;
//
//	// blocks filled in Jacobian ( cameras R,t and 3D points )
//	std::vector<double*> parameter_blocks;
//	std::vector<double> mutable_points;
//
//	// Data wrapper for refinement:
//	openMVG::Hash_Map<openMVG::IndexT, std::vector<double> > map_poses;
//
//	// Setup Poses data & subparametrization  -- UNCOMMENT FOR RUN OPENMVG
//	/*for (openMVG::sfm::Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose){
//	const openMVG::IndexT indexPose = itPose->first;
//
//	const openMVG::geometry::Pose3 & pose = itPose->second;
//	const openMVG::Mat3 R = pose.rotation();
//	const openMVG::Vec3 t = pose.translation();
//
//	double angleAxis[3];
//	ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
//	map_poses[indexPose].reserve(6); //angleAxis + translation
//	map_poses[indexPose].push_back(angleAxis[0]);
//	map_poses[indexPose].push_back(angleAxis[1]);
//	map_poses[indexPose].push_back(angleAxis[2]);
//	map_poses[indexPose].push_back(t(0));
//	map_poses[indexPose].push_back(t(1));
//	map_poses[indexPose].push_back(t(2));
//
//	double * parameter_block = &map_poses[indexPose][0];
//	problem.AddParameterBlock(parameter_block, 6);
//	parameter_blocks.push_back(parameter_block);
//	}*/
//
//	// ----------------------------------------------------------------------------------------------------
//	// TODO: covarince propagation for diferent camera types ( usually not used in practise )
//	// => if the intrinsic parameters differ from camera to camera the covariance computation fail
//
//	openMVG::Hash_Map<openMVG::IndexT, std::vector<double> > map_intrinsics;
//	// Setup Intrinsics data & subparametrization
//	for (const auto& itIntrinsic : sfm_data.GetIntrinsics()) {
//		const openMVG::IndexT idIntrinsics = itIntrinsic.first;
//
//		assert(isValid(itIntrinsic.second->getType()));
//		map_intrinsics[idIntrinsics] = itIntrinsic.second->getParams();
//
//		double * parameter_block = &map_intrinsics[idIntrinsics][0];
//		problem.AddParameterBlock(parameter_block, map_intrinsics[idIntrinsics].size());
//
//		// Focal length
//		if (itIntrinsic.second->initialFocalLengthPix() > 0) {
//			// If we have an initial guess, we only authorize a margin around this value.
//			assert(map_intrinsics[idIntrinsics].size() >= 1);
//			const unsigned int maxFocalErr = 0.2 * (MAXIMUM)(itIntrinsic.second->w(), itIntrinsic.second->h());
//			problem.SetParameterLowerBound(parameter_block, 0, (double)itIntrinsic.second->initialFocalLengthPix() - maxFocalErr);
//			problem.SetParameterUpperBound(parameter_block, 0, (double)itIntrinsic.second->initialFocalLengthPix() + maxFocalErr);
//		}
//		else {
//			// We don't have an initial guess, but we assume that we use a converging lens, so the focal length should be positive.
//			problem.SetParameterLowerBound(parameter_block, 0, 0.0);
//		}
//
//		// Optical center - don't refine the optical center
//		std::vector<int> vec_constant_params;
//		vec_constant_params.push_back(1);
//		vec_constant_params.push_back(2);
//		ceres::SubsetParameterization *subset_parameterization =
//			new ceres::SubsetParameterization(map_intrinsics[idIntrinsics].size(), vec_constant_params);
//		problem.SetParameterization(parameter_block, subset_parameterization);
//	}
//	// ----------------------------------------------------------------------------------------------------
//
//	// Set a LossFunction to be less penalized by false measurements
//	//  - set it to NULL if you don't want use a lossFunction.
//	ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(openMVG::Square(4.0));
//	// TODO: make the LOSS function and the parameter an option
//
//	// For all visibility add reprojections errors:
//	for (auto& landmarkIt : sfm_data.structure)
//	{
//		const openMVG::sfm::Observations & observations = landmarkIt.second.observations;
//		// Iterate over 2D observation associated to the 3D landmark
//		for (const auto& observationIt : observations)
//		{
//			// Build the residual block corresponding to the track observation:
//			const openMVG::sfm::View * view = sfm_data.views.at(observationIt.first).get();
//
//			// Each Residual block takes a point and a camera as input and outputs a 2
//			// dimensional residual. Internally, the cost function stores the observed
//			// image location and compares the reprojection against the observation.
//			ceres::CostFunction* cost_function = NULL;
//			//IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), observationIt.second.x);
//
//			/*if (cost_function)
//			problem.AddResidualBlock(cost_function,
//			p_LossFunction,
//			&map_intrinsics[view->id_intrinsic][0],
//			&map_poses[view->id_pose][0],
//			landmarkIt.second.X.data()); */  //Do we need to copy 3D point to avoid false motion, if failure ?
//		}
//		parameter_blocks.push_back(landmarkIt.second.X.data());
//		mutable_points.push_back(landmarkIt.second.X(0));
//		mutable_points.push_back(landmarkIt.second.X(1));
//		mutable_points.push_back(landmarkIt.second.X(2));
//	}
//
//	// Configure Jacobian engine 
//	double cost = 0.0;
//	ceres::Problem::EvaluateOptions evalOpt;
//	evalOpt.parameter_blocks = parameter_blocks;
//	evalOpt.num_threads = 8;
//	evalOpt.apply_loss_function = true;
//
//	// create Jacobain 
//	problem.Evaluate(evalOpt, &cost, NULL, NULL, &jacobian);
//
//	// Configure covarivnce engine ( find the indexes of the most distatnt points etc. )
//	setPts2Fix(opt, mutable_points.size() / 3, mutable_points.data());
//	opt._numCams = sfm_data.views.size();
//	opt._camParams = 6;
//	opt._numPoints = sfm_data.structure.size();
//	opt._numObs = jacobian.num_rows / 2;
//	opt._algorithm = 2;
//	opt._epsilon = 1e-10;
//	opt._lambda = -1;
//	opt._svdRemoveN = -1;
//	opt._maxIterTE = -1;
//}

bool OpenmvgIO::read(const string input_dir, Scene& scene) {
#ifdef USE_OPENMVG
	std::cout << "Loading a OpenMVG scene: " << input_dir << '\n';
	openMVG::sfm::SfM_Data sfm_data;
	loadSceneOpenMVG(input_file, sfm_data);
	scene._jacobian = ceres::CRSMatrix();
	scene._options = cov::Options();
	openmvgSfM2Jacobian(sfm_data, scene._jacobian, scene._options);    // work for just 9 params for camera representation
#else
	std::cerr << "Load OpenMVG library to enable it's input.";
	exit(1);
#endif
}

bool OpenmvgIO::write(const string output_dir, Scene& scene) {
	std::cerr << "OpenMVG output is not implemented yet.";
	exit(1);
}