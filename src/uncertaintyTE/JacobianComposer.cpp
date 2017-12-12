// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/* 
 * File:   JacobianComposer.cpp
 * Author: Michal Polic
 * 
 * Created on November 2, 2017, 4:41 PM
 */

#include "uncertainty.h"
#include "JacobianComposer.h"

#ifndef DBL_MIN
    #define DBL_MIN -1.7e+308
#endif

#define SIMPLE_PINHOLE std::string("SIMPLE_PINHOLE")
#define PINHOLE std::string("PINHOLE")
#define SIMPLE_RADIAL std::string("SIMPLE_RADIAL")
#define RADIAL std::string("RADIAL")


JacobianComposer::JacobianComposer() {}

JacobianComposer::~JacobianComposer() {}
 

inline double dist(Point3D &p1, Point3D &p2) {
	return sqrt((p1._X[0] - p2._X[0])*(p1._X[0] - p2._X[0]) + 
				(p1._X[1] - p2._X[1])*(p1._X[1] - p2._X[1]) + 
				(p1._X[2] - p2._X[2])*(p1._X[2] - p2._X[2]));
}

double dist(Point3D &p1, Point3D &p2, Point3D &p3) {
	return dist(p1, p2) + dist(p1, p3) + dist(p2, p3);
}

void JacobianComposer::findPts2Fix(cov::Options &opt, int N, double* points3Darr) {
	// refactor points into the structure Point3D
	std::vector<Point3D> pts3D;
	for (int i = 0; i < N; i++)
		pts3D.push_back(Point3D(i, points3Darr + (3*i)));

	// the simplest variant of RANSAC -> find triple of most distant points
	findPts2Fix(opt, N, pts3D);
}


void JacobianComposer::findPts2Fix(cov::Options &opt, int N, map<int, Point3D> &points3D){
    // Convert points to vector
    std::vector<Point3D> pts3D = std::vector<Point3D>();
    pts3D.reserve(N);
    for(auto const &p3D : points3D)
        pts3D.push_back(p3D.second);
    
    // the simplest variant of RANSAC -> find triple of most distant points
	findPts2Fix(opt, N, pts3D);
}

void JacobianComposer::findPts2Fix(cov::Options &opt, int N, std::vector<Point3D> &pts3D) {
	// the simplest variant of RANSAC -> find triple of most distant points
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
		double d = dist(pts3D[p1], pts3D[p2], pts3D[p3]);
		if (d > max_dist) {
			max_dist = d;
			opt._pts2fix = new int[3]{ p1, p2, p3 };
		}
		i++;
	}
	std::sort(opt._pts2fix, opt._pts2fix + 2);
}

void JacobianComposer::Scene2Problem(Scene &s, ceres::Problem* problem, std::string cam_model) {
    // Images
    for( auto &i : s._images ) {
        Image &img = i.second;
        Camera &c = s._cameras[img._cam_id];

		// Optional step for comparison -> sort the observations according point ids
		std::sort(img._point2D.begin(), img._point2D.end());

        // Observations ( don't contain the 2D points which don't correspond to any 3D point  )
		for (auto &p2D : img._point2D) {
			Point3D &p3D = s._points3D[p2D._point3D_id];

			// The loss function 
			ceres::LossFunction* loss_function = NULL;

			// Cost function & add residual   ( SIMPLE_PINHOLE, PINHOLE, SIMPLE_RADIAL, RADIAL )
			if (cam_model.compare(SIMPLE_PINHOLE) == 0) {
				ceres::CostFunction* cost_function = ceres::examples::SimplePinholeCameraReprojectionError::Create(p2D._xy[0] - c._uv[0], p2D._xy[1] - c._uv[1], c._f, c._r);
				problem->AddResidualBlock(cost_function, loss_function, img._aa, img._C, p3D._X);
			}
			else if (cam_model.compare(PINHOLE) == 0) {
				ceres::CostFunction* cost_function = ceres::examples::PinholeCameraReprojectionError::Create(p2D._xy[0] - c._uv[0], p2D._xy[1] - c._uv[1], c._r);
				problem->AddResidualBlock(cost_function, loss_function, img._aa, img._C, &(img._f), p3D._X);
			}
			else if (cam_model.compare(SIMPLE_RADIAL) == 0) {
				ceres::CostFunction* cost_function = ceres::examples::SimpleRadialCameraReprojectionError::Create(p2D._xy[0] - c._uv[0], p2D._xy[1] - c._uv[1]);
				problem->AddResidualBlock(cost_function, loss_function, img._aa, img._C, &(img._f), img._r, p3D._X);
			}
			else if (cam_model.compare(RADIAL) == 0) {
				ceres::CostFunction* cost_function = ceres::examples::RadialCameraReprojectionError::Create(p2D._xy[0] - c._uv[0], p2D._xy[1] - c._uv[1]);
				problem->AddResidualBlock(cost_function, loss_function, img._aa, img._C, &(img._f), img._r, p3D._X);
			}
		}
    }
}


void JacobianComposer::scene2Jacobian(std::string cam_model, std::string algorithm, Scene &s) {
    // Compose the problem
    ceres::Problem problem;
	Scene2Problem( s, &problem, cam_model);
    
    // Images & Cameras
    int n_imgs = 0, n_obs = 0;
    std::vector<double*> parameter_blocks;
    for(auto &i : s._images) {
        if ( i.second._point2D.size() == 0 ) continue;
        Image &img = i.second;
        Camera &c = s._cameras[img._cam_id];
        parameter_blocks.push_back(img._aa);
        parameter_blocks.push_back(img._C);
		
		// camera parameters 
		if (cam_model.compare(PINHOLE) == 0) {
			parameter_blocks.push_back(&(img._f));
		}
		else if (cam_model.compare(SIMPLE_RADIAL) == 0) {
			parameter_blocks.push_back(&(img._f));
			parameter_blocks.push_back(img._r);
		}
		else if (cam_model.compare(RADIAL) == 0) {
			parameter_blocks.push_back(&(img._f));
			parameter_blocks.push_back(img._r);
		}
        n_imgs++;
        n_obs += img._point2D.size();
    }
    
    // Points in 3D
    for(auto &p3D : s._points3D) 
        parameter_blocks.push_back( p3D.second._X );
    
    // Configure Jacobian engine 
    double cost = 0.0;
    ceres::Problem::EvaluateOptions evalOpt;
    evalOpt.parameter_blocks = parameter_blocks;
    evalOpt.apply_loss_function = true;

    // Build Jacobain 
	s._jacobian = ceres::CRSMatrix();
    problem.Evaluate(evalOpt, &cost, NULL, NULL, &(s._jacobian));

    // Configure the covariance engine ( find the indexes of the most distant points etc. )
	s._options = cov::Options();
    findPts2Fix(s._options, s._points3D.size(), s._points3D);
	s._options._numCams = n_imgs;

	if (cam_model.compare(SIMPLE_PINHOLE) == 0) {
		s._options._camParams = 6;
	}
	else if (cam_model.compare(PINHOLE) == 0) {
		s._options._camParams = 7;
	}
	else if (cam_model.compare(SIMPLE_RADIAL) == 0) {
		s._options._camParams = 8;
	}
	else if (cam_model.compare(RADIAL) == 0) {
		s._options._camParams = 9;
	}
	s._options._numPoints = s._points3D.size();
	s._options._numObs = n_obs;

    s._options._camParams = cov::EAlgorithm_stringToEnum(algorithm);

	s._options._epsilon = 1e-10;
	s._options._lambda = -1;
	s._options._svdRemoveN = 7;
	s._options._maxIterTE = -1;
}   
