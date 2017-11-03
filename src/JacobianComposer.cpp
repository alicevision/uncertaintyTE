/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   JacobianComposer.cpp
 * Author: policmic
 * 
 * Created on November 2, 2017, 4:41 PM
 */

#include "JacobianComposer.h"
#include "snavely_reprojection_error.h"

#ifndef DBL_MIN
    #define DBL_MIN -1.7e+308
#endif

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

void JacobianComposer::findPts2Fix(cov::Options &opt, int N, map<int, Point3D> &points3D){
    // Convert points to vector
    std::vector<Point3D> pts3D = std::vector<Point3D>();
    pts3D.reserve(N);
    for(auto const &p3D : points3D)
        pts3D.push_back(p3D.second);
    
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

void JacobianComposer::Scene2Problem(Scene &s, ceres::Problem* problem) {
    // Images
    for( auto &i : s._images ) {
        Image &img = i.second;
        Camera &c = s._cameras[img._cam_id];

        // Observations ( don't contain the 2D points which don't correspond to any 3D point  )
        for ( auto &p2D : img._point2D ){
            Point3D &p3D = s._points3D[p2D._point3D_id];
            
            // Cost function
            ceres::CostFunction* cost_function = ceres::examples::PinholeCameraReprojectionError::Create( p2D._xy[0], p2D._xy[1], c._f, c._r );
                    
            // The loss function
            ceres::LossFunction* loss_function = NULL;
            
            // Residual block
            problem->AddResidualBlock( cost_function, loss_function, img._aa, img._C, p3D._X );
        }
    }
}


void JacobianComposer::scene2Jacobian(std::string cam_model, Scene &s, ceres::CRSMatrix &jacobian, cov::Options &opt){
    // Compose the problem
    ceres::Problem problem;
    Scene2Problem( s, &problem );
    
    // Images & Cameras
    int n_imgs = 0, n_obs = 0;
    std::vector<double*> parameter_blocks;
    for(auto &i : s._images) {
        if ( i.second._point2D.size() == 0 ) continue;
        Image &img = i.second;
        Camera &c = s._cameras[img._cam_id];
        parameter_blocks.push_back(img._aa);
        parameter_blocks.push_back(img._C);
        n_imgs++;
        n_obs += img._point2D.size();
    }
    //parameter_blocks.push_back(&c._f);
    //parameter_blocks.push_back(c._r);
    
    // Points in 3D
    for(auto &p3D : s._points3D) 
        parameter_blocks.push_back( p3D.second._X );
    
    // Configure Jacobian engine 
    double cost = 0.0;
    ceres::Problem::EvaluateOptions evalOpt;
    evalOpt.parameter_blocks = parameter_blocks;
    evalOpt.num_threads = 8;
    evalOpt.apply_loss_function = true;

    // Build Jacobain 
    problem.Evaluate(evalOpt, &cost, NULL, NULL, &jacobian);

    // Configure the covariance engine ( find the indexes of the most distant points etc. )
    findPts2Fix(opt, s._points3D.size(), s._points3D);
    opt._numCams = n_imgs;
    opt._camParams = 6;         // TODO: different cameras
    opt._numPoints = s._points3D.size();
    opt._numObs = n_obs;
    opt._algorithm = 2;
    opt._epsilon = 1e-10;
    opt._lambda = -1;
    opt._svdRemoveN = -1;
    opt._maxIterTE = -1;    
}   
