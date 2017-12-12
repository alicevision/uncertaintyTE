// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: sameeragarwal@google.com (Sameer Agarwal)
//
// Templated struct implementing the camera model and residual
// computation for bundle adjustment used by Noah Snavely's Bundler
// SfM system. This is also the camera model/residual for the bundle
// adjustment problems in the BAL dataset. It is templated so that we
// can use Ceres's automatic differentiation to compute analytic
// jacobians.
//
// For details see: http://phototour.cs.washington.edu/bundler/
// and http://grail.cs.washington.edu/projects/bal/
#pragma once

#include "ceres/rotation.h"

namespace ceres {

template<typename T> inline
void AngleAxisRotatePoint(const T aa[3], const T pt[3], T result[3], T angle_scale[3]) {
T angle_axis[3] = { aa[0] * angle_scale[0], aa[1] * angle_scale[1], aa[2] * angle_scale[2] };
const T theta2 = DotProduct(angle_axis, angle_axis);
if (theta2 > T(std::numeric_limits<double>::epsilon())) {
	// Away from zero, use the rodriguez formula
	//
	//   result = pt costheta +
	//            (w x pt) * sintheta +
	//            w (w . pt) (1 - costheta)
	//
	// We want to be careful to only evaluate the square root if the
	// norm of the angle_axis vector is greater than zero. Otherwise
	// we get a division by zero.
	//
	const T theta = sqrt(theta2);
	const T costheta = cos(theta);
	const T sintheta = sin(theta);
	const T theta_inverse = 1.0 / theta;

	const T w[3] = { angle_axis[0] * theta_inverse,
		angle_axis[1] * theta_inverse,
		angle_axis[2] * theta_inverse };

	// Explicitly inlined evaluation of the cross product for
	// performance reasons.
	const T w_cross_pt[3] = { w[1] * pt[2] - w[2] * pt[1],
		w[2] * pt[0] - w[0] * pt[2],
		w[0] * pt[1] - w[1] * pt[0] };
	const T tmp = (w[0] * pt[0] + w[1] * pt[1] + w[2] * pt[2]) * (T(1.0) - costheta);

	result[0] = pt[0] * costheta + w_cross_pt[0] * sintheta + w[0] * tmp;
	result[1] = pt[1] * costheta + w_cross_pt[1] * sintheta + w[1] * tmp;
	result[2] = pt[2] * costheta + w_cross_pt[2] * sintheta + w[2] * tmp;
}
else {
	// Near zero, the first order Taylor approximation of the rotation
	// matrix R corresponding to a vector w and angle w is
	//
	//   R = I + hat(w) * sin(theta)
	//
	// But sintheta ~ theta and theta * w = angle_axis, which gives us
	//
	//  R = I + hat(w)
	//
	// and actually performing multiplication with the point pt, gives us
	// R * pt = pt + w x pt.
	//
	// Switching to the Taylor expansion near zero provides meaningful
	// derivatives when evaluated using Jets.
	//
	// Explicitly inlined evaluation of the cross product for
	// performance reasons.
	const T w_cross_pt[3] = { angle_axis[1] * pt[2] - angle_axis[2] * pt[1],
		angle_axis[2] * pt[0] - angle_axis[0] * pt[2],
		angle_axis[0] * pt[1] - angle_axis[1] * pt[0] };

	result[0] = pt[0] + w_cross_pt[0];
	result[1] = pt[1] + w_cross_pt[1];
	result[2] = pt[2] + w_cross_pt[2];
}
}

template<typename T> inline
void AngleAxisRotatePointRC(const T aa[3], const T pt[3], T result[3]) {
    const T theta2 = DotProduct(aa, aa);
    if (theta2 > T(std::numeric_limits<double>::epsilon())) {
            const T theta = sqrt(theta2);
            const T costheta = cos(theta);
            const T sintheta = sin(theta);
            const T theta_inverse = 1.0 / theta;

            const T w[3] = {aa[0] * theta_inverse,
                            aa[1] * theta_inverse,
                            aa[2] * theta_inverse };

            // Cross product
            const T w_cross_pt[3] = {   w[1] * pt[2] - w[2] * pt[1],
                                        w[2] * pt[0] - w[0] * pt[2],
                                        w[0] * pt[1] - w[1] * pt[0] };
            const T tmp = (w[0] * pt[0] + w[1] * pt[1] + w[2] * pt[2]) * (T(1.0) - costheta);

            result[0] = pt[0] * costheta + w_cross_pt[0] * sintheta + w[0] * tmp;
            result[1] = pt[1] * costheta + w_cross_pt[1] * sintheta + w[1] * tmp;
            result[2] = pt[2] * costheta + w_cross_pt[2] * sintheta + w[2] * tmp;
    }
    else {
            const T w_cross_pt[3] = {   aa[1] * pt[2] - aa[2] * pt[1],
                                        aa[2] * pt[0] - aa[0] * pt[2],
                                        aa[0] * pt[1] - aa[1] * pt[0] };
            result[0] = pt[0] + w_cross_pt[0];
            result[1] = pt[1] + w_cross_pt[1];
            result[2] = pt[2] + w_cross_pt[2];
    }
}



namespace examples {

    
struct SimplePinholeCameraReprojectionError {
	SimplePinholeCameraReprojectionError(const double observed_x, const double observed_y, const double f, const double *r) :
        _observed_x(observed_x), _observed_y(observed_y), _f(f), _r1(r[0]), _r2(r[1]) {}

    template <typename T>
    bool operator()(const T* const aa, const T* const C, const T* const point, T* residuals) const {
        T p[3];
        T pC[3];
        
        pC[0] = point[0] - C[0];
        pC[1] = point[1] - C[1];
        pC[2] = point[2] - C[2];
        
        ceres::AngleAxisRotatePointRC(aa, pC, p);
        T x = p[0] / p[2];
        T y = p[1] / p[2];
        
        T dist = x*x + y*y;
        T distortion = T(1.0) + T(_r1)*dist + T(_r2)*dist*dist;
        T predicted_x = T(_f) * x * distortion;
        T predicted_y = T(_f) * y * distortion;

        residuals[0] = predicted_x - T(_observed_x); 
        residuals[1] = predicted_y - T(_observed_y);
        return true;
    }

  // Factory to hide the construction of the CostFunction object from the client code.
  static ceres::CostFunction* Create(const double observed_x, const double observed_y, const double f, const double *r) {
    return (new ceres::AutoDiffCostFunction<SimplePinholeCameraReprojectionError, 2, 3, 3, 3>(new SimplePinholeCameraReprojectionError(observed_x, observed_y, f, r)));
  }
  double _observed_x;
  double _observed_y;
  double _f;
  double _r1;
  double _r2;
};


struct PinholeCameraReprojectionError {
	PinholeCameraReprojectionError(const double observed_x, const double observed_y, const double *r) :
		_observed_x(observed_x), _observed_y(observed_y), _r1(r[0]), _r2(r[1]) {}

	template <typename T>
	bool operator()(const T* const aa, const T* const C, 
		const T* const f, const T* const point, T* residuals) const {
		T p[3];
		T pC[3];

		pC[0] = point[0] - C[0];
		pC[1] = point[1] - C[1];
		pC[2] = point[2] - C[2];

		ceres::AngleAxisRotatePointRC(aa, pC, p);
		T x = p[0] / p[2];
		T y = p[1] / p[2];

		T dist = x*x + y*y;
		T distortion = T(1.0) + T(_r1)*dist + T(_r2)*dist*dist;
		T predicted_x = f[0] * x * distortion;
		T predicted_y = f[0] * y * distortion;

		residuals[0] = predicted_x - T(_observed_x);
		residuals[1] = predicted_y - T(_observed_y);
		return true;
	}

	// Factory to hide the construction of the CostFunction object from the client code.
	static ceres::CostFunction* Create(const double observed_x, const double observed_y, const double *r) {
		return (new ceres::AutoDiffCostFunction<PinholeCameraReprojectionError, 2, 3, 3, 1, 3>(new PinholeCameraReprojectionError(observed_x, observed_y, r)));
	}
	double _observed_x;
	double _observed_y;
	double _r1;
	double _r2;
};

    
struct SimpleRadialCameraReprojectionError {
	double _observed_x;
	double _observed_y;
	
	SimpleRadialCameraReprojectionError(const double observed_x, const double observed_y) : _observed_x(observed_x), _observed_y(observed_y) {}

	// The generic residula function
	template <typename T>
	bool operator()(const T* const aa,  const T* const C, 
					const T* const f,   const T* const r, 
					const T* const point, T* residuals) const {
		T p[3];
		T pC[3];

		pC[0] = point[0] - C[0];
		pC[1] = point[1] - C[1];
		pC[2] = point[2] - C[2];

		ceres::AngleAxisRotatePointRC(aa, pC, p);
		T x = p[0] / p[2];
		T y = p[1] / p[2];

		T dist = x*x + y*y;
		T distortion = T(1.0) + r[0]*dist;
		T predicted_x = f[0] * x * distortion;
		T predicted_y = f[0] * y * distortion;

		residuals[0] = predicted_x - T(_observed_x);
		residuals[1] = predicted_y - T(_observed_y);
		return true;
	}

	// Factory to hide the construction of the CostFunction object from the client code.
	static ceres::CostFunction* Create(const double observed_x, const double observed_y) {
		return (new ceres::AutoDiffCostFunction<SimpleRadialCameraReprojectionError, 2, 3, 3, 1, 1, 3>(
			new SimpleRadialCameraReprojectionError(observed_x, observed_y)));
	}
};
    

struct RadialCameraReprojectionError {
	double _observed_x;
	double _observed_y;

	RadialCameraReprojectionError(const double observed_x, const double observed_y) : _observed_x(observed_x), _observed_y(observed_y) {}

	// The generic residula function
	template <typename T>
	bool operator()(const T* const aa, const T* const C,
		const T* const f, const T* const r,
		const T* const point, T* residuals) const {
		T p[3];
		T pC[3];

		pC[0] = point[0] - C[0];
		pC[1] = point[1] - C[1];
		pC[2] = point[2] - C[2];

		ceres::AngleAxisRotatePointRC(aa, pC, p);
		T x = p[0] / p[2];
		T y = p[1] / p[2];

		T dist = x*x + y*y;
		T distortion = T(1.0) + r[0] * dist + r[1] * dist*dist;
		T predicted_x = f[0] * x * distortion;
		T predicted_y = f[0] * y * distortion;

		residuals[0] = predicted_x - T(_observed_x);
		residuals[1] = predicted_y - T(_observed_y);
		return true;
	}

	// Factory to hide the construction of the CostFunction object from the client code.
	static ceres::CostFunction* Create(const double observed_x, const double observed_y) {
		return (new ceres::AutoDiffCostFunction<RadialCameraReprojectionError, 2, 3, 3, 1, 2, 3>(
			new RadialCameraReprojectionError(observed_x, observed_y)));
	}
};

    
// Templated pinhole camera model for used with Ceres.  The camera is
// parameterized using 9 parameters: 3 for rotation, 3 for translation, 1 for
// focal length and 2 for radial distortion. The principal point is not modeled
// (i.e. it is assumed be located at the image center).  
struct SnavelyReprojectionError {
	SnavelyReprojectionError(double observed_x, double observed_y, double camera_w0, double camera_w1, double camera_w2, double camera_w3,
		double camera_w4, double camera_w5, double camera_w6, double camera_w7, double camera_w8, double point_w0, double point_w1, double point_w2) :
		observed_x(observed_x), observed_y(observed_y), camera_w0(camera_w0), camera_w1(camera_w1), camera_w2(camera_w2), camera_w3(camera_w3),
		camera_w4(camera_w4), camera_w5(camera_w5), camera_w6(camera_w6), camera_w7(camera_w7), camera_w8(camera_w8), point_w0(point_w0), point_w1(point_w1), point_w2(point_w2) {}

  template <typename T>
  bool operator()(const T* const camera, const T* const point, T* residuals) const {
	T p[3];
	T pC[3];
	T angle_axis_w[3] = {T(camera_w0),T(camera_w1),T(camera_w2)};

	// camera[3,4,5] are the camera center C
	pC[0] = (point[0] * T(point_w0)) - (camera[3] * T(camera_w3));
	pC[1] = (point[1] * T(point_w1)) - (camera[4] * T(camera_w4));
	pC[2] = (point[2] * T(point_w2)) - (camera[5] * T(camera_w5));

	// camera[0,1,2] are the angle-axis rotation.
    ceres::AngleAxisRotatePoint(camera, pC, p, angle_axis_w);

    const T& focal = (camera[6] * T(camera_w6));
    T xp = p[0] / p[2];
    T yp = p[1] / p[2];

    // Apply second and fourth order radial distortion.
    const T& l1 = (camera[7] * T(camera_w7));
    const T& l2 = (camera[8] * T(camera_w8));
	T r2 = xp*xp + yp*yp;
    T distortion = T(1.0) + l1*r2 + l2*r2*r2;

    // Compute final projected point position.
    T predicted_x = focal * xp * distortion;
    T predicted_y = focal * yp * distortion;

    // The error is the difference between the predicted and observed position.
	residuals[0] = predicted_x - T(observed_x); 
    residuals[1] = predicted_y - T(observed_y);
    return true;
  }

  // Factory to hide the construction of the CostFunction object from the client code.
  static ceres::CostFunction* Create(const double observed_x, const double observed_y, const double camera_w0, const double camera_w1, const double camera_w2, 
	  const double camera_w3, const double camera_w4, const double camera_w5, const double camera_w6, const double camera_w7, const double camera_w8,
	  const double point_w0, const double point_w1, const double point_w2) {
    return (new ceres::AutoDiffCostFunction<SnavelyReprojectionError, 2, 9, 3>(new SnavelyReprojectionError(observed_x, observed_y, 
		camera_w0, camera_w1, camera_w2, camera_w3, camera_w4, camera_w5, camera_w6, camera_w7, camera_w8, point_w0, point_w1, point_w2)));
  }
  double camera_w0;
  double camera_w1;
  double camera_w2;
  double camera_w3;
  double camera_w4;
  double camera_w5;
  double camera_w6;
  double camera_w7;
  double camera_w8;
  double point_w0;
  double point_w1;
  double point_w2;
  double observed_x;
  double observed_y;
};


// Fix also focal and radial 
struct FixOneCamCenterCoordinateFocalRadialReprojectionError {
	FixOneCamCenterCoordinateFocalRadialReprojectionError(double observed_x, double observed_y, double* cam) :
		observed_x(observed_x), observed_y(observed_y), cam_(cam) {}

	template <typename T>
	bool operator()(const T* const camera, const T* const point, T* residuals) const {
		T p[3];
		T pC[3];

		// camera[3,4,5] are the camera center C
		pC[0] = point[0] - T(cam_[3]);
		pC[1] = point[1] - camera[4];
		pC[2] = point[2] - camera[5];

		// camera[0,1,2] are the angle-axis rotation.
		ceres::AngleAxisRotatePoint(camera, pC, p);
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Apply second and fourth order radial distortion.
		T r2 = xp*xp + yp*yp;
		T distortion = T(1.0) + T(cam_[7]) * r2 + T(cam_[8]) * r2*r2;

		// Compute final projected point position.
		T predicted_x = T(cam_[6]) * xp * distortion;
		T predicted_y = T(cam_[6]) * yp * distortion;

		// The error is the difference between the predicted and observed position.
		residuals[0] = (predicted_x - T(observed_x)); // / sigma_reprojection;
		residuals[1] = predicted_y - T(observed_y);
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create(const double observed_x, const double observed_y, double* cam) {
		return (new ceres::AutoDiffCostFunction<FixOneCamCenterCoordinateFocalRadialReprojectionError, 2, 9, 3>(
			new FixOneCamCenterCoordinateFocalRadialReprojectionError(observed_x, observed_y, cam)));
	}

	double observed_x;
	double observed_y;
	double* cam_;
};



struct FixFocalReprojectionError {
	FixFocalReprojectionError(double observed_x, double observed_y, double* cam) :
		observed_x(observed_x), observed_y(observed_y), cam_(cam) {}

	template <typename T>
	bool operator()(const T* const camera, const T* const point, T* residuals) const {
		T p[3];
		T pC[3];

		// camera[3,4,5] are the camera center C
		pC[0] = point[0] - camera[3];
		pC[1] = point[1] - camera[4];
		pC[2] = point[2] - camera[5];

		// camera[0,1,2] are the angle-axis rotation.
		ceres::AngleAxisRotatePoint(camera, pC, p);
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Apply second and fourth order radial distortion.
		T r2 = xp*xp + yp*yp;
		T distortion = T(1.0) + camera[7] * r2 + camera[8] * r2*r2;

		// Compute final projected point position.
		T predicted_x = T(cam_[6]) * xp * distortion;
		T predicted_y = T(cam_[6]) * yp * distortion;

		// The error is the difference between the predicted and observed position.
		residuals[0] = (predicted_x - T(observed_x)); // / sigma_reprojection;
		residuals[1] = predicted_y - T(observed_y);
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create(const double observed_x, const double observed_y, double* cam) {
		return (new ceres::AutoDiffCostFunction<FixFocalReprojectionError, 2, 9, 3>(
			new FixFocalReprojectionError(observed_x, observed_y, cam)));
	}

	double observed_x;
	double observed_y;
	double* cam_;
};



// Templated pinhole camera model for used with Ceres.  The camera is
// parameterized using 6 parameters: 3 for rotation, 3 for translation and 
// fixed constant 1 for focal length and 2 for radial distortion. 
// The principal point is not modeled (i.e. it is assumed be located at the image center).  
struct FixFocalRadialReprojectionError {
	FixFocalRadialReprojectionError(double observed_x, double observed_y, double focal, 
		double radial1, double radial2) : observed_x(observed_x), observed_y(observed_y),
		focal_(focal), radial1_(radial1), radial2_(radial2){}

	template <typename T>
	bool operator()(const T* const camera, const T* const point, T* residuals) const {
		T p[3];
		T pC[3];

		// camera[3,4,5] are the camera center C
		pC[0] = point[0] - camera[3];
		pC[1] = point[1] - camera[4];
		pC[2] = point[2] - camera[5];

		// camera[0,1,2] are the angle-axis rotation.
		ceres::AngleAxisRotatePoint(camera, pC, p);
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Apply second and fourth order radial distortion.
		T r2 = xp*xp + yp*yp;
		T distortion = T(1.0) + T(radial1_)*r2 + T(radial2_)*r2*r2;

		// Compute final projected point position.
		T predicted_x = T(focal_) * xp * distortion;
		T predicted_y = T(focal_) * yp * distortion;

		// The error is the difference between the predicted and observed position.
		residuals[0] = (predicted_x - T(observed_x)); // / sigma_reprojection;
		residuals[1] = predicted_y - T(observed_y);
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create(const double observed_x,
		const double observed_y, double focal, double radial1, double radial2) {
		return (new ceres::AutoDiffCostFunction<FixFocalRadialReprojectionError, 2, 9, 3>(
			new FixFocalRadialReprojectionError(observed_x, observed_y, focal, radial1, radial2)));
	}

	double observed_x;
	double observed_y;
	double focal_;
	double radial1_;
	double radial2_;
};

struct FixCamsReprojectionError {
	FixCamsReprojectionError(double observed_x, double observed_y, double* cam) : 
		observed_x(observed_x), observed_y(observed_y), cam_(cam) {}

	template <typename T>
	bool operator()(const T* const point, T* residuals) const {
		T p[3];
		T pC[3];
		T camera[9];
		camera[0] = T(cam_[0]);
		camera[1] = T(cam_[1]);
		camera[2] = T(cam_[2]);
		camera[3] = T(cam_[3]);
		camera[4] = T(cam_[4]);
		camera[5] = T(cam_[5]);
		camera[6] = T(cam_[6]);
		camera[7] = T(cam_[7]);
		camera[8] = T(cam_[8]);

		// camera[3,4,5] are the camera center C
		pC[0] = point[0] - camera[3];
		pC[1] = point[1] - camera[4];
		pC[2] = point[2] - camera[5];

		// camera[0,1,2] are the angle-axis rotation.
		ceres::AngleAxisRotatePoint(camera, pC, p);
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Apply second and fourth order radial distortion.
		T r2 = xp*xp + yp*yp;
		T distortion = T(1.0) + camera[7]*r2 + camera[8]*r2*r2;

		// Compute final projected point position.
		T predicted_x = camera[6] * xp * distortion;
		T predicted_y = camera[6] * yp * distortion;

		// The error is the difference between the predicted and observed position.
		residuals[0] = (predicted_x - T(observed_x)); // / sigma_reprojection;
		residuals[1] = predicted_y - T(observed_y);
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create(const double observed_x, const double observed_y, double* cam) {
		return (new ceres::AutoDiffCostFunction<FixCamsReprojectionError, 2, 3>(
			new FixCamsReprojectionError(observed_x, observed_y, cam)));
	}

	double observed_x;
	double observed_y;
	double* cam_;
};


struct FixPtsReprojectionError {
	FixPtsReprojectionError(double observed_x, double observed_y, double* pt) :
		observed_x(observed_x), observed_y(observed_y), pt_(pt) {}

	template <typename T>
	bool operator()(const T* const camera, T* residuals) const {
		T p[3];
		T pC[3];

		// camera[3,4,5] are the camera center C
		pC[0] = T(pt_[0]) - camera[3];
		pC[1] = T(pt_[1]) - camera[4];
		pC[2] = T(pt_[2]) - camera[5];

		// camera[0,1,2] are the angle-axis rotation.
		ceres::AngleAxisRotatePoint(camera, pC, p);

		const T& focal = camera[6];
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Apply second and fourth order radial distortion.
		const T& l1 = camera[7];
		const T& l2 = camera[8];
		T r2 = xp*xp + yp*yp;
		T distortion = T(1.0) + l1*r2 + l2*r2*r2;

		// Compute final projected point position.
		T predicted_x = focal * xp * distortion;
		T predicted_y = focal * yp * distortion;

		// The error is the difference between the predicted and observed position.
		residuals[0] = (predicted_x - T(observed_x)); // / sigma_reprojection;
		residuals[1] = predicted_y - T(observed_y);
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create(const double observed_x,
		const double observed_y, double* pt) {
		return (new ceres::AutoDiffCostFunction<FixPtsReprojectionError, 2, 9>(
			new FixPtsReprojectionError(observed_x, observed_y, pt)));
	}

	double observed_x;
	double observed_y;
	double* pt_;
};

struct FixPtsFocalRadialReprojectionError {
	FixPtsFocalRadialReprojectionError(double observed_x, double observed_y, double* cam, double* pt) : 
		observed_x(observed_x), observed_y(observed_y), cam_(cam), pt_(pt) {}

	template <typename T>
	bool operator()(const T* const camera, T* residuals) const {
		T p[3];
		T pC[3];

		// camera[3,4,5] are the camera center C
		pC[0] = T(pt_[0]) - camera[3];
		pC[1] = T(pt_[1]) - camera[4];
		pC[2] = T(pt_[2]) - camera[5];

		// camera[0,1,2] are the angle-axis rotation.
		ceres::AngleAxisRotatePoint(camera, pC, p);

		const T& focal = T(cam_[6]);
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Apply second and fourth order radial distortion.
		const T& l1 = T(cam_[7]);
		const T& l2 = T(cam_[8]);
		T r2 = xp*xp + yp*yp;
		T distortion = T(1.0) + l1*r2 + l2*r2*r2;

		// Compute final projected point position.
		T predicted_x = focal * xp * distortion;
		T predicted_y = focal * yp * distortion;

		// The error is the difference between the predicted and observed position.
		residuals[0] = (predicted_x - T(observed_x)); // / sigma_reprojection;
		residuals[1] = predicted_y - T(observed_y);
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create(const double observed_x,
		const double observed_y, double* cam, double* pt) {
		return (new ceres::AutoDiffCostFunction<FixPtsFocalRadialReprojectionError, 2, 9>(
			new FixPtsFocalRadialReprojectionError(observed_x, observed_y, cam, pt)));
	}

	double observed_x;
	double observed_y;
	double* pt_;
	double* cam_;
};

// Templated pinhole camera model for used with Ceres.  The camera is
// parameterized using 10 parameters. 4 for rotation, 3 for
// translation, 1 for focal length and 2 for radial distortion. The
// principal point is not modeled (i.e. it is assumed be located at
// the image center).
struct SnavelyReprojectionErrorWithQuaternions {
  // (u, v): the position of the observation with respect to the image
  // center point.
  SnavelyReprojectionErrorWithQuaternions(double observed_x, double observed_y)
      : observed_x(observed_x), observed_y(observed_y) {}

  template <typename T>
  bool operator()(const T* const camera_rotation,
                  const T* const camera_translation_and_intrinsics,
                  const T* const point,
                  T* residuals) const {
    const T& focal = camera_translation_and_intrinsics[3];
    const T& l1 = camera_translation_and_intrinsics[4];
    const T& l2 = camera_translation_and_intrinsics[5];

    // Use a quaternion rotation that doesn't assume the quaternion is
    // normalized, since one of the ways to run the bundler is to let Ceres
    // optimize all 4 quaternion parameters unconstrained.
    T p[3];
    QuaternionRotatePoint(camera_rotation, point, p);

    p[0] += camera_translation_and_intrinsics[0];
    p[1] += camera_translation_and_intrinsics[1];
    p[2] += camera_translation_and_intrinsics[2];

    // Compute the center of distortion. The sign change comes from
    // the camera model that Noah Snavely's Bundler assumes, whereby
    // the camera coordinate system has a negative z axis.
    T xp = p[0] / p[2];
    T yp = p[1] / p[2];

    // Apply second and fourth order radial distortion.
    T r2 = xp*xp + yp*yp;
    T distortion = T(1.0) + r2  * (l1 + l2  * r2);

    // Compute final projected point position.
    T predicted_x = focal * distortion * xp;
    T predicted_y = focal * distortion * yp;

    // The error is the difference between the predicted and observed position.
    residuals[0] = predicted_x - T(observed_x);
    residuals[1] = predicted_y - T(observed_y);
    return true;
  }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create(const double observed_x,
                                     const double observed_y) {
    return (new ceres::AutoDiffCostFunction<
            SnavelyReprojectionErrorWithQuaternions, 2, 4, 6, 3>(
                new SnavelyReprojectionErrorWithQuaternions(observed_x,
                                                            observed_y)));
  }

  double observed_x;
  double observed_y;
};

}  // namespace examples
}  // namespace ceres

