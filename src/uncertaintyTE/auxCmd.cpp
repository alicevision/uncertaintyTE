// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/*
* File:   auxCmd.cpp
* Author: Michal Polic
*/

#include "auxCmd.h"
#include "ceres/rotation.h"

#ifdef USE_MATLAB
  #include <mex.h>
#endif

#include <random>
#include <limits>


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
	double max_dist = std::numeric_limits<double>::min();
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


