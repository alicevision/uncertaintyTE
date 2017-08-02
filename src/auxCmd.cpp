#ifndef AUXCMD_HEADER_INCLUDED
	#include "auxCmd.h"
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

void loadJacobian(std::ifstream& file, int algorithm, ceres::CRSMatrix& jacobian, cov::Options& options) {
	if (!file.good()) {
		std::cerr << "\nThe input file doesn't exist.\n";
		exit(1);
	}

	int numJ;
	int fixPts[3];
	options._algorithm = algorithm;
	file >> options._lambda >> options._numCams >> options._camParams >> options._numPoints >> options._numObs;
	file >> fixPts[0] >> fixPts[1] >> fixPts[2];
	if (fixPts[0] >= 0 && fixPts[1] >= 0 && fixPts[2] >= 0)
		options._pts2fix = new int[3]{ fixPts[0], fixPts[1], fixPts[2] };

	file >> jacobian.num_rows >> jacobian.num_cols >> numJ;
	std::vector<int> rows(jacobian.num_rows + 1);
	std::vector<int> cols(numJ);
	std::vector<double> values(numJ);

	for (int i = 0; i <= jacobian.num_rows; ++i)
		file >> rows[i];
	for (int i = 0; i < numJ; ++i)
		file >> cols[i];
	for (int i = 0; i < numJ; ++i)
		file >> values[i];

	jacobian.rows = rows;
	jacobian.cols = cols;
	jacobian.values = values;
}

void saveResults(std::string& process_file_name, const std::string& current_dir, cov::Options& options, cov::Statistic& statistic,
	int num_camera_covar_values, double* camUnc) {
	std::cout << "\nPrinting the results to file... ";
	replaceAll(process_file_name, "in", "out");
	replaceAll(process_file_name, ".jacob", ".cov");
	std::ofstream outfile(current_dir + process_file_name);
	outfile << "# ---- Covariance v0.1 ----\n";
	outfile << "# Number of cameras: " << options._numCams << "\n";
	outfile << "# Number of camera parameters: " << options._camParams << "\n";
	outfile << "# Number of points in 3D: " << options._numPoints << "\n";
	outfile << "# Number of observations: " << options._numObs << "\n";
	outfile << "# Used algorithm: " << algorihm2str(options._algorithm) << "\n";
	if (options._algorithm == TAYLOR_EXPANSION) {
		if (statistic.fixedPts != NULL)
			outfile << "# Fixed points: " << statistic.fixedPts[0] << ", " << statistic.fixedPts[1] << ", " << statistic.fixedPts[2] << "\n";
		outfile << "# Used lambda: " << statistic.lambda << "\n";
		outfile << "# Loading jacobian time: " << statistic.timeCreateJ << "s\n";
		outfile << "# Normalization of jacobian time: " << statistic.timeNormJ << "s\n";
		outfile << "# Compose information matrix time: " << statistic.timeMultiplyJJ << "s\n";
		outfile << "# Split infromation matrix time: " << statistic.timeSplitJJ << "s\n";
		outfile << "# Inverse V time: " << statistic.timeInvV << "s\n";
		outfile << "# Compose Z time: " << statistic.timeComposeZ << "s\n";
		outfile << "# Inverse Z time: " << statistic.timeInvZ << "s\n";
		outfile << "# Taylor expansion time: " << statistic.timeTE << "s\n";
		outfile << "# TE number of iterations: " << statistic.cycle_change.size() << "s\n";
		outfile << "# TE cycle change: ";
		for (int i = 0; i < statistic.cycle_change.size(); ++i)
			outfile << statistic.cycle_change.at(i) << " ";
		outfile << "\n";
	}
	outfile << "# Time of the algorithm: " << statistic.timeAll << "s\n";
	outfile << "# -----------------------\n";
	for (int i = 0; i < options._numCams; ++i) {
		for (int j = 0; j < num_camera_covar_values; ++j)
			outfile << camUnc[i * num_camera_covar_values + j] << " ";
		outfile << "\n";
	}
	for (int i = 0; i < options._numPoints; ++i) {
		for (int j = 0; j < 6; ++j)
			outfile << camUnc[i * 6 + j] << " ";
		outfile << "\n";
	};
	outfile.close();
	std::cout << "[done]\n";
}