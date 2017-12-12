
#include "IO.h"
#include "auxCmd.h"


bool saveResults(const std::string& filepath, cov::Options& options, cov::Statistic& statistic,
    int num_camera_covar_values, const double* camUnc, const double* ptsUnc)
{
	std::cout << "\nPrinting the results to file... ";
    std::ofstream outfile(filepath);
	outfile << "# ---- Covariance v0.1 ----\n";
	outfile << "# Number of cameras: " << options._numCams << "\n";
	outfile << "# Number of camera parameters: " << options._camParams << "\n";
	outfile << "# Number of points in 3D: " << options._numPoints << "\n";
	outfile << "# Number of observations: " << options._numObs << "\n";
    outfile << "# Used algorithm: " << EAlgorithm_enumToString(options._algorithm) << "\n";
    if (options._algorithm == cov::eAlgorithmSvdTaylorExpansion) {
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
		outfile << "# Point uncertainty time: " << statistic.timePtsUnc << "s\n";
		outfile << "# TE number of iterations: " << statistic.cycle_change.size() << "\n";
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
			outfile << ptsUnc[i * 6 + j] << " ";
		outfile << "\n";
	};
	outfile.close();
	std::cout << "[done]\n";
	return true;
}


bool IO::writeCov2File(const std::string& filepath, Scene& scene, cov::Statistic& statistic) {
    return saveResults(filepath, scene._options, statistic, scene._uncertainty._nbCovarianceValuePerCam, &scene._uncertainty._camerasUnc[0], &scene._uncertainty._pointsUnc[0]);
}
