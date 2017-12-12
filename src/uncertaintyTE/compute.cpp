// This file is part of the AliceVision project.
// This Source Code Form is subject to the terms of the Mozilla Public License,
// v. 2.0. If a copy of the MPL was not distributed with this file,
// You can obtain one at https://mozilla.org/MPL/2.0/.

/*
* File:   compute.cpp
* Author: Michal Polic
*/

#include "compute.h"
#include "auxCmd.h"
#include "ScaledSparseMatrix.h" 
#include "ScaledDenseMatrix.h" 
#include <Eigen/QR>
#include <Eigen/SparseCore>

#ifdef USE_MATLAB
    #include <mex.h>
#endif

using namespace std;
using namespace Eigen;

typedef ScaledSparseMatrix SSM;
typedef ScaledDenseMatrix SDM;
typedef std::unique_ptr<SSM> uSM;
typedef std::unique_ptr<SDM> uDM;
typedef SparseMatrix<float, RowMajor> SM;
typedef MatrixXf DM;

double timeDuration(const tp& from, const tp& to) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(to - from).count() * 10e-3;
}

tp t(const tp& s, const std::string& txt) {
	cout << " " << timeDuration(s, Clock::now()) << "s\n";
	cout << txt;
	return Clock::now();
}

tp t(const tp& s, const std::string& txt, double *time) {
	*time = timeDuration(s, Clock::now());
	cout << " " << (*time) << "s\n";
	cout << txt;
	return Clock::now();
}

void symmetrizeMatrix(int N, SDM* A) {
#pragma omp parallel for
	for (int i = 0; i < N; ++i) {
		for (int j = i; j < N; ++j)
			A->set(j, i, A->sval(i, j));
	}
}

int factorial(int n) {
	if (n > 1)
		return n * factorial(n - 1);
	else
		return 1;
}

void iZupdate(SDM *iZ, double coeff, SDM *iZadd) {
	// *iZ += k * iZadd;
	for (int i = 0; i < iZ->ncols(); ++i) {
		for (int j = 0; j < iZ->nrows(); ++j)
			iZ->set(j, i, (iZ->val(j, i) + (coeff * iZadd->c()) * iZadd->sval(j, i)));
	}
	iZ->set_c(1);
}

void composeZ(tp *s, cov::Options &options, cov::Statistic &statistic, SSM &J, double **diagRightScaleJ, SSM *Y, SDM *Z) {
	if(options._debug) J.printBlock2Matlab3("Jbs", 0, 0, J.nrows(), J.ncols());
	
	// Scale Jacobian from left by a vector
	double csJ = 1;
	J.scaleMat(RIGHT, diagRightScaleJ, &csJ);
	if(options._debug) J.printBlock2Matlab3("Jscale", 0, 0, J.nrows(), J.ncols());

	for (int i = 0; i < J.ncols(); ++i)
		(*diagRightScaleJ)[i] = 1 / (csJ * (*diagRightScaleJ)[i]);
	*s = t(*s, "Computing sM ... ", &(statistic.timeNormJ));

	// Compute cJJ * sJJ = sJ'sJ
	SSM Jt(J.trn());		// create transpose matix and save reference to it
	if(options._debug) Jt.printBlock2Matlab3("Jt", 0, 0, Jt.nrows(), Jt.ncols());

	SSM M(Jt * J);
	*s = t(*s, "Split sM -> sU, sW, sV ... ", &(statistic.timeMultiplyJJ));
	if(options._debug) M.printBlock2Matlab3("M", 0, 0, M.nrows(), M.ncols());

	// Split cJJ -> U, W, V
	SSM *U = new SSM(), *W = new SSM(), *V = new SSM(), *iV = new SSM();
	int camBlockSize = options._numCams * options._camParams;
	M.splitTo3Blocks(camBlockSize, camBlockSize, U, W, V);
	*s = t(*s, "Computing siV ... ", &(statistic.timeSplitJJ));
	if(options._debug) U->printBlock2Matlab3("U", 0, 0, U->nrows(), U->ncols());
	if(options._debug) W->printBlock2Matlab3("W", 0, 0, W->nrows(), W->ncols());
	if(options._debug) V->printBlock2Matlab3("V", 0, 0, V->nrows(), V->ncols());

	// Compute inverse of V: sV -> isV
	V->inv3x3blockSymmDiag(iV);
	delete V;
	*s = t(*s, "Computing sZ ... ", &(statistic.timeInvV));

	// Compute Z = Schur complement of V;  U,iV,W  -> Z
	*Y = std::move((*W) * (*iV));
	SSM sWt = std::move(W->trn());
	SSM YsWt = std::move((*Y) * sWt);
	SSM sZ = std::move((*U) - YsWt);
	*Z = std::move(SDM(sZ));
	*s = t(*s, "Computing scaled inverse of sZ ", &(statistic.timeComposeZ));

	delete U;
	delete W;
	delete iV;
}

// Return inverse of SDM **dZ is saved into the matrix dZ
void teInverse(tp *s, int N, cov::Options &options, cov::Statistic &statistic, SDM *iZ) {
	// Z -> Z + lambda I
	double lambda = options._lambda;
	if (lambda == -1)
		lambda = pow(10, -1.2653*log10(N) - 2.9415);
	for (int i = 0; i < N; ++i)
		iZ->set(i, i, (iZ->sval(i, i) + (lambda / iZ->c())));
	statistic.lambda = lambda;
	cout << "using lambda: " << lambda << " ... ";
#ifdef USE_MATLAB
	mexPrintf("using lambda: %e\n\n", lambda);
#endif

	// Z -> iZ
	iZ->inv();
	*s = t(*s, "Taylor expansion ... ", &(statistic.timeInvZ));
#ifdef USE_MATLAB
	mexPrintf("Taylor expansion ... ");
#endif

	// TE
	double old_change = DBL_MAX, k, change;
	SDM *iZorig = new SDM(*iZ);
	SDM iZadd = std::move((*iZ) * (*iZ));
	for (int i = 1; i < 20; ++i) {
		k = pow(lambda, i) / factorial(i - 1);
		change = abs(k * iZadd.absMax());
#ifdef USE_MATLAB
		mexPrintf("\n cykle %d, coeff: %e, change: %e ", i, k, change);
#endif
		cout << "\n>>> cykle " << i << ", coeff: " << k << ", change: " << change;
		//cout << "\n>>>> iZadd.absMax()" << iZadd.absMax();
		if (change < 1e-5 || change > old_change)
			break;
		old_change = change;
		statistic.cycle_change.push_back(change);

		iZupdate(iZ, k, &iZadd);
		iZadd *= (*iZorig);
	}
#ifdef USE_MATLAB
	mexPrintf("\n Taylor expansion finished.\n");
#endif
	cout << "\nTaylor expansion have been done in ";
	*s = t(*s, "Refactor solution to the output ... ", &(statistic.timeTE));
}

void removeScaleJ4Z(tp *s, const double *diagRightScaleJ, cov::Options &opt, SDM *iZ, double *camUnc) {
	if (camUnc == NULL) return;
	int l = -1;
	for (int i = 0; i < opt._numCams; ++i) {
		int st = i * opt._camParams;
		for (int j = st; j < st + opt._camParams; ++j) {
			for (int k = j; k < st + opt._camParams; ++k) {
				if (diagRightScaleJ != NULL)
					camUnc[++l] = (iZ->val(j, k) + iZ->val(k, j) / 2) * diagRightScaleJ[j] * diagRightScaleJ[k];
				else
					camUnc[++l] = (iZ->val(j, k) + iZ->val(k, j) / 2);
			}
		}
	}
	cout << " " << timeDuration(*s, Clock::now()) << "s\n";
}

void svdInverse(magma_int_t *info, int N, cov::Options &options, SSM &J, SDM *iZ) {
	SDM Z(J.trn() * J);
	double *sv = (double*)malloc(N * sizeof(double));
	double *U = (double*)malloc(N*N * sizeof(double));
	double *Vt = (double*)malloc(N*N * sizeof(double));
	assert(sv != NULL);
	assert(U != NULL);
	assert(Vt != NULL);
	int lwork, *iwork;
	double *hwork;

	// Use Intel MKL Lapack SVD   ( aprox. 3-4x faster then GPU variant working with double )
	switch (options._algorithm) {
    case cov::eAlgorithmSvdQrIteration:				// additional memory requirements: N*N + 3*N + 2*N*32 DOUBLE
		lwork = N*N + 3 * N + 2 * N * 32;
		hwork = (double*)malloc(lwork * sizeof(double));
		assert(hwork != NULL);
		lapackf77_dgesvd(lapack_vec_const(MagmaAllVec), lapack_vec_const(MagmaAllVec), &N, &N,
			Z.getMatPtr(), &N, sv, U, &N, Vt, &N, hwork, &lwork, info);
		free(hwork);
		break;

    case cov::eAlgorithmSvdDeviceAndconquer:		// additional memory requirements: 4*N*N + 7*N DOUBLE  + 8*N INT
		lwork = 4 * N*N + 7 * N;
		hwork = (double*)malloc(lwork * sizeof(double));
		assert(hwork != NULL);
		iwork = (int*)malloc(8 * N * sizeof(int));
		assert(iwork != NULL);
		lapackf77_dgesdd(lapack_vec_const(MagmaAllVec), &N, &N,
			Z.getMatPtr(), &N, sv, U, &N, Vt, &N, hwork, &lwork, iwork, info);
		free(hwork);
		free(iwork);
		break;
	}
	TESTING_CHECK(*info);

	// Combine all matrices back to the pseudo-inverse sdZ -> sdiM
	if (options._svdRemoveN != -1) {
		// U = U * diag(1/sv);   for values i < _svdRemoveN
#pragma omp parallel for 
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j)
				U[i*N + j] *= (i < (N - options._svdRemoveN) ? 1 / sv[i] : 0);  // standard use: options._svdRemoveN = 7
		}
	}
	else {
		// U = U * diag(1/sv);   for values sv(j) > eps
		double eps = options._epsilon;
		if (eps < 0)
			eps = 1e-10;
#pragma omp parallel for 
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j)
				U[i*N + j] *= (sv[i] > eps ? 1 / sv[i] : 0);
		}
	}

	// iUVW = U * Vt
	double alpha = 1, beta = 0;
	blasf77_dgemm(lapack_trans_const(MagmaNoTrans), lapack_trans_const(MagmaNoTrans), &N, &N, &N,
		&alpha, U, &N, Vt, &N, &beta, iZ->getMatPtr(), &N);
	free(sv);
	free(U);
	free(Vt);
}

void findICP(int numObs, int camParams, int numCams, int *h_Jcols, int **camsIds, int **ptsIds) {
	(*camsIds) = (int*)malloc(numObs * sizeof(int));
	assert((*camsIds) != NULL);
	(*ptsIds) = (int*)malloc(numObs * sizeof(int));
	assert((*ptsIds) != NULL);
	int step = 2 * (camParams + 3);
	int camOffset = camParams * numCams;
#pragma omp parallel for
	for (int i = 0; i < numObs; ++i) {
		(*camsIds)[i] = h_Jcols[i * step] / camParams;
		(*ptsIds)[i] = (h_Jcols[i * step + camParams] - camOffset) / 3;
	}
}

void exCSPts(int numObs, int numPoints, int *ptsIds, int *maxCams, int **csPts) {
	(*maxCams) = 0;
	(*csPts) = (int*)malloc((numPoints + 1) * sizeof(int));
	assert((*csPts) != NULL);
	memset((void*)(*csPts), 0, (numPoints + 1) * sizeof(int));
	int pid = 0, actCams;
	for (int i = 0; i < numObs; ++i) {
		if (pid != ptsIds[i]) {
			actCams = i - (*csPts)[pid];
			(*csPts)[++pid] = i;
			if (actCams >(*maxCams))
				(*maxCams) = actCams;
		}
	}
	(*csPts)[++pid] = numObs;
}

// Fix the reconstruction by fixing three points in Jacobian
void fixPts(tp *s, int *pts, cov::Options &opt, cov::Statistic &statistic, SSM *J) {
	if (pts == NULL) return;
	sort(pts, pts + 2);
	cout << "id: " << pts[0] << ", " << pts[1] << ", " << pts[2] << " ... ";
#ifdef USE_MATLAB
	mexPrintf("fixed points: %d, %d, %d\n", pts[0], pts[1], pts[2]);
#endif
	statistic.fixedPts = new int[3]{ pts[0], pts[1], pts[2] };

	// column ids to remove
	int remCols[9];
	auto A = J->get_sA();
	int camOffset = opt._numCams * opt._camParams;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)
			remCols[i * 3 + j] = camOffset + pts[i] * 3 + j;
	}

	// new ids 
	int *newIds = (int*)malloc(A->ncols * sizeof(int));
	assert(newIds != NULL);
	memset(newIds, 0, A->ncols * sizeof(int));
	for (int i = 0; i < 9; i++)
		newIds[remCols[i]] = -1;
	int off = 0;
	for (int i = 0; i < A->ncols; i++) {
		off += newIds[i];
		newIds[i] = i + off;
	}

	// decompose and compose jacobian
	vector<int> rows;
	vector<int> cols;
	vector<double> vals;
	vector<int> cols2;
	vector<double> vals2;
	vector<bool> use_cols;
	rows.push_back(0);
	int offset = 0;
	for (int i = 0; i < A->nrows; i++) {
		for (int j = A->row[i]; j < A->row[i + 1]; j++) {
			bool use = true;
			for (int k = 0; k < 9; k++) {
				if (A->col[j] == remCols[k])
					use = false;
			}
			if (use) {
				cols.push_back(newIds[A->col[j]]);
				vals.push_back(A->val[j]);
			}
		}
		rows.push_back(cols.size());
	}

	// remove redundant cameras
	int error_ofset = 0;
	for (int i = 1; i < rows.size(); i++) {
		bool use = true;
		if (rows[i] != (i * (opt._camParams + 3) - error_ofset) ) { // some camera without point corespondence 
			error_ofset += 3;
			use = false;
		}
		else {
			use_cols.push_back(use);
			use_cols.push_back(use);
			use_cols.push_back(use);
		}
		for (int j = 0; j < opt._camParams; j++)
			use_cols.push_back(use);
	}
	for (int i = 0; i < cols.size(); i++) {
		if (use_cols[i]) {
			cols2.push_back(cols[i]);
			vals2.push_back(vals[i]);
		}
	}
	rows.clear();
	int nr = cols2.size() / (opt._camParams + 3) + 1;
	rows.resize(nr);
	for (int i = 0; i < nr; i++)
		rows.at(i) = i * (opt._camParams + 3);

	free(newIds);
	J->set_sA(make_shared<CRS>( nr - 1, A->ncols - 9, rows.data(), cols2.data(), vals2.data()));
	*s = t(*s, "Computing sJJ ... ", &(statistic.timeFixJ));
}

void findCams2Point(cov::Options &options, SSM* J, std::vector< std::vector<int> > &pts2cams_ids) {
	int numPars = options._camParams + 3;
	int cams_offset = options._numCams * options._camParams;
	for (int i = 0; i < (J->nrows() / 2); i++) {
		int ptId = (J->col(i * 2 * numPars + options._camParams) - cams_offset) / 3;
		int camId = J->col(i * 2 * numPars) / options._camParams;
		std::vector<int> tmp = pts2cams_ids[ptId];
		tmp.push_back(camId);
		pts2cams_ids[ptId] = tmp;
	}
}

void composeCams2PointJacobian(cov::Options &opt, int point_id, std::vector<int> &cams, SSM *J, double *scale, Eigen::MatrixXd &A, Eigen::MatrixXd &B) {
	int ncams = cams.size();
	int numPars = opt._camParams + 3;
	int cams_offset = opt._numCams * opt._camParams;
	A = Eigen::MatrixXd(2 * ncams, 3);
	B = Eigen::MatrixXd::Zero(2 * ncams, (opt._camParams + 2) * ncams);
	Eigen::VectorXd sA(3);
	Eigen::VectorXd sB((opt._camParams + 2) * ncams);

	// Fill A,B
	int t = 0;
	for (int i = 0; i < (J->nrows() / 2); i++) {
		
		int k = i * 2 * numPars;
		if (point_id == ((J->col(k+opt._camParams) - cams_offset) / 3)) {
			
			// cameras and observations
			for (int j = 0; j < opt._camParams; j++) {
				B(t * 2 + 0, t*(opt._camParams + 2) + j) = J->val(k + j);
				B(t * 2 + 1, t*(opt._camParams + 2) + j) = J->val(k + j + numPars);
				sB(t*(opt._camParams + 2) + j) = 1/scale[J->col(k + j)];
			}
			B(t * 2 + 0, t*(opt._camParams + 2) + opt._camParams + 0) = -1;
			B(t * 2 + 1, t*(opt._camParams + 2) + opt._camParams + 1) = -1;
			sB(t*(opt._camParams + 2) + opt._camParams + 0) = 1;
			sB(t*(opt._camParams + 2) + opt._camParams + 1) = 1;

			// points
			for (int j = 0; j < 3; j++) {
				A(t * 2 + 0, j) = J->val(k + opt._camParams + j);
				A(t * 2 + 1, j) = J->val(k + opt._camParams + j + numPars);
				sA(j) = 1/scale[J->col(k + opt._camParams + j)];
			}
			t++;
		}
	}
	A *= sA.asDiagonal();
	B *= sB.asDiagonal();
}

void composeCamCovariances(cov::Options &opt, std::vector<int> &cams, SDM *iZ, const double *scale, Eigen::MatrixXd &Sigma) {
	int ncams = cams.size();
	int numCamU = opt._camParams + 2;
	int N = iZ->ncols();  // is N x N matrix
	Sigma = Eigen::MatrixXd::Zero(ncams * numCamU, ncams * numCamU);

	for (int i = 0; i < ncams; i++) {
		int sigOff = i * numCamU;
		int izOff = cams[i] * opt._camParams;

		// Copy submatrix iZ to submatrix Sigma
		for (int j = 0; j < opt._camParams; j++) {		// rows 
			for (int k = 0; k < opt._camParams; k++) {	// columns
				Sigma(sigOff + j, sigOff + k) = iZ->val(izOff + j, izOff + k) * scale[izOff + j] * scale[izOff + k];
			}
		}
		Sigma(sigOff + opt._camParams + 0, sigOff + opt._camParams + 0) = 1;
		Sigma(sigOff + opt._camParams + 1, sigOff + opt._camParams + 1) = 1;
	}
}

void fillPtUnc2OutArray(int ptId, Eigen::MatrixXd &Sigma, cov::Options &opt, double *out) {
	// fixed points offset
	for (int i = 0; i < 3; i++)
		if (ptId >= opt._pts2fix[i]){ ptId++; }

	int offset = ptId * 6;  // one point is represented by 6 vlaues 
	out[offset + 0] = Sigma(0, 0);
	out[offset + 1] = Sigma(0, 1);
	out[offset + 2] = Sigma(0, 2);
	out[offset + 3] = Sigma(1, 1);
	out[offset + 4] = Sigma(1, 2);
	out[offset + 5] = Sigma(2, 2);
}


void computeCovariances(cov::Options &options, cov::Statistic &statistic, ceres::CRSMatrix &jacobian, double *camUnc, double *ptsUnc) {
	if (camUnc == NULL || ptsUnc == NULL) return;
#if defined(_OPENMP)
	omp_set_num_threads(8);
#endif
	Eigen::setNbThreads(8);
	magma_queue_t queue;
	magma_int_t info = magma_init();
	TESTING_CHECK(info);
	magma_queue_create(info, &queue);
	magma_print_environment();

    cout << "\n------ " << EAlgorithm_enumToString(options._algorithm) << " ------\n";
	tp s = Clock::now(); tp s1 = s; cout << "Creating sJ ... ";
       
	// Create sparse matrix with separated scale coefficient 
	SSM *J = new SSM(jacobian.num_rows, jacobian.num_cols, jacobian.rows.data(), jacobian.cols.data(), jacobian.values.data());        
	if(options._debug) J->printBlock2Matlab3("J",0,0,J->nrows(),J->ncols());
        
        // Main algorithm
	int Ncams = options._numCams * options._camParams;
	int Npar = Ncams + options._numPoints * 3;
	double *diagRightScaleJ = NULL;
	SDM *iZ = NULL;
	switch (options._algorithm) {
    case cov::eAlgorithmSvdQrIteration:
		s = t(s, "Computing sJJ ... ", &(statistic.timeCreateJ));
		iZ = new SDM(Npar, Npar);
		svdInverse(&info, Npar, options, *J, iZ);
		removeScaleJ4Z(&s, diagRightScaleJ, options, iZ, camUnc);
		break;

    case cov::eAlgorithmSvdDeviceAndconquer:
		s = t(s, "Computing sJJ ... ", &(statistic.timeCreateJ));
		iZ = new SDM(Npar, Npar);
		svdInverse(&info, Npar, options, *J, iZ);
		removeScaleJ4Z(&s, diagRightScaleJ, options, iZ, camUnc);
		break;

    case cov::eAlgorithmSvdTaylorExpansion:
		s = t(s, "Fix pts sJ ... ", &(statistic.timeCreateJ));
		fixPts(&s, options._pts2fix, options, statistic, J);
		if(options._debug) J->printBlock2Matlab3("Jfix", 0, 0, J->nrows(), J->ncols());

		SSM *Y = new SSM();
		iZ = new SDM(Ncams, Ncams);
		composeZ(&s, options, statistic, *J, &diagRightScaleJ, Y, iZ);   // "iZ" contains Z	
		if(options._debug) iZ->printBlock2Matlab3("Z", 0, 0, iZ->nrows(), iZ->ncols());

		teInverse(&s, Ncams, options, statistic, iZ);		// "iZ" is inversed to iZ
		if(options._debug) iZ->printBlock2Matlab3("iZ", 0, 0, iZ->nrows(), iZ->ncols());

		removeScaleJ4Z(&s, diagRightScaleJ, options, iZ, camUnc);
		memset(ptsUnc, 0, 6 * (options._numPoints-3) * sizeof(double));
		delete Y;
		break;
	}

	// Covariances for points - propagation of implicit function
	std::vector< std::vector<int> > pts2cams_ids(options._numPoints-3, std::vector<int>());
	findCams2Point(options, J, pts2cams_ids);
#pragma omp parallel for
	for (int i = 0; i < pts2cams_ids.size(); i++) {
		Eigen::MatrixXd A, B, pinvA, Sigma_cam_u;
		composeCams2PointJacobian(options, i, pts2cams_ids[i], J, diagRightScaleJ, A, B);
		composeCamCovariances(options, pts2cams_ids[i], iZ, diagRightScaleJ, Sigma_cam_u);
		pinvA = A.completeOrthogonalDecomposition().pseudoInverse();
		Eigen::MatrixXd Sigma_pt = (pinvA * B) * Sigma_cam_u * (B.transpose() * pinvA.transpose());
		fillPtUnc2OutArray(i, Sigma_pt, options, ptsUnc);
	}
	statistic.timePtsUnc = timeDuration(s, Clock::now());
    cout << "\nPoints uncertainty computed in ... " << (statistic.timePtsUnc) << "s\n";

	free(diagRightScaleJ);
	delete iZ;
	delete J;

	statistic.timeAll = timeDuration(s1, Clock::now());
	cout << "\nAlgorithm done in ... " << (statistic.timeAll) << "s\n";

	magma_queue_destroy(queue);
	magma_finalize();
}
