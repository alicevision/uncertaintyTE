#include "compute.h"
#include "auxCmd.h"
#include "ScaledSparseMatrix.h" 
#include "ScaledDenseMatrix.h" 


using namespace std;
using namespace Eigen;

typedef ScaledSparseMatrix SSM;
typedef ScaledDenseMatrix SDM;
typedef std::unique_ptr<SSM> uSM;
typedef std::unique_ptr<SDM> uDM;
typedef SparseMatrix<float, RowMajor> SM;
typedef MatrixXf DM;

#ifdef _WIN32	
double timeDuration(tp from, tp to) {
	return std::chrono::duration_cast<std::chrono::nanoseconds>(to - from).count() * 1e-9;
}
#else
double timeDuration(tp from, tp to) {
	return 0;
}
#endif
 
tp t(tp s, char* txt) {
	cout << " " << timeDuration(s, Clock::now()) << "s\n";
	cout << txt;
	return Clock::now();
}

tp t(tp s, char* txt, double *time) {
	*time = timeDuration(s, Clock::now());
	cout << " " << (*time) << "s\n";
	cout << txt;
	return Clock::now();
}

std::string algorihm2str2(int alg) {
	switch (alg) {
	case 0: return "SVD_QR_ITERATION";
	case 1: return "SVD_DEVIDE_AND_CONQUER";
	case 2: return "TAYLOR_EXPANSION";
	default: return "not defined";
	}
}

void symmetrizeMatrix(int N, SDM* A) {
#pragma omp parallel for
	for (int i = 0; i < N; ++i) {
		for (int j = i; j < N; ++j)
			A->set(j, i, A->sval(i,j));
	}
}

int factorial(int n){
	if (n > 1)
		return n * factorial(n - 1);
	else
		return 1;
}

void iZupdate(SDM *iZ, double coeff, SDM *iZadd ) {
	// *iZ += k * iZadd;
	for (int i = 0; i < iZ->ncols(); ++i) {
		for (int j = 0; j < iZ->nrows(); ++j)
			iZ->set(j, i, (iZ->val(j,i) + (coeff * iZadd->c()) * iZadd->sval(j,i)));
	}
	iZ->set_c(1);
}

void composeZ(tp *s, cov::Options &options, cov::Statistic &statistic, SSM &J, double **diagRightScaleJ, SSM *Y, SDM *Z) {
	// Scale Jacobian from left by a vector
	double csJ = 1;
	J.scaleMat(RIGHT, diagRightScaleJ, &csJ);
	for (int i = 0; i < J.ncols(); ++i)
		(*diagRightScaleJ)[i] = 1 / (csJ * (*diagRightScaleJ)[i]);
	*s = t(*s, "Computing sM ... ", &(statistic.timeNormJ));

	// Compute cJJ * sJJ = sJ'sJ
	SSM Jt( J.trn() );		// create transpose matix and save reference to it
	SSM M(Jt * J);
	*s = t(*s, "Split sM -> sU, sW, sV ... ", &(statistic.timeMultiplyJJ));

	// Split cJJ -> U, W, V
	SSM *U = new SSM(), *W = new SSM(), *V = new SSM(), *iV = new SSM();
	int camBlockSize = options._numCams * options._camParams;
	M.splitTo3Blocks(camBlockSize, camBlockSize, U, W, V);
	*s = t(*s, "Computing siV ... ", &(statistic.timeSplitJJ));

	// Compute inverse of V: sV -> isV
	V->inv3x3blockSymmDiag(iV);
	delete V;
	*s = t(*s, "Computing sZ ... ", &(statistic.timeInvV));

	// Compute Z = Schur complement of V;  U,iV,W  -> Z
	*Y = std::move( (*W) * (*iV) );
	SSM sWt = std::move( W->trn() );
	SSM YsWt = std::move( (*Y) * sWt );
	SSM sZ = std::move( (*U) - YsWt );
	*Z = std::move( SDM(sZ) );
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
	SDM iZadd = std::move( (*iZ) * (*iZ) );
	for (int i = 1; i < 20; ++i) {
		k = pow(lambda,i) / factorial(i - 1);
		change = abs(k) * iZadd.absMax();
   #ifdef USE_MATLAB
		mexPrintf("\n cykle %d, coeff: %e, change: %e ",i,k,change);
   #endif
		cout << "\n>>> cykle " << i << ", coeff: " << k << ", change: " << change;
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
					camUnc[++l] = (iZ->val(j,k) + iZ->val(k, j) / 2) * diagRightScaleJ[j] * diagRightScaleJ[k];
				else
					camUnc[++l] = (iZ->val(j, k) + iZ->val(k, j) / 2);
			}	
		}
	}
	cout << " " << timeDuration(*s, Clock::now()) << "s\n";
}

void svdInverse(magma_int_t *info, int N, cov::Options &options, SSM &J, SDM *iZ) {
	SDM Z(J.trn() * J);	
	double *sv = (double*)malloc(N*sizeof(double));
	double *U = (double*)malloc(N*N*sizeof(double));
	double *Vt = (double*)malloc(N*N*sizeof(double));
	assert(sv != NULL);
	assert(U != NULL);
	assert(Vt != NULL);
	int lwork, *iwork;
	double *hwork;

	// Use Intel MKL Lapack SVD   ( aprox. 3-4x faster then GPU variant working with double )
	switch (options._algorithm) {
	case SVD_QR_ITERATION:				// additional memory requirements: N*N + 3*N + 2*N*32 DOUBLE
		lwork = N*N + 3 * N + 2 * N * 32;
		hwork = (double*)malloc(lwork * sizeof(double));
		assert(hwork != NULL);
		lapackf77_dgesvd(lapack_vec_const(MagmaAllVec), lapack_vec_const(MagmaAllVec), &N, &N,
			Z.getMatPtr(), &N, sv, U, &N, Vt, &N, hwork, &lwork, info);
		free(hwork);
		break;

	case SVD_DEVIDE_AND_CONQUER:		// additional memory requirements: 4*N*N + 7*N DOUBLE  + 8*N INT
		lwork = 4 * N*N + 7 * N;
		hwork = (double*)malloc(lwork * sizeof(double));
		assert(hwork != NULL);
		iwork = (int*)malloc(8 * N*sizeof(int));
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
	} else {
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
	(*csPts) = (int*)malloc((numPoints + 1)* sizeof(int));
	assert((*csPts) != NULL);
	memset((void*)(*csPts), 0, (numPoints + 1)* sizeof(int));
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
		for(int j = 0; j < 3; j++)
			remCols[i*3 + j] = camOffset + pts[i] * 3 + j;
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
	rows.push_back(0); 
	int offset = 0;
	for (int i = 0; i < A->nrows; i++) {
		for (int j = A->row[i]; j < A->row[i + 1]; j++) {
			bool use = true;
			for (int k = 0; k < 9; k++) {
				if (A->col[j] == remCols[k]) 
					use = false;
			}
			if (use){
				cols.push_back(newIds[A->col[j]]);
				vals.push_back(A->val[j]);
			}
		}
		rows.push_back(cols.size());
	}

	free(newIds);
	J->set_sA(make_shared<CRS>(A->nrows, A->ncols - 9, rows.data(), cols.data(), vals.data()));
	*s = t(*s, "Computing sJJ ... ", &(statistic.timeFixJ));
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

	cout << "\n------ " << algorihm2str2(options._algorithm) << " ------\n";
	tp s = Clock::now(); tp s1 = s; cout << "Creating sJ ... ";

	// Create sparse matrix with separated scale coefficient 
	SSM *J = new SSM( jacobian.num_rows, jacobian.num_cols, jacobian.rows.data(), jacobian.cols.data(), jacobian.values.data());

	// Main algorithm
	int Ncams = options._numCams * options._camParams;
	int Npar = Ncams + options._numPoints * 3;
	double *diagRightScaleJ = NULL; 
	SDM *iZ = NULL;
	switch (options._algorithm) {
		case SVD_QR_ITERATION:
			s = t(s, "Computing sJJ ... ", &(statistic.timeCreateJ));
			iZ = new SDM(Npar, Npar);
			svdInverse(&info, Npar, options, *J, iZ);
			removeScaleJ4Z(&s, diagRightScaleJ, options, iZ, camUnc);
			break;
	
		case SVD_DEVIDE_AND_CONQUER:
			s = t(s, "Computing sJJ ... ", &(statistic.timeCreateJ));
			iZ = new SDM(Npar, Npar);
			svdInverse(&info, Npar, options, *J, iZ);
			removeScaleJ4Z(&s, diagRightScaleJ, options, iZ, camUnc);
			break;
	
		case TAYLOR_EXPANSION:
			s = t(s, "Fix pts sJ ... ", &(statistic.timeCreateJ));
			fixPts(&s, options._pts2fix, options, statistic, J);
			SSM *Y = new SSM(); 
			iZ = new SDM(Ncams, Ncams);
			composeZ(&s, options, statistic, *J, &diagRightScaleJ, Y, iZ);   // "iZ" contains Z
			teInverse(&s, Ncams, options, statistic, iZ);		// "iZ" is inversed to iZ
			removeScaleJ4Z(&s, diagRightScaleJ, options, iZ, camUnc);
			memset(ptsUnc, 0, 6 * options._numPoints * sizeof(double));
			free(diagRightScaleJ);
			delete Y;
			break;
	}
	delete iZ;
	delete J;
	
	statistic.timeAll = timeDuration(s1, Clock::now());
	cout << "\nAlgorithm done in ... " << (statistic.timeAll) << "s\n";

	magma_queue_destroy(queue);
	magma_finalize();

//s = t(s, "Computing points uncertainty ... ");
//	cout << "Points uncertainty ...";
//	clock_t spts = clock();
//	int *camsIds, *ptsIds, *csPts, maxCams;
//	findICP(options._numObs, options->_camParams, options->_numCams, jacobian->_cols.data(), &camsIds, &ptsIds);
//	thrust::sort_by_key(ptsIds, ptsIds + options->_numObs, camsIds);
//	exCSPts(options->_numObs, options->_numPoints, ptsIds, &maxCams, &csPts);
//
////#pragma omp parallel for
//	for (int i = 0; i < options->_numPoints; ++i) {
//		int Ncams = csPts[i + 1] - csPts[i];
//		DM WVpt(Ncams*options->_camParams, 3);
//		DM iUVWpt(Ncams*options->_camParams, Ncams*options->_camParams);
//		DM Cpt(3, 3);
//		// W -> Wpt; iUVW -> iUVWpt
//		for (int j = 0; j < Ncams; ++j) {
//			int row = camsIds[csPts[i] + j] * options->_camParams;
//			WVpt.block(j*options->_camParams, 0, options->_camParams, 3) = Y.block(row, i * 3, options->_camParams, 3);
//			for (int k = 0; k < Ncams; ++k)
//				iUVWpt.block(j*options->_camParams, k*options->_camParams, options->_camParams, options->_camParams) =
//				iUVW.block(row, camsIds[csPts[i] + k] * options->_camParams, options->_camParams, options->_camParams);
//		}
//		// Dense multiplication of iUVW, iUVWpt -> list of covariance of points
//		Cpt = WVpt.transpose() * iUVWpt * WVpt;
//
//		// Unscale values by Jacobian coeficients and write them to the output
//		int ptScaleOffset = options->_camParams * options->_numCams + 3 * i;
//		ptsUnc[i * 6 + 0] = Cpt(0, 0) * jacobScale[ptScaleOffset] * jacobScale[ptScaleOffset];
//		ptsUnc[i * 6 + 1] = ((Cpt(0, 1) + Cpt(1, 0)) / 2) * jacobScale[ptScaleOffset] * jacobScale[ptScaleOffset + 1];
//		ptsUnc[i * 6 + 2] = ((Cpt(0, 2) + Cpt(2, 0)) / 2) * jacobScale[ptScaleOffset] * jacobScale[ptScaleOffset + 2];
//		ptsUnc[i * 6 + 3] = Cpt(1, 1) * jacobScale[ptScaleOffset + 1] * jacobScale[ptScaleOffset + 1];
//		ptsUnc[i * 6 + 4] = ((Cpt(1, 2) + Cpt(2, 1)) / 2) * jacobScale[ptScaleOffset + 1] * jacobScale[ptScaleOffset + 2];
//		ptsUnc[i * 6 + 5] = Cpt(2, 2) * jacobScale[ptScaleOffset + 2] * jacobScale[ptScaleOffset + 2];
//	}
//	cout << " " << timeDuration(spts, clock()) << "s\n";
//	cout << "Complete in time: " << timeDuration(s1, clock()) << "s";
}