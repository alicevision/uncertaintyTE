#include "ScaledDenseMatrix.h"


ScaledDenseMatrix::ScaledDenseMatrix() {}

ScaledDenseMatrix::ScaledDenseMatrix(int nrows, int ncols) : _nr(nrows), _nc(ncols) {
	_sA = (double*)malloc(_nr*_nc * sizeof(double));
	assert(_sA != NULL);
}

ScaledDenseMatrix::ScaledDenseMatrix(int nrows, int ncols, double* values) : _nr(nrows), _nc(ncols) {
	_sA = (double*)malloc(_nr*_nc * sizeof(double));
	assert(_sA != NULL);
#pragma omp parallel for
	for (int i = 0; i < _nc; i++) {
		for (int j = 0; j < _nr; j++)
			_sA[j + i*_nr] = values[i + j*_nc];
	}
}

ScaledDenseMatrix::ScaledDenseMatrix(const ScaledDenseMatrix &B) : _cA(B._cA), _nr(B._nr), _nc(B._nc) {
	_sA = (double*)malloc(_nr*_nc * sizeof(double));
	assert(_sA != NULL);
	memcpy(_sA, B._sA, _nr*_nc * sizeof(double));
}

ScaledDenseMatrix::ScaledDenseMatrix(ScaledSparseMatrix &A) {
	_nr = A.nrows();
	_nc = A.ncols();
	_cA = A.c();
	_sA = (double*)malloc(_nr*_nc * sizeof(double));		//magma_dmalloc_pinned(&_sA, _nr*_nc);
	assert(_sA != NULL);
	memset(_sA, 0, _nr*_nc * sizeof(double));
	//#pragma omp parallel for
	for (int i = 0; i < _nr; i++) {
		for (int j = A.row(i); j < A.row(i + 1); j++)
			_sA[i + A.col(j) *_nr ] = A.sval(j);
	}
}

ScaledDenseMatrix::ScaledDenseMatrix(ScaledSparseMatrix &&A) {
	_nr = A.nrows();
	_nc = A.ncols();
	_cA = A.c();
	_sA = (double*)malloc(_nr*_nc * sizeof(double));		//magma_dmalloc_pinned(&_sA, _nr*_nc);
	assert(_sA != NULL);
	memset(_sA, 0, _nr*_nc * sizeof(double));
	//#pragma omp parallel for
	for (int i = 0; i < _nr; i++) {
		for (int j = A.row(i); j < A.row(i + 1); j++)
			_sA[i + A.col(j) *_nr] = A.sval(j);
	}
}

ScaledDenseMatrix::~ScaledDenseMatrix() {
	free(_sA);
}

int ScaledDenseMatrix::nrows() {
	return _nr;
}
int ScaledDenseMatrix::ncols() {
	return _nc;
}
double ScaledDenseMatrix::c() const {
	return _cA;
}
void ScaledDenseMatrix::set_c(double c) {
	_cA = c;
}
double* ScaledDenseMatrix::getMatPtr() {
	return _sA;
}
double ScaledDenseMatrix::sval(int i) const {
	return _sA[i];
}
double ScaledDenseMatrix::sval(int row, int col) const {
	return _sA[row + col*_nr];
}
double ScaledDenseMatrix::val(int i) const {
	return (_cA * _sA[i]);
}
double ScaledDenseMatrix::val(int row, int col) const {
	return (_cA * _sA[row + col*_nr]);
}
void ScaledDenseMatrix::set(const int i, const double val) {
	_sA[i] = val;
}
void ScaledDenseMatrix::set(const int i, const int j, const double val) {
	_sA[i + j*_nr] = val;
}

void ScaledDenseMatrix::print() {
	std::cout << _cA << "*";
	magma_dprint(_nr, _nc, _sA, _nr);
}

#ifdef USE_MATLAB
void ScaledDenseMatrix::printBlock2Matlab(char* name, int row_from, int col_from, int row_to, int col_to) {
	std::cout << "\n\n" << name << " = [";
	for (int i = row_from; i < row_to; i++) {
		for (int j = col_from; j < col_to; j++)
			std::cout << val(i, j) << " ";
		std::cout << (i == (row_to - 1) ? "" : ";");
	}
	std::cout << "];\n\n\n";
}

void ScaledDenseMatrix::printBlock2Matlab2(char* name, int row_from, int col_from, int row_to, int col_to) {
	mexPrintf("\n\n %s = [",name);
	for (int i = row_from; i < row_to; i++) {
		for (int j = col_from; j < col_to; j++)
			mexPrintf("%e ", val(i, j));
		mexPrintf("%s ", (i == (row_to - 1) ? "" : "; ... \n"));
	}
	mexPrintf("];\n\n\n");
}

void ScaledDenseMatrix::printBlock2Matlab3(char* name, int row_from, int col_from, int row_to, int col_to) {
	std::ofstream file(std::string(name) + std::string(".txt"));
	file << name << " = [";
	for (int i = row_from; i < row_to; i++) {
		for (int j = col_from; j < col_to; j++)
			file << val(i, j) << " ";
		file << (i == (row_to - 1) ? "" : ";");
	}
	file << "];";
	file.close();
}
#endif

/*
Scale the matrix by diagonal matrix represented by A = cLR * sLR * sA
from left L or right R
*/
void ScaledDenseMatrix::scaleMat(int type, double **sLR, double *cLR) {
	double asum = 0;
	switch (type) {
	case LEFT:
		// auxiliary vectors
		*sLR = (double*)malloc(nrows() * sizeof(double));
		assert((*sLR) != NULL);
		memset(*sLR, 0, nrows() * sizeof(double));

		// encout scale 
		for (int i = 0; i < ncols(); ++i) {
			for (int j = 0; j < nrows(); ++j)
				(*sLR)[j] += _sA[i * nrows() + j] * _sA[i * nrows() + j];
		}
		for (int j = 0; j < nrows(); ++j) {
		//asum += (*sLR)[j];
			(*sLR)[j] = sqrt((*sLR)[j]);
		}

		// scale matrix
		for (int i = 0; i < ncols(); ++i) {
			for (int j = 0; j < nrows(); ++j)
				_sA[i * nrows() + j] *= 1 / (*sLR)[j];
		}
		//for (int j = 0; j < nrows(); ++j)
		//	(*sLR)[j] = (*sLR)[j] / asum;
		//*cLR = asum;
		*cLR = 1;

		break;
	case RIGHT:
		// auxiliary vectors
		(*sLR) = (double*)malloc(ncols() * sizeof(double));
		assert((*sLR) != NULL);
		memset((*sLR), 0, ncols() * sizeof(double));

		// encout scale 
		for (int i = 0; i < ncols(); ++i) {
			for (int j = 0; j < nrows(); ++j)
				(*sLR)[i] += _sA[i * nrows() + j] * _sA[i * nrows() + j];
			//asum += (*sLR)[i];
			(*sLR)[i] = sqrt((*sLR)[i]);
		}

		// scale matrix
		for (int i = 0; i < ncols(); ++i) {
			for (int j = 0; j < nrows(); ++j)
				_sA[i * nrows() + j] *= 1 / (*sLR)[i];
			//(*sLR)[i] = (*sLR)[i] / asum;
		}
		//*cLR = asum;
		*cLR = 1;
		break;
	}
}

// Find maximal coefficient 
double ScaledDenseMatrix::absMax() {
	double amax = DBL_MIN;
	for (int i = 0; i < _nc; ++i) {
		for (int j = 0; j < _nr; ++j) {
			if (amax < abs(_sA[i*_nr + j]))
				amax = abs(_sA[i*_nr + j]);
		}
	}
	return abs(_cA) * amax;
}

// Lower triangle is coppied to the upper triangle
void ScaledDenseMatrix::symmL2U() {
	if (_nr != _nc) {
		std::cerr << "Matrix can't be symmetrize because number of rows is not equal number of columns";
		exit(1);
	}
#pragma omp parallel for
	for (int i = 0; i < _nr; ++i) {
		for (int j = i + 1; j < _nr; ++j)
			_sA[j*_nr + i] = _sA[i*_nr + j];
	}
}

// General inverse: "Z" -> iZ
void ScaledDenseMatrix::inv() {
	// "Z" -> cL * sL * "sZ" * sR * cR
	double cL = 1, cR = 1;
	double *sL, *sR;
	scaleMat(LEFT, &sL, &cL);
	scaleMat(RIGHT, &sR, &cR);

	// LU decomposition
	magma_int_t     *ipiv, iunused[1], info;
	TESTING_CHECK(magma_imalloc_cpu(&ipiv, std::min(_nr, _nc)));
	lapackf77_dgetrf(&(_nr), &(_nc), _sA, &(_nr), ipiv, &info);
	if (info != 0)
		std::cerr << "Lapack LU decomposition error.";

	// init size of the work array
	double *work, unused[1] = { 0 }, tmp = 0;
	magma_int_t lwork = -1;
	lapackf77_dgetri(&(_nr), unused, &(_nr), iunused, &tmp, &lwork, &info);
	if (info != 0)
		std::cerr << "Lapack init work array size error.";
	lwork = static_cast<int>(tmp);
	TESTING_CHECK(magma_dmalloc_cpu(&work, static_cast<int>(lwork)));

	// inverse
	lapackf77_dgetri(&(_nr), _sA, &(_nr), ipiv, work, &lwork, &info);
	if (info != 0)
		std::cerr << "Lapack inverse error.";
	_cA = 1 / _cA;
	
	// inverse: cL -> icL, sL -> isL, .... cR -> icR, sR -> isR
	cL = 1 / cL;
	for (int i = 0; i < nrows(); ++i)
		sL[i] = 1 / sL[i];
	cR = 1 / cR;
	for (int i = 0; i < ncols(); ++i)
		sR[i] = 1 / sR[i];
	// compose: iZ
	_cA *= cL * cR;
	for (int i = 0; i < ncols(); ++i) {
		for (int j = 0; j < nrows(); ++j)
			_sA[i*nrows() + j] *= sL[i] * sR[j];
	}

	free(sL);
	free(sR);
	magma_free_cpu(ipiv);
	magma_free_cpu(work);
}

// Deep copy function
void ScaledDenseMatrix::copy(const ScaledDenseMatrix& other) {
	_nr = other._nr;
	_nc = other._nc;
	_cA = other._cA;
	_sA = (double*)malloc(_nr * _nc * sizeof(double));
	assert(_sA != NULL);
	memcpy(_sA, other._sA, _nr * _nc * sizeof(double));
}

// Add or subtract
void ScaledDenseMatrix::addSDM(int sign, ScaledDenseMatrix& B) {
	if (nrows() != B.nrows() || ncols() != B.ncols()) {
		std::cerr << "Addition or subtraction of matrices is not possible. Matrices have different number of rows or columns.";
		exit(1);
	}
	for (int i = 0; i < _nc; ++i) {
		for (int j = 0; j < _nr; ++j)
			_sA[i*_nr + j] = _cA * _sA[i*_nr + j] + sign * B._cA * B._sA[i*_nr + j];
	}
	_cA = 1;
}

// Deep copy assignment operator
ScaledDenseMatrix& ScaledDenseMatrix::operator=(const ScaledDenseMatrix &other) {
	//std::cout << "called copy operator";
	if (this != &other) {					// self-assignment check
		copy(other);
	}
	return *this;
}

// Multiplication operator
ScaledDenseMatrix& ScaledDenseMatrix::operator*=(ScaledDenseMatrix& B) {
	// Scale matrices
	double cl = 1, cr = 1;
	double *sL, *sR;
	scaleMat(LEFT, &sL, &cl);
	B.scaleMat(RIGHT, &sR, &cr);

	// Perform multiplication
	double alpha = 1, beta = 0;
	double *hC = (double*)malloc(_nr * B._nc * sizeof(double));
	assert(hC != NULL);
	blasf77_dgemm(lapack_trans_const(MagmaNoTrans),
		lapack_trans_const(MagmaNoTrans), &(_nr), &(B._nc), &(_nc),
		&alpha, _sA, &(_nr),
		B._sA, &(B._nr),
		&beta, hC, &(_nr));

	// Copy to this
	free(_sA);
	_sA = (double*)malloc(_nr * B._nc * sizeof(double));
	assert(_sA != NULL);
	memset(_sA, 0, _nr * B._nc * sizeof(double));
	memcpy(_sA, hC, _nr * B._nc * sizeof(double));

	// Unscale matrices
	double c = 1 / (cr * sL[0] * sR[0]);
	_cA = (_cA * cl * B._cA) / c;
	for (int i = 0; i < B._nc; ++i) {
		for (int j = 0; j < _nr; ++j)
			_sA[i*_nr + j] *= c * cr * sL[j] * sR[i];
	}
	for (int i = 0; i < B._nc; ++i) {
		for (int j = 0; j < B._nr; ++j)
			B._sA[i*B._nr + j] *= cr * sR[i];
	}
	free(sL);
	free(sR);
	free(hC);
	return *this;
}

// Multiplication with scalar
ScaledDenseMatrix& ScaledDenseMatrix::operator*=(double &B) {
	_cA *= B;
	return *this;
}

// Addition operator
ScaledDenseMatrix& ScaledDenseMatrix::operator+=(ScaledDenseMatrix& B) {
	addSDM(1, B);
	return *this;
}

// Subtraction operator
ScaledDenseMatrix& ScaledDenseMatrix::operator-=(ScaledDenseMatrix& B) {
	addSDM(-1, B);
	return *this;
}













/*
Multiply two matricies with scale
- create copy of lhs and chage it by call *=
*/
ScaledDenseMatrix operator*(ScaledDenseMatrix lhs, ScaledDenseMatrix& rhs) {
	return lhs *= rhs;
}

/*
Multiply by scalar
- create copy of lhs and chage it by call *=
*/
ScaledDenseMatrix operator*(double lhs, ScaledDenseMatrix& rhs) {
	return rhs *= lhs;
}
ScaledDenseMatrix operator*(ScaledDenseMatrix& lhs, double rhs) {
	return lhs *= rhs;
}

/*
Add two matrices with scale
- create copy of lhs and chage it by call +=
*/
ScaledDenseMatrix operator+(ScaledDenseMatrix lhs, ScaledDenseMatrix& rhs) {
	return lhs += rhs;
}

/*
Subtract two matrices with scale
- create copy of lhs and chage it by call -=
*/
ScaledDenseMatrix operator-(ScaledDenseMatrix lhs, ScaledDenseMatrix& rhs) {
	return lhs -= rhs;
}