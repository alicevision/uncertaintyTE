#ifndef COMPUTE_HEADER_INCLUDED
	#include "compute.h"
#endif
#ifndef SSM_HEADER_INCLUDED
	#define SSM_HEADER_INCLUDED
#endif

using namespace Eigen;

// Magma sparse matrix multiplied by scalar
magma_d_matrix operator*(double lhs, magma_d_matrix rhs) {
	for (int i = 0; i < rhs.nnz; i++)
		rhs.val[i] *= lhs;
	return rhs;
}

struct CRS {
	int nrows = 0;
	int ncols = 0;
	int nnz = 0;
	int *row = NULL;
	int *col = NULL;
	double *val = NULL;

	CRS() {};
	
	CRS(const int nr, const int nc) : nrows(nr), ncols(nc) {
		row = (int*)malloc((nrows + 1) * sizeof(int));
		assert(row != NULL);
		memset((void*)row, 0, nrows + 1 * sizeof(int));
	}

	CRS(const int nr, const int nc, const int nnz1) : nrows(nr), ncols(nc), nnz(nnz1) {
		row = (int*)malloc((nrows + 1) * sizeof(int));
		assert(row != NULL);
		col = (int*)malloc(nnz * sizeof(int));
		assert(col != NULL);
		val = (double*)malloc(nnz * sizeof(double));
		assert(val != NULL);
		memset((void*)row, 0, nrows + 1 * sizeof(int));
	}

	CRS(const int nr, const int nc, const int *rows, const int *cols, const double *values): nrows(nr), ncols(nc) {
		nnz = rows[nr];
		row = (int*)malloc((nrows + 1) * sizeof(int));
		assert(row != NULL);
		col = (int*)malloc(nnz * sizeof(int));
		assert(col != NULL);
		val = (double*)malloc(nnz * sizeof(double));
		assert(val != NULL);
		memcpy((void*)row, (void*)rows, (nrows + 1) * sizeof(int));
		memcpy((void*)col, (void*)cols, nnz * sizeof(int));
		memcpy((void*)val, (void*)values, nnz * sizeof(double));
	}

	~CRS() {
		free(row);
		free(col);
		free(val);
	}

	//CRS& operator=(CRS& B) {
	//	nrows = B.nrows;
	//	ncols = B.ncols;
	//	nnz = B.nnz;
	//	row = (int*)malloc((nrows + 1) * sizeof(int));
	//	assert(row == NULL);
	//	col = (int*)malloc(nnz * sizeof(int));
	//	assert(col == NULL);
	//	val = (double*)malloc(nnz * sizeof(double));
	//	assert(val == NULL);
	//	memcpy((void*)row, (void*)B.row, (nrows + 1) * sizeof(int));
	//	memcpy((void*)col, (void*)B.col, nnz * sizeof(int));
	//	memcpy((void*)val, (void*)B.val, nnz * sizeof(double));
	//	return *this;
	//}

	CRS& operator=(CRS&& B) {
		nrows = B.nrows;
		ncols = B.ncols;
		nnz = B.nnz;
		row = std::move(B.row);
		col = std::move(B.col);
		val = std::move(B.val);
		return *this;
	}
};

// Matrix scaled by scalar structure
class ScaledSparseMatrix {
private:	
	double _cA = 1;
	std::shared_ptr<CRS> _sA;

public:
	ScaledSparseMatrix(): _cA(1) {
		_sA = NULL;
	}

	ScaledSparseMatrix(const ScaledSparseMatrix& o): _cA(o._cA) {
		_sA = std::make_shared<CRS>(o._sA->nrows, o._sA->ncols, o._sA->row, o._sA->col, o._sA->val);
	}

	ScaledSparseMatrix(ScaledSparseMatrix&& o) : _cA(o._cA) {
		_sA = std::move(o.get_sA());
	}

	ScaledSparseMatrix(std::shared_ptr<CRS> &A) : _cA(1) {
		_sA = A;
	}

	ScaledSparseMatrix(const int nrows, const int ncols, const int nnz): _cA(1) {
		_sA = std::make_shared<CRS>(nrows, ncols, nnz);
	}

	ScaledSparseMatrix(const int nrows, const int ncols, const int *rows, const int *cols, const double *values) : _cA(1) {
		_sA = std::make_shared<CRS>(nrows, ncols, rows, cols, values);
	}


	std::shared_ptr<CRS> get_sA() {
		return _sA;
	}
	double get_cA() {
		return _cA;
	}
	void set_sA(std::shared_ptr<CRS> &A) {
		_sA = A;
	}
	void set_sA(std::shared_ptr<CRS> &&A) {
		_sA = std::move(A);
	}
	void set_cA(double c) {
		_cA = c;
	}
	int nrows() const {
		return _sA->nrows;
	}
	int ncols() const {
		return _sA->ncols;
	}
	int nnz() const {
		return _sA->nnz;
	}
	int col(const int i) const {
		return _sA->col[i];
	}
	int row(const int i) const{
		return _sA->row[i];
	}
	double c() const {
		return _cA;
	}
	double sval(int i) const {
		return _sA->val[i];
	}
	double val(int i) const{
		return (_cA * _sA->val[i]);
	}
	double val(int i, int j) const {
		for (int k = _sA->row[i]; k < _sA->row[i + 1]; ++k) {
			if (_sA->col[k] == j)
				return _cA * _sA->val[k];
		}
		return 0;
	}
	// Print submatrix up to then rows / cols 
	void print () const {
		std::cout << "\n X = [\n";
		for (int i = 0; i < min(10, nrows()); ++i) {
			for (int j = 0; j < min(10, ncols()); ++j) {
				std::cout << val(i, j) << " ";
			}
			std::cout << "\n";
		}
		std::cout << "]";
	}
	void printAll() const {
		std::cout << "\n X = [\n";
		for (int i = 0; i < nrows(); ++i) {
			for (int j = 0; j < ncols(); ++j) {
				std::cout << val(i, j) << " ";
			}
			std::cout << "\n";
		}
		std::cout << "]";
	}

	void printBlock2Matlab(char* name, int row_from, int col_from, int row_to, int col_to) {
		std::cout << "\n\n" << name << " = zeros(" << (row_to-row_from) << ", " << (col_to-col_from) << ");\n";
		for (int i = row_from; i < row_to; i++) {
			for (int j = row(i); j < row(i + 1); j++) {
				if (col(j) >= col_from && col(j) < col_to)
					std::cout << name << "(" << (i - row_from + 1) << "," << (col(j) - col_from + 1) << ") = " << val(j) << ";";
			}
		}
		std::cout << "\n\n\n";
	}

	void inv3x3blockSymmDiag(ScaledSparseMatrix *ssm) {
		// define output arrays
		ssm->set_sA(std::make_shared<CRS>(nrows(), ncols(), nnz()));
		auto csr = ssm->get_sA();

		// compute and save the inverse values
		memcpy(csr->row, _sA->row, (nrows() + 1) * sizeof(int));
		memcpy(csr->col, _sA->col, nnz() * sizeof(int));
		for (int i = 0; i < nnz(); i = i + 9) {
			double div = val(i)*val(i+4)*val(i+8) - val(i)*val(i+5)*val(i + 5) - val(i+1)*val(i+1)*val(i+8) + 2*val(i+1)*val(i+2)*val(i+5) - val(i+2)*val(i+2)*val(i+4);
			csr->val[i] = (val(i+4)*val(i+8) - val(i+5)*val(i+5)) / div;
			csr->val[i + 1] = (-val(i+1)*val(i+8) + val(i+2)*val(i+5)) / div;
			csr->val[i + 2] = (val(i+1)*val(i+5) - val(i+2)*val(i+4)) / div;
			csr->val[i + 3] = csr->val[i + 1];
			csr->val[i + 4] = (val(i)*val(i+8) - val(i+2)*val(i+2)) / div;
			csr->val[i + 5] = (-val(i)*val(i+5) + val(i+1)*val(i+2)) / div;
			csr->val[i + 6] = csr->val[i + 2];
			csr->val[i + 7] = csr->val[i + 5];
			csr->val[i + 8] = (val(i)*val(i+4) - val(i+1)*val(i+1)) / div;
		}
	}

	/*
	Split matrix to 3 submatrices defined by boundary (usefull for symmetric matrices)
	A = _sA(0..border_x-1, 0..border_y-1)
	B = _sA(0..border_x-1, border_y..end)
	C = _sA(border_x..end, border_y..end)
	*/
	void splitTo3Blocks(int border_x, int border_y, ScaledSparseMatrix *A, ScaledSparseMatrix *B, ScaledSparseMatrix *C) {
		int numA = 0, numB = 0, numC = 0;
		int N = nrows(), M = ncols();

		// count number of values in each matrix
		for (int i = 0; i < N; ++i) {	// rows
			for (int j = row(i); j < row(i + 1); ++j) {
				if (i < border_x && col(j) < border_y)
					numA++;
				if (i < border_x && col(j) >= border_y)
					numB++;
				if (i >= border_x && col(j) >= border_y)
					numC++;
			}
		}

		// create vectors for matrices 
		A->set_sA(std::make_shared<CRS>(border_x, border_y, numA));
		auto csrA = A->get_sA();
		B->set_sA(std::make_shared<CRS>(border_x, M - border_y, numB));
		auto csrB = B->get_sA();
		C->set_sA(std::make_shared<CRS>(N - border_x, M - border_y, numC));
		auto csrC = C->get_sA();

		// fill all vectors - rows, cols, and vals for each matrix
		bool changeA = false, changeB = false, changeC = false; 
		numA = 0, numB = 0, numC = 0;
		for (int i = 0; i < N; ++i) {	// rows
			for (int j = row(i); j < row(i + 1); ++j) {		// columns
				if (i < border_x && col(j) < border_y) {		// matrix A
					changeA = true;
					csrA->col[numA] = col(j);
					csrA->val[numA] = val(j);
					numA++;
				}
				if (i < border_x && col(j) >= border_y){		// matrix B
					changeB = true;
					csrB->col[numB] = col(j) - border_y;
					csrB->val[numB] = val(j);
					numB++;
				}
				if (i >= border_x && col(j) >= border_y){		// matrix C
					changeC = true;
					csrC->col[numC] = col(j) - border_y;
					csrC->val[numC] = val(j);
					numC++;
				}
			}
			if (changeA) {
				csrA->row[i + 1] = numA;
				changeA = false;
			}
			if (changeB) {
				csrB->row[i + 1] = numB;
				changeB = false;
			}
			if (changeC) {
				csrC->row[i + 1 - border_x] = numC;
				changeC = false;
			}
		}
	}

	/*
	Scale the matrix by diagonal matrix represented by A = cLR * sLR * sA   
	from left L or right R
	*/
	void scaleMat(int type, double **sLR, double *cLR) {
		double asum = 0;
		switch (type) {
		case LEFT:
			// auxiliary vectors
			(*sLR) = (double*)malloc(nrows()*sizeof(double));
			assert((*sLR) != NULL);
			memset(*sLR, 0, nrows()*sizeof(double));

			// encout scale 
			for (int i = 0; i < nrows(); i++) {
				for (int j = row(i); j < row(i + 1); j++)
					(*sLR)[i] += _sA->val[j] * _sA->val[j];	
				asum += (*sLR)[i];
				(*sLR)[i] = sqrt((*sLR)[i]);
			}

			// scale matrix
			for (int i = 0; i < nrows(); i++) {
				for (int j = row(i); j < row(i + 1); j++)
					_sA->val[j] *= 1 / (*sLR)[i];					
				(*sLR)[i] = (*sLR)[i] / asum;
			}
			*cLR = asum;

			break;
		case RIGHT:
			// auxiliary vectors
			(*sLR) = (double*)malloc(ncols()*sizeof(double));
			assert((*sLR) != NULL);
			memset(*sLR, 0, ncols()*sizeof(double));

			// encout scale 
			for (int i = 0; i < nnz(); i++)
				(*sLR)[col(i)] += _sA->val[i] * _sA->val[i];
			for (int i = 0; i < ncols(); i++) {
				asum += (*sLR)[i];
				(*sLR)[i] = sqrt((*sLR)[i]);
			}
			asum = sqrt(asum);

			// scale matrix
			for (int i = 0; i < nnz(); i++) 
				_sA->val[i] *= 1 / (*sLR)[col(i)];
			for (int i = 0; i < ncols(); i++)
				(*sLR)[i] = (*sLR)[i] / asum;
			*cLR = asum;
			break;
		}
	}

	/*
	Remove the diagonal scale of a matrix and apply the coefficient cLR at matrix coefficient
	*/
	void unscaleMat(int type, const double *sLR, const double cLR) {
		switch (type) {
		case LEFT:
			#pragma omp parallel for
			for (int i = 0; i < nrows(); i++) {
				for (int j = row(i); j < row(i + 1); j++)
					_sA->val[j] *= (cLR * sLR[i]);
			}
			break;
		case RIGHT:
			#pragma omp parallel for
			for (int i = 0; i < nnz(); ++i){
				_sA->val[i] *= (cLR * sLR[col(i)]);
			}
			break;
		}
	}

	// Transpose operator
	ScaledSparseMatrix trn() {
		SparseMatrix<double,RowMajor> eA = Map<SparseMatrix<double, RowMajor> >(nrows(), ncols(), nnz(), _sA->row, _sA->col, _sA->val, NULL);
		SparseMatrix<double,RowMajor> eAt = SparseMatrix<double,RowMajor>(eA.transpose());
		return ScaledSparseMatrix(ncols(), nrows(), eAt.outerIndexPtr(), eAt.innerIndexPtr(), eAt.valuePtr());
	}

	// Sparse multiplication using CPU
	void sparseMult(ScaledSparseMatrix& B) {
		int N = nrows();
		int M = B.ncols();
		
		ScaledSparseMatrix Bt( B.trn() );
		auto crs = std::make_shared<CRS>(N,M);
		std::vector<int> cols;
		std::vector<double> vals;

		for (int i = 0; i < N; ++i) {
			crs->row[i] = cols.size();

			for (int j = 0; j < M; ++j) {
				// multiply one sparse row of _sA with one sparse row of B._sA and save it to C
				int from1 = row(i);
				int to1 = row(i + 1);
				int from2 = Bt.row(j);
				int to2 = Bt.row(j + 1);
				int c1 = col(from1);
				int c2 = Bt.col(from2);
				double res = 0;

				// fast check
				if ((c1 > c2 && c1 > Bt.col(to2 - 1)) || (c1 < c2 && col(to1 - 1) < c2))
					continue;

				// regular process
				while (from1 != to1 && from2 != to2) {
					if (c1 > c2){
						c2 = Bt.col(++from2);

					}
					else {
						if (c1 < c2) {
							c1 = col(++from1);

						}
						else {							// col1 == col2
							res += sval(from1) * Bt.sval(from2);
							c1 = col(++from1);
							c2 = Bt.col(++from2);
						}
					}
				}

				if (res != 0.0) {
					cols.push_back(j);
					vals.push_back(res);
				}
			}
		}
		crs->row[N] = cols.size();

		// copy the solution to this
		crs->nnz = cols.size();
		crs->col = (int*)malloc( cols.size() * sizeof(int));
		assert(crs->col != NULL);
		crs->val = (double*)malloc( cols.size() * sizeof(double));
		assert(crs->val != NULL);
		memcpy(crs->col, cols.data(), cols.size() * sizeof(int));
		memcpy(crs->val, vals.data(), cols.size() * sizeof(double));
		set_sA(crs);

		// free the memmory
		cols.clear();
		vals.clear();
	}

	void sparseMultEigen(ScaledSparseMatrix& B) {
		auto sB = B.get_sA();
		SparseMatrix<double, RowMajor> eA = Map<SparseMatrix<double, RowMajor> >(nrows(), ncols(), nnz(), _sA->row, _sA->col, _sA->val, NULL);
		SparseMatrix<double, RowMajor> eB = Map<SparseMatrix<double, RowMajor> >(B.nrows(), B.ncols(), B.nnz(), sB->row, sB->col, sB->val, NULL);
		SparseMatrix<double, RowMajor> eC(eA * eB);
		set_sA(std::make_shared<CRS>(nrows(), B.ncols(), eC.outerIndexPtr(), eC.innerIndexPtr(), eC.valuePtr()));
	}

	void sparseMultMagma(ScaledSparseMatrix& B) {
		// init magma queue
		magma_queue_t queue;
		magma_int_t info = magma_init();
		magma_queue_create(info, &queue);
		magma_d_matrix hA = { Magma_CSR }, hB = { Magma_CSR }, hC = { Magma_CSR }, dA = { Magma_CSR }, dB = { Magma_CSR }, dC = { Magma_CSR };

		// init magma representation of the matrix
		magma_dcsrset(nrows(), ncols(), _sA->row, _sA->col, _sA->val, &hA, queue);
		magma_dcsrset(B.nrows(), B.ncols(), B.get_sA()->row, B.get_sA()->col, B.get_sA()->val, &hB, queue);

		// tranfer to GPU
		magma_dmtransfer(hA, &dA, Magma_CPU, Magma_DEV, queue);
		magma_dmtransfer(hB, &dB, Magma_CPU, Magma_DEV, queue);

		// sparse multiplication
		magma_dcuspmm(dA, dB, &dC, queue);

		// transfer back 
		magma_dmtransfer(dC, &hC, Magma_DEV, Magma_CPU, queue);
		
		// data to standard representation
		set_sA(std::make_shared<CRS>(nrows(), B.ncols(), hC.row, hC.col, hC.val));

		//free the memory
		magma_free_cpu(&hA);
		magma_free_cpu(&hB);
		magma_free_cpu(&hC);
		magma_dmfree(&dA, queue);
		magma_dmfree(&dB, queue);
		magma_dmfree(&dC, queue);
		magma_queue_destroy(queue);
		magma_finalize();
	}

	ScaledSparseMatrix& operator=(ScaledSparseMatrix& B) {
		_cA = B.c();
		_sA = B.get_sA();
		return *this;
	}

	ScaledSparseMatrix& operator=(ScaledSparseMatrix&& B) {
		_cA = B.c();
		_sA = std::move(B.get_sA());
		return *this;
	}
	
	// Multiplication operator
	ScaledSparseMatrix& operator*=(ScaledSparseMatrix& B) {
		// Scale matrices
		double cl = 1, cr = 1;
		double *sL = NULL;
		double *sR = NULL;
		scaleMat(LEFT, &sL, &cl);
		B.scaleMat(RIGHT, &sR, &cr);

		// Multiply two scaled matrices
		sparseMultEigen(B);
		//sparseMultMagma(B);

		// Unscale matrices and setup coefficients
		double c = 1 / (_sA->val[0] * sL[0] * sR[0]);
		unscaleMat(LEFT, sL, c);
		unscaleMat(RIGHT, sR, 1); 
		_cA = (_cA * B._cA * cl * cr) / c;
		B.set_cA(B.get_cA()*cr);
		B.unscaleMat(RIGHT, sR, 1);
		free(sL);
		free(sR);
		return *this;
	}

	// Addition operator
	ScaledSparseMatrix& operator+=(ScaledSparseMatrix& B) {
		addSSM(1, B);
		return *this;
	}

	// Subtraction operator
	ScaledSparseMatrix& operator-=(ScaledSparseMatrix& B) {
		addSSM(-1, B);
		return *this;
	}

	// Add or subtract
	void addSSM(int sign, ScaledSparseMatrix& B) {
		if (nrows() != B.nrows() || ncols() != B.ncols()) {
			std::cerr << "Addition or subtraction of matrices is not possible. Matrices have different number of rows or columns.";
			exit(1);
		}
		// count the number of values in output array
		int sumVal = 0;
		for (int i = 0; i < nrows(); ++i) {
			int col1 = row(i);
			int col2 = B.row(i);
			while (col1 < row(i + 1) || col2 < B.row(i + 1)) {
				// test if one of the indexes is at the and
				if (col1 == row(i + 1)) {
					sumVal += B.row(i + 1) - col2;
					break;
				}
				if (col2 == B.row(i + 1)) {
					sumVal += row(i + 1) - col1;
					break;
				}
				// we can compare the column indexes
				if (col(col1) == B.col(col2)) {
					col1++;
					col2++;
				}else if (col(col1) < B.col(col2)) {
					col1++;
				}else {
					col2++;
				}
				sumVal++;
			}
		}
		
		// create the arrays
		auto crs = std::make_shared<CRS>(nrows(), ncols(), sumVal);
		
		//fill the values
		sumVal = 0;
		for (int i = 0; i < nrows(); ++i) {
			int col1 = row(i);
			int col2 = B.row(i);
			while (col1 < row(i + 1) || col2 < B.row(i + 1)) {
				// test if one of the indexes is at the and
				if (col1 == row(i + 1)) {
					for (int j = col2; j < B.row(i + 1); ++j) {
						crs->col[sumVal] = B.col(j);
						crs->val[sumVal] = sign * B.val(j);
						sumVal++;
					}
					break;
				}
				if (col2 == B.row(i + 1)) {
					for (int j = col1; j < row(i + 1); ++j) {
						crs->col[sumVal] = col(j);
						crs->val[sumVal] = val(j);
						sumVal++;
					}
					sumVal += row(i + 1) - col1;
					break;
				}
				// we can compare the column indexes
				if (col(col1) == B.col(col2)) {
					crs->col[sumVal] = col(col1);
					crs->val[sumVal] = val(col1) + sign * B.val(col2);
					col1++;
					col2++;
				}
				else if (col(col1) < B.col(col2)) {
					crs->col[sumVal] = col(col1);
					crs->val[sumVal] = val(col1);
					col1++;
				}
				else {
					crs->col[sumVal] = B.col(col2);
					crs->val[sumVal] = sign * B.val(col2);
					col2++;
				}
				sumVal++;
			}
			crs->row[i + 1] = sumVal;
		}

		// fill the matrix
		_cA = 1;
		set_sA(crs);
	}

};

/*
Multiply two matricies with scale 
   - create copy of lhs and chage it by call *= 
*/
ScaledSparseMatrix operator*(ScaledSparseMatrix lhs, ScaledSparseMatrix& rhs){
	return lhs *= rhs;
}

/*
Add two matrices with scale
- create copy of lhs and chage it by call +=
*/
ScaledSparseMatrix operator+(ScaledSparseMatrix lhs, ScaledSparseMatrix& rhs) {
	return lhs += rhs;
}

/*
Subtract two matrices with scale
- create copy of lhs and chage it by call -=
*/
ScaledSparseMatrix operator-(ScaledSparseMatrix lhs, ScaledSparseMatrix& rhs) {
	return lhs -= rhs;
}