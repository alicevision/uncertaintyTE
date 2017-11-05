#pragma once
#include "compute.h"


// Magma sparse matrix multiplied by scalar
// magma_d_matrix operator*(double lhs, magma_d_matrix rhs);

struct CRS {
	int nrows = 0;
	int ncols = 0;
	int nnz = 0;
	int *row = NULL;
	int *col = NULL;
	double *val = NULL;

	CRS();

	CRS(const int nr, const int nc);

	CRS(const int nr, const int nc, const int nnz1);

	CRS(const int nr, const int nc, const int *rows, const int *cols, const double *values);

	~CRS();

	CRS& operator=(CRS&& B);
};

// Matrix scaled by scalar structure
class ScaledSparseMatrix {
private:
	double _cA = 1;
	std::shared_ptr<CRS> _sA;

public:
	ScaledSparseMatrix();

	ScaledSparseMatrix(const ScaledSparseMatrix& o);

	ScaledSparseMatrix(ScaledSparseMatrix&& o);

	ScaledSparseMatrix(std::shared_ptr<CRS> &A);

	ScaledSparseMatrix(const int nrows, const int ncols, const int nnz);

	ScaledSparseMatrix(const int nrows, const int ncols, const int *rows, const int *cols, const double *values);


	std::shared_ptr<CRS> get_sA();
	double get_cA();
	void set_sA(std::shared_ptr<CRS> &A);
	void set_sA(std::shared_ptr<CRS> &&A);
	void set_cA(double c);
	int nrows() const;
	int ncols() const;
	int nnz() const;
	int col(const int i) const;
	int row(const int i) const;
	double c() const;
	double sval(int i) const;
	double val(int i) const;
	double val(int i, int j) const;
	// Print submatrix up to then rows / cols 
	void print() const;
	void printAll() const;

	void printBlock2Matlab(char* name, int row_from, int col_from, int row_to, int col_to);
	void printBlock2Matlab2(char* name, int row_from, int col_from, int row_to, int col_to);
	void printBlock2Matlab3(char* name, int row_from, int col_from, int row_to, int col_to);

	void inv3x3blockSymmDiag(ScaledSparseMatrix *ssm);

	/*
	Split matrix to 3 submatrices defined by boundary (usefull for symmetric matrices)
	A = _sA(0..border_x-1, 0..border_y-1)
	B = _sA(0..border_x-1, border_y..end)
	C = _sA(border_x..end, border_y..end)
	*/
	void splitTo3Blocks(int border_x, int border_y, ScaledSparseMatrix *A, ScaledSparseMatrix *B, ScaledSparseMatrix *C);

	/*
	Scale the matrix by diagonal matrix represented by A = cLR * sLR * sA
	from left L or right R
	*/
	void scaleMat(int type, double **sLR, double *cLR);

	/*
	Remove the diagonal scale of a matrix and apply the coefficient cLR at matrix coefficient
	*/
	void unscaleMat(int type, const double *sLR, const double cLR);

	// Transpose operator
	ScaledSparseMatrix trn();

	// Sparse multiplication using CPU
	void sparseMult(ScaledSparseMatrix& B);

	void sparseMultEigen(ScaledSparseMatrix& B);

	//void sparseMultMagma(ScaledSparseMatrix& B);

	ScaledSparseMatrix& operator=(ScaledSparseMatrix& B);

	ScaledSparseMatrix& operator=(ScaledSparseMatrix&& B);

	// Multiplication operator
	ScaledSparseMatrix& operator*=(ScaledSparseMatrix& B);

	// Addition operator
	ScaledSparseMatrix& operator+=(ScaledSparseMatrix& B);

	// Subtraction operator
	ScaledSparseMatrix& operator-=(ScaledSparseMatrix& B);

	// Add or subtract
	void addSSM(int sign, ScaledSparseMatrix& B);
};

/*
Multiply two matricies with scale
   - create copy of lhs and chage it by call *=
*/
ScaledSparseMatrix operator*(ScaledSparseMatrix lhs, ScaledSparseMatrix& rhs);

/*
Add two matrices with scale
- create copy of lhs and chage it by call +=
*/
ScaledSparseMatrix operator+(ScaledSparseMatrix lhs, ScaledSparseMatrix& rhs);

/*
Subtract two matrices with scale
- create copy of lhs and chage it by call -=
*/
ScaledSparseMatrix operator-(ScaledSparseMatrix lhs, ScaledSparseMatrix& rhs);