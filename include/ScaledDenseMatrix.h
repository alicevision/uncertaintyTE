#pragma once
#include "compute.h"
#include "ScaledSparseMatrix.h" 
#include <float.h>


class ScaledDenseMatrix {
private:
	int _nr = 0;
	int _nc = 0;
	double _cA = 1;
	double *_sA = NULL;

public:
	ScaledDenseMatrix();

	ScaledDenseMatrix(int nrows, int ncols);

	ScaledDenseMatrix(int nrows, int ncols, double* values);

	ScaledDenseMatrix(const ScaledDenseMatrix &B);

	ScaledDenseMatrix(ScaledSparseMatrix &A);

	ScaledDenseMatrix(ScaledSparseMatrix &&A);

	~ScaledDenseMatrix();

	int nrows();
	int ncols();
	double c() const;
	void set_c(double c);
	double* getMatPtr();
	double sval(int i) const;
	double sval(int row, int col) const;
	double val(int i) const;
	double val(int row, int col) const;
	void set(const int i, const double val);
	void set(const int i, const int j, const double val);

	void print();

	void printBlock2Matlab(std::string name, int row_from, int col_from, int row_to, int col_to);
	void printBlock2Matlab2(std::string name, int row_from, int col_from, int row_to, int col_to);
	void printBlock2Matlab3(std::string name, int row_from, int col_from, int row_to, int col_to);

	/*
	Scale the matrix by diagonal matrix represented by A = cLR * sLR * sA
	from left L or right R
	*/
	void scaleMat(int type, double **sLR, double *cLR);

	// Find maximal coefficient 
	double absMax();

	// Lower triangle is coppied to the upper triangle
	void symmL2U();

	// General inverse: "Z" -> iZ
	void inv();

	// Deep copy function
	void copy(const ScaledDenseMatrix& other);

	// Add or subtract
	void addSDM(int sign, ScaledDenseMatrix& B);

	// Deep copy assignment operator
	ScaledDenseMatrix& operator=(const ScaledDenseMatrix &other);

	// Multiplication operator
	ScaledDenseMatrix& operator*=(ScaledDenseMatrix& B);

	// Multiplication with scalar
	ScaledDenseMatrix& operator*=(double &B);

	// Addition operator
	ScaledDenseMatrix& operator+=(ScaledDenseMatrix& B);

	// Subtraction operator
	ScaledDenseMatrix& operator-=(ScaledDenseMatrix& B);

};

/*
Multiply two matricies with scale
- create copy of lhs and chage it by call *=
*/
ScaledDenseMatrix operator*(ScaledDenseMatrix lhs, ScaledDenseMatrix& rhs);

/*
Multiply by scalar
- create copy of lhs and chage it by call *=
*/
ScaledDenseMatrix operator*(double lhs, ScaledDenseMatrix& rhs);
ScaledDenseMatrix operator*(ScaledDenseMatrix& lhs, double rhs);

/*
Add two matrices with scale
- create copy of lhs and chage it by call +=
*/
ScaledDenseMatrix operator+(ScaledDenseMatrix lhs, ScaledDenseMatrix& rhs);

/*
Subtract two matrices with scale
- create copy of lhs and chage it by call -=
*/
ScaledDenseMatrix operator-(ScaledDenseMatrix lhs, ScaledDenseMatrix& rhs);