#include <iostream>
#include <fstream>

#include "compute.h"
#include "ScaledSparseMatrix.h" 
#include "SparseMatrix/uncTest_SparseMatrix.h"

typedef ScaledSparseMatrix SSM;
#define EPS 1e-10

template<typename T>
bool checkArray(int N, T* arr1, T* arr2) {
	for (int i = 0; i < N; ++i) {
		if (arr1[i] != arr2[i])
			return false;
	}
	return true;
}

template<typename T>
bool checkArray(int N, double c1, T* arr1, double c2, T* arr2, double eps) {
	for (int i = 0; i < N; ++i) {
		double mean = std::abs((c1 * arr1[i] + c2 * arr2[i]) / 2);
		if (std::abs(c1 * arr1[i] - c2 * arr2[i]) >(mean * eps))
			return false;
	}
	return true;
}

bool loadJacobian(std::ifstream &file, int *nrows, int *ncols, int *nnz, int **rows, int **cols, double **vals) {
	if (file.is_open()){
		file >> (*nrows) >> (*ncols) >> (*nnz);
		(*rows) = (int*)malloc(((*nrows) + 1) * sizeof(int));
		(*cols) = (int*)malloc((*nnz) * sizeof(int));
		(*vals) = (double*)malloc((*nnz) * sizeof(double));
		for (int i = 0; i <= (*nrows); ++i)
			file >> (*rows)[i];
		for (int i = 0; i < (*nnz); ++i)
			file >> (*cols)[i];
		for (int i = 0; i < (*nnz); ++i)
			file >> (*vals)[i];
		file.close();
		return true;
	}
	return false;
}

bool Test_Multipl1() {
	// INPUT
	int nrows_in, ncols_in, nnz_in, *r_in = NULL, *c_in = NULL;
	double *v_in = NULL;
	std::ifstream file_in("d:/School/Research/Reconstruction3D/Software/Uncertainty/te_inversion/build/Debug/in/multipl_in_01.txt", std::ios_base::in);
	if (!loadJacobian(file_in, &nrows_in, &ncols_in, &nnz_in, &r_in, &c_in, &v_in))
		return false;
	SSM A(nrows_in, ncols_in, r_in, c_in, v_in);

	// OUTPUT - GT
	int nrows_out, ncols_out, nnz_out, *r_out = NULL, *c_out = NULL;
	double *v_out = NULL;
	std::ifstream file_out("d:/School/Research/Reconstruction3D/Software/Uncertainty/te_inversion/build/Debug/in/multipl_out_01.txt", std::ios_base::in);
	if (!loadJacobian(file_out, &nrows_out, &ncols_out, &nnz_out, &r_out, &c_out, &v_out))
		return false;
	SSM AA(nrows_out, ncols_out, r_out, c_out, v_out);

	// PROCESS
	SSM At(A.trn());
	SSM TEST_AA(At * A);

	// CHECK RESULT
	if (checkArray(AA.nrows() + 1, TEST_AA.get_sA()->row, AA.get_sA()->row) &&
		checkArray(AA.ncols(), TEST_AA.get_sA()->col, AA.get_sA()->col) &&
		checkArray(AA.nnz(), TEST_AA.get_cA() ,TEST_AA.get_sA()->val, AA.get_cA(), AA.get_sA()->val, EPS))
		return true;
	else
		return false;
}


bool Test_Multipl2() {
	// INPUT
	int nrows_in, ncols_in, nnz_in, *r_in = NULL, *c_in = NULL;
	double *v_in = NULL;
	std::ifstream file_in("d:/School/Research/Reconstruction3D/Software/Uncertainty/te_inversion/build/Debug/in/multipl_in_02.txt", std::ios_base::in);
	if (!loadJacobian(file_in, &nrows_in, &ncols_in, &nnz_in, &r_in, &c_in, &v_in))
		return false;
	SSM A(nrows_in, ncols_in, r_in, c_in, v_in);

	// OUTPUT - GT
	int nrows_out, ncols_out, nnz_out, *r_out = NULL, *c_out = NULL;
	double *v_out = NULL;
	std::ifstream file_out("d:/School/Research/Reconstruction3D/Software/Uncertainty/te_inversion/build/Debug/in/multipl_out_02.txt", std::ios_base::in);
	if (!loadJacobian(file_out, &nrows_out, &ncols_out, &nnz_out, &r_out, &c_out, &v_out))
		return false;
	SSM AA(nrows_out, ncols_out, r_out, c_out, v_out);

	// PROCESS
	SSM At(A.trn());
	SSM TEST_AA(At * A);

	// CHECK RESULT
	if (checkArray(AA.nrows() + 1, TEST_AA.get_sA()->row, AA.get_sA()->row) &&
		checkArray(AA.ncols(), TEST_AA.get_sA()->col, AA.get_sA()->col) &&
		checkArray(AA.nnz(), TEST_AA.get_cA(), TEST_AA.get_sA()->val, AA.get_cA(), AA.get_sA()->val, EPS))
		return true;
	else
		return false;
}


int main(int argc, char* argv[]) {
	std::cout << "Test_Multipl1: " << Test_Multipl1() << "\n";
	std::cout << "Test_Multipl2: " << Test_Multipl2() << "\n";
	exit(0);
}