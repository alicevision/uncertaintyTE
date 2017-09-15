#include "compute.h"
#include "ScaledDenseMatrix.h" 

typedef ScaledDenseMatrix SDM;
#define EPS_INV1 1e-7
#define EPS_INV2 1e-2

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


bool loadMatrix(std::ifstream &file, int *nrows, int *ncols, double **vals) {
	if (file.is_open()) {
		file >> (*nrows) >> (*ncols);
		int n = (*nrows) * (*ncols);
		(*vals) = (double*)malloc(n * sizeof(double));
		for (int i = 0; i < n; ++i)
			file >> (*vals)[i];
		file.close();
		return true;
	}
	return false;
}

bool Test_inversion1() {
	// INPUT
	int nrows_in, ncols_in;
	double *v_in = NULL;
	std::ifstream file_in("d:/School/Research/Reconstruction3D/Software/Uncertainty/te_inversion/build/Debug/in/inv_in_01.txt", std::ios_base::in);
	if (!loadMatrix(file_in, &nrows_in, &ncols_in, &v_in))
		return false;
	SDM A(nrows_in, ncols_in, v_in);

	// OUTPUT - GT
	int nrows_out, ncols_out;
	double *v_out = NULL;
	std::ifstream file_out("d:/School/Research/Reconstruction3D/Software/Uncertainty/te_inversion/build/Debug/in/inv_out_01.txt", std::ios_base::in);
	if (!loadMatrix(file_out, &nrows_out, &ncols_out, &v_out))
		return false;
	SDM iA(nrows_out, ncols_out, v_out);

	// PROCESS
	A.inv();

	// CHECK RESULT
	if (checkArray(nrows_out * ncols_out, iA.c(), iA.getMatPtr(), A.c(), A.getMatPtr(), EPS_INV1))
		return true;
	else
		return false;
}

bool Test_inversion2() {
	// INPUT
	int nrows_in, ncols_in;
	double *v_in = NULL;
	std::ifstream file_in("d:/School/Research/Reconstruction3D/Software/Uncertainty/te_inversion/build/Debug/in/inv_in_02.txt", std::ios_base::in);
	if (!loadMatrix(file_in, &nrows_in, &ncols_in, &v_in))
		return false;
	SDM A(nrows_in, ncols_in, v_in);

	// OUTPUT - GT
	int nrows_out, ncols_out;
	double *v_out = NULL;
	std::ifstream file_out("d:/School/Research/Reconstruction3D/Software/Uncertainty/te_inversion/build/Debug/in/inv_out_02.txt", std::ios_base::in);
	if (!loadMatrix(file_out, &nrows_out, &ncols_out, &v_out))
		return false;
	SDM iA(nrows_out, ncols_out, v_out);

	// PROCESS
	A.inv();

	// CHECK RESULT
	if (checkArray(nrows_out * ncols_out, iA.c(), iA.getMatPtr(), A.c(), A.getMatPtr(), EPS_INV2))
		return true;
	else
		return false;
}

int main(int argc, char* argv[]) {
	std::cout << "Test_inversion1: " << Test_inversion1() << "\n";
	std::cout << "Test_inversion2: " << Test_inversion2() << "\n";
	exit(0);
}