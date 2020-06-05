#include "util.h"
#include "core.h"

int main(int argc, char** argv) {
	int n = 0, nz = 0;
	size_t sn = 0;
	double* val = nullptr, *step_val = nullptr, *b = nullptr;
	int *nodes = nullptr, *row = nullptr, *col_index = nullptr;
	int *step_row = nullptr, *step_col_index = nullptr;
	const char* filename = "bcsstk01.mtx";
	const char* algo = "blas";
	double time = 0.;

	read_ccs(filename, n, nz, val, row, col_index);
	b = new double[n * 2]{ 0. };
	double* x = b + n;
	get_vector(n, b);
	double* y = new double[n];
	double* dense = new double[n*n]{ 0. };
	get_factor(n, val, row, col_index, dense);
	read_factor(n, nz, val, row, col_index);

	if (algo == "base") {
		double* val_t = new double[nz] {0.};
		int* row_t = new int[nz] {0};
		int* col_index_t = new int[n + 1]{0};
		transpose(n, nz, val, row, col_index, val_t, row_t, col_index_t);
		time += base_gauss_lower(n, val_t, row_t, col_index_t, x, b);
		#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			y[i] = x[i];
		}
		time += base_gauss_upper(n, val, row, col_index, y, x);
		delete[] val_t;
		delete[] row_t;
		delete[] col_index_t;
	}
	else if (algo == "custom") {
		get_supernodes(n, nz, val, row, col_index, nodes, sn, step_val, step_row, step_col_index);
		time += supernodal_lower(sn, nodes, x, val, row, col_index);
		#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			y[i] = x[i];
		}
		time += supernodal_upper(sn, nodes, y, val, row, col_index);
	}
	else if (algo == "blas") {
		get_supernodes(n, nz, val, row, col_index, nodes, sn, step_val, step_row, step_col_index);
		time += supernodal_blas_lower(n, sn, nodes, x, step_val, step_row, step_col_index);
		#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			y[i] = x[i];
		}
		time += supernodal_blas_upper(n, nz, sn, nodes, y, step_val, step_row, step_col_index);
	}
	else {
		std::cout << "\nUnknown type of algo. Exit";
		return -1;
	}
	delete[] val;
	delete[] row;
	delete[] col_index;
	read_ccs(filename, n, nz, val, row, col_index);
	get_vector(n, b);
	time = pardiso_solution(n, val, row, col_index, b, x);
	check_result(n, x, y);
	std::cout << "\nTime: " << time;

	delete[] b;
	delete[] y;
	delete[] val;
	delete[] row;
	delete[] col_index;
	delete[] step_col_index;
	delete[] step_row;
	delete[] step_val;
	delete[] dense;
	delete[] nodes;
	return 0;
}
