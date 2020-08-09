#include "util.h"
#include "core.h"
#include "string.h"


int main(int argc, char** argv) {
	int n = 0, nz = 0;
	size_t sn = 0;
	double* val = nullptr, *step_val = nullptr, *b = nullptr;
	int *nodes = nullptr, *row = nullptr, *col_index = nullptr;
	int *step_row = nullptr, *step_col_index = nullptr;
	const char* filename = "bcsstk01.mtx";
	const char* algo = "base";
	//const char* filename = argv[1];
	//const char* algo = argv[2];
	bool pardiso_compare = false;
	bool dtrsv_compare = true;
	double time = 0.;
	//n = 8;
	//get_random_L_bin(n, nz, sn, 0.45);

	//FILE *fp;
	//row = new int[nz];
	//val = new double[nz];
	//col_index = new int[n + 1];
	//nodes = new int[sn];
	//auto si = sizeof(int);
	//if ((fp = fopen("BIN_MAT.bin", "rb")) == NULL) {
	//	printf("Cannot open file.\n");
	//	exit(1);
	//}
	//fread(row, sizeof(int), nz, fp);
	//fread(col_index, sizeof(int), n + 1, fp);
	//fread(nodes, sizeof(int), sn, fp);

	//for (int i = 0; i < nz; ++i) {
	//	val[i] = 1 + rand() % 10;
	//}

	//tmp(n, nz, sn);

	read_ccs(filename, n, nz, val, row, col_index);
	b = new double[n * 2]{ 0. };
	double* x = b + n;
	get_vector(n, b);
	double* y = new double[n];
	double* dense = new double[n*n]{ 0. };
	get_factor(n, val, row, col_index, dense);
	read_factor(n, nz, val, row, col_index);

	if (strcmp(algo, "base") == 0) {
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
	else if (strcmp(algo, "custom") == 0) {
		get_supernodes(n, nz, val, row, col_index, nodes, sn, step_val, step_row, step_col_index);
		time += supernodal_lower(sn, nodes, x, val, row, col_index);
		#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			y[i] = x[i];
		}
		time += supernodal_upper(sn, nodes, y, val, row, col_index);
	}
	else if (strcmp(algo, "blas") == 0) {
		get_supernodes(n, nz, val, row, col_index, nodes, sn, step_val, step_row, step_col_index);
		time += supernodal_blas_lower(n, sn, nodes, x, step_val, step_row, step_col_index);
		#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			y[i] = x[i];
		}
		time += supernodal_blas_upper(n, nz, sn, nodes, y, step_val, step_row, step_col_index);
	}
	else if (strcmp(algo, "dtrsv") == 0) {
		get_vector(n, b);
		sparse_matrix_t A;
		mkl_sparse_d_create_csc(&A, SPARSE_INDEX_BASE_ZERO, n, n, col_index, col_index + 1, row, val);
		matrix_descr descr;
		descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
		descr.mode = SPARSE_FILL_MODE_LOWER;
		descr.diag = SPARSE_DIAG_NON_UNIT;
		double t1 = omp_get_wtime();
		mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, 1., A, descr, b, x);
		time += omp_get_wtime() - t1;
		std::cout << "Dtrsv time: " << time << "\n";
	}
	else {
		std::cout << "\nUnknown type of algo " << algo << ". Exit";
		return -1;
	}

	if (pardiso_compare) {
		delete[] val;
		delete[] row;
		delete[] col_index;
		read_ccs(filename, n, nz, val, row, col_index);
		get_vector(n, b);
		pardiso_solution(n, val, row, col_index, b, x);
		std::cout << time << "\n";
		check_result(n, x, y);
	}
	if (dtrsv_compare) {
		get_vector(n, b);
		sparse_matrix_t A;
		mkl_sparse_d_create_csc(&A, SPARSE_INDEX_BASE_ZERO, n, n, col_index, col_index + 1, row, val);
		matrix_descr descr;
		descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
		descr.mode = SPARSE_FILL_MODE_LOWER;
		descr.diag = SPARSE_DIAG_NON_UNIT;
		mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, 1., A, descr, b, x);
		std::cout << "\n" << time << "\n";
		check_result(n, x, y);
	}

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
