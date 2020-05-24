#include "util.h"
#include "core.h"


int main(int argc, char** argv) {
	int n = 0, nz = 0;
	size_t sn;
	int *nodes = nullptr, *row = nullptr, *col_index = nullptr;
	double* val = nullptr, *b = nullptr, *dense = nullptr;
	const char* filename = "mtx_nonpad.mtx";
	double time = 0.;

	read_ccs(filename, n, nz, val, row, col_index);

	//dense = new double[n*n]{ 0. };
	//get_factor(n, val, row, col_index, dense);

	//read_factor(n, nz, val, row, col_index);
	b = new double[n * 2]{ 0. };
	double* x = b + n;
	dense = new double[n*n]{ 0. };
	get_vector(n, b);
	get_supernodes(n, row, col_index, nodes, sn);
	time = supernodal_lower(sn, nodes, b, x, val, row, col_index);
	//time = supernodal_blas_lower(n, sn, nodes, b, x, val, row, col_index);

	//read_ccs(filename, n, nz, val, row, col_index);
	//b = new double[n * 2]{ 0. };
	//double* x = b + n;
	//dense = new double[1];
	//get_vector(n, b);
	//time = pardiso_solution(n, val, row, col_index, b, x);
	//read_factor(n, nz, val, row, col_index);

	double* val_t = new double[nz] {0};
	int* row_t = new int[nz] {0};
	int* col_index_t = new int[n + 1]{ 0 };
	transpose(n, nz, val, row, col_index, val_t, row_t, col_index_t);
	print_result(n, x, b);
	check_result(n, val_t, row_t, col_index_t, x, b);
	std::cout << "\nTIME " << time;

	delete[] b;
	delete[] val;
	delete[] row;
	delete[] col_index;
	delete[] val_t;
	delete[] row_t;
	delete[] col_index_t;
	delete[] dense;
	delete[] nodes;
	return 0;
}
