#include "util.h"

void read_ccs(const char* name, int& n, int& nz, double*& val, int*& row, int*& col_index) {
	std::ifstream file(name);
	while (file.peek() == '%') {
		file.ignore(2048, '\n');
	}

	file >> n >> n >> nz;
	col_index = new int[n + 1];
	col_index[0] = 0;
	val = new double[nz];
	row = new int[nz];
	int col_i, row_i, cnt = 0, cnt1 = 1;

	for (int i = 0; i < nz; ++i) {
		file >> row_i >> col_i >> val[i];
		row[i] = row_i - 1;
		if (col_i == cnt1) {
			cnt++;
			continue;
		}
		else {
			col_index[cnt1] = cnt;
			cnt++;
			cnt1++;
		}
	}
	col_index[n] = nz;
	for (int i = 0; i <= n; ++i) {
		std::cout << col_index[i] << " ";
	}
	std::cout << "\n";
}

void get_vector(int n, double* b) {
	//std::srand(unsigned(std::time(0)));
	#pragma omp parallel for num_threads(4)
	for (int i = 0; i < n; ++i) {
		//b[i] = rand() % 10 + 1;
		b[i] = i + 1.;
		b[n + i] = b[i];
	}
}

void print_result(int n, double* x, double* b) {
	std::cout.setf(std::ios::fixed);
	std::cout.precision();
	std::cout << "\n   x: \t\t\t   b:\n";
	for (int i = 0; i < n; ++i) {
		std::cout << x[i] << "\t=\t" << b[i] << "\n";
	}
}

void sparse_to_dense(int n, double* val, int* row, int* col, double* dense) {
	for (int i = 0; i < n; ++i) {
		for (int ind = col[i]; ind < col[i + 1]; ++ind) {
			dense[row[ind] * n + i] = val[ind];
		}
	}
}

void dense_to_sparse(int n, double* dense, double*& val_s, int*& row_s, int*& col_ind_s) {
	std::vector<int> row;
	std::vector<int> col;
	std::vector<double> val;
	int* j_counter = new int[n] {0};

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (dense[j*n + i] != 0.) {
				row.push_back(j);
				j_counter[i] += 1;
				val.push_back(dense[j*n + i]);
			}
		}
	}
	size_t nz = val.size();
	val_s = new double[nz];
	row_s = new int[nz];
	for (int i = 0; i < nz; ++i) {
		val_s[i] = val[i];
		row_s[i] = row[i];
		//std::cout << val_s[i] << "\t";
	}
	col_ind_s = new int[n + 1];
	col_ind_s[0] = 0;
	for (int i = 1; i < n + 1; ++i) {
		col_ind_s[i] = j_counter[i - 1] + col_ind_s[i - 1];
		//std::cout << col_ind_s[i] << "\t";
	}

	std::ofstream file;
	file.open("matrix.txt");
	file << n << "\n";
	file << nz << "\n";
	file << n+1 << "\n";
	for (int i = 0; i < nz; ++i) {
		file << val_s[i] << " ";
	}
	file << "\n";
	for (int i = 0; i < nz; ++i) {
		file << row_s[i] << " ";
	}
	file << "\n";
	for (int i = 0; i <= n; ++i) {
		file << col_ind_s[i] << " ";
	}
	file << "\n";
	file.close();
}

void get_factor(int n, double*& val, int*& row, int*& col_index, double*& dense) {
	sparse_to_dense(n, val, row, col_index, dense);
	//std::cout.setf(std::ios::fixed);
	//std::cout.precision();
	//for (int i = 0; i < n; ++i) {
	//	for (int j = 0; j <= i; ++j) {
	//		std::cout << dense[i * n + j] << "\t";
	//	}
	//	printf("\n");
	//}
	//printf("\n");

	char uplo = 'U';
	int info;
	dpotrf(&uplo, &n, dense, &n, &info);
	std::cout << "\ninfo: \t" << info << "\n";
	//for (int i = 0; i < n; ++i) {
	//	for (int j = 0; j <= i; ++j) {
	//		std::cout << dense[i * n + j] << "\t";
	//	}
	//	std::cout << "\n";
	//}
	//std::cout << "\n";

	delete[] val;
	delete[] row;
	delete[] col_index;
	val = nullptr;
	row = nullptr, col_index = nullptr;
	dense_to_sparse(n, dense, val, row, col_index);
	//printf("\n");
}

void read_factor(int& n, int& nz, double*& val, int*& row, int*& col_ind) {
	std::ifstream file("matrix.txt");
	if (file.is_open())
	{
		file >> n;
		int rows;
		file >> nz >> rows;
		val = new double[nz];
		row = new int[nz];
		col_ind = new int[n + 1];
		for (int i = 0; i < nz; ++i) {
			file >> val[i];
		}
		for (int i = 0; i < nz; ++i) {
			file >> row[i];
		}
		for (int i = 0; i < rows; ++i) {
			file >> col_ind[i];
		}
		file.close();
	}
}

void transpose(int n, int nz, double*& val, int*& row, int*& col_index,
	double*& val_t, int*& row_t, int*& col_index_t) {
	for (int i = 0; i < nz; ++i)
	{
		col_index_t[row[i] + 1]++;
	}
	int S = 0, tmp;
	for (int i = 1; i <= n; ++i) {
		tmp = col_index_t[i];
		col_index_t[i] = S;
		S += tmp;
	}
	for (int i = 0; i < n; i++)
	{
		int j1 = col_index[i]; int j2 = col_index[i + 1];
		int Col = i;
		for (int j = j1; j < j2; j++)
		{
			double V = val[j];
			int RIndex = row[j];
			int IIndex = col_index_t[RIndex + 1];
			val_t[IIndex] = V;
			row_t[IIndex] = Col;
			col_index_t[RIndex + 1]++;
		}
	}
	//std::cout << "\n";
	//for (int i = 0; i < nz; ++i) {
	//	std::cout << val_t[i] << " ";
	//}
	//std::cout << "\n";

	//for (int i = 0; i < nz; ++i) {
	//	std::cout << row_t[i] << " ";
	//}
	//std::cout << "\n";

	//for (int i = 0; i <= n; ++i) {
	//	std::cout << col_index_t[i] << " ";
	//}
	//std::cout << "\n";
}

void check_result(int n, double* val, int* col, int* row, double* x, double* b) {
	double resudual = 0.;
	double sum;
	std::cout << "\n1-Norm of Resudial vector: \n";
	for (int ir = 0; ir < n; ++ir) {
		sum = 0.;
		for (int j = row[ir]; j < row[ir + 1]; ++j) {
			sum += val[j] * x[col[j]];
		}
		resudual += abs(b[ir] - sum);
		//std::cout << b[ir] - sum << "\t";
	}
	std::cout.setf(std::ios::fixed);
	std::cout.precision(32);
	std::cout << resudual << "\n";
}
