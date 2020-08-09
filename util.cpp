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
	//for (int i = 0; i <= n; ++i) {
	//	std::cout << col_index[i] << " ";
	//}
	//std::cout << "\n";
}

void get_vector(int n, double* b) {
	#pragma omp parallel for num_threads(4)
	for (int i = 0; i < n; ++i) {
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

double get_random_L_bin(int n, int& nz, size_t& sn, double density) {
	int start_ind = 0, dim = 0, pattern_elems = 0, random_index = 0;
	srand(static_cast <unsigned> (47));
	//int max_node_size = n / 10 / log(n);
	int max_node_size = 3;


	std::vector<int> col_index_v;
	std::vector<int> supernodes_v;
	std::vector<int> pattern;
	col_index_v.push_back(0);

	supernodes_v.push_back(0);
	FILE *fp;
	auto si = sizeof(int);
	int tmp_ind = 0;
	if ((fp = fopen("BIN_MAT.bin", "wb")) == NULL) {
		printf("Cannot open file.");
		return 1;
	}
	double t1 = omp_get_wtime();

	while (start_ind < n) {
		dim = rand() % max_node_size + 1;
		if (dim >= n - start_ind) {
			dim = n - start_ind;
		}

		pattern_elems = (int)(density * (n - start_ind - dim));
		for (int el = 0; el < pattern_elems; ++el) {
			pattern.emplace_back(rand() % (n - start_ind - dim) + start_ind + dim);
		}
		std::sort(pattern.begin(), pattern.end());
		pattern.erase(std::unique(pattern.begin(), pattern.end()), pattern.end());

		for (int i = 0; i < dim; ++i) {
			for (int j = i; j < dim; ++j) {
				tmp_ind = start_ind + j;
				//row_v.push_back(start_ind + j);
				fwrite(&tmp_ind, sizeof(int), 1, fp);
			}
			for (int k = 0; k < pattern.size(); ++k) {
			//	//val_v.push_back(static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
			//	//row_v.push_back(pattern[k]);
				fwrite(&pattern[k], sizeof(int), 1, fp);

			}
			nz += dim - i + pattern.size();
			//nz += dim - i;
			col_index_v.push_back(nz);
		}
		start_ind += dim;
		pattern.clear();
		supernodes_v.emplace_back(dim + supernodes_v.back());
	}
	//std::cout << "SIZE: " << col_index_v.size() << '\n';
	for (int i = 0; i < col_index_v.size(); ++i) {
		fwrite(&col_index_v[i], si, 1, fp);
	}
	for (int i = 0; i < supernodes_v.size(); ++i) {
		fwrite(&supernodes_v[i], si, 1, fp);
	}
	//fwrite(&col_index_v, sizeof(std::vector<int>::value_type), col_index_v.size(), fp);
	//fwrite(&supernodes_v, sizeof(std::vector<int>::value_type), supernodes_v.size(), fp);
	sn = supernodes_v.size();
	fclose(fp);
	double t2 = omp_get_wtime();

	return t2 - t1;
}

double get_random_L(int n, int& nz, size_t& sn, double density, double*& val, int*& row, int*& col_index, int*& nodes) {
	int start_ind = 0, dim = 0, pattern_elems = 0, random_index = 0;
	srand(static_cast <unsigned> (47));

	std::vector<double> val_v;
	std::vector<int> row_v;
	std::vector<int> col_index_v;
	std::vector<int> supernodes_v;
	std::vector<int> pattern;
	int max_node_size = 100;

	double t1 = omp_get_wtime();
	col_index_v.push_back(0.);
	supernodes_v.push_back(0);

	while (start_ind < n) {
		dim = rand() % max_node_size + 1;
		if (dim >= n - start_ind) {
			dim = n - start_ind;
		}

		pattern_elems = (int)(density * (n - start_ind - dim));
		for (int el = 0; el < pattern_elems; ++el) {
			pattern.emplace_back(rand() % (n - start_ind - dim) + start_ind + dim);
		}
		std::sort(pattern.begin(), pattern.end());
		pattern.erase(std::unique(pattern.begin(), pattern.end()), pattern.end());

		for (int i = 0; i < dim; ++i) {
			for (int j = i; j < dim; ++j) {
				val_v.push_back(static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
				row_v.push_back(start_ind + j);
			}
			for (int k = 0; k < pattern.size(); ++k) {
				val_v.push_back(static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
				row_v.push_back(pattern[k]);
			}
			nz += dim - i + pattern.size();
			col_index_v.push_back(nz);
		}
		start_ind += dim;
		pattern.clear();
		supernodes_v.emplace_back(dim);
	}
	double t2 = omp_get_wtime() - t1;
	//double t3 = omp_get_wtime();
	//std::ofstream file;
	//file.open("random_factor_L.txt");
	//file << n << "\n";
	//file << nz << "\n";
	//file << supernodes_v.size() << "\n";
	//for (int i = 0; i < nz; ++i) {
	//	file << val_v[i] << " ";
	//}
	//file << "\n";
	//for (int i = 0; i < nz; ++i) {
	//	file << row_v[i] << " ";
	//}
	//file << "\n";
	//for (int i = 0; i <= n; ++i) {
	//	file << col_index_v[i] << " ";
	//}
	//file << "\n";
	//int sum = 0;
	//for (int i = 1; i < supernodes_v.size(); ++i) {
	//	sum += supernodes_v[i];
	//	file << sum << " ";
	//}
	//file << "\n";
	//file.close();
	//std::cout << "\n" << omp_get_wtime() - t3 << "\n";
	//val = new double[nz];
	//row = new int[nz];
	//col_index = new int[n + 1];
	//for (int i = 0; i < nz; ++i) {
	//	val[i] = val_v[i];
	//	row[i] = row_v[i];
	//}
	//sn = supernodes_v.size() + 1;
	//nodes = new int[sn];

	//double t3 = omp_get_wtime();
	//nodes[0] = 0;
	//for (int i = 1; i < sn; ++i) {
	//	nodes[i] = nodes[i - 1] + nodes[i - 1];
	//}
	//double t4 = omp_get_wtime() - t3;
	//delete[] unic_indices;
	return t2;
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
	//std::cout << "\ninfo: \t" << info << "\n";
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

double check_result(int n, double* x1, double* x2) {
	double sum = 0., norm = 0.;
	//std::cout << "\nL_2(sol - pardiso_sol): \n";
	for (int i = 0; i < n; ++i) {
		sum += pow(x1[i] - x2[i], 2);
		norm += x1[i] * x1[i];
	}
	//std::cout.setf(std::ios::fixed);
	//std::cout.precision(32);
	std::cout << sqrt(sum) / sqrt(norm);
	return sum;
}
