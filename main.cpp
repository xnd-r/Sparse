#include <iostream>
#include <vector>
#include <fstream>
//#include <ctime>
//#include <random>

void read_ccs(int& n, double*& val, int*& row, int*& col_index) {
	std::ifstream file("mtx.mtx");
	while (file.peek() == '%') {
		file.ignore(2048, '\n');
	}

	int nz;
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
}

void get_vector(int n, double* b) {
	//std::srand(unsigned(std::time(0)));
	for (int i = 0; i < n; ++i) {
		//b[i] = rand() % 10 + 1;
		b[i] = i + 1;
		b[n + i] = b[i];
	}
}

void get_supernodes(int n, int* row, int*col_index, int*& nodes, int& sn) {
	std::vector<int> supernodes;
	supernodes.push_back(0);
	int isn = 0;
	int dim = 1;
	int st, fn;
	do {
		for (int j = 0; j < dim; ++j) {
			if (row[col_index[isn + j] + dim - j] != isn + dim) {
				isn += dim;
				supernodes.push_back(dim);
				dim = 0;
				break;
			}
		}
		dim += 1;
	} while (isn < n);

	nodes = new int[supernodes.size()];
	nodes[0] = 0;
	sn = supernodes.size();
	for (int i = 1; i < sn; ++i) {
		nodes[i] = nodes[i - 1] + supernodes[i];
	}
}

void get_tr_part(int isn, int dim, double* b, double* x, double* val, int* col) {
	x[isn] = x[isn] / val[col[isn]];
	double sum;
	int cnt1, cnt2;
	for (int k = 1; k < dim; ++k) {
		sum = 0.;
		cnt1 = 0, cnt2 = k;
		for (int cnt = 0; cnt < k; ++cnt) {
			sum += val[col[isn + cnt1] + cnt2] * x[isn + cnt];
			cnt1++;
			cnt2--;
		}
		x[isn + k] = (x[isn + k] - sum) / val[col[isn + k]];
	}
}

void get_rect_part(int isn, int dim, double* b, double* x, double* val, int* row, int* col) {
	for (int i = isn; i < isn + dim; ++i) {
		for (int j = col[i] + dim - i + isn; j < col[i + 1]; ++j) {
			x[row[j]] -= x[i] * val[j];
		}
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


int main(int* argc, char** argv){
	int n; int sn;
	int* supernodes = nullptr;
	double* val = nullptr;
	int* row = nullptr, *col_index = nullptr;

	read_ccs(n, val, row, col_index);
	double* b = new double[n * 2]{ 0. };
	double* x = b + n;
	get_supernodes(n, row, col_index, supernodes, sn);
	get_vector(n, b);
	
	// supernodal algorithm
	for (int i = 0; i < sn - 1; ++i) {
		get_tr_part(supernodes[i], supernodes[i + 1] - supernodes[i], b, x, val, col_index);
		get_rect_part(supernodes[i], supernodes[i + 1] - supernodes[i], b, x, val, row, col_index);
	}

	//print_matrix(n, val, row, col_index);
	//check_result(n, val, col, row, x, b);
	print_result(n, x, b);
	
	delete[] b;
	delete[] val;
	delete[] row;
	delete[] col_index;
	return 0;
}
