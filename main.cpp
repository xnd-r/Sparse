#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <random>
#include <omp.h>
#include <fstream>
#include <string>


void get_matrix(int n, std::vector<int>& val, std::vector<int>& col, std::vector<int>& row,
	bool generate, bool to_file, bool from_file){
	if (generate) {
		std::srand(unsigned(std::time(0)));
		std::default_random_engine generator;
		std::uniform_real_distribution<float> distribution(0.f, 1.f);
		std::vector<int> j_conunter;
		int* i_conunter = new int[n] {0};
    
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < i; ++j) {
				if (distribution(generator) < .35f) {
					col.push_back(j);
					i_conunter[i] += 1;
					val.push_back(std::rand() % 10 + 1);
				}
			}
			col.push_back(i);
			i_conunter[i] += 1;
			val.push_back(std::rand() % 10 + 1);
		}
		row.push_back(0);
		for (int i = 1; i < n + 1; ++i) {
			row.push_back(i_conunter[i - 1] + row[i - 1]);
		}
		delete[] i_conunter;
		if (to_file) {
			std::ofstream file;
			//std::string name = "data\\";
			//name.append(std::to_string(n).append("_").append(std::to_string(val.size())).append(".txt"));
			//file.open(name);
			file.open("data\\matrix.txt");
			file << n << "\n";
			file << val.size() << "\n";
			file << row.size() << "\n";
			for (std::vector<int>::iterator it = val.begin(); it != val.end(); ++it) {
				file << *it << " ";
			}
			file << "\n";
			for (std::vector<int>::iterator it = col.begin(); it != col.end(); ++it) {
				file << *it << " ";
			}
			file << "\n";
			for (std::vector<int>::iterator it = row.begin(); it != row.end(); ++it) {
				file << *it << " ";
			}
			file << "\n";
			file.close();
		}
	}
	if (from_file) {
		std::ifstream file ("data\\matrix.txt");
		if (file.is_open())
		{
			file >> n;
			int nz, rows;
			file >> nz >> rows;
			int value;
			for (int i = 0; i < nz; ++i) {
				file >> value;
				val.push_back(value);
			}
			for (int i = 0; i < nz; ++i) {
				file >> value;
				col.push_back(value);
			}
			for (int i = 0; i < rows; ++i) {
				file >> value;
				row.push_back(value);
			}
			file.close();
		}
	}
}

void print_matrix(int n, std::vector<int>& val, std::vector<int>& col, std::vector<int>& row){
	std::cout << "Matrix:\n";
	for (int ir = 0; ir < row.size()-1; ++ir){
		int start = row[ir];
		int finish = row[ir+1];
		int j = 0;
		while(start< finish){
			while(j < n) {
				if (j != col[start]){
					std::cout << "0" << "\t";
					j += 1;
				}
				else{
					std::cout << val[start] << "\t";
					start += 1;
					j += 1;
					break;
				}
			}
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	//std::cout << "\n\nValues:\n";
	//for (std::vector<int>::iterator it = val.begin(); it != val.end(); ++it) {
	//	std::cout << *it << "\t";
	//}
	//std::cout << "\n\nColumns:\n";
	//for (std::vector<int>::iterator it = col.begin(); it != col.end(); ++it) {
	//	std::cout << *it << "\t";
	//}
	//std::cout << "\n\nRow indexes:\n";
	//for (std::vector<int>::iterator it = row.begin(); it != row.end(); ++it) {
	//	std::cout << *it << "\t";
	//}
	//std::cout << "\n\nMatrix dimension:\n" << n;
	//std::cout << "\n\nNon-zero elements:\n" << val.size() << "\n\n";
}


void get_vector(int n, double* b) {
	std::srand(unsigned(std::time(0)));
	for (int i = 0; i < n; ++i) {
		b[i] = rand() % 10 + 1;
	}
}

void base_gauss_reverse(int n, std::vector<int>& val, std::vector<int>& col,
	std::vector<int>& row, double* x, double* b, int nthreads) {
	get_vector(n, b);
	double sum;
	double t1 = omp_get_wtime();
	x[0] = b[0] / val[row[0]];
	for (int ir = 1; ir < n; ++ir) {
		sum = 0.;
		#pragma omp parallel for num_threads(nthreads)
		for (int j = row[ir]; j < row[ir + 1]-1; ++j) {
			sum += val[j] * x[col[j]];
		}
		x[ir] = (b[ir] - sum) / val[row[ir + 1]-1];
	}
	std::cout << "Time: \t" << omp_get_wtime() - t1 << "\tthreads " << nthreads << "\n\n";

}

void get_y1(std::vector<double>& tr_part, int isn, int dim, double* b, double* x) {
	x[isn] = b[isn] / tr_part[0];
	double sum;
	int cnt;
	for (int i = 1; i < dim; ++i) {
		sum = 0.;
		cnt = i * (i + 1) / 2;
		for (int j = 0; j < i; ++j) {
			sum += tr_part[cnt+j] * x[isn + j];
		}
		x[isn+i] = (b[isn + i] - sum) / tr_part[cnt + i];
	}
}
void get_y2(int n, std::vector<int>& val, std::vector<int>& col,
	std::vector<int>& row, double* x, double* b, int isn, int dim) {
	for (int i = isn + dim; i < n; ++i) {
		//x[i] = 0.;
		for (int j = row[i]; j < row[i + 1]; j++) {
			if (col[j] >= isn and col[j] < isn + dim) {
				b[i] -= val[j] * x[col[j]];
			}
		}
		//b[i] = b[i] - x[i];
	}
}

void supernode_gauss_reverse(int n, std::vector<int>& val, std::vector<int>& col,
	std::vector<int>& row, double* x, double* b, int nthreads) {
	std::vector<int> supernodes;
	std::vector<double> tr_part;
	get_vector(n, b);
	for (int i = 0; i < n; ++i) {
		std::cout << b[i] << "\t";
	}
	std::cout << "\n";
	int isn = 0; int dim = 1;
	while (isn < n-1) {
		tr_part.clear();
		tr_part.push_back(val[row[isn + 1] - 1]);
		while (true) {
			for (int j = row[isn + dim]; j < row[isn + dim + 1]; ++j) {
				if (col[row[isn + 1] - 1] == col[j]) {
					if (row[isn + dim + 1] - j == dim + 1) {
						dim += 1;
						for (int el = j; el < j + dim; ++el) {
							tr_part.push_back(val[el]);
						}
						if (dim + isn == n) {
							break;
						}
						continue;
					}
					else {
						break;
					}
				}
			}
			get_y1(tr_part, isn, dim, b, x);
			get_y2(n, val, col, row, x, b, isn, dim);
			//std::cout << "\nMatrix:\n";
			//for (std::vector<double>::iterator it = tr_part.begin(); it != tr_part.end(); ++it) {
			//	std::cout << *it << "\t";
			//}
			//std::cout << "\n";


			//std::cout << "***" << std::endl;
			isn += dim;

			supernodes.push_back(dim);
			dim = 1;
			break;
		}
	}
	if (isn == n - 1){
		supernodes.push_back(dim);
		x[n - 1] = b[n - 1] / val[val.size() - 1];
		//std::cout << "\nMatrix:\n" << val[val.size()-1] << "\n";

		// lma for last elem
	}
	//for (std::vector<int>::iterator it = supernodes.begin(); it != supernodes.end(); ++it) {
	//	std::cout << *it << "\t";
	//}
}
int main(int* argc, char** argv){
	int n = 6; // atoi(argv[1]);
	int nthreads = 1; // atoi(argv[2]);
	std::vector<int> val; 
	std::vector<int> col; 
	std::vector<int> row; 
	double* b = new double[n * 2]{ 0. };
	double* x = b + n;

	get_matrix(n, val, col, row, true, false, false);
	print_matrix(n, val, col, row);
	//base_gauss_reverse(n, val, col, row, x, b, nthreads);



	std::cout << "\n";
	supernode_gauss_reverse(n, val, col, row, x, b, nthreads);
	std::cout.setf(std::ios::fixed);
	//std::cout.precision();
	for (int i = 0; i < n; ++i) {
		std::cout << x[i] << "\t=\t" << b[i] << "\n";
	}
	delete[] b;
	return 0;
}
