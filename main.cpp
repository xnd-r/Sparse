#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <random>

void get_matrix(int n, std::vector<int>& val, std::vector<int>& col, std::vector<int>& row){
	std::srand(unsigned(std::time(0)));
	std::default_random_engine generator;
	std::uniform_real_distribution<float> distribution(0., 1.);
	std::vector<int> j_conunter;
	int* i_conunter = new int[n]{0};
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < i; ++j){
			if (distribution(generator) < 0.05f){
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
	for(int i = 1; i < n + 1; ++i){
		row.push_back(i_conunter[i-1] + row[i-1]);
	}
	delete[] i_conunter;
}

void print_matrix(int n, std::vector<int>& val, std::vector<int>& col, std::vector<int>& row){
	std::cout << "Generated Matrix:\n";
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
	std::cout << "\n\nValues:\n";
	for (std::vector<int>::iterator it = val.begin(); it != val.end(); ++it) {
		std::cout << *it << "\t";
	}
	std::cout << "\n\nColumns:\n";
	for (std::vector<int>::iterator it = col.begin(); it != col.end(); ++it) {
		std::cout << *it << "\t";
	}
	std::cout << "\n\nRow indexes:\n";
	for (std::vector<int>::iterator it = row.begin(); it != row.end(); ++it) {
		std::cout << *it << "\t";
	}
	std::cout << "\n\nMatrix dimension:\n" << n;
	std::cout << "\n\nNon-zero elements:\n" << val.size();
}

int main(int* argc, char** argv){
	int n = atoi(argv[1]);
	std::vector<int> val; 
	std::vector<int> col; 
	std::vector<int> row; 
	get_matrix(n, val, col, row);
	print_matrix(n, val, col, row);

    
	return 0;
}

