#pragma once
#include <fstream>
#include <iostream>
#include <vector>
#include <mkl.h>
#include <omp.h>
#include <cstdlib>
#include <algorithm>

void read_ccs(const char* name, int& n, int& nz, double*& val, int*& row, int*& col_index);
void get_vector(int n, double* b);
void print_result(int n, double* x, double* b);
void sparse_to_dense(int n, double* val, int* row, int* col, double* dense);
void dense_to_sparse(int n, double* dense, double*& val_s, int*& row_s, int*& col_ind_s);
void get_factor(int n, double*& val, int*& row, int*& col_index, double*& dense);
void read_factor(int& n, int& nz, double*& val, int*& row, int*& col_ind);
void transpose(int n, int nz, double*& val, int*& row, int*& col_index, double*& val_t, int*& row_t, int*& col_index_t);
double check_result(int n, double* x1, double* x2);
double get_random_L(int n, int& nz, size_t& sn, double density, double*& val, int*& row, int*& col_index, int*& supernodes);
double get_random_L_bin(int n, int& nz, size_t& sn, double density);
