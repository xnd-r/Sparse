#pragma once
#include <mkl.h>
#include <omp.h>
#include <vector>
#include <numeric>

double base_gauss_lower(int n, double* val, int* row, int* col_index, double* x, double* b);
double base_gauss_upper(int n, double* val, int* row, int* col_index, double* x, double* b);
double get_supernodes(int n, int* row, int*col_index, int*& nodes, size_t& sn);
double get_supernodes_new(int n, double* val, int* row, int*col_index, int*& nodes, size_t& sn);
void get_tr_part_lower(int isn, int dim, double* b, double* x, double* val, int* col);
void get_rect_part_lower(int isn, int dim, double* b, double* x, double* val, int* row, int* col);
double supernodal_lower(size_t sn, int* supernodes, double* b, double* x, double* val, int* row, int* col_index);
double supernodal_blas_lower(int n, size_t sn, int* supernodes, double* b, double* x, double* val, int* row, int* col_index);
double supernodal_blas_upper(int n, int nz, size_t sn, int* supernodes, double* b, double* x, double* val, int* row, int* col_index);
double pardiso_solution(int n, double* val, int* row, int*col_index, double* b, double* x);
