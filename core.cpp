#include "core.h"

double base_gauss_lower(int n, double* val, int* row, int* col_index, double* x, double* b) {
	double sum;
	double t1 = omp_get_wtime();
	x[0] = b[0] / val[row[0]];
	for (int ir = 1; ir < n; ++ir) {
		sum = 0.;
		for (int j = col_index[ir]; j < col_index[ir + 1] - 1; ++j) {
			sum += val[j] * x[row[j]];
		}
		x[ir] = (b[ir] - sum) / val[col_index[ir + 1] - 1];
	}
	return omp_get_wtime() - t1;
}

double base_gauss_upper(int n, double* val, int* row, int* col_index, double* x, double* b) {
	double sum;
	double t1 = omp_get_wtime();
	x[n - 1] = b[n - 1] / val[col_index[n - 1]];
	for (int ir = n - 2; ir >= 0; --ir) {
		sum = 0.;
		for (int j = col_index[ir] + 1; j < col_index[ir + 1]; ++j) {
			sum += val[j] * x[row[j]];
		}
		x[ir] = (b[ir] - sum) / val[col_index[ir]];
	}
	return omp_get_wtime() - t1;
}

void get_current_pattern(int i, int n, double* val, int* row, int* col_index, std::vector<int>& supernodes, int& isn, std::vector<int>& cur_pattern) {
	int amo;
	for (; i < n; ++i) {
		amo = col_index[isn + 1] - col_index[isn];
		if (amo == 1) {
			supernodes.push_back(1);
			isn += 1;
			continue;
		}
		else {
			for (int j = 0; j < amo; ++j) {
				if (val[col_index[isn] + j] != 0) {
					cur_pattern.push_back(row[col_index[isn] + j]);
				}
			}
			break;
		}
	}
	if (cur_pattern.size() > 0) {
		cur_pattern.erase(cur_pattern.begin());
	}
}

double get_supernodes(int n, int nz, double* val, int* row, int*col_index, int*& nodes, size_t& sn,
	double*& step_val, int*& step_row, int*& step_col_index) {
	int amo;
	int isn = 0; // start index of supernode
	int dim = 1;
	int cnt = 0;
	std::vector<int> supernodes;
	std::vector<int> cur_pattern;
	std::vector<int> next_pattern;
	step_col_index = new int[n + 1];
	std::vector<double> step_val_v(val, val + nz);
	std::vector<int> step_row_v(row, row + nz);

	for (int i = 0; i <= n; ++i) {
		step_col_index[i] = col_index[i];
	}
	int flag = 0;
	int next_elem;

	double t1 = omp_get_wtime();
	get_current_pattern(0, n, val, row, col_index, supernodes, isn, cur_pattern);
	for (int i = isn + 1; i < n; ++i) {
		amo = col_index[i + 1] - col_index[i] - 1;
		for (int j = 0; j < amo; ++j) {
			next_pattern.push_back(row[col_index[i] + 1 + j]);
		}
		if (cur_pattern.size() > 0) {
			next_elem = cur_pattern[0];
			if (next_elem == row[col_index[i]]) {
				cur_pattern.erase(cur_pattern.begin());
			}
			else {
				cur_pattern.clear();
				next_pattern.clear();
				supernodes.push_back(dim);
				dim = 1;
				get_current_pattern(i, n, val, row, col_index, supernodes, i, cur_pattern);
				continue;
			}
		}
		if (cur_pattern == next_pattern) {
			dim++;
			for (int d = 0; d < dim - 1; ++d) {
				step_val_v.insert(step_val_v.begin() + col_index[i] + cnt, 0.0);
				step_row_v.insert(step_row_v.begin() + col_index[i] + cnt, row[col_index[i]] - dim + 1 + d);
				for (int k = i + 1; k <= n; ++k) {
					step_col_index[k]++;
				}
				cnt++;
			}
			next_pattern.clear();
			if (cur_pattern.size() == 0) {
				supernodes.push_back(dim);
				dim = 1;
				i++;
				get_current_pattern(i, n, val, row, col_index, supernodes, i, cur_pattern);
				continue;
			}
			cur_pattern.erase(cur_pattern.begin());
		}
		else {
			cur_pattern.clear();
			next_pattern.clear();
			supernodes.push_back(dim);
			dim = 1;
			get_current_pattern(i, n, val, row, col_index, supernodes, i, cur_pattern);
		}
	}
	double t2 = omp_get_wtime() - t1;

	int nz_new = step_col_index[n];
	step_val = new double[nz_new];
	step_row = new int[nz_new];
	for (int i = 0; i < nz_new; ++i) {
		step_val[i] = step_val_v[i];
		step_row[i] = step_row_v[i];
	}
	sn = supernodes.size() + 1;
	nodes = new int[sn];

	double t3 = omp_get_wtime();
	nodes[0] = 0;
	for (int i = 1; i < sn; ++i) {
		nodes[i] = nodes[i - 1] + supernodes[i - 1];
	}
	double t4 = omp_get_wtime() - t3;
	return t2 + t4;
}

void get_tr_part_lower(int isn, int dim, double* x, double* val, int* col) {
	x[isn] = x[isn] / val[col[isn]];
	double sum;
	int cnt2;
	for (int k = 1; k < dim; ++k) {
		sum = 0.;
		cnt2 = k;
		for (int cnt = 0; cnt < k; ++cnt) {
			sum += val[col[isn + cnt] + cnt2] * x[isn + cnt];
			cnt2--;
		}
		x[isn + k] = (x[isn + k] - sum) / val[col[isn + k]];
	}
}

void get_rect_part_lower(int isn, int dim, double* x, double* val, int* row, int* col) {
	for (int i = isn; i < isn + dim; ++i) {
		for (int j = col[i] + dim - i + isn; j < col[i + 1]; ++j) {
			x[row[j]] -= x[i] * val[j];
		}
	}
}

double supernodal_lower(size_t sn, int* supernodes, double* x, double* val, int* row, int* col_index) {
	double t1 = omp_get_wtime();
	for (size_t i = 0; i < sn - 1; ++i) {
		get_tr_part_lower(supernodes[i], supernodes[i + 1] - supernodes[i], x, val, col_index);
		get_rect_part_lower(supernodes[i], supernodes[i + 1] - supernodes[i], x, val, row, col_index);
	}
	return omp_get_wtime() - t1;
}


void get_tr_part_upper(int isn, int dim, double* x, double* val, int* row, int* col) {
	x[isn + dim - 1] = x[isn + dim - 1] / val[col[isn + dim - 1]];
	double sum;
	int cnt = 0;
	for (int j = isn + dim - 2; j >= isn; --j) {
		sum = 0.;
		for (int k = col[j] + 1; k <= col[j] + 1 + cnt; ++k) {
			sum += val[k] * x[row[k]];
		}
		x[j] = (x[j] - sum) / val[col[j]];
		cnt++;
	}
}

void get_rect_part_upper(int isn, int dim, double* x, double* val, int* row, int* col) {
	double sum;
	int cnt = 0;
	double t1 = omp_get_wtime();
	for (int i = isn; i < isn + dim; ++i) {
		sum = 0.;
		for (int j = col[i] + dim - cnt; j < col[i + 1]; ++j) {
			sum += val[j] * x[row[j]];
		}
		x[i] -= sum;
		cnt++;
	}
}

double supernodal_upper(size_t sn, int* supernodes, double* x, double* val, int* row, int* col_index) {
	double t1 = omp_get_wtime();
	for (int i = (int)sn-2; i >= 0; --i) {
		get_rect_part_upper(supernodes[i], supernodes[i + 1] - supernodes[i], x, val, row, col_index);
		get_tr_part_upper(supernodes[i], supernodes[i + 1] - supernodes[i], x, val, row, col_index);

	}
	return omp_get_wtime() - t1;
}


double supernodal_blas_lower(int n, size_t sn, int* supernodes, double* x, double* val, int* row, int* col_index) {
	int tr_dim;
	int shift = 0;
	int rec_dim;
	int lda, ldb;
	double* c = new double[n] {0.};
	double t1 = omp_get_wtime();
	for (int i = 0; i < sn - 1; ++i) {
		tr_dim = supernodes[i + 1] - supernodes[i];
		lda = (col_index[supernodes[i + 1]] - col_index[supernodes[i]]) / tr_dim; // TODO: refactor
		cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
			tr_dim, 1, 1., val + shift, lda, x + supernodes[i], tr_dim);

		ldb = supernodes[i];
		rec_dim = col_index[supernodes[i] + 1] - col_index[supernodes[i]] - tr_dim;

		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
			rec_dim, 1, tr_dim, 1., val + shift + tr_dim, lda, x + supernodes[i], tr_dim, 1., c, lda);

		shift += (col_index[supernodes[i + 1]] - col_index[supernodes[i]]); // TODO: refactor

		for (int k = 0; k < rec_dim; ++k) {
			x[row[col_index[supernodes[i]] + tr_dim + k]] -= c[k];
			c[k] = 0;
		}
	}
	double t2 = omp_get_wtime() - t1;
	delete[] c;
	return t2;
}


double supernodal_blas_upper(int n, int nz, size_t sn, int* supernodes, double* x, double* val, int* row, int* col_index) {
	int tr_dim;
	int shift = 0;
	int rec_dim;
	int lda, ldb = 1;
	double* B = new double[n] {0.};
	double* c = new double[n] {0.};
	double t1 = omp_get_wtime();
	for (int i = (int)sn-2; i >= 0; --i) {
		tr_dim = supernodes[i + 1] - supernodes[i];
		shift = col_index[supernodes[i]];
		lda = (col_index[supernodes[i + 1]] - col_index[supernodes[i]]) / tr_dim; // TODO: refactor
		rec_dim = col_index[supernodes[i] + 1] - col_index[supernodes[i]] - tr_dim;

		for (int i = 0; i < rec_dim; ++i) {
			B[i] = x[row[shift + tr_dim + i]];
		}

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			tr_dim, 1, rec_dim, 1., val + shift + tr_dim, lda, B, 1, 1., c, 1);

		for (int k = 0; k < tr_dim; ++k) {
			x[supernodes[i] + k] -= c[k];
			c[k] = 0.;
		}
		cblas_dtrsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
			tr_dim, 1, 1., val + shift, lda, x + supernodes[i], ldb);
	}
	double t2 = omp_get_wtime() - t1;
	delete[] B;
	delete[] c;
	return t2;
}



double pardiso_solution(int n, double* val, int* row, int*col_index, double* b, double* x) {
	// Reference: https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-fortran/top/sparse-solver-routines/intel-mkl-pardiso-parallel-direct-sparse-solver-interface/pardiso.html
	MKL_INT mtype = 2;       /* Real and symmetric positive definite */
	MKL_INT nrhs = 1;     /* Number of right hand sides. */
	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures s*/
	void *pt[64];
	/* Pardiso control parameters. */
	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error, msglvl;
	/* Auxiliary variables. */
	MKL_INT i;
	double ddum;          /* Double dummy */
	MKL_INT idum;         /* Integer dummy. */
/* -------------------------------------*/
/* .. Setup Pardiso control parameters. */
/* -------------------------------------*/
	for (i = 0; i < 64; i++)
	{
		iparm[i] = 0;
	}
	iparm[0] = 1;         /* No solver default */
	iparm[1] = 2;         /* Fill-in reordering from METIS */
	iparm[3] = 0;         /* No iterative-direct algorithm */
	iparm[4] = 0;         /* No user fill-in reducing permutation */
	iparm[5] = 0;         /* Write solution into x */
	iparm[7] = 2;         /* Max numbers of iterative refinement steps */
	iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
	iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm[13] = 0;        /* Output: Number of perturbed pivots */
	iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;       /* Output: Mflops for LU factorization */
	iparm[19] = 0;        /* Output: Numbers of CG Iterations */
	iparm[34] = 1;        /* PARDISO use C-style indexing for ia and ja arrays */
	maxfct = 1;           /* Maximum number of numerical factorizations. */
	mnum = 1;             /* Which factorization to use. */
	msglvl = 1;           /* Print statistical information in file */
	error = 0;            /* Initialize error flag */
/* ----------------------------------------------------------------*/
/* .. Initialize the internal solver memory pointer. This is only  */
/*   necessary for the FIRST call of the PARDISO solver.           */
/* ----------------------------------------------------------------*/
	for (i = 0; i < 64; i++)
	{
		pt[i] = 0;
	}
	/* --------------------------------------------------------------------*/
	/* .. Reordering and Symbolic Factorization. This step also allocates  */
	/*    all memory that is necessary for the factorization.              */
	/* --------------------------------------------------------------------*/
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, val, col_index, row, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
	/* ----------------------------*/
	/* .. Numerical factorization. */
	/* ----------------------------*/
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, val, col_index, row, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	printf("\nFactorization completed ... ");

	/* -----------------------------------------------*/
	/* .. Back substitution and iterative refinement. */
	/* -----------------------------------------------*/
	// Reference: https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-fortran/top/sparse-solver-routines/intel-mkl-pardiso-parallel-direct-sparse-solver-interface.html#intel-mkl-pardiso-parallel-direct-sparse-solver-interface_GUID-70E6283E-35F8-4935-84DB-E4FBDD5B5A12
	phase = 33;
	iparm[7] = 1;         /* Max numbers of iterative refinement steps. */
	double t1 = omp_get_wtime();
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, val, col_index, row, &idum, &nrhs, iparm, &msglvl, b, x, &error);
	double t2 = omp_get_wtime();
	if (error != 0)
	{
		printf("\nERROR during solution: %d", error);
		exit(3);
	}
	printf("\nSolve completed ... ");
	//printf("\nThe solution of the system is: ");
	//for (i = 0; i < n; i++)
	//{
	//	printf("\n x [%d] = % f", i, x[i]);
	//}
	//printf("\n");
	/* --------------------------------------*/
	/* .. Termination and release of memory. */
	/* --------------------------------------*/
	phase = -1;           /* Release internal memory. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, &ddum, col_index, row, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	return t2 - t1;
}
