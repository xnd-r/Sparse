#include "core.h"

double get_supernodes(int n, int* row, int*col_index, int*& nodes, size_t& sn) {
	std::vector<int> supernodes;
	int isn = 0;
	int dim = 1;
	int f = 0;
	double t1 = omp_get_wtime();
	do {
		for (int j = 0; j < dim; ++j) {
			if ((row[col_index[isn + j] + dim - j] != isn + dim)
				|| ((col_index[isn + j] + dim - j) >= col_index[isn + dim])) {
				isn += dim;
				supernodes.push_back(dim);
				f = 1;
				break;
			}
		}
		if (f == 0) {
			dim++;
		}
		else {
			dim = 1;
			f = 0;
		}
	} while (isn + dim < n);
	supernodes.push_back(dim);
	double t2 = omp_get_wtime() - t1;

	sn = supernodes.size() + 1;
	nodes = new int[sn];

	double t3 = omp_get_wtime();
	nodes[0] = 0;
	for (int i = 1; i < sn; ++i) {
		nodes[i] = nodes[i - 1] + supernodes[i - 1];
	}
	double t4 = omp_get_wtime() - t3;

	return (t2 + t4);
}

void get_tr_part_lower(int isn, int dim, double* b, double* x, double* val, int* col) {
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

void get_rect_part_lower(int isn, int dim, double* b, double* x, double* val, int* row, int* col) {
	for (int i = isn; i < isn + dim; ++i) {
		for (int j = col[i] + dim - i + isn; j < col[i + 1]; ++j) {
			x[row[j]] -= x[i] * val[j];
		}
	}
}

double supernodal_lower(size_t sn, int* supernodes, double* b, double* x, double* val, int* row, int* col_index) {
	double t1 = omp_get_wtime();
	for (size_t i = 0; i < sn - 1; ++i) {
		get_tr_part_lower(supernodes[i], supernodes[i + 1] - supernodes[i], b, x, val, col_index);
		get_rect_part_lower(supernodes[i], supernodes[i + 1] - supernodes[i], b, x, val, row, col_index);
	}
	double t2 = omp_get_wtime();
	return (t2 - t1);
}

double supernodal_upper(size_t sn, int* supernodes, double* b, double* x, double* val, int* row, int* col_index) {
	double t1 = omp_get_wtime();
	for (size_t i = sn - 1; i > 0; --i) {
	
	}
	double t2 = omp_get_wtime();
	return (t2 - t1);
}


double supernodal_blas_lower(int n, size_t sn, int* supernodes, double* b, double* x, double* val, int* row, int* col_index) {
	int tr_dim;
	int shift = 0;
	int rec_dim;
	int lda, ldb;
	double* c = new double[n] {0.};
	double t1 = omp_get_wtime();
	for (int i = 0; i < sn - 1; ++i) {
		tr_dim = supernodes[i + 1] - supernodes[i];
		cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
			tr_dim, 1, 1., val + shift, tr_dim, x + supernodes[i], tr_dim);

		lda = (col_index[supernodes[i + 1]] - col_index[supernodes[i]]) / tr_dim;
		ldb = supernodes[i];
		rec_dim = lda - tr_dim;

		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
			rec_dim, 1, tr_dim, 1., val + shift + tr_dim, lda, x + supernodes[i], tr_dim, 1., c + supernodes[i], lda);

		shift += (col_index[supernodes[i + 1]] - col_index[supernodes[i]]);

		for (int k = 0; k < rec_dim; ++k) {
			x[row[col_index[supernodes[i]] + tr_dim + k]] -= c[k+ldb];
		}
	}
	double t2 = omp_get_wtime() - t1;
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
	phase = 331;
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
