#include "common.h"

#include <iostream>
#include <iomanip>
#include <random>
#include <cstring>
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> uni_rd(0.0, 1.0); // call uni_rd to get a random number in [0,1)

//#include <vector>
using namespace std;

//
void rd_wf(const my_int& dim, double* wf)
{
	for (my_int i = 0; i < dim; i++)
	{
		wf[i] = (uni_rd(gen) - 0.5);
	}
	// normalize
	cblas_dscal(dim, 1.0 / cblas_dnrm2(dim, wf, 1), wf, 1);
}

void rd_wf(const my_int& dim, my_cmplx* wf)
{
	my_cmplx II(0, 1);
	for (my_int i = 0; i < dim; i++)
	{
		wf[i] = (uni_rd(gen) - 0.5) + (uni_rd(gen) - 0.5) * II;
	}
	// normalize
	cblas_zdscal(dim, 1.0 / cblas_dznrm2(dim, wf, 1), wf, 1);
}

void rd_wf_product(const my_int& dim, double* wf)
{
	for (my_int i = 0; i < dim; i++)
	{
		wf[i] = 0;
	}
	my_int idx = uni_rd(gen) * dim;
	wf[idx] = 1.0;
	//cout << "norm of random product states"  << cblas_dnrm2(dim, wf, 1) << endl;
}

void uni_wf(const my_int& dim, double* wf)
{
	for (my_int i = 0; i < dim; i++)
	{
		wf[i] = 1;
	}
	// normalize
	cblas_dscal(dim, 1.0 / cblas_dnrm2(dim, wf, 1), wf, 1);
}

// bitwise operations
my_int Bits_Count1(const my_int& b)
{
	my_int a = b;
	my_int cnt = 0;
	while (a != 0)
	{
		++cnt;
		a &= (a - 1);
	}
	return cnt;
}

my_int Bits_CycleRight(const my_int& s, const my_int& n)
{
	//return (s >> 1) + ((s & 1) << (n-1));
	return (s >> 1) | (s << (n - 1));
}

my_int Bits_CycleLeft(const my_int& s, const my_int& n)
{
	//return ((s << 1) | (s >> (n - 1))) - ((s & (1 << (n - 1))) << 1);
	return ((s << 1) | (s >> (n - 1))) - ((s >> (n - 1)) << n);
}

my_int Bits_Reflect(const my_int& s, const my_int& n)
{
	my_int aux = 0;
	for (my_int i = 0; i < n; i++)
	{
		aux |= (aux | (((s >> i) & 1) << (n - 1 - i)));
	}
	return aux;
}

my_int Bits_Invert(const my_int& s, const my_int& n)
{
	my_int aux = (1 << n) - 1;
	return s ^ aux;
}

my_int Bits_HammingDis(const my_int& a, const my_int& b)
{
	return Bits_Count1(a ^ b);
}

void Bits_Decomposition(const my_int& s, const my_int& n, const my_int& size_A, my_int* sites_A, my_int& a, my_int& b)
{
	// 
	int size_B = n - size_A;
	my_int aux_B = 0;
	//
	a = 0;
	b = 0;
	for (my_int ia = 0; ia < size_A; ia++)
	{
		a += (((s >> sites_A[ia]) & 1) << ia);
		aux_B += (1 << sites_A[ia]);
	}
	//
	int ib = 0;
	aux_B = (~aux_B);
	for (my_int i = 0; i < n; i++)
	{
		if ((aux_B >> i) & 1)
		{
			b += ((s >> i) & 1) << ib;
			ib++;
		}
	}
}

void Bits_Print(const my_int& a, const my_int& n)
{
	cout << "  ";
	for (my_int ix = 0; ix < n; ix++) cout << ((a >> (n - ix - 1)) & 1);
	cout << "  ";
}

my_int nchoosek(const my_int& n, const my_int& _k)
{
	if (_k > n) return 0;
	my_int k = _k < (n - _k) ? _k : (n - _k);
	if (0 == k) return 1;
	if (1 == k) return n;
	double aux = 1.0;
	for (my_int i = 0; i < k; i++)
	{
		aux *= double(n - i) / double(k - i);
	}
	return (my_int)(aux + 1e-2);
}

void Debug_Check_Conj_dense(double* mat, my_int dim)
{
	double aux = 0;
	for (my_int i = 0; i < dim; i++)
	{
		for (my_int j = i + 1; j < dim; j++)
		{
			aux += abs(mat[i * dim + j] - mat[j * dim + i]);
		}
	}
	cout << "abs(H - H^T) = " << aux << endl;
}

// compute fermionic sign for 1 site operator  
double Cal_chargesign(const my_int& numk, const my_int& i) {
	my_int numk1 = numk >> (i + 1);
	my_int cnt = 0;
	while (numk1 != 0)
	{
		++cnt;
		numk1 &= (numk1 - 1);
	}
	double sign = cnt % 2 == 0 ? 1.0 : -1.0;
	return sign;
}

// compute fermionic sign for 2 sites
double Cal_chargesign(const my_int& numk, const my_int& i, const my_int& j) {
	double sign;
	if (j < i) {
		cout << "Error! In function Cal_chargesign(const my_int &numk, const my_int &i, const my_int &j) there must be i <= j";
		exit(-100);
	}
	else if (j <= i + 1) {
		sign = 1;
	}
	else {
		my_int cnt = 0;
		for (my_int site = i + 1; site < j; site++) {
			cnt += (numk >> site) & 1;
		}
		sign = cnt % 2 == 0 ? 1.0 : -1.0;
	}
	return sign;
}

double num_conj(double a)
{
	return a;
}
my_cmplx num_conj(my_cmplx a) {
	return std::conj(a);
}

double my_vec_dot(const my_int& dim, double* x, const my_int incx, double* y, const my_int incy)
{
	return cblas_ddot(dim, x, incx, y, incy);
}
my_cmplx my_vec_dot(const my_int& dim, my_cmplx* x, const my_int incx, my_cmplx* y, const my_int incy)
{
	my_cmplx aux;
	cblas_zdotc_sub(dim, x, incx, y, incy, &aux);
	return aux;
}

// parse input
void GetParaFromInput_int(const char* fname, const char* string_match, my_int& para) {
	FILE* f_in = fopen(fname, "r");
	//	ifstream f_in(fname);
	char testchar[40], line[80];
	sprintf(testchar, string_match);
	my_int len = strlen(testchar);
	while (strncmp(line, testchar, len) != 0) {
		fgets(line, 200, f_in);
		char* p = strtok(line, "=");
		stringstream ss;
		p = strtok(NULL, "=");
		ss << p;
		ss >> para;
	}
	//delete[]f_in;
}

void GetParaFromInput_real(const char* fname, const char* string_match, double& para) {
	FILE* f_in = fopen(fname, "r");
	//	ifstream f_in(fname);
	char testchar[40], line[80];
	sprintf(testchar, string_match);
	my_int len = strlen(testchar);
	while (strncmp(line, testchar, len) != 0) {
		fgets(line, 200, f_in);
		char* p = strtok(line, "=");
		stringstream ss;
		p = strtok(NULL, "=");
		ss << p;
		ss >> para;
	}
}

void GetParaFromInput_char(const char* fname, const char* string_match, char& para) {
	FILE* f_in = fopen(fname, "r");
	//	ifstream f_in(fname);
	char testchar[40], line[80];
	sprintf(testchar, string_match);
	my_int len = strlen(testchar);
	while (strncmp(line, testchar, len) != 0) {
		fgets(line, 200, f_in);
		char* p = strtok(line, "=");
		stringstream ss;
		p = strtok(NULL, "=");
		ss << p;
		ss >> para;
	}
	fclose(f_in);
	f_in = NULL;
}

void Vec_fwrite_double(const char* fname, double* data, const my_int& dsize)
{
	FILE* f_out;
	f_out = fopen(fname, "wb");
	fwrite(data, sizeof(double), dsize, f_out);
	fclose(f_out);
}

void Vec_fread_double(const char* fname, double* data, const my_int& dsize)
{
	FILE* f_in;
	f_in = fopen(fname, "rb");
	fread(data, sizeof(double), dsize, f_in);
	fclose(f_in);
}


// ===============================================================================
// Matrix evd (use mkl lapacke function) 
void DenseMatrixEigenSolver(my_int matrix_layout, char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w)
{
	int eigtype = 1;	// 0, dsyevd; 1, dsyev; 2 dsyevr;
	if (0 == eigtype)
	{
		LAPACKE_dsyev(matrix_layout, jobz, uplo, n, a, lda, w);
	}
	if (1 == eigtype)
	{
		LAPACKE_dsyevd(matrix_layout, jobz, uplo, n, a, lda, w);
	}
	if (2 == eigtype)
	{
		char range = 'A'; // all eigenvalues
		lapack_int* isuppz = new lapack_int[2 * n];
		double abstol = 0;
		LAPACKE_dsyevr(matrix_layout, jobz, range, uplo,
			n, a, lda, NULL, NULL, NULL, NULL,
			abstol, &n, w, a, n, isuppz);
	}
}

void DenseMatrixEigenSolver_FInterface(my_int matrix_layout, char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w)
{
	my_int eigtype = 0;	// 0, dsyevd; 1, dsyev; 2 dsyevr;
	lapack_int info;
	cout << "DenseMatrixEigenSolver_FInterface(my_int matrix_layout, char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w):" << endl;
	cout << "n = " << n << endl << endl;

	// dsyev
	if (0 == eigtype)
	{
		lapack_int lwork;
		double wkopt;
		double* work;
		lwork = -1;
		dsyev(&jobz, &uplo, &n, a, &lda, w, &wkopt, &lwork, &info);
		lwork = (my_int)wkopt;
		cout << "lwork = " << lwork << endl;
		work = (double*)malloc(lwork * sizeof(double));
		// Solve eigenproblem
		dsyev(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
		free(work);
	}
	// dsyevd
	if (1 == eigtype)
	{
		lapack_int lwork, liwork;
		lapack_int* iwork;
		double* work;
		lwork = 2 * n * n + 6 * n + 1;
		work = (double*)malloc(lwork * sizeof(double));
		liwork = 5 * n + 3;
		iwork = (lapack_int*)malloc(liwork * sizeof(my_int));
		// Solve eigenproblem
		dsyevd(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info);
		free(work);
		free(iwork);
	}

	// dsyevr 
	if (2 == eigtype)
	{
		lapack_int il, iu, ldz, info, lwork, liwork;
		ldz = lda;
		double abstol, vl, vu;
		lapack_int iwkopt;
		lapack_int* iwork;
		double wkopt;
		double* work;
		lapack_int* isuppz = new lapack_int[n];

		abstol = -1;
		il = 1;
		iu = n;
		lwork = -1;
		liwork = -1;
		dsyevr("Vectors", "Indices", "Upper", &n, a, &lda, &vl, &vu, &il, &iu,
			&abstol, &n, w, a, &ldz, isuppz, &wkopt, &lwork, &iwkopt, &liwork,
			&info);
		lwork = (my_int)wkopt;
		work = (double*)malloc(lwork * sizeof(double));
		liwork = iwkopt;
		iwork = (lapack_int*)malloc(liwork * sizeof(my_int));
		cout << "lwork = " << lwork << ", liwork = " << liwork << endl;
		// Solve eigenproblem
		dsyevr("Vectors", "Indices", "Upper", &n, a, &lda, &vl, &vu, &il, &iu,
			&abstol, &n, w, a, &ldz, isuppz, work, &lwork, iwork, &liwork,
			&info);

		free(iwork);
		free(work);
		delete[]isuppz;
	}

	/* Check for convergence */
	if (info > 0) {
		printf("The algorithm failed to compute eigenvalues.\n");
		exit(1);
	}
}

void DenseMatrixEigenSolver(my_int matrix_layout, char jobz, char uplo, lapack_int n, my_cmplx* a, lapack_int lda, double* w)
{
	my_int eigtype = 1;	// 0, dsyevd; 1, dsyev; 2 dsyevr;
	if (0 == eigtype)
	{
		LAPACKE_zheev(matrix_layout, jobz, uplo, n, a, lda, w);
	}
	if (1 == eigtype)
	{
		LAPACKE_zheevd(matrix_layout, jobz, uplo, n, a, lda, w);
	}
	if (2 == eigtype)
	{
		char range = 'A'; // all eigenvalues
		lapack_int* isuppz = new lapack_int[2 * n];
		double abstol = 0;
		LAPACKE_zheevr(matrix_layout, jobz, range, uplo,
			n, a, lda, NULL, NULL, NULL, NULL,
			abstol, &n, w, a, n, isuppz);
		free(isuppz);
	}
}

void DenseMatrixEigenSolver_FInterface(my_int matrix_layout, char jobz, char uplo, lapack_int n, my_cmplx* a, lapack_int lda, double* w)
{
	my_int eigtype = 0;	// 0, dsyevd; 1, dsyev; 2 dsyevr;
	lapack_int info;
	// dsyev
	cout << "DenseMatrixEigenSolver_FInterface(my_int matrix_layout, char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w):" << endl;
	cout << "n = " << n << endl;
	if (0 == eigtype)
	{
		lapack_int lwork;
		my_cmplx wkopt;
		my_cmplx* work;
		/* rwork dimension should be at least max(1,3*n-2) */
		double* rwork = new double[3 * n - 2];
		lwork = -1;
		zheev(&jobz, &uplo, &n, a, &lda, w, &wkopt, &lwork, rwork, &info);
		lwork = (my_int)wkopt.real();
		cout << "lwork = " << lwork;
		work = (my_cmplx*)malloc(lwork * sizeof(my_cmplx));
		// Solve eigenproblem
		zheev(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &info);
		delete[]rwork;
		free((void*)work);
	}
	// dsyevd
	if (1 == eigtype)
	{
		lapack_int lwork, lrwork, liwork;
		lapack_int iwkopt;
		lapack_int* iwork;
		double rwkopt;
		double* rwork;
		my_cmplx wkopt;
		my_cmplx* work;
		/* Query and allocate the optimal workspace */
		lwork = -1;
		lrwork = -1;
		liwork = -1;
		zheevd("Vectors", "Lower", &n, a, &lda, w, &wkopt, &lwork, &rwkopt,
			&lrwork, &iwkopt, &liwork, &info);
		lwork = (my_int)wkopt.real();
		work = (my_cmplx*)malloc(lwork * sizeof(my_cmplx));
		lrwork = (my_int)rwkopt;
		rwork = (double*)malloc(lrwork * sizeof(double));
		liwork = iwkopt;
		iwork = (my_int*)malloc(liwork * sizeof(my_int));
		/*
		lwork = 2 * n * n + 6 * n + 1;
		work = (my_cmplx*)malloc(lwork * sizeof(my_cmplx));
		liwork = 5 * n + 3;
		iwork = (lapack_int*)malloc(liwork * sizeof(my_int));
		*/
		// Solve eigenproblem
		zheevd(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
		/* Free workspace */
		free((void*)iwork);
		free((void*)rwork);
		free((void*)work);
	}
	/* Check for convergence */
	if (info > 0) {
		printf("The algorithm failed to compute eigenvalues.\n");
		exit(1);
	}
}

/*void MatrixSvd(my_int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, double* a, lapack_int lda, double* s, double* u, lapack_int ldu, double* vt, lapack_int ldvt) {
	double* superb = new double[m < n ? m : n];
	LAPACKE_dgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
}*/


void RK4_onestep(const my_int& dim, complex<double>* vals, my_int* cols, my_int* PointerBE, complex<double>* wf1, complex<double>* wf0, complex<double>& tau)
{
	sparse_matrix_t A;
	mkl_sparse_z_create_csr(&A, SPARSE_INDEX_BASE_ZERO, dim, dim, PointerBE, PointerBE + 1, cols, vals);
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_HERMITIAN;
	//descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	descrA.mode = SPARSE_FILL_MODE_UPPER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_optimize(A);

	complex<double> ii(0, 1.0);
	complex<double>* k1, * k2, * k3, * k4;
	k1 = new complex<double>[dim];
	k2 = new complex<double>[dim];
	k3 = new complex<double>[dim];
	k4 = new complex<double>[dim];

	// k1 = tau*H wf0 
	mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, tau, A, descrA, wf0, 0, k1);
	//mkl_zcsrmv(&trans, &dim, &dim, &tau, matdescra, vals, cols, PointerBE, PointerBE + 1, wf0, &num_zero, k1);
	// k2 = tau*H (wf0 + 0.5*k1) = k1 + 0.5*tau*H*k1 
	cblas_zcopy(dim, k1, 1, k2, 1);
	mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 0.5 * tau, A, descrA, k1, 1, k2);
	//mkl_zcsrmv(&trans, &dim, &dim, &tau_half, matdescra, vals, cols, PointerBE, PointerBE + 1, k1, &num_one, k2);
	// k3 = tau*H (wf0 + 0.5*k2) = k1 + 0.5*tau*H*k2
	cblas_zcopy(dim, k1, 1, k3, 1);
	mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 0.5 * tau, A, descrA, k2, 1, k3);
	//mkl_zcsrmv(&trans, &dim, &dim, &tau_half, matdescra, vals, cols, PointerBE, PointerBE + 1, k2, &num_one, k3);
	// k4 = -ii*H (wf0 + k3)     = k1 + tau*H*k3
	cblas_zcopy(dim, k1, 1, k4, 1);
	mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, tau, A, descrA, k3, 1, k4);
	//mkl_zcsrmv(&trans, &dim, &dim, &tau, matdescra, vals, cols, PointerBE, PointerBE + 1, k3, &num_one, k4);
	// wf1 = wf0 + (1/6)*(k1 + 2*k2 + 2*k3 + k4) ?
	// wf1 = (1/6)*(k1 + 2*k2 + 2*k3 + k4)
	complex<double> coeff = 1.0 / 6.0;
	for (my_int i = 0; i < dim; i++) {
		wf1[i] = wf0[i] + coeff * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
	}
	// normalization of wf1
	double res = cblas_dznrm2(dim, wf1, 1);
	for (my_int i = 0; i < dim; i++) {
		wf1[i] = wf1[i] / res;
	}
	//	cout << "check norm: " << cblas_dznrm2(Dim, wf1, 1) << endl;
	delete[]k1;
	delete[]k2;
	delete[]k3;
	delete[]k4;
}

void Lanczos(const my_int& dim, double* vals, my_int* cols, my_int* PointerBE, const my_int& M, char job_vec, double* alpha, double* l_vecs, double* l_eigvecs)
{
	sparse_matrix_t A;
	mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, dim, dim, PointerBE, PointerBE + 1, cols, vals);
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	descrA.mode = SPARSE_FILL_MODE_UPPER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_optimize(A);
	// 
	my_int steps = M;
	// Lanczos matrix elements
	double* beta = new double[steps];		// beta[0] is unused

	// if we do not need lanzcos vecters, only 3 vectors are used
	if ('N' == job_vec)
	{
		double* phi0 = new double[dim];
		double* phi1 = new double[dim];
		double* phi2 = new double[dim];
		// copy initial vector form 1st dim elements form l_vecs 
		cblas_dcopy(dim, l_vecs, 1, phi0, 1);

		// alpha[0], beta[1] -- beta[0] is not required
		mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, phi0, 0, phi1);
		alpha[0] = cblas_ddot(dim, phi0, 1, phi1, 1);
		cblas_daxpy(dim, -alpha[0], phi0, 1, phi1, 1);
		beta[1] = cblas_dnrm2(dim, phi1, 1);
		cblas_dscal(dim, 1.0 / beta[1], phi1, 1);

		// alpha[m], beta[m+1], m = [0,steps-2]
		for (my_int m = 1; m < steps - 1; m++)
		{
			// alpha[m], m = [0,steps-2]
			mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, phi1, 0, phi2);
			alpha[m] = cblas_ddot(dim, phi1, 1, phi2, 1);
			// beta[m+1], m = [0,steps-2]
			cblas_daxpy(dim, -alpha[m], phi1, 1, phi2, 1);
			cblas_daxpy(dim, -beta[m], phi0, 1, phi2, 1);
			beta[m + 1] = cblas_dnrm2(dim, phi2, 1);
			cblas_dscal(dim, 1.0 / beta[m + 1], phi2, 1);
			cblas_dcopy(dim, phi1, 1, phi0, 1);
			cblas_dcopy(dim, phi2, 1, phi1, 1);
		}
		// alpha[m], m = steps-1
		int m = steps - 1;
		mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, phi1, 0, phi2);
		alpha[m] = cblas_ddot(dim, phi1, 1, phi2, 1);

		LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'N', steps, alpha, &beta[0] + 1, NULL, steps);

		delete[]phi0;
		delete[]phi1;
		delete[]phi2;
	}

	// lanzcos vecters are required 
	if ('V' == job_vec)
	{
		// alpha[0], beta[1] -- beta[0] is not required
		mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, &l_vecs[0], 0, &l_vecs[dim]);
		alpha[0] = cblas_ddot(dim, &l_vecs[0], 1, &l_vecs[dim], 1);
		cblas_daxpy(dim, -alpha[0], &l_vecs[0], 1, &l_vecs[dim], 1);
		beta[1] = cblas_dnrm2(dim, &l_vecs[dim], 1);
		cblas_dscal(dim, 1.0 / beta[1], &l_vecs[dim], 1);

		// alpha[0], beta[1] -- beta[0] is not required
		for (my_int m = 1; m < steps - 1; m++)
		{
			// alpha[m], m = [0,steps-2]
			mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, &l_vecs[m * dim], 0, &l_vecs[(m + 1) * dim]);
			alpha[m] = cblas_ddot(dim, &l_vecs[m * dim], 1, &l_vecs[(m + 1) * dim], 1);
			// beta[m+1], m = [0,steps-2]
			cblas_daxpy(dim, -alpha[m], &l_vecs[m * dim], 1, &l_vecs[(m + 1) * dim], 1);
			cblas_daxpy(dim, -beta[m], &l_vecs[(m - 1) * dim], 1, &l_vecs[(m + 1) * dim], 1);
			beta[m + 1] = cblas_dnrm2(dim, &l_vecs[(m + 1) * dim], 1);
			cblas_dscal(dim, 1.0 / beta[m + 1], &l_vecs[(m + 1) * dim], 1);
		}
		// alpha[m], m = steps-1
		my_int m = steps - 1;
		double* tmpvec = new double[dim];
		mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, &l_vecs[m * dim], 0, tmpvec);
		alpha[m] = cblas_ddot(dim, &l_vecs[m * dim], 1, tmpvec, 1);
		delete[]tmpvec;

		LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'I', steps, alpha, &beta[0] + 1, l_eigvecs, steps);

		// test wf0
		/*
		double* wf0 = new double[dim];
		double* wf1 = new double[dim];
		for (my_int i = 0; i < dim; i++)
		{
			wf0[i] = my_vec_dot(M, &l_eigvecs[0], M, &l_vecs[i], dim);
		}
		double aux0 = cblas_dnrm2(dim, wf0, 1);
		cout << "aux0: " << aux0 << endl;
		cblas_dscal(dim, 1.0 / aux0, wf0, 1);
		mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, wf0, 0, wf1);
		cout << "<wf0|H|wf0> = " << setprecision(14) << cblas_ddot(dim, wf0, 1, wf1, 1) << endl;
		delete[]wf0;
		delete[]wf1;
		*/
	}
	delete[]beta;
	mkl_sparse_destroy(A);
}

void Lanczos(const my_int& dim, my_cmplx* vals, my_int* cols, my_int* PointerBE, const my_int& M, char job_vec,
	double* alpha, my_cmplx* l_vecs, double* l_eigvecs)
{
	sparse_matrix_t A;
	mkl_sparse_z_create_csr(&A, SPARSE_INDEX_BASE_ZERO, dim, dim, PointerBE, PointerBE + 1, cols, vals);
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_HERMITIAN;
	descrA.mode = SPARSE_FILL_MODE_UPPER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_optimize(A);
	// 
	my_int steps = M;
	// Lanczos matrix elements
	double* beta = new double[steps];		// beta[0] is unused

	// if we do not need lanzcos vecters, only 3 vectors are used
	if ('N' == job_vec)
	{
		my_cmplx* phi0 = new my_cmplx[dim];
		my_cmplx* phi1 = new my_cmplx[dim];
		my_cmplx* phi2 = new my_cmplx[dim];
		// copy initial vector form 1st dim elements form l_vecs 
		cblas_zcopy(dim, l_vecs, 1, phi0, 1);

		// alpha[0], beta[1] -- beta[0] is not required
		mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, phi0, 0, phi1);
		my_cmplx aux;
		cblas_zdotc_sub(dim, phi0, 1, phi1, 1, &aux);
		alpha[0] = aux.real();
		my_cmplx minus_aux = -aux;
		cblas_zaxpy(dim, &minus_aux, phi0, 1, phi1, 1);
		beta[1] = cblas_dznrm2(dim, phi1, 1);
		cblas_zdscal(dim, 1.0 / beta[1], phi1, 1);

		// alpha[m], beta[m+1], m = [0,steps-2]
		for (my_int m = 1; m < steps - 1; m++)
		{
			// alpha[m], m = [0,steps-2]
			mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, phi1, 0, phi2);
			cblas_zdotc_sub(dim, phi1, 1, phi2, 1, &aux);
			alpha[m] = aux.real();
			// beta[m+1], m = [0,steps-2]
			minus_aux = -aux;
			cblas_zaxpy(dim, &minus_aux, phi1, 1, phi2, 1);
			my_cmplx minus_beta = -beta[m];
			cblas_zaxpy(dim, &minus_beta, phi0, 1, phi2, 1);
			beta[m + 1] = cblas_dznrm2(dim, phi2, 1);
			cblas_zdscal(dim, 1.0 / beta[m + 1], phi2, 1);
			cblas_zcopy(dim, phi1, 1, phi0, 1);
			cblas_zcopy(dim, phi2, 1, phi1, 1);
		}
		// alpha[m], m = steps-1
		my_int m = steps - 1;
		mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, phi1, 0, phi2);
		cblas_zdotc_sub(dim, phi1, 1, phi2, 1, &aux);
		alpha[m] = aux.real();

		LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'N', steps, alpha, &beta[0] + 1, NULL, steps);

		delete[]phi0;
		delete[]phi1;
		delete[]phi2;
	}

	// lanzcos vecters are required 
	//my_cmplx* alpha_cmplx = new my_cmplx[dim];
	if ('V' == job_vec)
	{
		// alpha[0], beta[1] -- beta[0] is not required
		my_cmplx aux;
		mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, &l_vecs[0], 0, &l_vecs[dim]);
		cblas_zdotc_sub(dim, &l_vecs[0], 1, &l_vecs[dim], 1, &aux);
		alpha[0] = aux.real();
		my_cmplx minus_aux = -aux;
		cblas_zaxpy(dim, &minus_aux, &l_vecs[0], 1, &l_vecs[dim], 1);
		beta[1] = cblas_dznrm2(dim, &l_vecs[dim], 1);
		cblas_zdscal(dim, 1.0 / beta[1], &l_vecs[dim], 1);

		// alpha[m], beta[m+1], m = [0,steps-2]
		for (my_int m = 1; m < steps - 1; m++)
		{
			// alpha[m], m = [0,steps-2]
			mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, &l_vecs[m * dim], 0, &l_vecs[(m + 1) * dim]);
			cblas_zdotc_sub(dim, &l_vecs[m * dim], 1, &l_vecs[(m + 1) * dim], 1, &aux);
			alpha[m] = aux.real();
			// alpha[m+1], m = [0,steps-2]
			minus_aux = -aux;
			cblas_zaxpy(dim, &minus_aux, &l_vecs[m * dim], 1, &l_vecs[(m + 1) * dim], 1);
			my_cmplx minus_beta = -beta[m];
			cblas_zaxpy(dim, &minus_beta, &l_vecs[(m - 1) * dim], 1, &l_vecs[(m + 1) * dim], 1);
			beta[m + 1] = cblas_dznrm2(dim, &l_vecs[(m + 1) * dim], 1);
			cblas_zdscal(dim, 1.0 / beta[m + 1], &l_vecs[(m + 1) * dim], 1);
		}
		// alpha[m], m = steps-1
		my_int m = steps - 1;
		my_cmplx* tmpvec = new my_cmplx[dim];
		mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, &l_vecs[m * dim], 0, tmpvec);
		cblas_zdotc_sub(dim, &l_vecs[m * dim], 1, tmpvec, 1, &aux);
		delete[]tmpvec;
		alpha[m] = aux.real();

		//LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'I', steps, alpha, &beta[0] + 1, l_eigvecs, steps);
		//for (my_int m = 0; m < steps; m++) alpha[m] = alpha_cmplx[m].real();
		//LAPACKE_zsteqr(LAPACK_ROW_MAJOR, 'I', steps, alpha, &beta[0] + 1, l_eigvecs, steps);
		LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'I', steps, alpha, &beta[0] + 1, l_eigvecs, steps);

		// test wf0
		/*
		my_cmplx* wf0 = new my_cmplx[dim];
		my_cmplx* wf1 = new my_cmplx[dim];
		for (my_int i = 0; i < dim; i++)
		{
			//wf0[i] = my_vec_dot(M, &l_eigvecs[0], M, &l_vecs[i], dim);
			wf0[i] = 0;
			for (my_int m = 0; m < M; m++)
			{
				wf0[i] += l_eigvecs[m * M] * l_vecs[m * dim + i];
			}
		}
		double aux0 = cblas_dznrm2(dim, wf0, 1);
		cout << "aux0: " << endl;
		cblas_zdscal(dim, 1.0 / aux0, wf0, 1);
		mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, wf0, 0, wf1);
		cout << "<wf0|H|wf0> = " << setprecision(14) << my_vec_dot(dim, wf0, 1, wf1, 1) << endl;
		delete[]wf0;
		delete[]wf1;
		*/
	}

	delete[]beta;
	mkl_sparse_destroy(A);
}


double VN_entropy(const my_int& dim, double* p) {
	double aux = 0;
	for (my_int i = 0; i < dim; i++) {
		if (p[i] > 1e-32) {
			aux += -p[i] * log(p[i]);
		}
	}
	return aux;
}


// extra vector operations 
// y = a*x
void myvec_ax(const my_int& dim, const my_cmplx& a, my_cmplx* x, my_cmplx* y)
{
	for (my_int i = 0; i < dim; i++)
	{
		y[i] = a * x[i];
	}
}
void myvec_ax(const my_int& dim, const double& a, double* x, double* y)
{
	for (my_int i = 0; i < dim; i++)
	{
		y[i] = a * x[i];
	}
}
// a = real(x)
void myvec_real(const my_int& dim, my_cmplx* x, double* a)
{
	for (my_int i = 0; i < dim; i++)
	{
		a[i] = x[i].real();
	}
}
// a = imag(x)
void myvec_imag(const my_int& dim, my_cmplx* x, double* a)
{
	for (my_int i = 0; i < dim; i++)
	{
		a[i] = x[i].imag();
	}
}
// a = x + i*y, x and y are real vectors
void myvec_xpiy(const my_int& dim, double* x, double* y, my_cmplx* a)
{
	my_cmplx II(0, 1);
	for (my_int i = 0; i < dim; i++)
	{
		a[i] = x[i] + II * y[i];
	}
}


void MatrixSvd(int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, double* a, lapack_int lda, double* s, double* u, lapack_int ldu, double* vt, lapack_int ldvt)
{
	double* superb = new double[min(m, n)];
	LAPACKE_dgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
}

void MatrixSvd(int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, my_cmplx* a, lapack_int lda, double* s, my_cmplx* u, lapack_int ldu, my_cmplx* vt, lapack_int ldvt)
{
	double* superb = new double[min(m, n)];
	LAPACKE_zgesvd(matrix_layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
}