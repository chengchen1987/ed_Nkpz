#include <iostream>
#include <iomanip> 
#include "sparsemat.h"

using namespace std;

// sparse matrix operations in CSR format 

// <wf0|H|wf0>
double SparseMat_Inner(SparseMat<double>& smat, double* wf0)
{
	// 
	sparse_matrix_t A;
	mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, smat.mat_dim, smat.mat_dim, smat.Pointer_BE, smat.Pointer_BE + 1, smat.cols, smat.vals);
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	descrA.mode = SPARSE_FILL_MODE_UPPER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_optimize(A);
	// 
	double* wf1 = new double[smat.mat_dim];
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, wf0, 0, wf1);
	double aux = cblas_ddot(smat.mat_dim, wf0, 1, wf1, 1);
	delete[]wf1;
	mkl_sparse_destroy(A);
	//
	return aux;
}

// <wf0|H|wf0>
double SparseMat_Inner(SparseMat<double>& smat, my_cmplx* wf0)
{
	// 
	sparse_matrix_t A;
	mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, smat.mat_dim, smat.mat_dim, smat.Pointer_BE, smat.Pointer_BE + 1, smat.cols, smat.vals);
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	descrA.mode = SPARSE_FILL_MODE_UPPER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_optimize(A);
	// 
	double* wf1_real = new double[smat.mat_dim];
	myvec_real(smat.mat_dim, wf0, wf1_real);
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, wf1_real, 0, wf1_real);
	double* wf1_imag = new double[smat.mat_dim];
	myvec_imag(smat.mat_dim, wf0, wf1_imag);
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, wf1_imag, 0, wf1_imag);
	my_cmplx* wf1 = new my_cmplx[smat.mat_dim];
	myvec_xpiy(smat.mat_dim, wf1_real, wf1_imag, wf1);

	my_cmplx aux;
	cblas_zdotc_sub(smat.mat_dim, wf0, 1, wf1, 1, &aux);
	if (aux.imag() > 1e-10)
		cout << "Error in : cblas_zdotc_sub() <- SparseMat_Inner()!" << endl;

	delete[]wf1;
	delete[]wf1_real;
	delete[]wf1_imag;
	mkl_sparse_destroy(A);
	//
	return aux.real();
}

// real Hamiltonian, real tau: imaginary time evolution 
void RK4_onestep(const SparseMat<double>& smat, const double& tau, double* wf0, double* wf1)
{
	// initialize sparse matrix 
	sparse_matrix_t A;
	mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, smat.mat_dim, smat.mat_dim, smat.Pointer_BE, smat.Pointer_BE + 1, smat.cols, smat.vals);
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	descrA.mode = SPARSE_FILL_MODE_UPPER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_optimize(A);

	// auxiliary vectors
	double* k1_real = new double[smat.mat_dim];
	double* k2_real = new double[smat.mat_dim];
	double* k3_real = new double[smat.mat_dim];
	double* k4_real = new double[smat.mat_dim];
	// k1 = tau*H*wf0 
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, tau, A, descrA, wf0, 0, k1_real);
	// k2 = tau*H (wf0 + 0.5*k1) = k1 + 0.5*tau*H*k1 
	cblas_dcopy(smat.mat_dim, k1_real, 1, k2_real, 1);
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 0.5 * tau, A, descrA, k1_real, 1, k2_real);
	// k3 = tau*H (wf0 + 0.5*k2) = k1 + 0.5*tau*H*k2
	cblas_dcopy(smat.mat_dim, k1_real, 1, k3_real, 1);
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 0.5 * tau, A, descrA, k2_real, 1, k3_real);
	// k4 = -ii*H (wf0 + k3) = k1 + tau*H*k3
	cblas_dcopy(smat.mat_dim, k1_real, 1, k4_real, 1);
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, tau, A, descrA, k3_real, 1, k4_real);
	// wf1 = wf0 + (1/6)*(k1 + 2*k2 + 2*k3 + k4) ?
	// wf1 = (1/6)*(k1 + 2*k2 + 2*k3 + k4)
	double coeff = 1.0 / 6.0;
	for (my_int i = 0; i < smat.mat_dim; i++) {
		wf1[i] = wf0[i] + coeff * (
			(k1_real[i] + 2.0 * k2_real[i] + 2.0 * k3_real[i] + k4_real[i])
			);
	}
	// normalization of wf1
	double res = cblas_dnrm2(smat.mat_dim, wf1, 1);
	myvec_ax(smat.mat_dim, 1 / res, wf1, wf1);
	//	cout << "check norm: " << cblas_dznrm2(Dim, wf1, 1) << endl;
	delete[]k1_real;
	delete[]k2_real;
	delete[]k3_real;
	delete[]k4_real;
	//
	mkl_sparse_destroy(A);
}

void RK4_onestep(const SparseMat<double>& smat, const my_cmplx& tau, my_cmplx* wf0, my_cmplx* wf1)
{
	my_cmplx II(0, 1);

	// initialize sparse matrix 
	sparse_matrix_t A;
	mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, smat.mat_dim, smat.mat_dim, smat.Pointer_BE, smat.Pointer_BE + 1, smat.cols, smat.vals);
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	descrA.mode = SPARSE_FILL_MODE_UPPER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_optimize(A);

	// auxiliary vectors
	double* k1_real = new double[smat.mat_dim];
	double* k1_imag = new double[smat.mat_dim];
	double* k2_real = new double[smat.mat_dim];
	double* k2_imag = new double[smat.mat_dim];
	double* k3_real = new double[smat.mat_dim];
	double* k3_imag = new double[smat.mat_dim];
	double* k4_real = new double[smat.mat_dim];
	double* k4_imag = new double[smat.mat_dim];
	// k1 = tau*H*wf0 
	my_cmplx* wf_aux = new my_cmplx[smat.mat_dim];
	for (my_int i = 0; i < smat.mat_dim; i++) wf_aux[i] = 0;
	myvec_ax(smat.mat_dim, tau, wf0, wf_aux);
	myvec_real(smat.mat_dim, wf_aux, k1_real);
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, k1_real, 0, k1_real);
	myvec_imag(smat.mat_dim, wf_aux, k1_imag);
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, k1_imag, 0, k1_imag);
	// k2 = tau*H (wf0 + 0.5*k1) = k1 + 0.5*tau*H*k1 
	myvec_xpiy(smat.mat_dim, k1_real, k1_imag, wf_aux);
	myvec_ax(smat.mat_dim, 0.5 * tau, wf_aux, wf_aux);
	myvec_real(smat.mat_dim, wf_aux, k2_real);
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, k2_real, 0, k2_real);
	myvec_imag(smat.mat_dim, wf_aux, k2_imag);
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, k2_imag, 0, k2_imag);
	vdAdd(smat.mat_dim, k1_real, k2_real, k2_real);
	vdAdd(smat.mat_dim, k1_imag, k2_imag, k2_imag);
	// k3 = tau*H (wf0 + 0.5*k2) = k1 + 0.5*tau*H*k2
	myvec_xpiy(smat.mat_dim, k2_real, k2_imag, wf_aux);
	myvec_ax(smat.mat_dim, 0.5 * tau, wf_aux, wf_aux);
	myvec_real(smat.mat_dim, wf_aux, k3_real);
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, k3_real, 0, k3_real);
	myvec_imag(smat.mat_dim, wf_aux, k3_imag);
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, k3_imag, 0, k3_imag);
	vdAdd(smat.mat_dim, k1_real, k3_real, k3_real);
	vdAdd(smat.mat_dim, k1_imag, k3_imag, k3_imag);
	// k4 = tau*H (wf0 + k3) = k1 + tau*H*k3
	myvec_xpiy(smat.mat_dim, k3_real, k3_imag, wf_aux);
	myvec_ax(smat.mat_dim, tau, wf_aux, wf_aux);
	myvec_real(smat.mat_dim, wf_aux, k4_real);
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, k4_real, 0, k4_real);
	myvec_imag(smat.mat_dim, wf_aux, k4_imag);
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, k4_imag, 0, k4_imag);
	vdAdd(smat.mat_dim, k1_real, k4_real, k4_real);
	vdAdd(smat.mat_dim, k1_imag, k4_imag, k4_imag);

	// wf1 = wf0 + (1/6)*(k1 + 2*k2 + 2*k3 + k4) 
	double coeff = 1.0 / 6.0;
	for (my_int i = 0; i < smat.mat_dim; i++) {
		wf1[i] = wf0[i] + coeff * (
			(k1_real[i] + 2.0 * k2_real[i] + 2.0 * k3_real[i] + k4_real[i]) +
			(k1_imag[i] + 2.0 * k2_imag[i] + 2.0 * k3_imag[i] + k4_imag[i]) * II
			);
	}
	// normalization of wf1
	double res = cblas_dznrm2(smat.mat_dim, wf1, 1);
	myvec_ax(smat.mat_dim, 1 / res, wf1, wf1);
	//cout << "check norm: " << cblas_dznrm2(smat.mat_dim, wf1, 1) << endl;

	delete[]k1_real;
	delete[]k1_imag;
	delete[]k2_real;
	delete[]k2_imag;
	delete[]k3_real;
	delete[]k3_imag;
	delete[]k4_real;
	delete[]k4_imag;
	delete[]wf_aux;
	//
	mkl_sparse_destroy(A);
}

// real Hamiltonian, real tau: imaginary time evolution 
void RK4_onestep(const SparseMat<my_cmplx>& smat, const my_cmplx& tau, my_cmplx* wf0, my_cmplx* wf1)
{
	// initialize sparse matrix 
	sparse_matrix_t A;
	mkl_sparse_z_create_csr(&A, SPARSE_INDEX_BASE_ZERO, smat.mat_dim, smat.mat_dim, smat.Pointer_BE, smat.Pointer_BE + 1, smat.cols, smat.vals);
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_HERMITIAN;;
	descrA.mode = SPARSE_FILL_MODE_UPPER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_optimize(A);

	// auxiliary vectors
	my_cmplx* k1 = new my_cmplx[smat.mat_dim];
	my_cmplx* k2 = new my_cmplx[smat.mat_dim];
	my_cmplx* k3 = new my_cmplx[smat.mat_dim];
	my_cmplx* k4 = new my_cmplx[smat.mat_dim];
	// k1 = tau*H*wf0 
	mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, tau, A, descrA, wf0, 0, k1);
	// k2 = tau*H (wf0 + 0.5*k1) = k1 + 0.5*tau*H*k1 
	cblas_zcopy(smat.mat_dim, k1, 1, k2, 1);
	mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 0.5 * tau, A, descrA, k1, 1, k2);
	// k3 = tau*H (wf0 + 0.5*k2) = k1 + 0.5*tau*H*k2
	cblas_zcopy(smat.mat_dim, k1, 1, k3, 1);
	mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 0.5 * tau, A, descrA, k2, 1, k3);
	// k4 = -ii*H (wf0 + k3) = k1 + tau*H*k3
	cblas_zcopy(smat.mat_dim, k1, 1, k4, 1);
	mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, tau, A, descrA, k3, 1, k4);
	// wf1 = wf0 + (1/6)*(k1 + 2*k2 + 2*k3 + k4) ?
	// wf1 = (1/6)*(k1 + 2*k2 + 2*k3 + k4)
	double coeff = 1.0 / 6.0;
	for (my_int i = 0; i < smat.mat_dim; i++) {
		wf1[i] = wf0[i] + coeff * (
			(k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
			);
	}
	// normalization of wf1
	double res = cblas_dznrm2(smat.mat_dim, wf1, 1);
	myvec_ax(smat.mat_dim, 1 / res, wf1, wf1);
	//	cout << "check norm: " << cblas_dznrm2(Dim, wf1, 1) << endl;
	delete[]k1;
	delete[]k2;
	delete[]k3;
	delete[]k4;
	//
	mkl_sparse_destroy(A);
}

// one-step of Lanczos time evolution
// |wf1> = exp(-1i*H*dt)|wf0>
void Lanczos_evo_onestep(const SparseMat<my_cmplx>& smat, const my_int& M, const my_cmplx& dt, my_cmplx* wf0, my_cmplx* wf1)
{
	my_cmplx II(0, 1);

	my_int dim = smat.mat_dim;
	for (my_int i = 0; i < dim; i++) wf1[i] = 0;

	double* l_eigvecs = new double[M * M];			// eigenvectors in Lanczos basis (real even for complex wavefunctions)
	double* alpha = new double[M];
	my_cmplx* l_vecs = new my_cmplx[M * dim];
	// initial state as the 1st row l_vecs
	cblas_zcopy(dim, wf0, 1, l_vecs, 1);

	Lanczos(dim, smat.vals, smat.cols, smat.Pointer_BE, M, 'V', alpha, l_vecs, l_eigvecs);

	my_cmplx* auxvec = new my_cmplx[M];
	// auxvec = exp(-1i*dt*alpha)*U^+*V^+*|Psi(t)> = exp(-1i*dt*alpha)*U^+*(1,0,0,...)^T
	for (my_int i = 0; i < M; i++)
	{
		auxvec[i] = l_eigvecs[i] * exp(-II * alpha[i] * dt);
	}

	// auxvec = U*exp(-1i*dt*alpha)*U^+*(1,0,0,...)^T = U*auxvec
	double* auxvec_real = new double[M];
	myvec_real(M, auxvec, auxvec_real);
	double* auxvec_imag = new double[M];
	myvec_imag(M, auxvec, auxvec_imag);
	double* auxvec_real_1 = new double[M];
	double* auxvec_imag_1 = new double[M];
	cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1, l_eigvecs, M, auxvec_real, 1, 0, auxvec_real_1, 1);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, M, M, 1, l_eigvecs, M, auxvec_imag, 1, 0, auxvec_imag_1, 1);

	myvec_xpiy(M, auxvec_real_1, auxvec_imag_1, auxvec);
	delete[]auxvec_real;
	delete[]auxvec_imag;
	delete[]auxvec_real_1;
	delete[]auxvec_imag_1;

	// wf1 = V*U*exp(-1i*dt*alpha)*U^+*(1,0,0,...)^T = (l_vecs^T)*auxvec
	my_cmplx cmplx_1(1, 0);
	my_cmplx cmplx_0(0, 0);
	cblas_zgemv(CblasRowMajor, CblasTrans, M, dim, &cmplx_1, l_vecs, dim, auxvec, 1, &cmplx_0, wf1, 1);

	cout << "check norm: " << cblas_dznrm2(dim, wf1, 1) << endl;
	cblas_zcopy(dim, wf1, 1, wf0, 1);

	delete[]auxvec;

	delete[]l_eigvecs;
	delete[]l_vecs;
	delete[]alpha;
}