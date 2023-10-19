#pragma once
#ifndef COMMON_H
#define COMMON_H

#define _CRT_SECURE_NO_DEPRECATE
#define PI 3.14159265358979323846

#include <cstdint>
#include <complex>

#define MKL_Complex16 std::complex<double>
#include <mkl.h>
#define my_int MKL_INT
#define my_cmplx std::complex<double>

void rd_wf(const my_int& dim, double* wf);
void rd_wf(const my_int& dim, my_cmplx* wf);
void uni_wf(const my_int& dim, double* wf);
void rd_wf_product(const my_int& dim, double* wf);

// some extra bitwise operations
my_int Bits_Count1(const my_int& a);
my_int Bits_CycleLeft(const my_int& s, const my_int &n);
my_int Bits_CycleRight(const my_int& s, const my_int& n);
my_int Bits_Reflect(const my_int& s, const my_int& n);
my_int Bits_Invert(const my_int&s, const my_int &n);
my_int Bits_HammingDis(const my_int& a, const my_int& b);
void Bits_Print(const my_int& a, const my_int& n);
void Bits_Decomposition(const my_int& s, const my_int& n, const my_int& size_A, my_int* sites_A, my_int& a, my_int& b);
// 
my_int nchoosek(const my_int& k, const my_int& n);

void Vec_fwrite_double(const char* fname, double* data, const my_int& dsize);
void Vec_fread_double(const char* fname, double* data, const my_int& dsize);
//
void Debug_Check_Conj_dense(double* mat, my_int dim);

double Cal_chargesign(const my_int& s, const my_int& i);	// minus sign is considered for fermion exchange 
double Cal_chargesign(const my_int& s, const my_int& i, const my_int& j);
double num_conj(double a);
my_cmplx num_conj(my_cmplx a);
double my_vec_dot(const my_int &n, double*x, const my_int incx, double*y, const my_int incy);
my_cmplx my_vec_dot(const my_int &n, my_cmplx* x, const my_int incx, my_cmplx* y, const my_int incy);

// parse input 
void GetParaFromInput_int(const char* fname, const char* string_match, my_int& para);
void GetParaFromInput_real(const char* fname, const char* string_match, double& para);
void GetParaFromInput_char(const char* fname, const char* string_match, char& para);

//void Lanczos_timestep(const my_int& Dim, double* vals, my_int* cols, my_int* PointerBE, const my_int& M, const double &dt, double *wft);

//void RK4_onestep(const my_int& Dim, my_cmplx* vals, my_int* cols, my_int* PointerBE, my_cmplx* wf1, my_cmplx* wf0, my_cmplx& tau);

void DenseMatrixEigenSolver(my_int matrix_layout, char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w);
void DenseMatrixEigenSolver_FInterface(my_int matrix_layout, char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w);

void DenseMatrixEigenSolver(my_int matrix_layout, char jobz, char uplo, lapack_int n, my_cmplx* a, lapack_int lda, double* w);
void DenseMatrixEigenSolver_FInterface(my_int matrix_layout, char jobz, char uplo, lapack_int n, my_cmplx* a, lapack_int lda, double* w);

void MatrixSvd(int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, double* a, lapack_int lda, double* s, double* u, lapack_int ldu, double* vt, lapack_int ldvt);
void MatrixSvd(int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, my_cmplx* a, lapack_int lda, double* s, my_cmplx* u, lapack_int ldu, my_cmplx* vt, lapack_int ldvt);

// 
double VN_entropy(const my_int& dim, double* p);


double my_vec_dot(const my_int& n, double* x, const my_int incx, double* y, const my_int incy);
my_cmplx my_vec_dot(const my_int& n, my_cmplx* x, const my_int incx, my_cmplx* y, const my_int incy);


void Lanczos(const my_int& Dim, double* vals, my_int* cols, my_int* PointerBE, const my_int& M, char job_vec,
	double* l_vals, double* l_vecs, double* l_eigvecs); // job_vec: ('N')'I' (no) Lanzcos vectors 
//void Lanczos(const my_int& Dim, my_cmplx* vals, my_int* cols, my_int* PointerBE, const my_int& M, char job_vec,
//	double* l_vals, my_cmplx* l_vecs, my_cmplx* l_eigvecs); // job_vec: ('N')'I' (no) Lanzcos vectors 
void Lanczos(const my_int& dim, my_cmplx* vals, my_int* cols, my_int* PointerBE, const my_int& M, char job_vec,
	double* alpha, my_cmplx* l_vecs, double* l_eigvecs);

// extra vector operations 
// y = a*x
void myvec_ax(const my_int& dim, const double& a, double* x, double* y);
void myvec_ax(const my_int& dim, const my_cmplx& a, my_cmplx* x, my_cmplx* y);
// a = real(x)
void myvec_real(const my_int& dim, my_cmplx* x, double* a);
// a = imag(x)
void myvec_imag(const my_int& dim, my_cmplx* x, double* a);
// a = x + i*y, x and y are real vectors
void myvec_xpiy(const my_int& dim, double* x, double* y, my_cmplx* a);
#endif // !COMMON_H