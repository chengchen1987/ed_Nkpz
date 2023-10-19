#pragma once

#include "common.h"

template <class Type>
class SparseMat {
public:
	SparseMat();
	~SparseMat();
	void SparseMat_Init(const my_int& m_dim, const my_int& m_len);
	void SparseMat_Copy(SparseMat& mat_b, my_int& m_dim, const my_int& m_len);
	void SparseMat_Copy_real2complex(SparseMat<double>& mat_b, my_int& m_dim, const my_int& m_len);
	void SparseMat_Write(char* fname);
	void SparseMat_Read(char* fname);
	void SparseMat_Clear();
	// size of matrix, Number of the nonzero elements (length of data)
	my_int mat_dim;
	my_int mat_nnz;
	// data, compressed sparse row format
	Type* vals;
	my_int* cols;
	my_int* Pointer_BE;
	// debug, get dense form 
	void SparseMat_ToDense();
	void SparseMat_sparsemv(double* wf, double* wf1);
	Type* denseform;
};

template <class Type>
SparseMat<Type>::SparseMat() {
}
template <class Type>
SparseMat<Type>::~SparseMat() {
}

template <class Type>
void SparseMat<Type>::SparseMat_Init(const my_int& m_dim, const my_int& m_len) {
	mat_dim = m_dim;
	mat_nnz = m_len;
	vals = new Type[mat_nnz];
	cols = new my_int[mat_nnz];
	Pointer_BE = new my_int[mat_dim + 1];
}

template <class Type>
void SparseMat<Type>::SparseMat_Copy(SparseMat& mat_b, my_int& m_dim, const my_int& m_len) {
	for (my_int i = 0; i < m_len; i++) {
		vals[i] = mat_b.vals[i];
		cols[i] = mat_b.cols[i];
	}
	for (my_int i = 0; i < m_dim + 1; i++) {
		Pointer_BE[i] = mat_b.Pointer_BE[i];
	}
}

template <class Type>
void SparseMat<Type>::SparseMat_Copy_real2complex(SparseMat<double>& mat_b, my_int& m_dim, const my_int& m_len) {
	for (my_int i = 0; i < m_len; i++) {
		vals[i] = mat_b.vals[i];
		cols[i] = mat_b.cols[i];
	}
	for (my_int i = 0; i < m_dim + 1; i++) {
		Pointer_BE[i] = mat_b.Pointer_BE[i];
	}
}

template <class Type>
void SparseMat<Type>::SparseMat_Clear() {
	delete[] vals;
	delete[] cols;
	delete[] Pointer_BE;
}

template <class Type>
void SparseMat<Type>::SparseMat_ToDense() {
	for (my_int i = 0; i < mat_dim * mat_dim; i++) {
		denseform[i] = 0;
	}
	for (my_int i = 0; i < mat_dim; i++) {
		for (my_int j = Pointer_BE[i]; j < Pointer_BE[i + 1]; j++) {
			denseform[i * mat_dim + cols[j]] += vals[j];
		}
	}
	// display
	/*
	for (my_int i = 0; i < mat_dim; i++) {
		for (my_int j = 0; j < mat_dim; j++) {
			std::cout << std::setw(10) << std::setprecision(4) << denseform[i * mat_dim + j];
		}
		std::cout << std::endl;
	}*/
}

// write and read ////////////////////////////////////////////////////
template <class Type>
void SparseMat<Type>::SparseMat_Write(char* fname) {
	FILE* fp = fopen(fname, "w");
	fwrite(&mat_dim, sizeof(my_int), 1, fp);
	fwrite(&mat_nnz, sizeof(my_int), 1, fp);
	fwrite(&vals, sizeof(Type), mat_nnz, fp);
	fwrite(&cols, sizeof(my_int), mat_nnz, fp);
	fwrite(&Pointer_BE, sizeof(my_int), mat_dim + 1, fp);
	fclose(fp);
}

template <class Type>
void SparseMat<Type>::SparseMat_Read(char* fname) {
	FILE* fr = fopen(fname, "w");
	fread(&mat_dim, sizeof(my_int), 1, fr);
	fread(&mat_nnz, sizeof(my_int), 1, fr);
	vals = new Type[mat_nnz];
	cols = new my_int[mat_nnz];
	Pointer_BE = new my_int[mat_dim + 1];
	fread(&vals, sizeof(Type), mat_nnz, fr);
	fread(&cols, sizeof(my_int), mat_nnz, fr);
	fread(&Pointer_BE, sizeof(my_int), mat_dim + 1, fr);
	fclose(fr);
}

template <class Type>
void SparseMat<Type>::SparseMat_sparsemv(double* wf, double* wf1)
{
	/*
	my_int* rows_start = new my_int[mat_dim];
	my_int* rows_end = new my_int[mat_dim];
	for (my_int i = 0; i < mat_dim; i++)
	{
		rows_start[i] = Pointer_BE[i];
		rows_end[i] = Pointer_BE[i+1];
	}

	sparse_matrix_t A;
	//mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, mat_dim, mat_dim, Pointer_BE, Pointer_BE + one, cols, vals);
	mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, mat_dim, mat_dim, rows_start, rows_end, cols, vals);
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	descrA.mode = SPARSE_FILL_MODE_UPPER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_optimize(A);
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, wf, 0, wf1);
	mkl_sparse_destroy(A);

	delete[]rows_start;
	delete[]rows_end;
	*/
	for (my_int i = 0; i < mat_dim; i++)
	{
		wf1[i] = 0;
	}
	for (my_int r = 0; r < mat_dim; r++)
	{
		for (my_int i = Pointer_BE[r]; i < Pointer_BE[r + 1]; i++)
		{
			wf1[cols[i]] += vals[i] * wf[r];
			if (cols[i] != r)
			{
				wf1[r] += vals[i] * wf[cols[i]];
			}
		}
	}

}

// sparse matrix operations in CSR format 
double SparseMat_Inner(SparseMat<double>& smat, double* wf0);
double SparseMat_Inner(SparseMat<double>& smat, my_cmplx* wf0);

// real Hamiltonian, real tau: imaginary time evolution, exp(-beta*H) 
void RK4_onestep(const SparseMat<double>& smat, const double& tau, double* wf0, double* wf1);
// real Hamiltonian, complex tau: real time evolution, exp(-i*t*H) 
void RK4_onestep(const SparseMat<double>& smat, const my_cmplx& tau, my_cmplx* wf0, my_cmplx* wf1);
// complex Hamiltonian, complex tau: real/imaginary time evolution 
void RK4_onestep(const SparseMat<my_cmplx>& smat, const my_cmplx& tau, my_cmplx* wf0, my_cmplx* wf1);

// one-step of Lanczos time evolution
// |wf1> = exp(-1i*H*dt)|wf0>
void Lanczos_evo_onestep(const SparseMat<my_cmplx>& smat, const my_int& M, const my_cmplx& dt, my_cmplx* wf0, my_cmplx* wf1);