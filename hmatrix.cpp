#include <cstdlib>
#include <cstdio>
#include <iostream> 
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <vector>
#include <tuple>
#include <algorithm>
using namespace std;

#include "basis.h" 
#include "hmatrix.h"
#include "sparsemat.h"

void Ham_HardCoreBoson<double>::DenseMat_Check_Hermitian()
{
	cout << "Check_Hermitian... " << endl;
	for (my_int i = 0; i < Dim; i++)
	{
		for (my_int j = 0; j < Dim; j++)
		{
			if (1e-10 < (DMat[i * Dim + j] - DMat[j * Dim + i]))
			{
				cout << "Hamtrix is not hermitian! " << i << " " << j << endl;
				cout << "DMat[i * Dim + j] = " << DMat[i * Dim + j] << ", DMat[j * Dim + i] = " << DMat[j * Dim + i] << endl;
			}
		}
	}

}

template<>
void Ham_HardCoreBoson<double>::Calc_Roag()
{
	// gaps
	double* gaps = new double[Dim - 1];
	for (unsigned long long i = 0; i < Dim - 1; i++) {
		gaps[i] = spec[i + 1] - spec[i];
	}
	// r0
	double* r0 = new double[Dim - 2];
	for (unsigned long long i = 0; i < Dim - 2; i++) {
		r0[i] = gaps[i + 1] / gaps[i];
	}
	// ratio of adjacent energy gaps 
	roag = new double[Dim - 2];
	for (unsigned long long i = 0; i < Dim - 2; i++) {
		roag[i] = r0[i] < 1 ? r0[i] : 1 / r0[i];
	}
	delete[]gaps;
	delete[]r0;

	my_int rstart = (Dim - 2) / 4;
	my_int rend = (Dim - 2) * 3 / 4;
	my_int rlength = rend - rstart;
	double mean_roag_mid = cblas_dasum(rlength, roag + rstart, 1) / rlength;
	cout << "mean_roag_mid = " << mean_roag_mid << endl;
}

template<>
void Ham_HardCoreBoson<double>::DenseMat_Get_psi_t(const my_int& p, const double& t, double* wf_r, double* wf_i)
{
	// evolution start from a product state
	// |wf> = exp(-iHt)|wf0>; wf0[p] = 1, wf0[q!=p] = 0 for all qs;
	// |wf> = U^{-1}*exp(-i*t*spec)*U*|wf0>; U*|wf0> = U[p,:]
	double* aux_r = new double[Dim];
	double* aux_i = new double[Dim];
	// exp(-i*t*spec)*U*|wf0>, U|p> = DMat[:,p]
	for (my_int s = 0; s < Dim; s++) {
		aux_r[s] = cos(-t * spec[s]) * DMat[s * Dim + p];;
		aux_i[s] = sin(-t * spec[s]) * DMat[s * Dim + p];;
	}
	// |wf> = U^{-1}*exp(-i*t*spec)*U*|wf0>
	cblas_dgemv(CblasRowMajor, CblasTrans, Dim, Dim, 1, DMat, Dim, aux_r, 1, 0, wf_r, 1);
	cblas_dgemv(CblasRowMajor, CblasTrans, Dim, Dim, 1, DMat, Dim, aux_i, 1, 0, wf_i, 1);
}

template<>
void Ham_HardCoreBoson<double>::DenseMat_Get_psi_t(my_int* inds, const my_int& ns, const double& t, double* wf_r, double* wf_i)
{
	// evolution start from a product state
	// |wf> = exp(-iHt)|wf0>; wf0[p] = 1, wf0[q!=p] = 0 for all qs;
	// |wf> = U^{-1}*exp(-i*t*spec)*U*|wf0>; U*|wf0> = U[p,:]
	double* aux_r = new double[Dim];
	double* aux_i = new double[Dim];
	// exp(-i*t*spec)*U*|wf0>, U|p> = DMat[:,p]
	for (my_int s = 0; s < Dim; s++) {
		aux_r[s] = 0;
		aux_i[s] = 0;
		for (my_int ip = 0; ip < ns; ip++)
		{
			aux_r[s] += cos(-t * spec[s]) * DMat[s * Dim + inds[ip]] / sqrt(ns);
			aux_i[s] += sin(-t * spec[s]) * DMat[s * Dim + inds[ip]] / sqrt(ns);
		}
	}
	// |wf> = U^{-1}*exp(-i*t*spec)*U*|wf0>
	cblas_dgemv(CblasRowMajor, CblasTrans, Dim, Dim, 1, DMat, Dim, aux_r, 1, 0, wf_r, 1);
	cblas_dgemv(CblasRowMajor, CblasTrans, Dim, Dim, 1, DMat, Dim, aux_i, 1, 0, wf_i, 1);
}

template<>
void Ham_HardCoreBoson<double>::DenseMat_Get_psi_t(double* wf0, const double& t, double* wf_r, double* wf_i)
{
	// evolution start from a general real wf0
	// |wf> = exp(-iHt)|wf0>; 
	//      = U^{-1}*exp(-i*t*spec)*U*|wf0>

	// U*|wf0>
	cblas_dgemv(CblasRowMajor, CblasNoTrans, Dim, Dim, 1, DMat, Dim, wf0, 1, 0, wf_r, 1);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, Dim, Dim, 1, DMat, Dim, wf0, 1, 0, wf_i, 1);
	// exp(-i*t*spec)*U*|wf0>
	for (my_int s = 0; s < Dim; s++) {
		wf_r[s] *= cos(-t * spec[s]);
		wf_i[s] *= sin(-t * spec[s]);
	}
	// |wf> = U^{-1}*exp(-i*t*spec)*U*|wf0>
	cblas_dgemv(CblasRowMajor, CblasTrans, Dim, Dim, 1, DMat, Dim, wf_r, 1, 0, wf_r, 1);
	cblas_dgemv(CblasRowMajor, CblasTrans, Dim, Dim, 1, DMat, Dim, wf_i, 1, 0, wf_i, 1);
}