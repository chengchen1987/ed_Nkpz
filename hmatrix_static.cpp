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

template <>
void Ham_HardCoreBoson<double>::Cal_H_ni_s(my_int i, double* ni_s)
{
	for (my_int s = 0; s < Dim; s++)
		ni_s[s] = ((basis->get_state(s) >> i) & 1);
}

template <>
void Ham_HardCoreBoson<double>::Cal_H_ninj_s(my_int i, my_int j, double* ninj_s)
{
	if (i == j)
		for (my_int s = 0; s < Dim; s++)
			ninj_s[s] = ((basis->get_state(s) >> i) & 1);
	if (i != j)
		for (my_int s = 0; s < Dim; s++)
		{
			ninj_s[s] = ((basis->get_state(s) >> i) & 1) * ((basis->get_state(s) >> j) & 1);
		}
}

template <>
void Ham_HardCoreBoson<double>::Cal_H_bibj_Smat(my_int i, my_int j, SparseMat<double>* _smat)
{
	// calculate all nonzero elements (maybe in a unprefered order), and then make a csr format
	my_int L = LatticeSize;
	double one = 1;

	// count nonzero elements 
	my_int index = 0;
	for (my_int s = 0; s < Dim; s++)
	{
		my_int state_s = basis->get_state(s);
		// set the first elements (with is diagonal) as nonzero
		index++;
		// (0 == r): no extra nonzero element
		if (i != j)
		{
			if ((((state_s >> i) & 1)) ^ (((state_s >> j) & 1)))
			{
				my_int state_s1 = state_s ^ (1 << i) ^ (1 << j);
				my_int s1 = basis->get_index(state_s1);
				if (-1 != s1)
				{
					index++;
				}
			}
		}
	}
	my_int nnz = index;

	// construct sparse matrix
	_smat->SparseMat_Init(Dim, nnz);
	index = 0;
	for (my_int s = 0; s < Dim; s++)
	{
		my_int state_s = basis->get_state(s);

		// i = 0; set the first elements (with is diagonal) as nonzero
		_smat->vals[index] = 0;
		_smat->cols[index] = s;
		_smat->Pointer_BE[s] = index;
		index++;
		if (i != j)
		{
			if ((((state_s >> i) & 1)) ^ (((state_s >> j) & 1)))
			{
				my_int state_s1 = state_s ^ (1 << i) ^ (1 << j);
				my_int s1 = basis->get_index(state_s1);
				if (-1 != s1)
				{

					_smat->vals[index] = -1;
					_smat->cols[index] = s1;
					index++;
				}
			}
		}
		_smat->Pointer_BE[s + 1] = index;
	}
}

template<>
void Ham_HardCoreBoson<double>::DenseMat_Static()
{
	my_int L = LatticeSize;

	// quntities from eigenvalues 
	char fspec[80];
	sprintf(fspec, "spec.bin");
	Vec_fwrite_double(fspec, spec, Dim);

	sprintf(fspec, "Hdiag.bin");
	Vec_fwrite_double(fspec, H_diag, Dim);

	Calc_Roag();
	sprintf(fspec, "roag.bin");
	Vec_fwrite_double(fspec, roag, Dim - 2);

	// quntities from eigenvectors
	double* ipr = new double[Dim];
	double* SE = new double[Dim];
	double* HI = new double[Dim];
	double* EE = new double[Dim];
	//double* nnr = new double[Dim * (L / 2 + 1)];

	double* wf = new double[Dim];
	double* wf_sq = new double[Dim];

	my_int time0 = time(0);
	for (my_int s = 0; s < Dim; s++)
	{
		Fetch_EigenVec(s, wf);
		vdSqr(Dim, wf, wf_sq);
		//nk0[s] = cblas_ddot(Dim, wf_sq, 1, nk0_p, 1);
		HI[s] = cblas_ddot(Dim, wf_sq, 1, H_diag, 1);
		ipr[s] = Cal_IPR_wfsq(wf_sq);
		SE[s] = Cal_SE_wfsq(wf_sq);
		EE[s] = Cal_halfchain_EE_wf(wf);
		//Cal_Hkpz_nnr_wfsq(wf_sq, &nnr[s * (L / 2 + 1)]);
	}
	my_int time1 = time(0);
	sprintf(fspec, "eig_HI.bin");
	Vec_fwrite_double(fspec, HI, Dim);
	sprintf(fspec, "eig_ipr.bin");
	Vec_fwrite_double(fspec, ipr, Dim);
	sprintf(fspec, "eig_se.bin");
	Vec_fwrite_double(fspec, SE, Dim);
	sprintf(fspec, "eig_EE.bin");
	Vec_fwrite_double(fspec, EE, Dim);

	delete[]ipr;
	delete[]SE;
	delete[]HI;
	delete[]EE;


	double* eig_ni = new double[Dim * L];
	double* eig_ninj = new double[Dim * L * L];
	double* ni_s = new double[Dim * L];
	for (my_int i = 0; i < L; i++) { Cal_H_ni_s(i, &ni_s[i * Dim]); }
	double* ninj_s = new double[Dim * L * L];
	for (my_int i = 0; i < L; i++)
		for (my_int j = 0; j < L; j++) { Cal_H_ninj_s(i, j, &ninj_s[(i * L + j) * Dim]); }

	for (my_int s = 0; s < Dim; s++)
	{
		Fetch_EigenVec(s, wf);
		vdSqr(Dim, wf, wf_sq);
		for (my_int i = 0; i < L; i++)
		{
			eig_ni[s * L + i] = cblas_ddot(Dim, wf_sq, 1, &ni_s[i * Dim], 1);
		}

		for (my_int i = 0; i < L; i++)
			for (my_int j = 0; j < L; j++)
			{
				eig_ninj[s * L + i * L + j] = cblas_ddot(Dim, wf_sq, 1, &ninj_s[(i * L + j) * Dim], 1);
			}
	}
	sprintf(fspec, "eig_ni.bin");
	Vec_fwrite_double(fspec, eig_ni, Dim * L);
	sprintf(fspec, "eig_ninj.bin");
	Vec_fwrite_double(fspec, eig_ninj, Dim * L * L);

	delete[]ni_s;
	delete[]ninj_s;
	delete[]eig_ni;
	delete[]eig_ninj;

	delete[]wf;
	delete[]wf_sq;
	//delete[]nnr;

	my_int time2 = time(0);

	cout << "Time others: " << time1 - time0 << ", time correlations: " << time2 - time1 << endl;
}