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
void Ham_HardCoreBoson<double>::Cal_Hkpz_nnr_s(my_int r, double* nnr_s)
{
	my_int L = LatticeSize;

	for (my_int s = 0; s < Dim; s++) {
		nnr_s[s] = 0;

		my_int state_s;
		state_s = basis_kpz->get_state(s);

		// check coeff_n; state_s can be the same with different sigmas, consider only once here
		/*
		my_int na;
		if (s > 0 && basis_kpz->get_state(s - 1) == state_s) { continue; }
		else if (s < Dim - 1 && basis_kpz->get_state(s + 1) == state_s) { na = 2; }
		else { na = 1; }

		for (my_int i = 0; i < L; i++)
		{
			my_int j = (i + r) % L;
			nnr_s[s] += ((basis_kpz->get_state(s) >> i) & 1) * ((basis_kpz->get_state(s) >> j) & 1) * na * na;
		}*/

		if (0 == r)
		{
			for (my_int i = 0; i < L; i++)
			{
				nnr_s[s] += ((basis_kpz->get_state(s) >> i) & 1);
			}
		}
		else
		{
			for (my_int i = 0; i < L; i++)
			{
				my_int j = (i + r) % L;
				nnr_s[s] += ((basis_kpz->get_state(s) >> i) & 1) * ((basis_kpz->get_state(s) >> j) & 1);
			}
		}
	}
}

template <>
void Ham_HardCoreBoson<double>::Cal_Hkpz_nnr_wfsq(double* wfsq, double* nnr)
{
	my_int L = LatticeSize;
	double one = 1;

	double* wf1 = new double[Dim];

	for (my_int r = 0; r <= L / 2; r++) nnr[r] = 0;

	for (my_int r = 0; r <= L / 2; r++)
	{
		for (my_int s = 0; s < Dim; s++) wf1[s] = 0;
		for (my_int s = 0; s < Dim; s++)
		{
			my_int state_s = basis_kpz->get_state(s);
			my_int na;
			if (s > 0 && basis_kpz->get_state(s - 1) == state_s) { continue; }
			else if (s < Dim - 1 && basis_kpz->get_state(s + 1) == state_s) { na = 2; }
			else { na = 1; }

			for (my_int i = 0; i < L; i++)
			{
				my_int j = (i + r) % L;
				// j = i
				if (0 == r && (state_s >> i) & 1)
				{
					for (my_int sa = s; sa < s + na; sa++)
					{
						nnr[r] += wfsq[s];
					}
				}
				else if ((((state_s >> i) & 1)) & (((state_s >> j) & 1)))
				{
					for (my_int sa = s; sa < s + na; sa++)
					{
						nnr[r] += wfsq[s];
					}
				}
			}
		}
	}
	delete[]wf1;
}

template <>
double Ham_HardCoreBoson<double>::Cal_Ob_smat(double* wf, SparseMat<double>* _smat)
{
	/*
	sparse_matrix_t A;
	sparse_status_t status = mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, Dim, Dim,
		_smat->Pointer_BE, _smat->Pointer_BE + 1, _smat->cols, _smat->vals);
	if (0 != status) {
		cout << "mkl_sparse_d_create_csr, sparse_status_t: " << status << endl;
	}
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	descrA.mode = SPARSE_FILL_MODE_UPPER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_optimize(A);

	cout << "c^dagger_i c_j matrix:" << endl;
	for (my_int r = 0; r < _smat->mat_dim; r++)
	{
		cout << "r = " << r << endl;
		for (my_int i = _smat->Pointer_BE[r]; i < _smat->Pointer_BE[r + 1]; i++)
		{
			cout << "c = " << _smat->cols[i] << "	, val = " << _smat->vals[i] << endl;
		}
	}

	double* wf1 = new double[Dim];
	sparse_status_t status1 = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, wf, 0, wf1);
	//if (0 != status1) {
	cout << "mkl_sparse_d_mv, sparse_status_t: " << status1 << endl;
	//}
	double aux = cblas_ddot(Dim, wf, 1, wf1, 1);

	double* wf2 = new double[Dim];
	_smat->denseform = new double[Dim * Dim];
	_smat->SparseMat_ToDense();
	cblas_dsymv(CblasRowMajor, CblasUpper, Dim, 1, _smat->denseform, Dim, wf, 1, 0, wf2, 1);

	double* wf3 = new double[Dim];
	_smat->SparseMat_sparsemv(wf, wf3);

	cout << "wf: " << endl;
	for (my_int r = 0; r < _smat->mat_dim; r++)
	{
		cout << "r = " << r << ", wf[r] = " << wf[r] << ", wf1[r] = " << wf1[r] << ", wf2[r] = " << wf2[r] << ", wf3[r] = " << wf3[r] << endl;
	}

	cout << "cblas_ddot(Dim, wf, 1, wf1, 1): " << aux << endl;
	//cout << "cblas_ddot(Dim, wf, 1, wf, 1): " << cblas_ddot(Dim, wf, 1, wf, 1) << endl;
	cout << "cblas_ddot(Dim, wf, 1, wf2, 1): " << cblas_ddot(Dim, wf, 1, wf2, 1) << endl;
	cout << "cblas_ddot(Dim, wf, 1, wf3, 1): " << cblas_ddot(Dim, wf, 1, wf3, 1) << endl;

	mkl_sparse_destroy(A);
	delete[]wf1;

	delete[]wf2;
	delete[]_smat->denseform;
	*/
	double* wf3 = new double[Dim];
	_smat->SparseMat_sparsemv(wf, wf3);
	double aux =  cblas_ddot(Dim, wf, 1, wf3, 1);
	delete[]wf3;
	return aux;
}

template <>
double Ham_HardCoreBoson<double>::Cal_Ob_smat(double* wf_r, double* wf_i, SparseMat<double>* _smat)
{
	sparse_matrix_t A;
	sparse_status_t status = mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, Dim, Dim, _smat->Pointer_BE, _smat->Pointer_BE + 1, _smat->cols, _smat->vals);
	if (0 != status) {
		cout << "mkl_sparse_d_create_csr, sparse_status_t: " << status << endl;
	}
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	descrA.mode = SPARSE_FILL_MODE_UPPER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	mkl_sparse_optimize(A);

	// smat|r+i> = smat|r> + smat|i>
	// <r+i|smat|r+i> = <r+i|(smat|r> + smat|i>) 
	//                = real: <r|smat|r> + <i|smat|i> 
	//				  + imag: <r|smat|i> + <i|smat|r> : this term is zero for physical observables
	double* wf_r1 = new double[Dim];
	double* wf_i1 = new double[Dim];

	sparse_status_t status1 = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, wf_r, 0, wf_r1);
	if (0 != status1) { cout << "mkl_sparse_d_mv, sparse_status_t: " << status1 << endl; }
	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descrA, wf_i, 0, wf_i1);
	double aux = cblas_ddot(Dim, wf_r, 1, wf_r1, 1) + cblas_ddot(Dim, wf_i, 1, wf_i1, 1);

	mkl_sparse_destroy(A);
	delete[]wf_r1;
	delete[]wf_i1;

	return aux;
}

template <>
void Ham_HardCoreBoson<double>::Cal_Hkpz_bbr_Smat(my_int r, SparseMat<double>* _smat)
{
	// calculate all nonzero elements (maybe in a unprefered order), and then make a csr format

	my_int L = LatticeSize;
	double one = 1;

	// count nonzero elements 
	my_int index = 0;
	for (my_int s = 0; s < Dim; s++)
	{
		my_int state_s = basis_kpz->get_state(s);
		my_int na;
		if (s > 0 && basis_kpz->get_state(s - 1) == state_s) { continue; }
		else if (s < Dim - 1 && basis_kpz->get_state(s + 1) == state_s) { na = 2; }
		else { na = 1; }

		for (my_int sa = s; sa < s + na; sa++)
		{
			// set the first elements (which is diagonal) as nonzero
			index++;
			for (my_int i = 0; i < L; i++)
			{
				my_int j = (i + r) % L;
				my_int site_i, site_j;
				if (i < j) {
					site_i = i;
					site_j = j;
				}
				else
				{
					site_i = j;
					site_j = i;
				}
				// (0 == r): no extra nonzero element
				if (0 < r)
				{
					//if ((((state_s >> i) & 1)) ^ (((state_s >> j) & 1)))
					if (((state_s >> site_i) & 1) && (0 == ((state_s >> site_j) & 1)))
					{
						my_int state_s1b = state_s ^ (1 << i) ^ (1 << j);
						my_int state_s1;
						my_int nb, lb, qb, gb;
						basis_kpz->Representative(state_s1b, state_s1, lb, qb, gb);
						my_int s1 = basis_kpz->get_index(state_s1);
						if (-1 != s1)
						{
							if (s1 > 1 && basis_kpz->get_state(s1 - 1) == state_s1) { nb = 2; s1 -= 1; }
							else if (s1 < Dim - 1 && basis_kpz->get_state(s1 + 1) == state_s1) { nb = 2; }
							else { nb = 1; }
							my_int aux = -1;
							for (my_int sb = s1; sb < s1 + nb; sb++)
							{
								// if (sa == sb): no extra nonzero element
								if (sa < sb)
								{
									index++;
								}
							}
						}
					}
				}
			}
		}
	}
	my_int nnz = index;
	std::cout << "Cal_Hkpz_bbr_Smat(), r = " << r << ", nnz = " << nnz << endl;

	// construct sparse matrix
	_smat->SparseMat_Init(Dim, nnz);
	index = 0;
	for (my_int s = 0; s < Dim; s++)
	{
		my_int state_s = basis_kpz->get_state(s);
		my_int na;
		if (s > 0 && basis_kpz->get_state(s - 1) == state_s) { continue; }
		else if (s < Dim - 1 && basis_kpz->get_state(s + 1) == state_s) { na = 2; }
		else { na = 1; }

		for (my_int sa = s; sa < s + na; sa++)
		{
			// i = 0; set the first elements (with is diagonal) as nonzero
			_smat->vals[index] = 0;
			_smat->cols[index] = sa;
			_smat->Pointer_BE[s] = index;
			index++;
			for (my_int i = 0; i < L; i++)
			{
				my_int j = (i + r) % L;
				// j = i
				if (0 == r)
				{
					_smat->vals[index - 1] += ((state_s >> i) & 1) ? 1 : 0;
				}
				// j > i
				else
				{
					my_int site_i, site_j;
					if (i < j) {
						site_i = i;
						site_j = j;
					}
					else
					{
						site_i = j;
						site_j = i;
					}
					if (((state_s >> site_i) & 1) && (0 == ((state_s >> site_j) & 1)))
						//if ((((state_s >> i) & 1)) ^ (((state_s >> j) & 1)))
					{
						my_int state_s1b = state_s ^ (1 << i) ^ (1 << j);
						my_int state_s1;
						my_int nb, lb, qb, gb;
						basis_kpz->Representative(state_s1b, state_s1, lb, qb, gb);
						my_int s1 = basis_kpz->get_index(state_s1);
						if (-1 != s1)
						{
							if (s1 > 1 && basis_kpz->get_state(s1 - 1) == state_s1) { nb = 2; s1 -= 1; }
							else if (s1 < Dim - 1 && basis_kpz->get_state(s1 + 1) == state_s1) { nb = 2; }
							else { nb = 1; }
							my_int aux = -1;
							for (my_int sb = s1; sb < s1 + nb; sb++)
							{
								if (sa == sb)
								{
									_smat->vals[_smat->Pointer_BE[sa]] += Cal_Hkpz_element(sa, sb, lb, qb, gb, one);
								}
								else if (sa < sb)
								{
									_smat->vals[index] = Cal_Hkpz_element(sa, sb, lb, qb, gb, one);
									_smat->cols[index] = sb;
									index++;
								}
							}
						}
					}
				}
			}
			_smat->Pointer_BE[sa + 1] = index;
		}
	}
	std::cout << "Cal_Hkpz_bbr_Smat(), r = " << r << ", nnz = " << nnz << endl;
	/*
	for (my_int i = 0; i < nnz; i++)
	{
		cout << _smat->cols[i] << " " << _smat->vals[i] << endl;;
	}
	for (my_int i = 0; i < Dim; i++)
	{
		cout << _smat->Pointer_BE[i] << " ";
	}
	*/
}



template <>
void Ham_HardCoreBoson<double>::DenseMat_Static_kpz()
{
	my_int L = LatticeSize;

	my_int k = basis_kpz->get_k();
	my_int p = basis_kpz->get_p();
	my_int z = basis_kpz->get_z();

	// quntities from eigenvalues 
	char fspec[80];
	sprintf(fspec, "spec_k%d_p%d_z%d.bin", k, p, z);
	Vec_fwrite_double(fspec, spec, Dim);

	sprintf(fspec, "Hdiag_k%d_p%d_z%d.bin", k, p, z);
	Vec_fwrite_double(fspec, H_diag, Dim);

	Calc_Roag();
	sprintf(fspec, "roag_k%d_p%d_z%d.bin", k, p, z);
	Vec_fwrite_double(fspec, roag, Dim - 2);

	// quntities from eigenvectors

	//double* nk0 = new double[Dim];
	double* ipr = new double[Dim];
	double* SE = new double[Dim];
	double* HI = new double[Dim];
	double* nnr = new double[Dim * (L / 2 + 1)];
	double* nnr_s = new double[Dim * (L / 2 + 1)];
	for (my_int r = 0; r < L / 2 + 1; r++) { Cal_Hkpz_nnr_s(r, &nnr_s[r * Dim]); }

	double* wf = new double[Dim];
	double* wf_sq = new double[Dim];

	my_int time0 = time(0);
	for (my_int s = 0; s < Dim; s++)
	{
		Fetch_EigenVec(s, wf);
		vdSqr(Dim, wf, wf_sq);
		HI[s] = cblas_ddot(Dim, wf_sq, 1, H_diag, 1);
		ipr[s] = Cal_IPR_wfsq(wf_sq);
		SE[s] = Cal_SE_wfsq(wf_sq);
		for (my_int r = 0; r < L / 2 + 1; r++)
		{
			nnr[s * (L / 2 + 1) + r] = cblas_ddot(Dim, wf_sq, 1, &nnr_s[r * Dim], 1);
		}
		//Cal_Hkpz_nnr_wfsq(wf_sq, &nnr[s * (L / 2 + 1)]);
	}
	my_int time1 = time(0);
	sprintf(fspec, "eig_HI_k%d_p%d_z%d.bin", k, p, z);
	Vec_fwrite_double(fspec, HI, Dim);
	sprintf(fspec, "eig_ipr_k%d_p%d_z%d.bin", k, p, z);
	Vec_fwrite_double(fspec, ipr, Dim);
	sprintf(fspec, "eig_se_k%d_p%d_z%d.bin", k, p, z);
	Vec_fwrite_double(fspec, SE, Dim);
	sprintf(fspec, "eig_nnr_k%d_p%d_z%d.bin", k, p, z);
	Vec_fwrite_double(fspec, nnr, Dim * (L / 2 + 1));

	delete[]ipr;
	delete[]SE;
	delete[]HI;
	delete[]nnr;
	delete[]wf;
	delete[]wf_sq;
	std::cout << "Time for diagonal Obs: " << time1 - time0 << endl;

	// this part does not cost too much time
	//if (0 == k && 1 == p && 1 == z)
	if (1)
	{
		double* bbr = new double[Dim * (L / 2 + 1)];
		for (my_int r = 0; r < L / 2 + 1; r++)
		{
			SparseMat<double> smat_bbr;
			Cal_Hkpz_bbr_Smat(r, &smat_bbr);
			for (my_int s = 0; s < Dim; s++)
			{
				double* wfs = new double[Dim];
				Fetch_EigenVec(s, wfs);

				bbr[s * (L / 2 + 1) + r] = Cal_Ob_smat(wfs, &smat_bbr);
				delete[]wfs;
			}
			smat_bbr.SparseMat_Clear();
		}
		sprintf(fspec, "eig_bbr_k%d_p%d_z%d.bin", k, p, z);
		Vec_fwrite_double(fspec, bbr, Dim * (L / 2 + 1));
		delete[]bbr;

		my_int time2 = time(0);

		std::cout << "Time for offdiagonal Obs (b_i^\dagger b_{r+i}): " << time2 - time1 << endl;
	}

	/*
	// save all eigenvectors, in reordered basis (sorted by Hdiag)
	vector <pair<double, my_int> > Fock_E_n;  // stores E_fock, and index of this Fock state
	for (my_int s = 0; s < Dim; s++)
	{
		pair<double, my_int> aux = make_pair(H_diag[s], s);
		Fock_E_n.push_back(aux);
	}
	std::sort(Fock_E_n.begin(), Fock_E_n.end());

	double* mat = new double[Dim * Dim];
	for (my_int a = 0; a < Dim; a++)
	{
		my_int a1 = Fock_E_n[a].second;
		for (my_int b = 0; b < Dim; b++)
		{
			my_int b1 = Fock_E_n[b].second;
			mat[a * Dim + b] = DMat[a1 * Dim + b1];
		}
	}

	char fh[80];
	sprintf(fh, "EigMat_reorder_k%d_p%d_z%d.bin", k, p, z);
	Vec_fwrite_double(fh, mat, Dim * Dim);
	delete[]mat;

	sprintf(fh, "EigMat_k%d_p%d_z%d.bin", k, p, z);
	Vec_fwrite_double(fh, DMat, Dim * Dim);
	*/
}

template <>
void Ham_HardCoreBoson<double>::SparseMat_build_bibj_kpz(my_int i, my_int j, SparseMat<double>* _smat)
{
	// build sparse matrices for b_i^\dagger b_j
			// count nonzero elements 
	my_int nnz = 0;
	for (my_int s = 0; s < Dim; s++)
	{
		my_int state_s = basis_kpz->get_state(s);
		int na;
		if (s > 0 && basis_kpz->get_state(s - 1) == state_s) { continue; }
		else if (s < Dim - 1 && basis_kpz->get_state(s + 1) == state_s) { na = 2; }
		else { na = 1; }

		if ((((state_s >> i) & 1)) ^ (((state_s >> j) & 1)))
		{
			my_int state_s1b = state_s ^ (1 << i) ^ (1 << j);
			my_int state_s1;
			my_int nb, lb, qb, gb;
			basis_kpz->Representative(state_s1b, state_s1, lb, qb, gb);
			my_int s1 = basis_kpz->get_index(state_s1);
			if (-1 != s1)
			{
				if (s1 > 1 && basis_kpz->get_state(s1 - 1) == state_s1) { nb = 2; s1 -= 1; }
				else if (s1 < Dim - 1 && basis_kpz->get_state(s1 + 1) == state_s1) { nb = 2; }
				else { nb = 1; }
				for (my_int sa = s; sa < s + na; sa++)
					for (my_int sb = s1; sb < s1 + nb; sb++)
					{
						if (sa <= sb)
						{
							nnz++;
						}
					}
			}
		}
	}
	if (0 != nnz)
	{
		_smat->SparseMat_Init(Dim, nnz);

		// build sparse matrix
		double one = 1;
		double minus_1 = -1;
		my_int index = 0;
		for (my_int s = 0; s < Dim; s++)
		{
			my_int state_s = basis_kpz->get_state(s);

			my_int na;
			if (s > 0 && basis_kpz->get_state(s - 1) == state_s) { continue; }
			else if (s < Dim - 1 && basis_kpz->get_state(s + 1) == state_s) { na = 2; }
			else { na = 1; }

			if ((((state_s >> i) & 1)) ^ (((state_s >> j) & 1)))
			{
				my_int state_s1b = state_s ^ (1 << i) ^ (1 << j);
				my_int state_s1;
				my_int nb, lb, qb, gb;
				basis_kpz->Representative(state_s1b, state_s1, lb, qb, gb);
				my_int s1 = basis_kpz->get_index(state_s1);
				if (-1 != s1)
				{
					if (s1 > 1 && basis_kpz->get_state(s1 - 1) == state_s1) { nb = 2; s1 -= 1; }
					else if (s1 < Dim - 1 && basis_kpz->get_state(s1 + 1) == state_s1) { nb = 2; }
					else { nb = 1; }
					for (my_int sa = s; sa < s + na; sa++)
					{
						my_int aux = -1;
						for (my_int sb = s1; sb < s1 + nb; sb++)
						{
							if (sa <= sb)
							{
								_smat->vals[index] = Cal_Hkpz_element(sa, sb, lb, qb, gb, one);
								_smat->cols[index] = sb;
								aux++;

								if (0 == aux) { _smat->Pointer_BE[sa] = index; }
								index++;
							}
						}
					}
				}
			}
		}
	}
}

template <>
void Ham_HardCoreBoson<double>::Cal_Hkpz_bbr_wf(double* wf, double* br)
{
	my_int L = LatticeSize;
	double one = 1;

	double* wf1 = new double[Dim];

	for (my_int r = 0; r <= L / 2; r++)
	{
		for (my_int s = 0; s < Dim; s++) wf1[s] = 0;
		for (my_int s = 0; s < Dim; s++)
		{
			my_int state_s = basis_kpz->get_state(s);
			my_int na;
			if (s > 0 && basis_kpz->get_state(s - 1) == state_s) { continue; }
			else if (s < Dim - 1 && basis_kpz->get_state(s + 1) == state_s) { na = 2; }
			else { na = 1; }

			for (my_int i = 0; i < L; i++)
			{
				my_int j = (i + r) % L;
				// j = i
				if (0 == r && (state_s >> i) & 1)
				{
					for (my_int sa = s; sa < s + na; sa++)
					{
						wf1[sa] += wf[sa];
					}
				}
				// j > i
				else if ((((state_s >> i) & 1)) ^ (((state_s >> j) & 1)))
				{
					my_int state_s1b = state_s ^ (1 << i) ^ (1 << j);
					my_int state_s1;
					my_int nb, lb, qb, gb;
					basis_kpz->Representative(state_s1b, state_s1, lb, qb, gb);
					my_int s1 = basis_kpz->get_index(state_s1);
					if (-1 != s1)
					{
						if (s1 > 1 && basis_kpz->get_state(s1 - 1) == state_s1) { nb = 2; s1 -= 1; }
						else if (s1 < Dim - 1 && basis_kpz->get_state(s1 + 1) == state_s1) { nb = 2; }
						else { nb = 1; }
						for (my_int sa = s; sa < s + na; sa++)
						{
							my_int aux = -1;
							for (my_int sb = s1; sb < s1 + nb; sb++)
							{
								wf1[sb] += Cal_Hkpz_element(sa, sb, lb, qb, gb, one) * wf[sa];
							}
						}
					}
				}
			}
		}
		br[r] = cblas_ddot(Dim, wf, 1, wf1, 1);
	}
	delete[]wf1;
}

template <>
void Ham_HardCoreBoson<double>::Cal_Hkpz_bbr_wf(my_cmplx* wf, double* br)
{
	my_int L = LatticeSize;
	double one = 1;

	my_cmplx* wf1 = new my_cmplx[Dim];

	for (my_int r = 0; r <= L / 2; r++)
	{
		for (my_int s = 0; s < Dim; s++) wf1[s] = 0;
		for (my_int s = 0; s < Dim; s++)
		{
			my_int state_s = basis_kpz->get_state(s);
			my_int na;
			if (s > 0 && basis_kpz->get_state(s - 1) == state_s) { continue; }
			else if (s < Dim - 1 && basis_kpz->get_state(s + 1) == state_s) { na = 2; }
			else { na = 1; }

			for (my_int i = 0; i < L; i++)
			{
				my_int j = (i + r) % L;
				// j = i
				if (0 == r && (state_s >> i) & 1)
				{
					for (my_int sa = s; sa < s + na; sa++)
					{
						wf1[sa] += wf[sa];
					}
				}
				// j > i
				else if ((((state_s >> i) & 1)) ^ (((state_s >> j) & 1)))
				{
					my_int state_s1b = state_s ^ (1 << i) ^ (1 << j);
					my_int state_s1;
					my_int nb, lb, qb, gb;
					basis_kpz->Representative(state_s1b, state_s1, lb, qb, gb);
					my_int s1 = basis_kpz->get_index(state_s1);
					if (-1 != s1)
					{
						if (s1 > 1 && basis_kpz->get_state(s1 - 1) == state_s1) { nb = 2; s1 -= 1; }
						else if (s1 < Dim - 1 && basis_kpz->get_state(s1 + 1) == state_s1) { nb = 2; }
						else { nb = 1; }
						for (my_int sa = s; sa < s + na; sa++)
						{
							my_int aux = -1;
							for (my_int sb = s1; sb < s1 + nb; sb++)
							{
								wf1[sb] += Cal_Hkpz_element(sa, sb, lb, qb, gb, one) * wf[sa];
							}
						}
					}
				}
			}
		}
		//br[r] = cblas_ddot(Dim, wf, 1, wf1, 1);
		my_cmplx aux = my_vec_dot(Dim, wf, 1, wf1, 1);
		br[r] = aux.real();
	}
	delete[]wf1;
}