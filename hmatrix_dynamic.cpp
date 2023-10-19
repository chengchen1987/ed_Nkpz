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
#include <ctime>
using namespace std;

#include "basis.h" 
#include "hmatrix.h"
#include "sparsemat.h"

template <>
void Ham_HardCoreBoson<double>::Cal_PIr_wfsq(double* wfsq, double* HamD_pk, double* PI_r)
{
	for (my_int r = 0; r < LatticeSize + 1; r++)
	{
		PI_r[r] = 0;
		for (my_int k = 0; k < Dim; k++)
		{
			if (1e-10 > abs(r - HamD_pk[k]))
			{
				PI_r[r] += wfsq[k];
			}
		}
	}
}

template <>
void Ham_HardCoreBoson<double>::DenseMat_Evo_form_s(my_int id, my_int ind)
{
	my_int L = LatticeSize;
	my_int nt = Params.evo_time_steps;

	// compute Fidelity
	double* Fid_t = new double[nt];
	double* ipr_t = new double[nt];
	double* PI_r_t = new double[nt];
	double* EE_t = new double[nt];

	double* ni_t = new double[nt * L];
	double* ni_s = new double[Dim * L];
	for (my_int i = 0; i < L; i++) { Cal_H_ni_s(i, &ni_s[i * Dim]); }

	double* ninj_t = new double[nt * L * L];
	double* ninj_s = new double[Dim * L * L];
	for (my_int i = 0; i < L; i++)
		for (my_int j = 0; j < L; j++) { Cal_H_ninj_s(i, j, &ninj_s[(i * L + j) * Dim]); }

	double* bibj_t = new double[nt * L * L];
	SparseMat<double>* smat_bibj = new SparseMat<double>[L * L];
	for (my_int i = 0; i < L; i++)
		for (my_int j = 0; j < L; j++) { Cal_H_bibj_Smat(i, j, &smat_bibj[i * L + j]); }

	for (my_int it = 0; it < nt; it++) {
		double t = Params.time_vec[it];
		double* wf_r = new double[Dim];
		double* wf_i = new double[Dim];
		double* wf_sq = new double[Dim];
		DenseMat_Get_psi_t(ind, t, wf_r, wf_i);
		for (my_int s = 0; s < Dim; s++) { wf_sq[s] = wf_r[s] * wf_r[s] + wf_i[s] * wf_i[s]; }

		// compute Fidelity, ipr, nnr, bbr
		// Fidelity
		Fid_t[it] = wf_r[ind] * wf_r[ind] + wf_i[ind] * wf_i[ind];
		// ipr
		ipr_t[it] = Cal_IPR_wfsq(wf_sq);
		// EE
		my_cmplx ii(0, 1);
		my_cmplx* wft = new my_cmplx[Dim];
		for (my_int s = 0; s < Dim; s++) { wft[s] = wf_r[s] + ii * wf_i[s]; }
		my_int* sites = new my_int[L / 2];
		for (my_int i = 0; i < L / 2; i++) { sites[i] = i; }
		EE_t[it] = Cal_EE_wf(wft, L / 2, sites);
		// ni
		for (my_int i = 0; i < L; i++)
		{
			ni_t[it * L + i] = cblas_ddot(Dim, wf_sq, 1, &ni_s[i * Dim], 1);
		}
		// ninj
		for (my_int i = 0; i < L; i++)
		{
			ninj_t[it * L * L + i * L + i] = ni_t[it * L + i];
			for (my_int j = i + 1; j < L; j++)
			{
				ninj_t[it * L * L + i * L + j] = cblas_ddot(Dim, wf_sq, 1, &ninj_s[(i * L + j) * Dim], 1);
				ninj_t[it * L * L + j * L + i] = ninj_t[it * L * L + i * L + j];
			}
		}
		// bibj
		for (my_int i = 0; i < L; i++)
		{
			bibj_t[it * L * L + i * L + i] = ni_t[it * L + i];
			for (my_int j = i + 1; j < L; j++)
			{
				bibj_t[it * L * L + i * L + j] = Cal_Ob_smat(wf_r, wf_i, &smat_bibj[i * L + j]);
				bibj_t[it * L * L + j * L + i] = bibj_t[it * L * L + i * L + j];
			}
		}

		delete[]wf_sq;
		delete[]wf_r;
		delete[]wf_i;
	}

	char fout[80];
	sprintf(fout, "evo_ipr_id%d.bin", id);
	Vec_fwrite_double(fout, ipr_t, nt);
	sprintf(fout, "evo_Fid_id%d.bin", id);
	Vec_fwrite_double(fout, Fid_t, nt);
	sprintf(fout, "evo_EE_id%d.bin", id);
	Vec_fwrite_double(fout, EE_t, nt);

	sprintf(fout, "evo_ni_id%d.bin", id);
	Vec_fwrite_double(fout, ni_t, nt * L);
	sprintf(fout, "evo_ninj_id%d.bin", id);
	Vec_fwrite_double(fout, ninj_t, nt * L * L);
	sprintf(fout, "evo_bibj_id%d.bin", id);
	Vec_fwrite_double(fout, bibj_t, nt * L * L);

	delete[]ipr_t;
	delete[]Fid_t;
	delete[]EE_t;
	delete[]ni_t;
	delete[]ninj_t;
	delete[]bibj_t;

	for (my_int r = 0; r < L * L; r++)
	{
		smat_bibj[r].SparseMat_Clear();
	}
	delete[]ninj_s;
	delete[]ni_s;
}

template <>
void Ham_HardCoreBoson<double>::DenseMat_Evo_form_inds(my_int id, const my_int& ns, my_int* inds)
{
	double* wf0 = new double[Dim];
	for (my_int i = 0; i < Dim; i++) { wf0[i] = 0; }
	for (my_int i = 0; i < ns; i++) { wf0[inds[i]] = 1.0 / sqrt(ns); }

	my_int L = LatticeSize;
	my_int nt = Params.evo_time_steps;

	// compute Fidelity
	double* Fid_t = new double[nt];
	double* ipr_t = new double[nt];
	double* PI_r_t = new double[nt];
	double* EE_t = new double[nt];

	double* ni_t = new double[nt * L];
	double* ni_s = new double[Dim * L];
	for (my_int i = 0; i < L; i++) { Cal_H_ni_s(i, &ni_s[i * Dim]); }

	double* ninj_t = new double[nt * L * L];
	double* ninj_s = new double[Dim * L * L];
	for (my_int i = 0; i < L; i++)
		for (my_int j = 0; j < L; j++) { Cal_H_ninj_s(i, j, &ninj_s[(i * L + j) * Dim]); }

	double* bibj_t = new double[nt * L * L];
	SparseMat<double>* smat_bibj = new SparseMat<double>[L * L];
	for (my_int i = 0; i < L; i++)
		for (my_int j = 0; j < L; j++) { Cal_H_bibj_Smat(i, j, &smat_bibj[i * L + j]); }

	for (my_int it = 0; it < nt; it++) {
		double t = Params.time_vec[it];
		double* wf_r = new double[Dim];
		double* wf_i = new double[Dim];
		double* wf_sq = new double[Dim];
		DenseMat_Get_psi_t(inds, ns, t, wf_r, wf_i);

		for (my_int s = 0; s < Dim; s++) { wf_sq[s] = wf_r[s] * wf_r[s] + wf_i[s] * wf_i[s]; }

		// compute Fidelity, ipr, nnr, bbr
		// Fidelity
		Fid_t[it] = cblas_ddot(Dim, wf0, 1, wf_sq, 1) * sqrt(ns);
		// ipr
		ipr_t[it] = Cal_IPR_wfsq(wf_sq);
		// EE
		my_cmplx ii(0, 1);
		my_cmplx* wft = new my_cmplx[Dim];
		for (my_int s = 0; s < Dim; s++) { wft[s] = wf_r[s] + ii * wf_i[s]; }
		my_int* sites = new my_int[L / 2];
		for (my_int i = 0; i < L / 2; i++) { sites[i] = i; }
		EE_t[it] = Cal_EE_wf(wft, L / 2, sites);
		// ni
		for (my_int i = 0; i < L; i++)
		{
			ni_t[it * L + i] = cblas_ddot(Dim, wf_sq, 1, &ni_s[i * Dim], 1);
		}
		// ninj
		for (my_int i = 0; i < L; i++)
		{
			ninj_t[it * L * L + i * L + i] = ni_t[it * L + i];
			for (my_int j = i + 1; j < L; j++)
			{
				ninj_t[it * L * L + i * L + j] = cblas_ddot(Dim, wf_sq, 1, &ninj_s[(i * L + j) * Dim], 1);
				ninj_t[it * L * L + j * L + i] = ninj_t[it * L * L + i * L + j];
			}
		}
		// bibj
		for (my_int i = 0; i < L; i++)
		{
			bibj_t[it * L * L + i * L + i] = ni_t[it * L + i];
			for (my_int j = i + 1; j < L; j++)
			{
				bibj_t[it * L * L + i * L + j] = Cal_Ob_smat(wf_r, wf_i, &smat_bibj[i * L + j]);
				bibj_t[it * L * L + j * L + i] = bibj_t[it * L * L + i * L + j];
			}
		}

		delete[]wf_sq;
		delete[]wf_r;
		delete[]wf_i;
	}

	char fout[80];
	sprintf(fout, "evo_ipr_id%d_vec.bin", id);
	Vec_fwrite_double(fout, ipr_t, nt);
	sprintf(fout, "evo_Fid_id%d_vec.bin", id);
	Vec_fwrite_double(fout, Fid_t, nt);
	sprintf(fout, "evo_EE_id%d_vec.bin", id);
	Vec_fwrite_double(fout, EE_t, nt);

	sprintf(fout, "evo_ni_id%d_vec.bin", id);
	Vec_fwrite_double(fout, ni_t, nt * L);
	sprintf(fout, "evo_ninj_id%d_vec.bin", id);
	Vec_fwrite_double(fout, ninj_t, nt * L * L);
	sprintf(fout, "evo_bibj_id%d_vec.bin", id);
	Vec_fwrite_double(fout, bibj_t, nt * L * L);

	delete[]ipr_t;
	delete[]Fid_t;
	delete[]EE_t;
	delete[]ni_t;
	delete[]ninj_t;
	delete[]bibj_t;

	for (my_int r = 0; r < L * L; r++)
	{
		smat_bibj[r].SparseMat_Clear();
	}
	delete[]ninj_s;
	delete[]ni_s;
}


template <>
void Ham_HardCoreBoson<double>::DenseMat_Evo_form_wf0(my_int id, double* wf0)
{
	my_int L = LatticeSize;
	my_int nt = Params.evo_time_steps;

	// compute Fidelity
	double* Fid_t = new double[nt];
	double* ipr_t = new double[nt];
	double* PI_r_t = new double[nt];
	double* EE_t = new double[nt];

	double* ni_t = new double[nt * L];
	double* ni_s = new double[Dim * L];
	for (my_int i = 0; i < L; i++) { Cal_H_ni_s(i, &ni_s[i * Dim]); }

	double* ninj_t = new double[nt * L * L];
	double* ninj_s = new double[Dim * L * L];
	for (my_int i = 0; i < L; i++)
		for (my_int j = 0; j < L; j++) { Cal_H_ninj_s(i, j, &ninj_s[(i * L + j) * Dim]); }

	double* bibj_t = new double[nt * L * L];
	SparseMat<double>* smat_bibj = new SparseMat<double>[L * L];
	for (my_int i = 0; i < L; i++)
		for (my_int j = 0; j < L; j++) { Cal_H_bibj_Smat(i, j, &smat_bibj[i * L + j]); }


	for (my_int it = 0; it < nt; it++) {
		double t = Params.time_vec[it];
		double* wf_r = new double[Dim];
		double* wf_i = new double[Dim];
		double* wf_sq = new double[Dim];
		//DenseMat_Get_psi_t(ind, t, wf_r, wf_i);
		DenseMat_Get_psi_t(wf0, t, wf_r, wf_i);
		for (my_int s = 0; s < Dim; s++) { wf_sq[s] = wf_r[s] * wf_r[s] + wf_i[s] * wf_i[s]; }

		// compute Fidelity, ipr, nnr, bbr
		// Fidelity
		Fid_t[it] = Fid_t[it] = cblas_ddot(Dim, wf0, 1, wf_sq, 1);
		// ipr
		ipr_t[it] = Cal_IPR_wfsq(wf_sq);
		// ni
		for (my_int i = 0; i < L; i++)
		{
			ni_t[it * L + i] = cblas_ddot(Dim, wf_sq, 1, &ni_s[i * Dim], 1);
		}
		// EE
		my_cmplx ii(0, 1);
		my_cmplx* wft = new my_cmplx[Dim];
		for (my_int s = 0; s < Dim; s++) { wft[s] = wf_r[s] + ii * wf_i[s]; }
		my_int* sites = new my_int[L / 2];
		for (my_int i = 0; i < L / 2; i++) { sites[i] = i; }
		EE_t[it] = Cal_EE_wf(wft, L / 2, sites);
		delete[]wft;
		// ninj
		for (my_int i = 0; i < L; i++)
		{
			ninj_t[it * L * L + i * L + i] = ni_t[it * L + i];
			for (my_int j = i + 1; j < L; j++)
			{
				ninj_t[it * L * L + i * L + j] = cblas_ddot(Dim, wf_sq, 1, &ninj_s[(i * L + j) * Dim], 1);
				ninj_t[it * L * L + j * L + i] = ninj_t[it * L * L + i * L + j];
			}
		}
		// bibj
		for (my_int i = 0; i < L; i++)
		{
			bibj_t[it * L * L + i * L + i] = ni_t[it * L + i];
			for (my_int j = i + 1; j < L; j++)
			{
				bibj_t[it * L * L + i * L + j] = Cal_Ob_smat(wf_r, wf_i, &smat_bibj[i * L + j]);
				bibj_t[it * L * L + j * L + i] = bibj_t[it * L * L + i * L + j];
			}
		}

		delete[]wf_sq;
		delete[]wf_r;
		delete[]wf_i;
	}

	char fout[80];
	sprintf(fout, "evo_ipr_id%d_span.bin", id);
	Vec_fwrite_double(fout, ipr_t, nt);
	sprintf(fout, "evo_Fid_id%d_span.bin", id);
	Vec_fwrite_double(fout, Fid_t, nt);
	sprintf(fout, "evo_EE_id%d_span.bin", id);
	Vec_fwrite_double(fout, EE_t, nt);

	sprintf(fout, "evo_ni_id%d_span.bin", id);
	Vec_fwrite_double(fout, ni_t, nt * L);
	sprintf(fout, "evo_ninj_id%d_span.bin", id);
	Vec_fwrite_double(fout, ninj_t, nt * L * L);
	sprintf(fout, "evo_bibj_id%d_span.bin", id);
	Vec_fwrite_double(fout, bibj_t, nt * L * L);

	delete[]ipr_t;
	delete[]Fid_t;
	delete[]EE_t;
	delete[]ni_t;
	delete[]ninj_t;
	delete[]bibj_t;

	for (my_int r = 0; r < L * L; r++)
	{
		smat_bibj[r].SparseMat_Clear();
	}
	delete[]ninj_s;
	delete[]ni_s;
}

template<>
void Ham_HardCoreBoson<double>::DenseMat_Dynamic()
{
	my_int L = LatticeSize;

	Params.GetTimeVec();
	char ft[80];
	sprintf(ft, "evo_tvec.bin");
	Vec_fwrite_double(ft, Params.time_vec, Params.evo_time_steps);

	vector<my_int> ind_vec = Get_IniStates();
	my_int n_ini = ind_vec.size();

	char fout[80];
	sprintf(fout, "evo_iniinfo.dat");
	ofstream ofini(fout);
	for (my_int ip = 0; ip < n_ini; ip++)
	{
		ofini << ip << "  " << ind_vec[ip] << "  " << basis->get_state(ind_vec[ip])
			<< "  " << H_diag[ind_vec[ip]] << "  " << (H_diag[ind_vec[ip]] - spec[0]) / (spec[Dim - 1] - spec[0]) << endl;
	}
	ofini.close();

	// evolution of initial states
	//#pragma omp parallel for schedule(dynamic)
	cout << "Evolutions: " << endl;
	for (my_int ip = 0; ip < n_ini; ip++)
	{
		my_int ind = ind_vec[ip];
		cout << ip << "	" << ind_vec[ip] << "	 " << basis->get_state(ind_vec[ip])
			<< "  " << H_diag[ind_vec[ip]] << "  " << (H_diag[ind_vec[ip]] - spec[0]) / (spec[Dim - 1] - spec[0]) << endl;
		DenseMat_Evo_form_s(ip, ind);

		// long-time limit average: O(t->infty) -> \sum_a c_a^2 <a|O|a>
		// <a|O|a> is computed in DenseMat_Static and printed to file
		// here we print c_a^2, and get O(t->infty) by postscripts
		double* Csq_p_alpha = new double[Dim];	// <p|alpha>, overlaps of initial and eigen states
		Cal_Csq_s_alpha(ind, Csq_p_alpha);
		char fname[80];
		sprintf(fname, "evo_csq_id%d.bin", ip);
		Vec_fwrite_double(fname, Csq_p_alpha, Dim);
		delete[]Csq_p_alpha;
	}

	// evolution of the combined states (spanned by kpz symetry operations)
	my_int _k, _p, _z;
	_k = 0; _p = 1; _z = 1;
	Basis_kpz _basis_kpz(L, basis->get_nup(), _k, _p, _z);
	//#pragma omp parallel for schedule(dynamic)
	for (my_int ip = 0; ip < n_ini; ip++)
	{
		my_int ind = ind_vec[ip];
		my_int s0 = basis->get_state(ind);

		vector<my_int> inds = Get_Spanned_States(s0);
		my_int ns = inds.size();

		DenseMat_Evo_form_inds(ip, ns, &inds[0]);
	}
}