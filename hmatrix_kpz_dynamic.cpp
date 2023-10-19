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
void Ham_HardCoreBoson<double>::DenseMat_Evo_form_s_kpz(my_int id, my_int ind)
{
	my_int L = LatticeSize;
	my_int nt = Params.evo_time_steps;
	// 
	my_int k = basis_kpz->get_k();
	my_int p = basis_kpz->get_p();
	my_int z = basis_kpz->get_z();

	// compute Fidelity, bbr, nnr
	double* Fid_t = new double[nt];
	double* ipr_t = new double[nt];

	double* nnr_t = new double[nt * (L / 2 + 1)];
	double* nnr_s = new double[Dim * (L / 2 + 1)];
	for (my_int r = 0; r < L / 2 + 1; r++) { Cal_Hkpz_nnr_s(r, &nnr_s[r * Dim]); }

	double* bbr_t = new double[nt * (L / 2 + 1)];
	SparseMat<double>* smat_bbr = new SparseMat<double>[L / 2 + 1];
	for (my_int r = 0; r < L / 2 + 1; r++)
	{
		Cal_Hkpz_bbr_Smat(r, &smat_bbr[r]);
	}

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
		// nnr
		for (my_int r = 0; r < L / 2 + 1; r++)
		{
			nnr_t[it * (L / 2 + 1) + r] = cblas_ddot(Dim, wf_sq, 1, &nnr_s[r * Dim], 1);
		}
		//Cal_Hkpz_nnr_wfsq(wf_sq, &nnr_t[it * (L / 2 + 1)]);
		// bbr
		for (my_int r = 0; r < L / 2 + 1; r++)
		{
			bbr_t[it * (L / 2 + 1) + r] = Cal_Ob_smat(wf_r, wf_i, &smat_bbr[r]);
		}
		delete[]wf_sq;
		delete[]wf_r;
		delete[]wf_i;
	}

	char fout[80];
	sprintf(fout, "evo_ipr_k%d_p%d_z%d_id%d.bin", k, p, z, id);
	Vec_fwrite_double(fout, ipr_t, nt);
	sprintf(fout, "evo_Fid_k%d_p%d_z%d_id%d.bin", k, p, z, id);
	Vec_fwrite_double(fout, Fid_t, nt);

	sprintf(fout, "evo_nnr_k%d_p%d_z%d_id%d.bin", k, p, z, id);
	Vec_fwrite_double(fout, nnr_t, nt * (L / 2 + 1));
	sprintf(fout, "evo_bbr_k%d_p%d_z%d_id%d.bin", k, p, z, id);
	Vec_fwrite_double(fout, bbr_t, nt * (L / 2 + 1));

	delete[]ipr_t;
	delete[]Fid_t;
	delete[]nnr_t;
	delete[]nnr_s;
	delete[]bbr_t;
	for (my_int r = 0; r < L / 2 + 1; r++)
	{
		smat_bbr[r].SparseMat_Clear();
	}
	delete[]smat_bbr;
}

template <>
void Ham_HardCoreBoson<double>::DenseMat_Dynamic_kpz()
{
	my_int L = LatticeSize;
	// 
	my_int k = basis_kpz->get_k();
	my_int p = basis_kpz->get_p();
	my_int z = basis_kpz->get_z();

	Params.GetTimeVec();
	char ft[80];
	sprintf(ft, "evo_tvec_k%d_p%d_z%d.bin", k, p, z);
	Vec_fwrite_double(ft, Params.time_vec, Params.evo_time_steps);
	
	vector<my_int> ind_vec = Get_IniStates_kpz();
	my_int n_ini = ind_vec.size();
	// trans to its represatatives
	char fout[80];
	sprintf(fout, "evo_iniinfo_k%d_p%d_z%d.dat", k, p, z);
	ofstream ofini(fout);
	for (my_int ip = 0; ip < n_ini; ip++)
	{
		ofini << ip << "	" << ind_vec[ip] << "	 " << basis_kpz->get_state(ind_vec[ip])
			<< "  " << H_diag[ind_vec[ip]] << "  " << (H_diag[ind_vec[ip]] - spec[0]) / (spec[Dim - 1] - spec[0]) << endl;
	}
	ofini.close();
	

#pragma omp parallel for schedule(dynamic)
	for (my_int ip = 0; ip < n_ini; ip++)
	{
		my_int ind = ind_vec[ip];
		cout << ip << "	" << ind_vec[ip] << "	 " << basis_kpz->get_state(ind_vec[ip])
			<< "  " << H_diag[ind_vec[ip]] << "  " << (H_diag[ind_vec[ip]] - spec[0]) / (spec[Dim - 1] - spec[0]) << endl;
		DenseMat_Evo_form_s_kpz(ip, ind);
	}

	ind_vec.clear();
}