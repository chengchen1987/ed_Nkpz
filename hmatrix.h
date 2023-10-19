#pragma once

#ifndef    HMATRIX_H
#define    HMATRIX_H

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
#include <random>
using namespace std;

#include "basis.h"
#include "sparsemat.h"

class Parameters
{
public:
	// Basis
	my_int LatticeSize, N_up;
	// Hamiltonian parameters 
	double t, V, alpha;
	my_int bc;

	// choose a model, do it later
	void SetDefault();
	void Read_input();

	// evolution
//	double evo_time_min;
	double evo_time_max;
	my_int evo_time_steps;
	my_int evo_time_IncStyle;	// 0/1 for linear/exponetial time increase
	my_int evo_nini;
	double* time_vec;
	void GetTimeVec();

	my_int evo_type; // 0/1 for Lanczos/RK4 evolution
	my_int evo_Lanczos_M;

	char evo_ni;
	char evo_ninj;
	char evo_IPR;
	char evo_SE;
	char evo_PI_r;
	char evo_EE;
	char evo_Pmax;
	// entanglement entropy <- information of subsystem:
	my_int EE_size_A;
	my_int* EE_sites_A;
};

class Dyn_DataStruc
{
public:
	void Initialize(Parameters& para, const my_int& _dim, const my_int& _prt_ind);
	void ReleaseSpace(Parameters& para);
	void PrintDynResults(Parameters& para);

	my_int prt_ind;

	my_int dim;

	double* ni_t;
	double* ninj_t;
	double* IPR_t;
	double* PI_r_t;
	double* HamD_k; // used to accelerate calculating Pi_r_t
	double* SE_t;
	double* EE_t;
};

template <class dtype>
class Ham_HardCoreBoson {
private:
	Basis* basis;
	Basis_k* basis_k;
	Basis_kp* basis_kp;
	Basis_kpz* basis_kpz;
public:
	~Ham_HardCoreBoson();
	Ham_HardCoreBoson(Basis* _basis);
	Ham_HardCoreBoson(Basis_k* _basis);
	Ham_HardCoreBoson(Basis_kp* _basis);
	Ham_HardCoreBoson(Basis_kpz* _basis);

	Parameters Params;

	my_int LatticeSize;

	my_int Dim;
	// evolution parameters
	double dt;
	my_int nt;
	// lattice 
	my_int N_hops;
	std::vector<std::tuple<my_int, my_int, double> > hops;
	void tV_dipolar_MakeLatticeHops(std::vector<std::tuple<my_int, my_int, double> >& hops);
	my_int N_ints;
	std::vector<std::tuple<my_int, my_int, double> > ints;
	void tV_dipolar_MakeLatticeInts(std::vector<std::tuple<my_int, my_int, double> >& ints);

	// dense matrix ==================================================== 
	void DenseMat_Display();
	void DenseMat_Check_Hermitian();
	void DenseMat_Diag_Compare_Disorder_LRI();
	// 
	void DenseMat_Build();
	void DenseMat_Build_k();
	dtype Cal_Hkp_element(my_int s, my_int s1, my_int l, my_int q, double& hj);
	void DenseMat_Build_kp();
	dtype Cal_Hkpz_element(my_int s, my_int s1, my_int l, my_int q, my_int g, double& hj);
	void DenseMat_Build_kpz();

	// sparse matrix
	void SparseMat_build_bibj_kpz(my_int i, my_int j, SparseMat<double>* _smat);

	dtype* DMat;
	double* H_diag;
	double* spec;
	double* roag;

	void DenseMat_Eig();
	void Calc_Roag();

	void DenseMat_Static();
	void DenseMat_Dynamic();
	void DenseMat_Evo_form_s(my_int id, my_int s);
	void DenseMat_Evo_form_inds(my_int id, const my_int& ns, my_int* inds);
	void DenseMat_Evo_form_wf0(my_int id, double* wf0);
	void Cal_Csq_s_alpha(const my_int& s, double* Csq_s_alpha);

	vector<my_int> Get_IniStates();
	vector<my_int> Get_Spanned_States(const my_int& s0);
	vector<my_int> Get_IniStates_kpz();

	void DenseMat_Static_kpz();
	void DenseMat_Dynamic_kpz();
	void DenseMat_Evo_form_s_kpz(my_int id, my_int s);
	void Cal_H_ni_s(my_int i, double* ni_s);
	void Cal_H_ninj_s(my_int i, my_int j, double* ninj_s);
	void Cal_H_bibj_Smat(my_int i, my_int j, SparseMat<double>* _smat);

	void DenseMat_Get_psi_t(const my_int& p, const double& t, double* wf_r, double* wf_i);
	void DenseMat_Get_psi_t(my_int* inds, const my_int& ns, const double& t, double* wf_r, double* wf_i);
	void DenseMat_Get_psi_t(double* wf0, const double& t, double* wf_r, double* wf_i);
	// hmaitrx, obs from wavefunctions 
	void Fetch_EigenVec(const my_int& s, dtype* wf);

	//void Cal_ni_wfsq(double* wfsq, double* ni);
	//void Cal_ninj_wfsq(double* wfsq, double* ninj);
	void Cal_PIr_wfsq(double* wfsq, double* HamD_pk, double* PI_r);
	double Cal_SE_wfsq(double* wfsq);
	double Cal_IPR_wfsq(double* wfsq);
	//double Cal_EE_wf(double* wf, const my_int& size_A, my_int& sites_A);
	double Cal_EE_wf(my_cmplx* wf, const my_int& size_A, my_int* sites_A);
	double Cal_halfchain_EE_wf(double* wf);
	double Cal_halfchain_EE_wf(my_cmplx* wf);

	void Cal_Hkpz_nnr_s(my_int r, double* nnr_s);
	void Cal_Hkpz_nnr_wfsq(double* wfsq, double* nnr);
	void Cal_Hkpz_bbr_wf(double* wf, double* bbr);
	void Cal_Hkpz_bbr_wf(my_cmplx* wf, double* bbr);
	void Cal_Hkpz_bbr_Smat(my_int r, SparseMat<double>* _smat);
	double Cal_Ob_smat(double* wf, SparseMat<double>* _smat);
	double Cal_Ob_smat(double* wf_r, double* wf_i, SparseMat<double>* _smat);

	//double Cal_EE_wf(double* wf, const my_int& size_A, my_int* sites_A);
	//double Cal_EE_wf(my_cmplx* wf, const my_int& size_A, my_int* sites_A);
	void Cal_product_spectrum(double* pspec);
	my_int Get_TargetFock_left(const double& Target_E, vector <pair<double, my_int> >& Fock_E_n);
	my_int Get_TargetFock_right(const double& Target_E, vector <pair<double, my_int> >& Fock_E_n);
};

template<class dtype>
Ham_HardCoreBoson<dtype>::Ham_HardCoreBoson(Basis* _basis) :
	basis(_basis),
	LatticeSize(basis->get_L()),
	Dim(basis->get_Dim())
{
	cout << "Dim = " << Dim << endl;

	Params.SetDefault();
	Params.Read_input();

	tV_dipolar_MakeLatticeInts(ints);
	tV_dipolar_MakeLatticeHops(hops);
}

template<class dtype>
Ham_HardCoreBoson<dtype>::~Ham_HardCoreBoson()
{
	delete[]DMat;
	delete[]spec;
	delete[]roag;
	delete[]H_diag;
}

template<class dtype>
void Ham_HardCoreBoson<dtype>::tV_dipolar_MakeLatticeHops(std::vector<std::tuple<my_int, my_int, double> >& hops)
{
	// note: keep i < j
	N_hops = Params.bc ? LatticeSize - 1 : LatticeSize;
	hops.reserve(N_hops);
	cout << "N_hops = " << N_hops << endl;
	//
	for (my_int i = 0; i < LatticeSize - 1; i++)
	{
		my_int j = i + 1;
		auto bar = std::make_tuple(i, j, -Params.t);
		hops.push_back(bar);
	}
	if (0 == Params.bc)
	{
		auto bar = std::make_tuple(0, LatticeSize - 1, -Params.t);
		hops.push_back(bar);
	}
	/*
	cout << "hops:" << endl;
	for (my_int ib = 0; ib < N_hops; ib++)
	{
		cout << setw(6) << std::get<0>(hops[ib]) << setw(6) << std::get<1>(hops[ib]) << setw(24) << std::get<2>(hops[ib]) << endl;
	}
	cout << endl;
	*/
}

template<class dtype>
void Ham_HardCoreBoson<dtype>::tV_dipolar_MakeLatticeInts(std::vector<std::tuple<my_int, my_int, double> >& my_ints)
{
	N_ints = LatticeSize * (LatticeSize - 1) / 2;
	ints.reserve(N_ints);
	cout << "N_ints = " << N_ints << endl;
	for (my_int i = 0; i < LatticeSize; i++)
		for (my_int j = i + 1; j < LatticeSize; j++)
		{
			my_int r = Params.bc ? j - i : min(j - i, i + LatticeSize - j);
			auto bar = std::make_tuple(i, j, Params.V * pow(r, -Params.alpha));
			ints.push_back(bar);
		}
	/*
	cout << "ints:" << endl;
	for (my_int ib = 0; ib < N_ints; ib++)
	{
		cout << setw(6) << std::get<0>(my_ints[ib]) << setw(6) << std::get<1>(my_ints[ib]) << setw(24) << std::get<2>(my_ints[ib]) << endl;
	}
	cout << endl;
	*/
}

template<class dtype>
void Ham_HardCoreBoson<dtype>::Cal_product_spectrum(double* pspec)
{
}

template<class dtype>
void Ham_HardCoreBoson<dtype>::DenseMat_Diag_Compare_Disorder_LRI()
{
	// LRI
	double* HI_LRI = new double[Dim];
	for (my_int s = 0; s < Dim; s++) HI_LRI[s] = 0;

	for (my_int s = 0; s < Dim; s++) 
	{
		my_int state_s;
		state_s = basis->get_state(s);
		for (my_int ib = 0; ib < N_ints; ib++)
		{
			my_int i = std::get<0>(ints[ib]);
			my_int j = std::get<1>(ints[ib]);
			my_int s_i = ((state_s >> i) & 1);
			my_int s_j = ((state_s >> j) & 1);
			if (s_i & s_j) HI_LRI[s] += std::get<2>(ints[ib]);
		}
	}
	char ft[80];
	sprintf(ft, "HI_V.bin");
	Vec_fwrite_double(ft, HI_LRI, Dim);
	delete[]HI_LRI;

	// disorder
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> ran_double_01(0, 1);
	my_int samples = 100;
	for (my_int isam = 0; isam < samples; isam++)
	{
		// onsite disorder 
		double* dis = new double[LatticeSize];
		for (my_int i = 0; i < LatticeSize; i++)
		{
			dis[i] = ran_double_01(gen);
		}

		double* HI_disorder = new double[Dim];
		for (my_int s = 0; s < Dim; s++) HI_disorder[s] = 0;

		for (my_int s = 0; s < Dim; s++) 
		{
			my_int state_s;
			state_s = basis->get_state(s);
			for (my_int i = 0; i < LatticeSize; i++)
			{
				HI_disorder[s] += dis[i] * ((state_s >> i) & 1);
			}
		}
		char ft[80];
		sprintf(ft, "HI_dis_r%d.bin", isam);
		Vec_fwrite_double(ft, HI_disorder, Dim);
		delete[]HI_disorder;
		delete[]dis;
	}
}

template<class dtype>
void Ham_HardCoreBoson<dtype>::DenseMat_Build()
{
	my_int matlen = Dim * Dim;
	cout << "Dim*Dim = " << matlen << endl;
	cout << "Estimated memory cost of the dense matrix: " << matlen * 8 / 1e9 << " GB" << endl;

	DMat = new double[matlen];
	for (my_int s = 0; s < matlen; s++) DMat[s] = 0;
	cout << "Allocate space for DMat..." << endl;
	H_diag = new double[Dim];
	//omp_set_num_threads(NUMTHREADS);
 //   #pragma omp parallel for schedule(dynamic)
	/*-------------------------------------------------------------------
	  H = -t*sum_{i < j}[b_i^+ b_j^- + b_i^- b_j^+]
	  + V * sum_{i < j} (i,r) r^(-3) n_i n_{i+r}
	  --------------------------------------------------------------------*/
	for (my_int s = 0; s < Dim; s++) {
		my_int state_s;
		state_s = basis->get_state(s);

		// diagoanl term 
		// V * sum_{i < j} (i,r) r^(-3) n_i n_{i+r}
		my_int diag_index = s + s * Dim;
		for (my_int ib = 0; ib < N_ints; ib++)
		{
			my_int i = std::get<0>(ints[ib]);
			my_int j = std::get<1>(ints[ib]);
			my_int s_i = ((state_s >> i) & 1);
			my_int s_j = ((state_s >> j) & 1);
			if (s_i & s_j) DMat[diag_index] += std::get<2>(ints[ib]);
		}
		H_diag[s] = DMat[diag_index];
		// offdiagonal term
		// -t*sum_{i < j}[b_i^+ b_j^- + b_i^- b_j^+]
		for (my_int ib = 0; ib < N_hops; ib++)
		{
			my_int i = std::get<0>(hops[ib]);
			my_int j = std::get<1>(hops[ib]);
			my_int s_i = ((state_s >> i) & 1);
			my_int s_j = ((state_s >> j) & 1);
			if ((s_i) && (!s_j))
			{
				my_int state_s1 = state_s ^ (1 << i) ^ (1 << j);
				my_int s1 = basis->get_index(state_s1);
				DMat[s * Dim + s1] += std::get<2>(hops[ib]);
				DMat[s1 * Dim + s] += std::get<2>(hops[ib]);
			}
		}
	}

	//DenseMat_Display();
}


template <class dtype>
Ham_HardCoreBoson<dtype>::Ham_HardCoreBoson(Basis_k* _basis) :
	basis_k(_basis),
	LatticeSize(basis_k->get_L()),
	Dim(basis_k->get_Dim())
{
	cout << "k = " << basis_k->get_k() << endl;
	cout << "Dim = " << Dim << endl;

	Params.SetDefault();
	Params.Read_input();

	if (1 == Params.bc)
	{
		cout << "bc must be 0 for momentum states! Set bc to 0 for PBC." << endl;
		Params.bc = 0;
	}
	//
	tV_dipolar_MakeLatticeInts(ints);
	tV_dipolar_MakeLatticeHops(hops);
}

template <class dtype>
Ham_HardCoreBoson<dtype>::Ham_HardCoreBoson(Basis_kp* _basis) :
	basis_kp(_basis),
	LatticeSize(basis_kp->get_L()),
	Dim(basis_kp->get_Dim())
{
	cout << "k = " << basis_kp->get_k() << endl;
	cout << "p = " << basis_kp->get_p() << endl;
	cout << "Dim = " << Dim << endl;

	Params.SetDefault();
	Params.Read_input();

	if (1 == Params.bc)
	{
		cout << "bc must be 0 for momentum states! Set bc to 0 for PBC." << endl;
		Params.bc = 0;
	}
	//
	tV_dipolar_MakeLatticeInts(ints);
	tV_dipolar_MakeLatticeHops(hops);
}

template <class dtype>
Ham_HardCoreBoson<dtype>::Ham_HardCoreBoson(Basis_kpz* _basis) :
	basis_kpz(_basis),
	LatticeSize(basis_kpz->get_L()),
	Dim(basis_kpz->get_Dim())
{
	cout << "k = " << basis_kpz->get_k() << endl;
	cout << "p = " << basis_kpz->get_p() << endl;
	cout << "z = " << basis_kpz->get_z() << endl;
	cout << "Dim = " << Dim << endl;

	Params.SetDefault();
	Params.Read_input();

	if (1 == Params.bc)
	{
		cout << "bc must be 0 for momentum states! Set bc to 0 for PBC." << endl;
		Params.bc = 0;
	}
	//
	tV_dipolar_MakeLatticeInts(ints);
	tV_dipolar_MakeLatticeHops(hops);
}

template <class dtype>
void Ham_HardCoreBoson<dtype>::DenseMat_Display()
{
	cout << "Dense matrix of H:" << endl;
	for (my_int i = 0; i < Dim; i++)
	{
		for (my_int j = 0; j < Dim; j++)
		{
			cout << setw(8) << setprecision(4) << DMat[i * Dim + j];
		}
		cout << endl;
	}
}


template <class dtype>
void Ham_HardCoreBoson<dtype>::DenseMat_Eig()
{
	spec = new double[Dim];
	// 'V'/'N' with/without eigenvectors
	char job_vec = 'V';
	//DenseMatrixEigenSolver(LAPACK_ROW_MAJOR, job_vec, 'U', Dim, DMat, Dim, spec);
	DenseMatrixEigenSolver_FInterface(LAPACK_ROW_MAJOR, job_vec, 'U', Dim, DMat, Dim, spec);

	cout << "spec[0] of DMat: " << setprecision(14) << spec[0] << endl;
}

template <class dtype>
my_int Ham_HardCoreBoson<dtype>::Get_TargetFock_left(const double& Target_E, vector <pair<double, my_int> >& Fock_E_n) {
	for (my_int p = 0; p < Dim; p++) {
		if (Fock_E_n[p].first > Target_E) {
			return p;
		}
	}
	return -1;
}

template <class dtype>
my_int Ham_HardCoreBoson<dtype>::Get_TargetFock_right(const double& Target_E, vector <pair<double, my_int> >& Fock_E_n) {
	for (my_int p = Dim - 1; p > 0; p--) {
		if (Fock_E_n[p].first < Target_E) {
			return p;
		}
	}
	return -1;
}

template <class dtype>
void Ham_HardCoreBoson<dtype>::Cal_Csq_s_alpha(const my_int& s, double* Csq_s_alpha) {
	//vdMul(Dim, &DMat[s * Dim], &DMat[s * Dim], Csq_s_alpha);
	for (my_int a = 0; a < Dim; a++)
	{
		Csq_s_alpha[a] = DMat[a * Dim + s] * DMat[a * Dim + s];
	}
}

#include "hmatrix_wf_obs.h"
#include "hmatrix_utils.h"

#endif