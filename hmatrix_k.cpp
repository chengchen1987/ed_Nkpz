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

void Ham_HardCoreBoson<double>::DenseMat_Build_k()
{
	my_int matlen = (my_int)Dim * (my_int)Dim;
	cout << "Dim*Dim = " << matlen << endl;
	cout << "Estimated memory cost of the dense matrix: " << matlen * 8 / 1e9 << " GB" << endl;

	double kk = basis_k->get_k() * 2 * PI / LatticeSize;

	DMat = new double[matlen];
	for (my_int s = 0; s < matlen; s++) DMat[s] = 0;
	//omp_set_num_threads(NUMTHREADS);
 //   #pragma omp parallel for schedule(dynamic)
	/*-------------------------------------------------------------------
	  H = -t*sum_{i < j}[b_i^+ b_j^- + b_i^- b_j^+]
	  + V * sum_{i < j} (i,r) r^(-3) n_i n_{i+r}
	  --------------------------------------------------------------------*/
	for (my_int s = 0; s < Dim; s++) {
		my_int state_s;
		state_s = basis_k->get_state(s);
		double Ra = basis_k->get_Rs(s);

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
		// offdiagonal term
		// -t*sum_{i < j}[b_i^+ b_j^- + b_i^- b_j^+]
		for (my_int ib = 0; ib < N_hops; ib++)
		{
			my_int i = std::get<0>(hops[ib]);
			my_int j = std::get<1>(hops[ib]);
			my_int s_i = ((state_s >> i) & 1);
			my_int s_j = ((state_s >> j) & 1);
			if ((s_i) ^ (s_j))
			//if ((s_i) && (!s_j))
			{
				my_int state_s1b = state_s ^ (1 << i) ^ (1 << j);
				my_int state_s1;
				my_int lb;
				basis_k->Representative(state_s1b, state_s1, lb);
				double Rb = basis_k->CheckState_Ra(state_s1);
				if (-1 != Rb) {
					my_int s1 = basis_k->get_index(state_s1);
					DMat[s * Dim + s1] += std::get<2>(hops[ib]) * sqrt(Ra / Rb) * cos(kk * lb);
					/*
					cout << s << " ";
					Bits_Print(state_s,LatticeSize);
					cout << s1 << " ";
					Bits_Print(state_s1b,LatticeSize);
					cout << "i = " << i << ", j = " << j << endl;
					cout << std::get<2>(hops[ib]) * sqrt(Ra / Rb) * cos(2 * PI * Params.k * lb / LatticeSize) << " " << std::get<2>(hops[ib]) << " " << sqrt(Ra / Rb) << " " << cos(2 * PI * Params.k * lb / LatticeSize) << endl;
					*/
				}
				//cout << "x" << endl;
			}
		}
	}

	/*
	cout << "Dense matrix of H:" << endl;
	for (my_int i = 0; i < Dim; i++)
	{
		for (my_int j = 0; j < Dim; j++)
		{
			cout << setw(8) << setprecision(4) << DMat[i * Dim + j];
		}
		cout << endl;
	}
	*/
}

void Ham_HardCoreBoson<my_cmplx>::DenseMat_Build_k()
{
	my_cmplx II(0, 1);

	my_int matlen = (my_int)Dim * (my_int)Dim;
	cout << "Dim*Dim = " << matlen << endl;
	cout << "Estimated memory cost of the dense matrix: " << matlen * 8 / 1e9 << " GB" << endl;

	double kk = basis_k->get_k() * 2 * PI / LatticeSize;

	DMat = new my_cmplx[matlen];
	for (my_int s = 0; s < matlen; s++) DMat[s] = 0;
	//omp_set_num_threads(NUMTHREADS);
 //   #pragma omp parallel for schedule(dynamic)
	/*-------------------------------------------------------------------
	  H = -t*sum_{i < j}[b_i^+ b_j^- + b_i^- b_j^+]
	  + V * sum_{i < j} (i,r) r^(-3) n_i n_{i+r}
	  --------------------------------------------------------------------*/
	//double kk = Params.k
	for (my_int s = 0; s < Dim; s++) {
		my_int state_s;
		state_s = basis_k->get_state(s);
		double Ra = basis_k->get_Rs(s);

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
		// offdiagonal term
		// -t*sum_{i < j}[b_i^+ b_j^- + b_i^- b_j^+]
		for (my_int ib = 0; ib < N_hops; ib++)
		{
			my_int i = std::get<0>(hops[ib]);
			my_int j = std::get<1>(hops[ib]);
			my_int s_i = ((state_s >> i) & 1);
			my_int s_j = ((state_s >> j) & 1);
			if ((s_i) ^ (s_j))
				//if ((s_i) && (!s_j))
			{
				my_int state_s1b = state_s ^ (1 << i) ^ (1 << j);
				my_int state_s1;
				my_int lb;
				basis_k->Representative(state_s1b, state_s1, lb);
				double Rb = basis_k->CheckState_Ra(state_s1);
				if (-1 != Rb) {
					my_int s1 = basis_k->get_index(state_s1);
					DMat[s * Dim + s1] += std::get<2>(hops[ib]) * sqrt(Ra / Rb) * cos(kk * lb);
					DMat[s * Dim + s1] += std::get<2>(hops[ib]) * sqrt(Ra / Rb) * sin(kk * lb) * II;
					/*
					cout << s << " ";
					Bits_Print(state_s, LatticeSize);
					cout << s1 << " ";
					Bits_Print(state_s1b, LatticeSize);
					cout << "i = " << i << ", j = " << j << endl;
					cout << std::get<2>(hops[ib]) * sqrt(Ra / Rb) * cos(2 * PI * Params.k * lb / LatticeSize) << " " << std::get<2>(hops[ib]) << " " << sqrt(Ra / Rb) << " " << cos(2 * PI * Params.k * lb / LatticeSize) << endl;
					*/
				}
				//cout << "x" << endl;
			}
		}
	}
}