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

//template <>
//double Ham_HardCoreBoson<my_cmplx>::Cal_Hkpz_element(my_int s, my_int s1, my_int l, my_int q, my_int g, double& hj)
//{}

template <>
double Ham_HardCoreBoson<double>::Cal_Hkpz_element(my_int s, my_int s1, my_int l, my_int q, my_int g, double& hj)
{
	my_int L = LatticeSize;
	double kk = basis_kpz->get_k() * 2 * PI / L;
	my_int p = basis_kpz->get_p();
	my_int z = basis_kpz->get_z();

	//int mp, mz, mpz;
	my_int mp = basis_kpz->get_mp(s);
	my_int mz = basis_kpz->get_mz(s);
	my_int mpz = basis_kpz->get_mpz(s);
	my_int sigma = basis_kpz->get_sigma(s);
	double Na = basis_kpz->get_Norm(sigma, mp, mz, mpz) * 2 * L * L / basis_kpz->get_Rs(s);

	my_int m1p = basis_kpz->get_mp(s1);
	my_int m1z = basis_kpz->get_mz(s1);
	my_int m1pz = basis_kpz->get_mpz(s1);
	my_int sigma1 = basis_kpz->get_sigma(s1);
	double Nb = basis_kpz->get_Norm(sigma1, m1p, m1z, m1pz) * 2 * L * L / basis_kpz->get_Rs(s1);;

	//if ((4 == s && 5 == s1) || (5 == s && 4 == s1))
	if (0)
	{
		cout << "z = " << z << endl;
		Bits_Print(basis_kpz->get_state(s), L);
		cout << " For s: sigma = " << sigma << ", mp = " << mp << ", mz = " << mz << ", mpz = " << mpz << endl;
		Bits_Print(basis_kpz->get_state(s1), L);
		cout << " For s1: sigma = " << sigma1 << ", mp = " << m1p << ", mz = " << m1z << ", mpz = " << m1pz << ", l = " << l << ", q = " << q << ", g = " << g << endl;
	}
	double coeff;
	coeff = (0 == q || sigma == p) ? 1 : -1;
	coeff *= (1 == z || 0 == g) ? 1 : -1;
	coeff *= hj * sqrt(Nb / Na);

	if (sigma1 == sigma)
	{
		// case 1: mp = mz = mpz = -1
		// case 3: mp = -1, mz != -1, mpz = -1
		if ((-1 == m1p && -1 == m1z && -1 == m1pz) || (-1 == m1p) && (-1 != m1z) && (-1 == m1pz))
		{
			return cos(kk * l) * coeff;

		}
		// case 2: mp != -1, mz = -1, mpz = -1
		// case 5: mp != -1, mz != -1, mpz != -1
		if ((-1 != m1p && -1 == m1z && -1 == m1pz) || (-1 != m1p && -1 != m1z && -1 != m1pz))
		{
			return (cos(kk * l) + sigma * p * cos(kk * (l - m1p))) / (1 + sigma * p * cos(kk * m1p)) * coeff;
		}
		// case 4: mp = -1, mz == -1, mpz != -1
		if ((-1 == m1p) && (-1 == m1z) && (-1 != m1pz))
		{
			return (cos(kk * l) + sigma * p * z * cos(kk * (l - m1pz))) / (1 + sigma * p * z * cos(kk * m1pz)) * coeff;
		}
	}

	if (sigma1 != sigma)
	{
		// case 1: mp = mz = mpz = -1
		// case 3: mp = -1, mz != -1, mpz = -1
		if ((-1 == m1p && -1 == m1z && -1 == m1pz) || (-1 == m1p) && (-1 != m1z) && (-1 == m1pz))
		{
			return sigma1 * sin(kk * l) * coeff;
		}
		// case 2: mp != -1, mz = -1, mpz = -1
		// case 5: mp != -1, mz != -1, mpz != -1
		if ((-1 != m1p && -1 == m1z && -1 == m1pz) || (-1 != m1p && -1 != m1z && -1 != m1pz))
		{
			return (sigma1 * sin(kk * l) + p * sin(kk * (l - m1p))) / (1 - sigma * p * cos(kk * m1p)) * coeff;
		}
		// case 4: mp = -1, mz == -1, mpz != -1, Normalization constants: 1 + sigma*p*z*cos(k*mpz)
		if ((-1 == m1p) && (-1 == m1z) && (-1 != m1pz))
		{
			return (sigma1 * sin(kk * l) + p * z * sin(kk * (l - m1pz))) / (1 - sigma * p * z * cos(kk * m1pz)) * coeff;
		}
	}
}

template <>
void Ham_HardCoreBoson<double>::DenseMat_Build_kpz()
{
	my_int matlen = (my_int)Dim * (my_int)Dim;
	cout << "Dim*Dim = " << matlen << endl;
	cout << "Estimated memory cost of the dense matrix: " << matlen * 8 / 1e9 << " GB" << endl;

	DMat = new double[matlen];
	cout << "Allocate space for DMat..." << endl;
	for (my_int s = 0; s < matlen; s++) DMat[s] = 0;
	H_diag = new double[Dim];
	//omp_set_num_threads(NUMTHREADS);
 //   #pragma omp parallel for schedule(dynamic)
	/*-------------------------------------------------------------------
	  H = -t*sum_{i < j}[b_i^+ b_j^- + b_i^- b_j^+]
	  + V * sum_{i < j} (i,r) r^(-3) n_i n_{i+r}
	  --------------------------------------------------------------------*/
	for (my_int s = 0; s < Dim; s++) {
		my_int state_s;
		state_s = basis_kpz->get_state(s);

		// check coeff_n; state_s can be the same with different sigmas, consider only once here
		my_int na;
		if (s > 0 && basis_kpz->get_state(s - 1) == state_s) { continue; }
		else if (s < Dim - 1 && basis_kpz->get_state(s + 1) == state_s) { na = 2; }
		else { na = 1; }
		//
		double Ra = basis_kpz->get_Rs(s);

		// diagoanl term 
		// V * sum_{i < j} (i,r) r^(-3) n_i n_{i+r}
		my_int diag_index = s + s * Dim;
		for (my_int ib = 0; ib < N_ints; ib++)
		{
			my_int i = std::get<0>(ints[ib]);
			my_int j = std::get<1>(ints[ib]);
			my_int s_i = ((state_s >> i) & 1);
			my_int s_j = ((state_s >> j) & 1);
 
			if (s_i & s_j)
			{
				DMat[diag_index] += std::get<2>(ints[ib]);
				H_diag[s] = DMat[diag_index];
				if (2 == na)
				{
					DMat[(s + 1) + (s + 1) * Dim] = DMat[diag_index];
					H_diag[s + 1] = DMat[diag_index];
				}

			}
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

					///*
					//if ((4 == s && 5 == s1) || (5 == s && 4 == s1))
					if (0)
					{
						cout << "-------------------------------------" << endl;
						cout << i << " " << j << endl;
						cout << "s = " << s << ", State[s]: ";
						Bits_Print(state_s, LatticeSize);
						cout << endl;
						cout << "s1 = " << s1 << ", State[s1b]: ";
						Bits_Print(state_s1b, LatticeSize);
						cout << endl;
					}
					//*/
					for (my_int sa = s; sa < s + na; sa++)
						for (my_int sb = s1; sb < s1 + nb; sb++)
						{
							DMat[sa * Dim + sb] += Cal_Hkpz_element(sa, sb, lb, qb, gb, std::get<2>(hops[ib]));
							///*
							//if ((4 == s && 5 == s1) || (5 == s && 4 == s1))
							if (0)
							{
								cout << "sa = " << sa
									<< ", sb = " << sb
									<< ", DMat[sa * Dim + sb] += " << Cal_Hkpz_element(sa, sb, lb, qb, gb, std::get<2>(hops[ib]))
									<< endl;
							}
							//*/
						}
				}
			}
		}
	}

	//DenseMat_Display();
}