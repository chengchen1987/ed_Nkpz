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

template <>
double Ham_HardCoreBoson<double>::Cal_Hkp_element(my_int s, my_int s1, my_int l, my_int q, double& hj)
{
	double aux;

	my_int L = LatticeSize;
	double kk = basis_kp->get_k() * 2 * PI / L;
	my_int p = basis_kp->get_p();
	my_int m = basis_kp->get_ms(s);
	my_int m1 = basis_kp->get_ms(s1);
	my_int sigma = basis_kp->get_sigma(s);
	my_int sigma1 = basis_kp->get_sigma(s1);
	double Na = -1 == m ? L * L / basis_kp->get_Rs(s) : L * L / basis_kp->get_Rs(s) * (1 + sigma * p * cos(kk * m));
	double Nb = -1 == m1 ? L * L / basis_kp->get_Rs(s1) : L * L / basis_kp->get_Rs(s1) * (1 + sigma1 * p * cos(kk * m1));

	//cout << " m1 = " << m1 << ", l = " << l << ", q = " << q << endl;
	//if ((0 == s && 1 == s1) || (1 == s1 && 0 == s))
	if (0)
	{
		Bits_Print(basis_kp->get_state(s), L);
		cout << " For s: sigma = " << sigma << ", m = " << m << endl;
		Bits_Print(basis_kp->get_state(s1), L);
		cout << " For s1: sigma = " << sigma1 << ", m = " << m1 << ", l = " << l << ", q = " << q << endl;
	}

	double coeff;
	coeff = (0 == q || sigma == p) ? 1 : -1;
	coeff *= hj * sqrt(Nb / Na);

	if (sigma1 == sigma)
	{
		if (-1 == m1)
		{
			aux = cos(kk * l) * coeff;
		}
		else
		{
			aux = (cos(kk * l) + sigma * p * cos(kk * (l - m1))) / (1 + sigma * p * cos(kk * m1)) * coeff;
		}
	}
	else if (sigma1 == -sigma)
	{
		if (-1 == m1)
		{
			aux = sigma1 * sin(kk * l) * coeff;
		}
		else
		{
			aux = (sigma1 * sin(kk * l) + p * sin(kk * (l - m1))) / (1 - sigma * p * cos(kk * m1)) * coeff;
		}
	}
	return aux;
}

template <>
void Ham_HardCoreBoson<double>::DenseMat_Build_kp()
{
	my_int matlen = (my_int)Dim * (my_int)Dim;
	cout << "Dim*Dim = " << matlen << endl;
	cout << "Estimated memory cost of the dense matrix: " << matlen * 8 / 1e9 << " GB" << endl;

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
		state_s = basis_kp->get_state(s);

		// check coeff_n; state_s can be the same with different sigmas, consider only once here
		my_int na;
		if (s > 0 && basis_kp->get_state(s - 1) == state_s)
			//break;
			continue;
		else if (s < Dim - 1 && basis_kp->get_state(s + 1) == state_s)
			na = 2;
		else
			na = 1;
		//
		double Ra = basis_kp->get_Rs(s);

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
				if (2 == na)
				{
					DMat[(s + 1) + (s + 1) * Dim] = DMat[diag_index];
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
				my_int nb, lb, qb;
				basis_kp->Representative(state_s1b, state_s1, lb, qb);

				my_int s1 = basis_kp->get_index(state_s1);
				if (-1 != s1)
				{
					if (s1 > 1 && basis_kp->get_state(s1 - 1) == state_s1)
					{
						nb = 2; s1 -= 1;
					}
					else if (s1 < Dim - 1 && basis_kp->get_state(s1 + 1) == state_s1)
					{
						nb = 2;
					}
					else
					{
						nb = 1;
					}

					//if ((0 == s && 1 == s1) || (1 == s1 && 0 == s))
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
					for (my_int sa = s; sa < s + na; sa++)
						for (my_int sb = s1; sb < s1 + nb; sb++)
						{
							DMat[sa * Dim + sb] += Cal_Hkp_element(sa, sb, lb, qb, std::get<2>(hops[ib]));
							//if ((0 == s && 1 == s1) || (0 == s1 && 1 == s))
							if (0)
							{
								cout << "sa = " << sa
									<< ", sb = " << sb
									<< ", DMat[sa * Dim + sb] += " << Cal_Hkp_element(sa, sb, lb, qb, std::get<2>(hops[ib]))
									<< endl;
							}
						}
				}
			}
		}
	}

	//DenseMat_Display();
	//DenseMat_Check_Hermitian();
}

