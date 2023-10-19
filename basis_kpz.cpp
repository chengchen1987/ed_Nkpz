#include <cstdlib>
#include <cstdio>
#include <iostream> 
#include <fstream>
#include <iomanip>
//#include <cmath>
//#include <cstring>
//#include <vector>
using namespace std;
#include "basis.h" 

Basis_kpz::Basis_kpz(my_int _L, my_int _nup, my_int _k, my_int _p, my_int _z) :
	L(_L),
	nup(_nup),
	k(_k),
	p(_p),
	z(_z)
{
	double kk = 2 * PI * k / L;
	Index = new my_int[1 << L];
	for (my_int i = 0; i < (1 << L); i++) Index[i] = -1;

	my_int* aux_state = new my_int[nchoosek(L, nup)];
	my_int* aux_state_r = new my_int[nchoosek(L, nup)];
	my_int* aux_state_sigma = new my_int[nchoosek(L, nup)];
	my_int* aux_state_mp = new my_int[nchoosek(L, nup)];
	my_int* aux_state_mz = new my_int[nchoosek(L, nup)];
	my_int* aux_state_mpz = new my_int[nchoosek(L, nup)];

	my_int j = 0;
	my_int totalcharge;
	my_int aux_R, aux_mp, aux_mz, aux_mpz;
	for (my_int s = 0; s < (1 << L); s++) {
		totalcharge = Bits_Count1(s);
		if (totalcharge == nup) {
			CheckState_R_mpz(s, aux_R, aux_mp, aux_mz, aux_mpz);
			// choose a or both sigma
			my_int N_sigma = (0 == k || (L == k * 2)) ? 1 : 2;
			my_int sigma = 1;
			for (my_int i_sigma = 0; i_sigma < N_sigma; i_sigma++)
			{
				my_int aux_r = aux_R;
				if (0 == get_Norm(sigma, aux_mp, aux_mz, aux_mpz))
				{
					aux_r = -1;
				}
				if (-1 == sigma && 0 != get_Norm(1, aux_mp, aux_mz, aux_mpz) 
					&& (-1 != aux_mp || -1 != aux_mpz))
				{
					aux_r = -1;
				}

				if (-1 != aux_r)
				{
					Index[s] = j;
					aux_state[j] = s;
					aux_state_r[j] = aux_r;
					aux_state_mp[j] = aux_mp;
					aux_state_mz[j] = aux_mz;
					aux_state_mpz[j] = aux_mpz;
					aux_state_sigma[j] = sigma;

					j++;
					/*
					cout << "j: " << j - 1 << ", R = " << aux_R << ", mp = " << aux_mp << ", mz = " << aux_mz << ", mpz = " << aux_mpz 
						<< ", sigma = " << sigma << ", Norm = " << get_Norm(sigma, aux_mp, aux_mz, aux_mpz) << endl;
					Bits_Print(s, L);
					cout  << endl;
					*/
				}
				sigma -= 2;
			}
		}
	}
	// 
	Dim = j;
	State = new my_int[Dim];
	State_R = new my_int[Dim];
	State_mp = new my_int[Dim];
	State_mz = new my_int[Dim];
	State_mpz = new my_int[Dim];
	State_sigma = new my_int[Dim];
	for (my_int s = 0; s < Dim; s++)
	{
		State[s] = aux_state[s];
		State_R[s] = aux_state_r[s];
		State_mp[s] = aux_state_mp[s];
		State_mz[s] = aux_state_mz[s];
		State_mpz[s] = aux_state_mpz[s];
		State_sigma[s] = aux_state_sigma[s];
	}
	delete[]aux_state;
	delete[]aux_state_r;
	delete[]aux_state_mp;
	delete[]aux_state_mz;
	delete[]aux_state_mpz;
	delete[]aux_state_sigma;

	cout << "Dim = " << Dim << endl;
	cout << "Estimated memory cost of Basis: " << (Dim * 4 * 6 + (1 << L) * 4) / 1e9 << " GB" << endl;
}

Basis_kpz::~Basis_kpz()
{
	delete[] State;
	delete[] Index;

	delete[] State_R;
	delete[] State_mp;
	delete[] State_mz;
	delete[] State_mpz;
	delete[] State_sigma;
}

// for translational invariant Hamiltonian
// check/get Ra
void Basis_kpz::CheckState_R_mpz(my_int s, my_int& r, my_int& mp, my_int& mz, my_int& mpz)
{
	r = -1; mp = -1; mz = -1; mpz = -1;
	// check minimun R for T^R|s> = |s>
	my_int t = s;
	//cout << "CheckState_R_mpz, s = "; Bits_Print(t, L); cout << endl;
	//cout << "Cherk T^R|s> = |s>" << endl;
	for (my_int i = 1; i <= L; i++)
	{
		//Bits_Print(t, L); cout << endl;
		t = Bits_CycleLeft(t, L);
		if (t < s) { return; }
		if (t == s)
		{
			if (0 == (k % (L / i))) { r = i; }
			break;
		}
	}
	// check minimum m_p for T^m_p*P|s> = |s>
	t = Bits_Reflect(t, L);
	//cout << "Check T^m_p*P|s> = |s>" << endl;
	for (my_int i = 0; i < r; i++)
	{
		//Bits_Print(t, L); cout << endl;
		if (t < s) { r = -1; return; }
		if (t == s) { mp = i; break; }
		t = Bits_CycleLeft(t, L);
	}
	// check minimum m_z for T^m_z*Z|s> = |s>
	t = Bits_Invert(s, L);
	//cout << "Check T^m_z*Z|s> = |s>" << endl;
	for (my_int i = 0; i < r; i++)
	{
		//Bits_Print(t, L); cout << endl;
		if (t < s) { r = -1; return; }
		if (t == s) { mz = i; break; }
		t = Bits_CycleLeft(t, L);
	}
	// check minimum m_pz for T^m_pz*PZ|s> = |s>
	t = Bits_Reflect(Bits_Invert(s, L), L);
	//cout << "Check T^m_pz*PZ|s> = |s>" << endl;
	for (my_int i = 0; i < r; i++)
	{
		//Bits_Print(t, L); cout << endl;
		if (t < s) { r = -1; return; }
		if (t == s) { mpz = i; break; }
		t = Bits_CycleLeft(t, L);
	}
	//if (-1 != r) { cout << " mp = " << mp << ", mz = " << mz << ", mpz = " << mpz << endl; }
}

void Basis_kpz::Representative(const my_int& s, my_int& sa, my_int& la, my_int& qa, my_int& ga)
{
	sa = s; la = 0; qa = 0; ga = 0;
	my_int auxt = s;
	for (my_int i = 1; i < L; i++)
	{
		auxt = Bits_CycleLeft(auxt, L);
		if (auxt < sa)
		{
			sa = auxt;
			la = i;
		}
	}
	//if (0 != la) return;
	// consider parity
	auxt = Bits_Reflect(s, L);
	for (my_int i = 0; i < L; i++)
	{
		if (auxt < sa)
		{
			sa = auxt;
			la = i;
			qa = 1;
		}
		auxt = Bits_CycleLeft(auxt, L);
		
	}
	//if (0 != la) return;
	// consider spin-inversion
	auxt = Bits_Invert(s, L);
	for (my_int i = 0; i < L; i++)
	{
		if (auxt < sa)
		{
			sa = auxt;
			la = i;
			qa = 0;
			ga = 1;
		}
		auxt = Bits_CycleLeft(auxt, L);
		
	}
	//if (0 != la) return;
	// consider parity and spin-inversion
	auxt = Bits_Reflect(Bits_Invert(s, L), L);
	for (my_int i = 0; i < L; i++)
	{
		if (auxt < sa)
		{
			sa = auxt;
			la = i;
			qa = 1;
			ga = 1;
		}
		auxt = Bits_CycleLeft(auxt, L);
	}
}

double Basis_kpz::get_Norm(my_int sigma, my_int mp, my_int mz, my_int mpz)
{
	double kk = 2 * PI * k / L;
	// case 1: mp = mz = mpz = -1, Normalization constants: 1
	if ((-1 == mp) && (-1 == mz) && (-1 == mpz))
	{
		return 1;
	}
	// case 2: mp != -1, mz = -1, mpz = -1, Normalization constants: 1 + sigma*p*cos(k*mp)
	else if ((-1 != mp) && (-1 == mz) && (-1 == mpz))
	{
		return 1 + sigma * p * cos(kk * mp);
	}
	// case 3: mp = -1, mz != -1, mpz = -1, Normalization constants: 1 + z*cos(k*mz)
	if ((-1 == mp) && (-1 != mz) && (-1 == mpz))
	{
		return 1 + z * cos(kk * mz);
	}
	// case 4: mp = -1, mz == -1, mpz != -1, Normalization constants: 1 + sigma*p*z*cos(k*mpz)
	if ((-1 == mp) && (-1 == mz) && (-1 != mpz))
	{
		return 1 + sigma * p * z * cos(kk * mpz);
	}
	// case 5: mp != -1, mz != -1, mpz != -1, Normalization constants: (1 + sigma*p*cos(k*mp))*(1 + z*cos(k*mz))
	if ((-1 != mp) && (-1 != mz) && (-1 != mpz))
	{
		return (1 + sigma * p * cos(kk * mp)) * (1 + z * cos(kk * mz));
	}
}

my_int Basis_kpz::get_state(const my_int& s)
{
	return State[s];
}

my_int Basis_kpz::get_Rs(const my_int& s)
{
	return State_R[s];
}

my_int Basis_kpz::get_mp(const my_int& s)
{
	return State_mp[s];
}

my_int Basis_kpz::get_mz(const my_int& s)
{
	return State_mz[s];
}

my_int Basis_kpz::get_mpz(const my_int& s)
{
	return State_mpz[s];
}

my_int Basis_kpz::get_sigma(const my_int& s)
{
	return State_sigma[s];
}

my_int Basis_kpz::get_index(const my_int& s)
{
	return Index[s];
}