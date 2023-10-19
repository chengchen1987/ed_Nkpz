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

Basis_kp::Basis_kp(my_int _L, my_int _nup, my_int _k, my_int _p) :
	L(_L),
	nup(_nup),
	k(_k),
	p(_p)
{
	Index = new my_int[1 << L];
	for (my_int i = 0; i < (1 << L); i++) Index[i] = -1;

	my_int* aux_state = new my_int[nchoosek(L, nup)];
	my_int* aux_state_r = new my_int[nchoosek(L, nup)];
	my_int* aux_state_sigma = new my_int[nchoosek(L, nup)];
	my_int* aux_state_m = new my_int[nchoosek(L, nup)];

	my_int j = 0;
	my_int totalcharge;
	my_int aux_R, aux_m;
	for (my_int s = 0; s < (1 << L); s++) {
		totalcharge = Bits_Count1(s);
		if (totalcharge == nup) {
			CheckState_R_m(s, aux_R, aux_m);
			// choose a or both sigma
			my_int N_sigma = 2;
			if (0 == k || (L == k * 2))
			{
				N_sigma = 1;
			}
			my_int sigma = 1;
			for (my_int i_sigma = 0; i_sigma < N_sigma; i_sigma++)
			{
				my_int aux_r = aux_R;
				if (-1 != aux_m)
				{
					if (0 == 1 + sigma * p * cos(k * aux_m * 2 * PI / L))
						aux_r = -1;
					if (-1 == sigma && 0 != (1 - sigma * p * cos(k * aux_m * 2 * PI / L)))
						aux_r = -1;
				}
				if (-1 != aux_r)
				{
					Index[s] = j;
					aux_state[j] = s;
					aux_state_r[j] = aux_r;
					aux_state_m[j] = aux_m;
					aux_state_sigma[j] = sigma;

					j++;

					/*
					cout << "j: " << j - 1 << ", R = " << aux_R << ", m = " << aux_m << ", sigma = " << sigma << endl;
					Bits_Print(s, L);
					cout << endl;
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
	State_m = new my_int[Dim];
	State_sigma = new my_int[Dim];
	for (my_int s = 0; s < Dim; s++)
	{
		State[s] = aux_state[s];
		State_R[s] = aux_state_r[s];
		State_m[s] = aux_state_m[s];
		State_sigma[s] = aux_state_sigma[s];
	}
	delete[]aux_state;
	delete[]aux_state_r;
	delete[]aux_state_m;
	delete[]aux_state_sigma;
}

Basis_kp::~Basis_kp()
{
	delete[] State;
	delete[] Index;

	delete[] State_R;
	delete[] State_m;
	delete[] State_sigma;
}

// for translational invariant Hamiltonian
// check/get Ra
void Basis_kp::CheckState_R_m(my_int s, my_int& r, my_int& m)
{
	r = -1; m = -1;
	my_int t = s;
	//cout << endl; cout << "i = " << 0 << " " << s << " "; Bits_Print(t, L); cout << endl;
	for (my_int i = 1; i <= L; i++)
	{
		t = Bits_CycleLeft(t, L);
		//cout << "i = " << i << " " << t << " "; Bits_Print(t, L + 1); cout << endl;
		if (t < s)
		{
			return;
		}
		if (t == s)
		{
			if (0 == (k % (L / i)))
			{
				r = i;
			}
			break;
		}
	}
	// 
	t = Bits_Reflect(s, L);

	/*
	cout << "Bits_Reflect()! " << endl;
	Bits_Print(s, L);
	cout << endl;
	Bits_Print(t, L);
	cout << endl;
	*/

	for (my_int i = 0; i < r; i++)
	{
		//cout << "i = " << i << " " << t << " "; Bits_Print(t, L); cout << endl;
		if (t < s)
		{
			r = -1; return;
		}
		if (t == s)
		{
			m = i; return;
		}
		t = Bits_CycleLeft(t, L);
	}
}

void Basis_kp::Representative(const my_int& s, my_int& sa, my_int& la, my_int& qa)
{
	sa = s; la = 0; qa = 0;
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
}

my_int Basis_kp::get_state(const my_int& s)
{
	return State[s];
}

my_int Basis_kp::get_Rs(const my_int& s)
{
	return State_R[s];
}

my_int Basis_kp::get_ms(const my_int& s)
{
	return State_m[s];
}

my_int Basis_kp::get_sigma(const my_int& s)
{
	return State_sigma[s];
}

my_int Basis_kp::get_index(const my_int& s)
{
	return Index[s];
}