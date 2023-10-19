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

Basis_k::Basis_k(my_int _L, my_int _nup, my_int _k) :
	L(_L),
	nup(_nup),
	k(_k)
{
	Index = new my_int[1 << L];
	for (my_int i = 0; i < (1 << L); i++) Index[i] = -1;
	my_int* aux_state = new my_int[nchoosek(L, nup)];
	my_int* aux_state_r = new my_int[nchoosek(L, nup)];
	my_int j = 0;
	int totalcharge;
	for (my_int s = 0; s < (1 << L); s++) {
		totalcharge = Bits_Count1(s);
		if (totalcharge == nup) {
			my_int auxR = CheckState_Ra(s);
			//Bits_Print(s, L);
			//cout << "Ra: " << auxR;
			if (-1 != auxR)
			{
				Index[s] = j;
				aux_state[j] = s;
				aux_state_r[j] = auxR;
				j++;
			}
		}
	}
	// 
	Dim = j;
	State = new my_int[Dim];
	State_R = new my_int[Dim];
	for (my_int s = 0; s < Dim; s++)
	{
		State[s] = aux_state[s];
		State_R[s] = aux_state_r[s];
	}
	delete[]aux_state;
	delete[]aux_state_r;
}

Basis_k::~Basis_k()
{
	delete[] State;
	delete[] Index;
}

// for translational invariant Hamiltonian
// check/get Ra
my_int Basis_k::CheckState_Ra(my_int s)
{
	my_int t = s;
	for (my_int i = 1; i <= L; i++)
	{
		t = Bits_CycleLeft(t, L);
		if (t < s)
		{
			return -1;
		}
		if (t == s)
		{
			if (0 == (k % (L / i))) return i;
			else return -1;
		}
	}
}

void Basis_k::Representative(const my_int& s, my_int& sa, my_int& la)
{
	sa = s;
	la = 0;

	my_int auxt = s;

	for (my_int i = 1; i < L; i++)
	{
		auxt = Bits_CycleLeft(auxt, L);
		if (auxt < sa)
		{
			sa = auxt;
			//la = min(i,L-i);
			la = i;
		}
	}
}

my_int Basis_k::get_state(const my_int& s)
{
	return State[s];
}

my_int Basis_k::get_Rs(const my_int& s)
{
	return State_R[s];
}

my_int Basis_k::get_index(const my_int& s)
{
	return Index[s];
}