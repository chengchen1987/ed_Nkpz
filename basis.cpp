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

Basis::Basis(my_int _L, my_int _nup) :
	L(_L),
	nup(_nup),
	Dim(nchoosek(L, nup))
{
	State = new my_int[Dim];
	Index = new my_int[1 << L];
	for (my_int i = 0; i < (1 << L); i++) Index[i] = -1;

	my_int j = 0;
	int totalcharge;
	for (my_int s = 0; s < (1 << L); s++) {
		totalcharge = Bits_Count1(s);
		if (totalcharge == nup) {
			State[j] = s;
			Index[s] = j;
			j++;
		}
	}
}



Basis::~Basis()
{
	delete[] State;
	delete[] Index;
}


my_int Basis::get_state(const my_int& s)
{
	return State[s];
}

my_int Basis::get_index(const my_int& s)
{
	return Index[s];
}

/*
my_int Basis::get_index(const my_int& s)
{
	if (s < 0) return -1;
	my_int bmin = 0, bmax = Dim;
	my_int b;
	while (bmin <= bmax) {
		b = bmin + (bmax - bmin) / 2;
		my_int aux = State[b];
		if (s == aux) {
			return b;
		}
		else if (s > aux) {
			bmin = b + 1;
		}
		else if (s < aux) {
			bmax = b - 1;
		}
	}
	return -1;
}
*/

/*
void Basis::PrintStates()
{
	ofstream ofd("states");
	for (my_int i = 0; i < Dim; i++) {
		ofd << setw(4) << i << setw(8) << State[i];
		ofd << "  ";
		for (my_int j = 0; j < L; j++) {
			ofd << ((State[i] >> (L - 1 - j)) & 1);
		}
		ofd << endl;
	}
	ofd.close();
}
*/

