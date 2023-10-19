#pragma once

#include "common.h"

class Basis
{
private:
	my_int L;
	my_int nup;

	my_int Dim;
	my_int* State;
	my_int* Index;

public:
	~Basis();
	// Basis with conserved nup
	Basis(my_int _L, my_int _nup);

	my_int get_L(void) { return L; }
	my_int get_nup(void) { return nup; }

	my_int get_Dim(void) { return Dim; }

	my_int get_state(const my_int& index);
	my_int get_index(const my_int& state);
};

// momentum states, conserved k
class Basis_k
{
private:
	my_int L;
	my_int nup;
	my_int k;

	my_int Dim;
	my_int* State;
	my_int* State_R;
	my_int* Index;
public:
	~Basis_k();
	// Basis with conserved nup, k
	Basis_k(my_int _L, my_int _nup, my_int _k);

	my_int CheckState_Ra(my_int s);
	void Representative(const my_int& s, my_int& sa, my_int& la);

	my_int get_L(void) { return L; }
	my_int get_nup(void) { return nup; }
	my_int get_k(void) { return k; }

	my_int get_Dim(void) { return Dim; }

	my_int get_state(const my_int& index);
	my_int get_Rs(const my_int& index);
	my_int get_index(const my_int& state);
	// 
};

// semi-momentum states + parity, conserved k,p
class Basis_kp
{
private:
	my_int L;
	my_int nup;
	my_int k, p;

	my_int Dim;
	my_int* State;
	my_int* Index;

	my_int* State_R;
	my_int* State_sigma;
	my_int* State_m;
public:
	~Basis_kp();
	// Basis with conserved nup, k
	Basis_kp(my_int _L, my_int _nup, my_int _k, my_int _p);

	void CheckState_R_m(my_int s, my_int &r, my_int &m);
	void Representative(const my_int& s, my_int& sa, my_int& la, my_int&qa);

	my_int get_L(void) { return L; }
	my_int get_nup(void) { return nup; }
	my_int get_k(void) { return k; }
	my_int get_p(void) { return p; }
	
	my_int get_Dim(void) { return Dim; }

	my_int get_state(const my_int& index);
	my_int get_Rs(const my_int& index);
	my_int get_ms(const my_int& index);
	my_int get_sigma(const my_int& index);

	my_int get_index(const my_int& state);
	// 
};

// semi-momentum states + parity + spin-inversion, conserved k,p,z for only the mz=0 sector
class Basis_kpz
{
private:
	my_int L;
	my_int nup;
	my_int k, p, z;

	my_int Dim;
	my_int* State;
	my_int* Index;

	my_int* State_R;
	my_int* State_sigma;
	my_int* State_mp;
	my_int* State_mz;
	my_int* State_mpz;
public:
	~Basis_kpz();
	// Basis with conserved nup, k
	Basis_kpz(my_int _L, my_int _nup, my_int _k, my_int _p, my_int _z);

	void CheckState_R_mpz(my_int s, my_int& r, my_int& mp, my_int &mz, my_int &mpz);
	void Representative(const my_int& s, my_int& sa, my_int& la, my_int& qa, my_int &ga);
	double get_Norm(my_int sigma, my_int mp, my_int mz, my_int mpz);

	my_int get_L(void) { return L; }
	my_int get_nup(void) { return nup; }
	my_int get_k(void) { return k; }
	my_int get_p(void) { return p; }
	my_int get_z(void) { return z; }

	my_int get_Dim(void) { return Dim; }

	my_int get_state(const my_int& index);
	my_int get_Rs(const my_int& index);
	my_int get_mp(const my_int& index);
	my_int get_mz(const my_int& index);
	my_int get_mpz(const my_int& index);
	my_int get_sigma(const my_int& index);

	my_int get_index(const my_int& state);
	// 
};