#include "hmatrix.h"
using namespace std;
#include <omp.h>

// initialize default values
void Parameters::SetDefault()
{	
	//
	t = 1;
	V = 1;
	alpha = 3;
	bc = 0;

	// evolution
	evo_time_max = 100;
	evo_time_steps = 100;
	evo_time_IncStyle = 0;	// 0/1 for linear/exponetial time increase

	evo_type = 0;
	evo_Lanczos_M = 10;

	evo_ni = 'y';
	evo_ninj = 'n';
	evo_IPR = 'n';
	evo_PI_r = 'n';
	evo_SE = 'n';
	evo_EE = 'n';
	evo_Pmax = 'n';
}

void Parameters::Read_input()
{
	GetParaFromInput_real("input.in", "t", t);
	GetParaFromInput_real("input.in", "V", V);
	GetParaFromInput_real("input.in", "alpha", alpha);
	GetParaFromInput_int("input.in", "bc", bc);
	
	cout << endl << "t-V Hamiltonian parameters:" << endl;
	cout << "t = " << t << endl;
	cout << "V = " << V << endl;
	cout << "bc = " << bc << ", 0/1 for obc/pbc" << endl;

	/*
	GetParaFromInput_real("input.in", "evo_time_max", evo_time_max);
	GetParaFromInput_int("input.in", "evo_time_steps", evo_time_steps);
	GetParaFromInput_int("input.in", "evo_time_IncStyles", evo_time_IncStyle);

	GetParaFromInput_int("input.in", "evo_type", evo_type);
	GetParaFromInput_int("input.in", "evo_Lanczos_M", evo_Lanczos_M);

	GetParaFromInput_char("input.in", "evo_ni", evo_ni);
	GetParaFromInput_char("input.in", "evo_ninj", evo_ninj);
	GetParaFromInput_char("input.in", "evo_IPR", evo_IPR);
	GetParaFromInput_char("input.in", "evo_PI_r", evo_PI_r);
	GetParaFromInput_char("input.in", "evo_SE", evo_SE);
	GetParaFromInput_char("input.in", "evo_EE", evo_EE);
	GetParaFromInput_char("input.in", "evo_Pmax", evo_Pmax);
	//
	GetTimeVec();
	*/

	/*
	// entanglement entropy <- information of subsystem:
	if ('y' == evo_EE)
	{
		ifstream ifin("input_EE_SubSysInfo.in", ios::in);
		if (ifin)
		{
			ifin >> EE_size_A;
			EE_sites_A = new my_int[EE_size_A];
			for (my_int i = 0; i < EE_size_A; i++)
			{
				ifin >> EE_sites_A[i];
			}
			ifin.close();
		}
		else
		{
			cout << "Subsystem size = 0! No entanglement entropy computed!" << endl;
			EE_size_A = 0;
		}
	}
	*/
}

void Parameters::GetTimeVec()
{
	GetParaFromInput_real("input.in", "evo_time_max", evo_time_max);
	GetParaFromInput_int("input.in", "evo_time_steps", evo_time_steps);
	GetParaFromInput_int("input.in", "evo_time_IncStyles", evo_time_IncStyle);
	GetParaFromInput_int("input.in", "evo_nini", evo_nini);

	my_int nt = evo_time_steps;
	time_vec = new double[nt];
	// linear time increase
	if (0 == evo_time_IncStyle)
	{
		for (my_int i = 0; i < nt; i++)
		{
			time_vec[i] = i * evo_time_max / nt;
		}
	}

	else if (1 == evo_time_IncStyle)
	{
		double tmaxlog = log10(evo_time_max);
		double tminlog = log10(1e-3);
		time_vec[0] = 0;
		for (my_int i = 1; i < nt; i++)
		{
			time_vec[i] = pow(10, tminlog + i * (tmaxlog - tminlog) / nt);
		}
	}
}

void Dyn_DataStruc::Initialize(Parameters& para, const my_int& _dim, const my_int& _prt_ind)
{
	dim = _dim;
	prt_ind = _prt_ind;

	// necessaraties for evo
	my_int t_len = para.evo_time_steps;
	//GetTimeVec(para);

	if ('y' == para.evo_ni)
		ni_t = new double[t_len * para.LatticeSize];
	if ('y' == para.evo_ninj)
		ninj_t = new double[t_len * para.LatticeSize * para.LatticeSize];
	if ('y' == para.evo_IPR)
		IPR_t = new double[t_len];
	if ('y' == para.evo_SE)
		SE_t = new double[t_len];
	if ('y' == para.evo_EE)
		EE_t = new double[t_len];
	if ('y' == para.evo_PI_r)
	{
		PI_r_t = new double[t_len * (para.LatticeSize + 1)];
		HamD_k = new double[_dim];
	}
}

void Dyn_DataStruc::ReleaseSpace(Parameters& para)
{
	// necessaraties for evo
	if ('y' == para.evo_IPR)
		delete[]IPR_t;
	if ('y' == para.evo_SE)
		delete[]SE_t;
	if ('y' == para.evo_ni)
		delete[]ni_t;
	if ('y' == para.evo_ninj)
		delete[]ninj_t;
	if ('y' == para.evo_EE)
		delete[]EE_t;
	if ('y' == para.evo_PI_r)
	{
		delete[]PI_r_t;
		delete[]HamD_k;
	}
}

void::Dyn_DataStruc::PrintDynResults(Parameters& para) {

	my_int t_len = para.evo_time_steps;
	// real-time evolution results 
	if ('y' == para.evo_ni)
	{
		char fname[80];
		sprintf(fname, "evo_ni_ind%d.bin", prt_ind);
		Vec_fwrite_double(fname, ni_t, t_len * para.LatticeSize);
	}

	if ('y' == para.evo_ninj)
	{
		char fname[80];
		sprintf(fname, "evo_ninj_ind%d.bin", prt_ind);
		Vec_fwrite_double(fname, ninj_t, t_len * para.LatticeSize * para.LatticeSize);
	}

	if ('y' == para.evo_PI_r)
	{
		char fname[80];
		sprintf(fname, "evo_PIr_ind%d.bin", prt_ind);
		Vec_fwrite_double(fname, PI_r_t, t_len * (para.LatticeSize + 1));
	}

	if ('y' == para.evo_IPR)
	{
		char fname[80];
		sprintf(fname, "evo_IPR_ind%d.bin", prt_ind);
		Vec_fwrite_double(fname, IPR_t, t_len);
	}

	if ('y' == para.evo_SE)
	{
		char fname[80];
		sprintf(fname, "evo_SE_ind%d.bin", prt_ind);
		Vec_fwrite_double(fname, SE_t, t_len);
	}

	if ('y' == para.evo_EE)
	{
		char fname[80];
		sprintf(fname, "evo_EE_ind%d.bin", prt_ind);
		Vec_fwrite_double(fname, EE_t, t_len);
	}
}
