#pragma once

template<class dtype>
double Ham_HardCoreBoson<dtype>::Cal_SE_wfsq(double* wfsq)
{
	return VN_entropy(Dim, wfsq);
}

template<class dtype>
double  Ham_HardCoreBoson<dtype>::Cal_IPR_wfsq(double* wfsq)
{
	double IPR = 1.0 / cblas_dnrm2(Dim, wfsq, 1);	// dnrm2 = sqrt(sum(wfsq^2))
	IPR *= IPR;
	return IPR;
}

template<class dtype>
void Ham_HardCoreBoson<dtype>::Fetch_EigenVec(const my_int& s, dtype* wf)
{
	for (my_int r = 0; r < Dim; r++)
	{
		wf[r] = DMat[s * Dim + r];
	}
}

/*
template<class dtype>
void Ham_HardCoreBoson<dtype>::Cal_ninj_wfsq(double* wfsq, double* ninj)
{
	for (my_int i = 0; i < LatticeSize * LatticeSize; i++) ninj[i] = 0;
	for (my_int s = 0; s < Dim; s++) {
		my_int state_s = basis_kpz->get_state(s);
		int na;
		if (s > 0 && basis_kpz->get_state(s - 1) == state_s) { continue; }
		else if (s < Dim - 1 && basis_kpz->get_state(s + 1) == state_s) { na = 2; }
		else { na = 1; }
		for (my_int i = 0; i < LatticeSize; i++)
		{
			for (my_int sa = s; sa < s + na; sa++)
			{
				ninj[i * LatticeSize + i] += ((state_s >> i) & 1) ? wfsq[sa] : 0;
			}
			for (my_int j = i + 1; j < LatticeSize; j++)
			{
				for (my_int sa = s; sa < s + na; sa++)
				{
					ninj[i * LatticeSize + j] += ((state_s >> i) & 1) & ((state_s >> j) & 1) ? wfsq[sa] : 0;
				}
			}
		}
	}
	for (my_int i = 0; i < LatticeSize; i++)
		for (my_int j = i + 1; j < LatticeSize; j++)
		{
			ninj[j * LatticeSize + i] = ninj[i * LatticeSize + j];
		}
}
*/

/*
template<class dtype>
double Ham_HardCoreBoson<dtype>::Cal_EE_wf(double* wf, const my_int& size_A, my_int* sites_A)
{
	my_int size_B = LatticeSize - size_A;
	my_int dim_A = 1 << size_A;
	my_int dim_B = 1 << size_B;
	my_int dim_2d = dim_A * dim_B;
	double* wf_2d = new double[dim_2d];
	for (my_int k = 0; k < Dim; k++)
	{
		my_int s = basis->get_state(k);
		my_int ind_a, ind_b;
		bits_decomposition(s, LatticeSize, size_A, sites_A, ind_a, ind_b);
		//
		wf_2d[ind_a + ind_b * dim_A] = wf[k];
	}
	double* lambda = new double[max(dim_A, dim_B)];
	MatrixSvd(LAPACK_ROW_MAJOR, 'N', 'N', dim_B, dim_A, wf_2d, dim_A, lambda, NULL, dim_A, NULL, dim_B);
	// singular value is the sqrt of eigenvalues of reduced density martrix
	vdSqr(dim_A, lambda, lambda);
	double ee = VN_entropy(dim_A, lambda);
	delete[]wf_2d;
	delete[]lambda;
	return ee;
}*/

template<class dtype>
double Ham_HardCoreBoson<dtype>::Cal_EE_wf(my_cmplx* wf, const my_int& size_A, my_int* sites_A)
{
	my_int size_B = LatticeSize - size_A;
	my_int dim_A = 1 << size_A;
	my_int dim_B = 1 << size_B;
	my_int dim_2d = dim_A * dim_B;
	my_cmplx* wf_2d = new my_cmplx[dim_2d];
	for (my_int k = 0; k < Dim; k++)
	{
		my_int s = basis->get_state(k);
		my_int ind_a, ind_b;
		Bits_Decomposition(s, LatticeSize, size_A, sites_A, ind_a, ind_b);
		//
		wf_2d[ind_a + ind_b * dim_A] = wf[k];
	}
	double* lambda = new double[max(dim_A, dim_B)];
	MatrixSvd(LAPACK_ROW_MAJOR, 'N', 'N', dim_B, dim_A, wf_2d, dim_A, lambda, NULL, dim_A, NULL, dim_B);

	// singular value is the sqrt of eigenvalues of reduced density martrix
	vdSqr(dim_A, lambda, lambda);
	double ee = VN_entropy(dim_A, lambda);
	delete[]wf_2d;
	delete[]lambda;
	return ee;
}

template<class dtype>
double Ham_HardCoreBoson<dtype>::Cal_halfchain_EE_wf(my_cmplx* wf)
{
	my_int size_A = LatticeSize / 2;
	my_int size_B = LatticeSize - size_A;
	my_int dim_A = 1 << size_A;
	my_int dim_B = 1 << size_B;
	my_int dim_2d = dim_A * dim_B;
	my_cmplx* wf_2d = new my_cmplx[dim_2d];
	for (my_int s = 0; s < dim_2d; s++) { wf_2d[s] = 0; }
	for (my_int s = 0; s < Dim; s++)
	{
		wf_2d[basis->get_state(s)] = wf[s];
	}
	double* lambda = new double[max(dim_A, dim_B)];
	MatrixSvd(LAPACK_ROW_MAJOR, 'N', 'N', dim_B, dim_A, wf_2d, dim_A, lambda, NULL, dim_A, NULL, dim_B);

	// singular value is the sqrt of eigenvalues of reduced density martrix
	vdSqr(dim_A, lambda, lambda);
	double ee = VN_entropy(dim_A, lambda);
	delete[]wf_2d;
	delete[]lambda;
	return ee;
}

template<class dtype>
double Ham_HardCoreBoson<dtype>::Cal_halfchain_EE_wf(double* wf)
{
	my_int size_A = LatticeSize / 2;
	my_int size_B = LatticeSize - size_A;
	my_int dim_A = 1 << size_A;
	my_int dim_B = 1 << size_B;
	my_int dim_2d = dim_A * dim_B;
	double* wf_2d = new double[dim_2d];
	for (my_int s = 0; s < dim_2d; s++) { wf_2d[s] = 0; }
	for (my_int s = 0; s < Dim; s++)
	{
		wf_2d[basis->get_state(s)] = wf[s];
	}
	double* lambda = new double[max(dim_A, dim_B)];
	MatrixSvd(LAPACK_ROW_MAJOR, 'N', 'N', dim_B, dim_A, wf_2d, dim_A, lambda, NULL, dim_A, NULL, dim_B);

	// singular value is the sqrt of eigenvalues of reduced density martrix
	vdSqr(dim_A, lambda, lambda);
	double ee = VN_entropy(dim_A, lambda);
	delete[]wf_2d;
	delete[]lambda;
	return ee;
}