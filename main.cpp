#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <ctime>
#include <cmath>
using namespace std;
#include "common.h"
#include "basis.h"
#include "hmatrix.h"

int main() {

	cout << "sizeof(my_int)" << sizeof(my_int) << endl;
	cout << "sizeof(lapack_int)" << sizeof(lapack_int) << endl;

	cout << "Exact diagonalization for HardCore Boson t-V model with long-range ionteractions(V)..." << endl;
	{
		// read lattice informationi from input files
		my_int LatticeSize, nup;
		my_int bc, k, p;
		GetParaFromInput_int("input.in", "LatticeSize", LatticeSize);
		GetParaFromInput_int("input.in", "nup", nup);
		GetParaFromInput_int("input.in", "bc", bc);
		//GetParaFromInput_real("input.in", "alpha", alpha)
		
		// a test on diagonal elements 
		if (0)
		{
			cout << "===================================================================" << endl;
			cout << "Dense matrix diagonalization with conserved nup (number of particles)" << endl;

			cout << "LatticeSize: " << LatticeSize << endl;
			cout << "No. of up-spins: " << nup << endl;
			// generate basis
			my_int time0 = time(0);
			Basis basis(LatticeSize, nup);
			my_int time1 = time(0);
			cout << "Time(s) to Build Basis: " << time1 - time0 << endl << endl;

			Ham_HardCoreBoson<double> ham(&basis);
			ham.DenseMat_Diag_Compare_Disorder_LRI();
			my_int time2 = time(0);
			cout << "Time(s) to Get Diag-elements: " << time2 - time1 << endl << endl;
		}
		 
		// dense matrix algorithm, conserved nup (number of particles)
		if (0)
		{
			cout << "===================================================================" << endl;
			cout << "Dense matrix diagonalization with conserved nup (number of particles)" << endl;

			cout << "LatticeSize: " << LatticeSize << endl;
			cout << "No. of up-spins: " << nup << endl;
			// generate basis
			my_int time0 = time(0);
			Basis basis(LatticeSize, nup);
			my_int time1 = time(0);
			cout << "Time(s) to Build Basis: " << time1 - time0 << endl << endl;

			Ham_HardCoreBoson<double> ham(&basis);

			// dense 
			ham.DenseMat_Build();
			my_int time2 = time(0);
			cout << "Time(s) to Build Dense Matrix: " << time2 - time1 << endl << endl;

			//double* HMat = new double[ham.Dim * ham.Dim];
			//for (my_int i = 0; i < ham.Dim * ham.Dim; i++) HMat[i] = ham.DMat[i];
			ham.DenseMat_Eig();
			// eigenvector test, stored in rows 
			/*double* wf = new double[ham.Dim];
			cblas_dgemv(CblasRowMajor, CblasNoTrans, ham.Dim, ham.Dim, 1, HMat, ham.Dim, &ham.DMat[0], 1, 0, wf, 1);
			double aux = cblas_ddot(ham.Dim, &ham.DMat[0], 1, wf, 1);
			cout << "<wf0|H|wf0> = " << aux << endl;*/
			my_int time3 = time(0);
			cout << "Time(s) to Diagonalize Dense Matrix: " << time3 - time2 << endl << endl;

			ham.DenseMat_Static();
			my_int time4 = time(0);
			cout << "Time(s) for static observables: " << time4 - time3 << endl << endl;

			ham.DenseMat_Dynamic();
			my_int time5 = time(0);
			cout << "Time(s) for dynamic observables: " << time5 - time4 << endl << endl;
			
			//Vec_fwrite_double("spec.bin", ham.spec, ham.Dim);
		}

		// dense matrix algorithm, conserved nup (number of particles) and k (momentum)
		if (1)
		{
			//GetParaFromInput_int("input.in", "k", k);

			cout << "===================================================================" << endl;
			cout << "Dense matrix diagonalization with conserved nup (number of particles) and k (momentum)" << endl;

			for (my_int k = 0; k < LatticeSize; k++)
				//for (my_int k = 1; k < 2; k++)
			{
				cout << "-----------------------------------------------------------------" << endl;
				cout << "LatticeSize: " << LatticeSize << endl;
				cout << "No. of up-spins: " << nup << endl;
				cout << "k: " << k << endl;
				// generate basis
				my_int time0 = time(0);
				Basis_k basis_k(LatticeSize, nup, k);
				my_int time1 = time(0);
				cout << "Time(s) to Build Basis: " << time1 - time0 << endl << endl;

				if (0 == k || (0 == LatticeSize % 2 && LatticeSize / 2 == k))
				{
					Ham_HardCoreBoson<double> ham(&basis_k);
					// dense 
					ham.DenseMat_Build_k();
					ham.DenseMat_Eig();

					char fspec[80];
					sprintf(fspec, "spec_k%d.bin", k);
					Vec_fwrite_double(fspec, ham.spec, ham.Dim);
				}
				else {
					Ham_HardCoreBoson<complex<double> > ham(&basis_k);
					// dense 
					ham.DenseMat_Build_k();
					ham.DenseMat_Eig();
					char fspec[80];
					sprintf(fspec, "spec_k%d.bin", k);
					Vec_fwrite_double(fspec, ham.spec, ham.Dim);
				}
			}

			//ham.DenseMat_Static();
		}

		// dense matrix algorithm, conserved nup (number of particles) and k (momentum) and p (parity)
		if (0)
		{
			//GetParaFromInput_int("input.in", "k", k);
			//GetParaFromInput_int("input.in", "p", p);
			cout << "===================================================================" << endl;
			cout << "Dense matrix diagonalization with conserved nup (number of particles) and k (momentum) and p (parity)" << endl;

			for (my_int k = 0; k < LatticeSize / 2 + 1; k++)
				//for (my_int k = 3; k < 4; k++)
				//for (my_int k = 1; k < 2; k++)
			{
				my_int np = (0 == k) || (k * 2 == LatticeSize) ? 2 : 1;

				for (my_int ip = 0; ip < np; ip++)
					//for (my_int ip = 0; ip < 1; ip++)
				{
					my_int p = 1 - 2 * ip;
					cout << "-----------------------------------------------------------------" << endl;
					cout << "LatticeSize: " << LatticeSize << endl;
					cout << "No. of up-spins: " << nup << endl;
					cout << "k: " << k << endl;
					cout << "p: " << p << endl;
					// generate basis
					my_int time0 = time(0);
					Basis_kp basis_kp(LatticeSize, nup, k, p);
					my_int time1 = time(0);
					cout << "Time(s) to Build Basis: " << time1 - time0 << endl << endl;

					Ham_HardCoreBoson<double> ham(&basis_kp);
					// dense 
					ham.DenseMat_Build_kp();
					ham.DenseMat_Eig();

					char fspec[80];
					sprintf(fspec, "spec_k%d_p%d.bin", k, p);
					Vec_fwrite_double(fspec, ham.spec, ham.Dim);
					//ham.DenseMat_Static();

					ham.Calc_Roag();
					sprintf(fspec, "roag_k%d_p%d.bin", k, p);
					Vec_fwrite_double(fspec, ham.roag, ham.Dim - 2);
				}
			}
		}

		// dense matrix algorithm, conserved nup (number of particles) and k (momentum) and p (parity) and z (spin-inversion/particle-hole)
		if (0)
		{
			//GetParaFromInput_int("input.in", "k", k);
			//GetParaFromInput_int("input.in", "p", p);
			cout << "===================================================================" << endl;
			cout << "Dense matrix diagonalization with conserved nup (number of particles) and k (momentum) and p (parity) and z (spin-inversion/particle-hole)" << endl;

			if (nup * 2 != LatticeSize)
			{
				cout << "Error!!! No spin-iversion/particle-hole sysmetry in this sector!";
				exit(1);
			}

			//for (my_int k = 0; k < LatticeSize/2 + 1; k++)
			//for (my_int k = 5; k < 6; k++)
			for (my_int ik = 0; ik < 2; ik++)
			{
				my_int k;
				if (0 == ik) { k = 0; }
				else if (1 == ik) { k = LatticeSize / 2; }

				my_int np = (0 == k) || (k * 2 == LatticeSize) ? 2 : 1;

				//for (my_int ip = 0; ip < np; ip++)
				for (my_int ip = 0; ip < 2; ip++)
				{
					my_int p = 1 - 2 * ip;
					for (my_int iz = 0; iz < 2; iz++)
					//for (my_int iz = 0; iz < 1; iz++)
					{
						my_int z = 1 - 2 * iz;
						cout << "-----------------------------------------------------------------" << endl;
						cout << "LatticeSize: " << LatticeSize << endl;
						cout << "No. of up-spins: " << nup << endl;
						cout << "k: " << k << endl;
						cout << "p: " << p << endl;
						cout << "z: " << z << endl;
						// generate basis
						my_int time0 = time(0);
						Basis_kpz basis_kpz(LatticeSize, nup, k, p, z);
						my_int time1 = time(0);
						cout << "Time(s) to Build Basis: " << time1 - time0 << endl << endl;

						Ham_HardCoreBoson<double> ham(&basis_kpz);
						// dense 
						ham.DenseMat_Build_kpz();
						my_int time2 = time(0);
						cout << "Time(s) to Build Dense Matrix: " << time2 - time1 << endl << endl;

						ham.DenseMat_Eig();
						my_int time3 = time(0);
						cout << "Time(s) to Diagonalize Dense Matrix: " << time3 - time2 << endl << endl;

						ham.DenseMat_Static_kpz();
						my_int time4 = time(0);
						cout << "Time(s) for static observables: " << time4 - time3 << endl << endl;

						if (0 == k && 1 == p && 1 == z)
						{
							//ham.DenseMat_Dynamic_kpz();
							my_int time5 = time(0);
							cout << "Time(s) for dynamic observables: " << time5 - time4 << endl << endl;
						}
					}
				}
			}
		}

	}
	return 0;
}