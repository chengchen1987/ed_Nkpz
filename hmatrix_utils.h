#pragma once


template <class dtype>
vector<my_int> Ham_HardCoreBoson<dtype>::Get_IniStates()
{
	vector<my_int> ind_vec;
	// read iniinfo if "evo_iniinfo.dat" exists
	ifstream fini("evo_iniinfo.dat", ios::in);
	if (fini)
	{
		cout << "Read initial states from evo_iniinfo.dat" << endl;
		// count lines 
		string line;
		my_int nl = 0;
		while (!fini.eof())
		{
			getline(fini, line);
			if (0 != line.size()) nl++;
		}
		fini.close();
		// read initial informations 
		double aux;
		ifstream fini("evo_iniinfo.dat", ios::in);
		for (my_int i = 0; i < nl; i++)
		{
			my_int s0;
			fini >> aux >> aux >> s0 >> aux >> aux;
			ind_vec.push_back(basis->get_index(s0));
		}
		fini.close();
	}
	else {

		// choose random initial produuct states according to its energy
		// sort product states byy energies
		vector <pair<double, my_int> > Fock_E_n;  // stores E_fock, and index of this Fock state
		for (my_int s = 0; s < Dim; s++)
		{
			pair<double, my_int> aux = make_pair(H_diag[s], s);
			Fock_E_n.push_back(aux);
		}
		std::sort(Fock_E_n.begin(), Fock_E_n.end());

		my_int ind_0, ind_1;
		// product-state spectrum in fixed energy density window
		if (0)
		{
			double epsi = 0.5;
			double E_min = spec[0];
			double E_max = spec[Dim - 1];
			double* sorted_Fock_epsi = new double[Dim];
			for (my_int s = 0; s < Dim; s++) {
				sorted_Fock_epsi[s] = (Fock_E_n[s].first - E_min) / (E_max - E_min);
			}
			ind_0 = Get_TargetFock_left((epsi - 0.02) * (E_max - E_min) + E_min, Fock_E_n);
			ind_1 = Get_TargetFock_right((epsi + 0.02) * (E_max - E_min) + E_min, Fock_E_n);
		}

		// product-state spectrum in fixed index range
		if (1)
		{
			srand(time(0));
			ind_0 = Dim / 4;
			ind_1 = Dim / 4 * 3;
		}

		// choose random initial states
		my_int nini = Params.evo_nini < (ind_1 - ind_0) ? Params.evo_nini : ind_1 - ind_0;
		vector<my_int> aux_vec;
		aux_vec.resize(ind_1 - ind_0);
		for (my_int i = 0; i < ind_1 - ind_0; i++) { aux_vec[i] = ind_0 + i; }
		std::random_shuffle(aux_vec.begin(), aux_vec.end());

		for (my_int i = 0; i < nini; i++) ind_vec.push_back(Fock_E_n[aux_vec[i]].second);
		aux_vec.clear();
	}
	return ind_vec;
}

template <class dtype>
vector<my_int> Ham_HardCoreBoson<dtype>::Get_IniStates_kpz()
{
	vector<my_int> ind_vec;
	// read iniinfo if "evo_iniinfo.dat" exists
	ifstream fini("evo_iniinfo.dat", ios::in);
	if (fini)
	{
		cout << "Read initial states from evo_iniinfo.dat" << endl;
		// count lines 
		string line;
		my_int nl = 0;
		while (!fini.eof())
		{
			getline(fini, line);
			if (0 != line.size()) nl++;
		}
		fini.close();
		// read initial informations 
		double aux;
		ifstream fini("evo_iniinfo.dat", ios::in);
		for (my_int i = 0; i < nl; i++)
		{
			my_int s0;
			fini >> aux >> aux >> s0 >> aux >> aux;
			my_int s, la, qa, ga;
			basis_kpz->Representative(s0, s, la, qa, ga);
			ind_vec.push_back(basis_kpz->get_index(s));
		}
		fini.close();
	}
	else {

		// choose random initial produuct states according to its energy
		// sort product states byy energies
		vector <pair<double, my_int> > Fock_E_n;  // stores E_fock, and index of this Fock state
		for (my_int s = 0; s < Dim; s++)
		{
			pair<double, my_int> aux = make_pair(H_diag[s], s);
			Fock_E_n.push_back(aux);
		}
		std::sort(Fock_E_n.begin(), Fock_E_n.end());

		my_int ind_0, ind_1;
		// product-state spectrum in fixed energy density window
		if (0)
		{
			double epsi = 0.5;
			double E_min = spec[0];
			double E_max = spec[Dim - 1];
			double* sorted_Fock_epsi = new double[Dim];
			for (my_int s = 0; s < Dim; s++) {
				sorted_Fock_epsi[s] = (Fock_E_n[s].first - E_min) / (E_max - E_min);
			}
			ind_0 = Get_TargetFock_left((epsi - 0.02) * (E_max - E_min) + E_min, Fock_E_n);
			ind_1 = Get_TargetFock_right((epsi + 0.02) * (E_max - E_min) + E_min, Fock_E_n);
		}

		// product-state spectrum in fixed index range
		if (1)
		{
			srand(time(0));
			ind_0 = Dim / 4;
			ind_1 = Dim / 4 * 3;
		}

		// choose random initial states
		my_int nini = Params.evo_nini < (ind_1 - ind_0) ? Params.evo_nini : ind_1 - ind_0;
		vector<my_int> aux_vec;
		aux_vec.resize(ind_1 - ind_0);
		for (my_int i = 0; i < ind_1 - ind_0; i++) { aux_vec[i] = ind_0 + i; }
		std::random_shuffle(aux_vec.begin(), aux_vec.end());

		for (my_int i = 0; i < nini; i++) ind_vec.push_back(Fock_E_n[aux_vec[i]].second);
		aux_vec.clear();
	}
	return ind_vec;
}

template <class dtype>
vector<my_int> Ham_HardCoreBoson<dtype>::Get_Spanned_States(const my_int& s0)
{
	vector<my_int> ind_vec;
	my_int L = LatticeSize;
	// ind vec stores all states that construct the a kpz state, k = 0, p = 1, z = 1
	my_int ns = L * 2 * 2;
	vector<double> coeff_vec;
	double coeff = 1;

	my_int s_t = s0;
	for (my_int r = 0; r < L; r++)
	{
		// T^r
		if (0 < r) { s_t = Bits_CycleLeft(s_t, L); }
		my_int check = 0;
		for (my_int ix = 0; ix < ind_vec.size(); ix++)
		{
			if (basis->get_index(s_t) == ind_vec[ix])
			{
				check++;
				coeff_vec[ix] += 1;
				break;
			}
		}
		if (0 == check)
		{
			ind_vec.push_back(basis->get_index(s_t));
			coeff_vec.push_back(coeff);
		}
		// T^r * Z
		my_int s_tz = Bits_Invert(s_t, L);
		check = 0;
		for (my_int ix = 0; ix < ind_vec.size(); ix++)
		{
			if (basis->get_index(s_tz) == ind_vec[ix])
			{
				check++;
				coeff_vec[ix] += 1;
				break;
			}
		}
		if (0 == check)
		{
			ind_vec.push_back(basis->get_index(s_tz));
			coeff_vec.push_back(coeff);
		}
		// T^r * P
		my_int s_tp = Bits_Reflect(s_t, L);
		check = 0;
		for (my_int ix = 0; ix < ind_vec.size(); ix++)
		{
			if (basis->get_index(s_tp) == ind_vec[ix])
			{
				check++;
				coeff_vec[ix] += 1;
				break;
			}
		}
		if (0 == check)
		{
			ind_vec.push_back(basis->get_index(s_tp));
			coeff_vec.push_back(coeff);
		}
		// T^r * P * Z
		my_int s_tpz = Bits_Reflect(s_tz, L);
		check = 0;
		for (my_int ix = 0; ix < ind_vec.size(); ix++)
		{
			if (basis->get_index(s_tpz) == ind_vec[ix])
			{
				check++;
				coeff_vec[ix] += 1;
				break;
			}
		}
		if (0 == check)
		{
			ind_vec.push_back(basis->get_index(s_tpz));
			coeff_vec.push_back(coeff);
		}
	}
	ns = ind_vec.size();

	cout << "spanned states: " << endl;
	for (my_int i = 0; i < ns; i++)
	{
		Bits_Print(basis->get_state(ind_vec[i]), L);
		cout << "  " << basis->get_state(ind_vec[i]) << " coeff: " << coeff_vec[i] << endl;
	}
	coeff_vec.clear();

	return ind_vec;
}