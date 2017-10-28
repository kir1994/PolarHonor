#include <fstream>
#include <iostream>
#include <complex>
#include <ctime>
#include <algorithm>

#include "../HonorCup/Config.h"

normal_distribution<double> noiseD(0., TARGET_SIGMA);
uniform_int_distribution<int> distr(0., CONSTELLATION_SIZE - 1);

void PrintStat(const vector < complex < double > >& Constellation)
{
	cout << "Power: ";
	cerr << avgPower(Constellation) << "\t";
	cout << endl;
	double m = HUGE_VAL;
	for (auto&& elem : Constellation)
		m = min(m, abs(elem));
	cout << "Min abs: ";
	cerr << m << "\t";
	cout << endl << "Capacity at " << TARGET_SNR << "dB: ";
	cerr << CapacityApprox(Constellation, TARGET_SIGMA, NUM_OF_ITERATIONS) << "\t";
	cout << endl << "SNR for " << TARGET_CAPACITY << " bits/symbol: ";
	cerr << FindSNR(Constellation);
}

void main(int argc, char ** argv)
{
	ifstream ifs(argv[1]);
	double eps = 0.0001;
	string mode = "";
	if (argc > 2)
	{
		eps = stod(argv[2]);
		if (argc > 3)
			mode = argv[3];
	}
	gen.seed(CAPACITY_SEED);
	for (unsigned i = 0; i < 2 * NUM_OF_ITERATIONS; ++i) {
		noiseTable[i] = noiseD(gen);
	}
	for (unsigned i = 0; i < NUM_OF_ITERATIONS; ++i) {
		uniformTable[i] = distr(gen);
	}


	if (!ifs)
		cout << "Cannot open file\n";
	else
	{
		vector < complex < double > > Constellation(CONSTELLATION_SIZE);
		for (auto&& elem : Constellation)
			ifs >> elem._Val[0] >> elem._Val[1];
		cout << "Original constellation:\n";
		PrintStat(Constellation);

		cout << endl << "\nMove point with eps " << eps << ": \n";
		auto tmpConst = Constellation;
		Normalize(tmpConst, MOVE_POINT, eps);
		PrintStat(tmpConst);

		if (mode == "--dumppoint")
			PrintConstellation(tmpConst);
/*
		cout << endl << "\nMove constellation with eps " << eps << ": \n";
		tmpConst = Constellation;
		Normalize(tmpConst, MOVE_CONSTELLATION, eps);
		PrintStat(tmpConst);
		if (mode == "--dumpconst")
			PrintConstellation(tmpConst);*/
	}
}