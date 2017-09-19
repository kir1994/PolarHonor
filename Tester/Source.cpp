#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>

using namespace std;

uniform_real_distribution<double> distr;
auto gen = mt19937_64(time(nullptr));

const double TARGET_CAPACITY = 3.6;
const unsigned NUM_OF_ITERATIONS = 100000;
const double EPS = 1e-4;
const unsigned CONSTELLATION_SIZE = 16;

double avgPower(const vector<complex<double> >& Constellation)
{
	double res = 0;

	for (auto elem : Constellation)
		res += pow(abs(elem), 2) / Constellation.size();

	return res;
}



double CapacityApprox(const vector<complex<double> >& Constellation, double Sigma, unsigned NumOfIterations)
{
	uniform_int_distribution<int> distr(0., Constellation.size() - 1);
	normal_distribution<double> noiseD(0., Sigma);
	vector<double> probs(Constellation.size());

	double globalSum = 0;

	for (unsigned i = 0; i < NumOfIterations; ++i)
	{
		int elem = distr(gen);
		complex<double> noise = { noiseD(gen),  noiseD(gen) };
		complex<double> transmitted = Constellation[elem] + noise;


		double sum = 0;

		for (unsigned j = 0; j < Constellation.size(); ++j)
		{
			probs[j] = exp(-(pow(abs(transmitted - Constellation[j]), 2) - pow(abs(noise), 2)) / (2 * Sigma * Sigma));
			sum += probs[j];
		}
		for (unsigned j = 0; j < Constellation.size(); ++j)
			probs[j] /= sum;

		double mutualInf = 0.;

		for (unsigned j = 0; j < Constellation.size(); ++j)
			if (probs[j] != 0)
				mutualInf += probs[j] * log2(probs[j]);
		globalSum += mutualInf;
	}

	return globalSum / NumOfIterations + log2(Constellation.size());
}

double FindSNR(const vector<complex<double> >& Constellation)
{
	double leftSNR = 10.;
	double rightSNR = 25;

	while (fabs(rightSNR - leftSNR) > EPS)
	{
		double SNR = (rightSNR + leftSNR) / 2;
		double Sigma = sqrt(1. / (2 * pow(10., SNR / 10)));
		double C = CapacityApprox(Constellation, Sigma, NUM_OF_ITERATIONS);

		if (C - TARGET_CAPACITY > EPS)
			rightSNR = SNR;
		else if (TARGET_CAPACITY - C > EPS)
			leftSNR = SNR;
		else
			return SNR;
	}
	return leftSNR;
}


void main(int argc, char ** argv)
{
	ifstream ifs(argv[1]);

	if (!ifs)
		cout << "Cannot open file\n";
	else
	{
		vector < complex < double > > Constellation(CONSTELLATION_SIZE);
		for (auto&& elem : Constellation)
			ifs >> elem._Val[0] >> elem._Val[1];
		cout << "Power: " << avgPower(Constellation) << endl;
		double m = 100.;
		for (auto&& elem : Constellation)
			m = min(m, abs(elem));
		cout << "Min abs: " << m << endl;
		cout << FindSNR(Constellation);
	}
}