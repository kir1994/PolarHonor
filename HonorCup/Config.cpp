#include <algorithm>
#include <iostream>
#include <functional>
#include <ctime>
#include <numeric>

#include "Config.h"

uniform_real_distribution<double> mutation;
normal_distribution<double> constNoise(0., 0.02);
normal_distribution<double> mutationNoise(0., 0.5);


normal_distribution<double> randConst(0., TARGET_SIGMA);

mt19937_64 gen = mt19937_64(time(nullptr));

vector<double> noiseTable(NUM_OF_ITERATIONS * 2);
vector<unsigned> uniformTable(NUM_OF_ITERATIONS);
uniform_real_distribution<double> distrR;
double probs[CONSTELLATION_SIZE];

vector<complex<double> > qam64 =
{ { -7,-7 },{ -7,-5 },{ -7,-3 },{ -7,-1 },{ -7,1 },{ -7,3 },{ -7,5 },{ -7,7 },{ -5,-7 },{ -5,-5 },{ -5,-3 },{ -5,-1 },{ -5,1 },{ -5,3 },{ -5,5 },{ -5,7 },{ -3,-7 },{ -3,-5 },{ -3,-3 },{ -3,-1 },{ -3,1 },{ -3,3 },{ -3,5 },{ -3,7 },{ -1,-7 },{ -1,-5 },{ -1,-3 },{ -1,-1 },{ -1,1 },{ -1,3 },{ -1,5 },{ -1,7 },{ 1,-7 },{ 1,-5 },{ 1,-3 },{ 1,-1 },{ 1,1 },{ 1,3 },{ 1,5 },{ 1,7 },{ 3,-7 },{ 3,-5 },{ 3,-3 },{ 3,-1 },{ 3,1 },{ 3,3 },{ 3,5 },{ 3,7 },{ 5,-7 },{ 5,-5 },{ 5,-3 },{ 5,-1 },{ 5,1 },{ 5,3 },{ 5,5 },{ 5,7 },{ 7,-7 },{ 7,-5 },{ 7,-3 },{ 7,-1 },{ 7,1 },{ 7,3 },{ 7,5 },{ 7,7 }};

vector<complex<double>> opt64 = 
{
	{ -0.250063, 1.03656 },{ 0.616321, 1.10861 },{ -0.252301, 0.328461 },{ -0.716828, 0.425236 },{ -0.582291, -0.784889 },{ -0.473162, -1.23552 },{ 0.163801, -0.0200623 },{ -0.702416, -0.0332262 },{ 0.341103, -0.117253 },{ 0.581277, 0.226812 },{ -0.468035, 0.053457 },{ 0.434349, -0.380458 },{ 0.202955, 1.04805 },{ 1.29581, 0.570222 },{ 1.07562, -1.15493 },{ -0.69784, -0.423914 },{ -0.0747872, -1.38111 },{ 0.843034, -0.277044 },{ -0.106171, -0.632955 },{ 1.1348, -0.192018 },{ 0.0144637, -0.263539 },{ 0.049037, 0.276871 },{ -0.526875, 1.39922 },{ 0.319336, 0.423293 },{ -0.138425, -0.944274 },{ 0.0586696, 0.480996 },{ -1.31638, -0.872515 },{ 0.189676, -0.377574 },{ 0.365309, 0.20641 },{ -0.92854, -0.992313 },{ 0.927081, -0.657673 },{ -0.472756, 0.320681 },{ -0.990197, 0.2213 },{ -1.39549, 0.376179 },{ -0.403325, -0.316222 },{ 1.10054, 0.978196 },{ 0.74698, -0.910935 },{ -0.278611, -0.0357268 },{ 0.868738, 0.653061 },{ -1.40977, -0.121128 },{ 0.0561273, 0.777762 },{ 0.640633, -0.116556 },{ 0.0563964, 1.42328 },{ 0.558984, -0.555844 },{ -0.804661, 1.0082 },{ -0.553244, 0.697041 },{ 1.41616, -0.518016 },{ -0.982134, -0.122245 },{ -1.07509, -0.536251 },{ 0.087592, -0.663674 },{ -0.365676, -0.52431 },{ -1.0687, 0.664589 },{ 0.402726, -0.75533 },{ 0.0678338, 0.0902357 },{ 0.646528, 0.511533 },{ 0.393079, 0.648807 },{ 0.521547, -1.34625 },{ 1.25941, 0.190665 },{ -0.140353, -0.111067 },{ 0.603982, 1.48502 },{ 0.24352, -1.04594 },{ 0.890294, 0.162086 },{ 0.487003, 0.853546 },{ -0.21764, 0.673205 } 
};


vector<complex<double>> opt2_64 = 
{
	{ -0.245758, 1.0134 },{ 0.703374, 1.16075 },{ -0.199034, 0.356605 },{ -0.725156, 0.480727 },{ -0.632054, -0.802358 },{ -0.507564, -1.24854 },{ 0.168489, -0.0655306 },{ -0.706528, -0.0687197 },{ 0.364601, -0.120937 },{ 0.688274, 0.170619 },{ -0.502382, 0.124408 },{ 0.452829, -0.394196 },{ 0.214005, 1.0454 },{ 1.26711, 0.555059 },{ 1.02046, -1.11093 },{ -0.632725, -0.368837 },{ -0.0712263, -1.40697 },{ 0.859338, -0.310685 },{ -0.113349, -0.540222 },{ 1.05981, -0.134143 },{ 0.00409272, -0.229995 },{ 0.0907329, 0.143148 },{ -0.562824, 1.34602 },{ 0.259541, 0.386315 },{ -0.214315, -0.934153 },{ -0.00427375, 0.484131 },{ -1.29384, -0.70881 },{ 0.201766, -0.394334 },{ 0.377331, 0.198263 },{ -0.965051, -1.05636 },{ 0.996955, -0.679078 },{ -0.435519, 0.376977 },{ -0.891441, 0.225431 },{ -1.32948, 0.306052 },{ -0.388175, -0.302448 },{ 1.07848, 0.988834 },{ 0.614565, -0.959355 },{ -0.31673, -0.0111726 },{ 0.838835, 0.625998 },{ -1.42492, -0.171445 },{ 0.00745251, 0.78068 },{ 0.638838, -0.113789 },{ -0.0810909, 1.38746 },{ 0.651755, -0.564092 },{ -1.03501, 1.04751 },{ -0.648456, 0.83559 },{ 1.39784, -0.444476 },{ -1.01303, -0.106056 },{ -0.976447, -0.497699 },{ 0.0802188, -0.687493 },{ -0.369111, -0.622045 },{ -1.07591, 0.625475 },{ 0.323827, -0.708584 },{ -0.103529, 0.115942 },{ 0.588362, 0.460224 },{ 0.322188, 0.629099 },{ 0.468233, -1.37371 },{ 1.41008, 0.0858935 },{ -0.192074, -0.248627 },{ 0.375961, 1.43902 },{ 0.171278, -1.04894 },{ 0.977562, 0.22119 },{ 0.488196, 0.851128 },{ -0.353955, 0.687761 } 
};


double avgPower(const vector<complex<double> >& Constellation)
{
	double res = 0;

	for (auto elem : Constellation)
		res += pow(abs(elem), 2) / Constellation.size();

	return res;
}


double CapacityApprox(const vector<complex<double> >& Constellation, double Sigma, unsigned NumOfIterations)
{
	double globalSum = 0;
	double N0 = 2 * Sigma * Sigma;

	for (unsigned i = 0; i < NumOfIterations; ++i)
	{
		int elem = uniformTable[i];
		complex<double> noise = { noiseTable[i << 1],  noiseTable[(i << 1) + 1] };
		complex<double> transmitted = Constellation[elem] + noise;

		double sum = 0;
		double an = abs(noise);
		an *= an;

		for (unsigned j = 0; j < CONSTELLATION_SIZE; ++j)
		{
			double c = abs(transmitted - Constellation[j]);
			c *= c;
			probs[j] = exp((an - c) / N0);
			sum += probs[j];
		}
		for (unsigned j = 0; j < CONSTELLATION_SIZE; ++j)
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


void Normalize(vector<complex<double> >& Constellation, handle_bad_points mode, double eps)
{
	if (mode != DO_NOTHING)
		Normalize(Constellation, DO_NOTHING);
	vector<complex<double> >::iterator curMin = Constellation.begin();
	double m = HUGE_VAL;
	for (auto&& elem = Constellation.begin(); elem != Constellation.end(); ++elem)
	{
		double tmp = abs(*elem);
		if (tmp < m)
		{
			curMin = elem; 
			m = tmp;
		}
	}

	double P = avgPower(Constellation);
	double sP = sqrt(P);

	if (m / sP < SMALL_MODULE)
	{
		double scale = (SMALL_MODULE + eps) / m;
		switch (mode)
		{
		case MOVE_POINT: 
			*curMin *= scale; 
			P = avgPower(Constellation);
			sP = sqrt(P);
			break;
		case MOVE_CONSTELLATION: 
		case DO_NOTHING:
		default: ;
		}
	}

	for (auto&& elem : Constellation)
		elem /= sP;
}





vector<complex<double> > GenerateRandomConstellation(unsigned ConstellationSize)
{
	vector<complex<double> > res(ConstellationSize);

	for (auto&& elem : res)
		elem = { randConst(gen) * SCALING_FACTOR,  randConst(gen) * SCALING_FACTOR };
	Normalize(res);
	return res;
}

vector<complex<double> > GenerateNoisyConstellation(unsigned ConstellationSize, unsigned Ind)
{
	if (Ind % 6 >= 2)
		return GenerateRandomConstellation(ConstellationSize);
	vector<complex<double> > res;
	switch (Ind % 6) {
	case 0:
		res = opt2_64;
		break;
	case 1:
		res = opt64;
		break;
	//case 2:
	//	res = circle16;
	//	break;
	//case 3:
	//	res = apsk12_4;
	//	break;
	//case 4:
	//	res = qam16;
	//	break;
	//case 5:
	//	res = best16;
	//	break;
	default:
		res = qam64;
	}

	for (auto&& elem : res)
	{
		elem.real(elem.real() + constNoise(gen));
		elem.imag(elem.imag() + constNoise(gen));
	}
	Normalize(res, MOVE_POINT);
	return res;
}

bool operator <(const vector<complex<double> >& a, const vector<complex<double> >& b)
{
	return a[0].real() < b[0].real();
}

volatile bool F = true;


void handler(int sig)
{
	F = false;
}

vector<complex<double> > ReconstructConstellation(const vector<complex<double> >& Constellation)
{
	vector<complex<double> > res;
	res.reserve(Constellation.size() * 4);

	for (auto&& elem : Constellation)
	{
		res.push_back({ elem.real(),  elem.imag() });
		res.push_back({ elem.real(), -elem.imag() });
		res.push_back({ -elem.real(),  elem.imag() });
		res.push_back({ -elem.real(), -elem.imag() });
	}

	return res;
}

vector<complex<double> > DEConstellationSearch(
	unsigned ConstellationSize, bool UseBest,
	unsigned PopulationSize, double CrossoverProbability, double AmplificationFactor,
	unsigned NumOfGenerations)
{
	vector<pair<double, vector<complex<double> > > > population(PopulationSize);

	for (unsigned i = 0; i < PopulationSize; ++i)
	{
		population[i].second = GenerateNoisyConstellation(ConstellationSize, i);
		population[i].first = CapacityApprox(population[i].second, TARGET_SIGMA, NUM_OF_ITERATIONS);
	}

	sort(population.begin(), population.end(), std::greater<pair<double, vector<complex<double> > > >());


	vector<unsigned> permutation(PopulationSize);
	iota(permutation.begin(), permutation.end(), 0);

	vector<complex<double> > tempConstellation(ConstellationSize);
	double curBest = population.begin()->first;
	unsigned cnt = 0;
	bool flag = false;
	for (unsigned i = 0; F && i < NumOfGenerations; ++i)
	{
		flag = false;
		for (unsigned k = 0; k < PopulationSize; ++k)
		{
			auto&& pos = find(permutation.begin(), permutation.end(), k);
			swap(*permutation.rbegin(), *pos);
			shuffle(permutation.begin(), prev(permutation.end()), gen);
			for (unsigned j = 0; j < ConstellationSize; ++j)
			{
				if (mutation(gen) > CrossoverProbability)
					tempConstellation[j] = population[k].second[j];
				else
				{
					if (!UseBest)
					{
						tempConstellation[j] = population[permutation[0]].second[j] +
							AmplificationFactor * (population[permutation[1]].second[j] - population[permutation[2]].second[j]);
					}
					else
					{
						tempConstellation[j] = population[0].second[j] +
							AmplificationFactor *
							(population[permutation[0]].second[j] + population[permutation[1]].second[j]
								- population[permutation[2]].second[j] - population[permutation[3]].second[j]);
					}
				}
			}
			Normalize(tempConstellation, MOVE_POINT);
			double tempMetric = CapacityApprox(tempConstellation, TARGET_SIGMA, NUM_OF_ITERATIONS);

			if (tempMetric > population[k].first)
			{
				population[k] = { tempMetric, tempConstellation };
				flag = true;
			}
		}
		sort(population.begin(), population.end(), greater<pair<double, vector<complex<double> > > >());
		if (i % 1 == 0)
			cout << "Generation " << i << ": best metric " << population.begin()->first << endl;
		if (flag == false)
		{
			cout << "Generation " << i << ": convergence\n";
			break;
		}
	}
	return population.begin()->second;
}


void PrintConstellation(const vector<complex<double> >& Constellation)
{
	cout << endl;
	for (auto&& elem : Constellation)
		cout << elem.real() << ' ' << elem.imag() << ' ';
	cout << endl;
}


//#include <dlib/optimization.h>
////using namespace dlib;
//
//typedef dlib::matrix<double, 0, 1> column_vector;
//
//
//double CapacityApproxDlib(const column_vector& constellation)
//{
//	vector<complex<double>> Constellation(CONSTELLATION_SIZE);
//	for (unsigned i = 0; i < CONSTELLATION_SIZE; ++i)
//		Constellation[i] = { constellation(2 * i), constellation(2 * i + 1)};
//	Normalize(Constellation, MOVE_POINT);
//
//	double globalSum = 0;
//	double N0 = 2 * TARGET_SIGMA * TARGET_SIGMA;
//
//	for (unsigned i = 0; i < NUM_OF_ITERATIONS; ++i)
//	{
//		int elem = uniformTable[i];
//		complex<double> noise = { noiseTable[i << 1],  noiseTable[(i << 1) + 1] };
//		complex<double> transmitted = Constellation[elem] +noise;
//
//		double sum = 0;
//		double an = abs(noise);
//		an *= an;
//
//		for (unsigned j = 0; j < CONSTELLATION_SIZE; ++j)
//		{
//			double c = abs(transmitted - Constellation[j]);
//			c *= c;
//			probs[j] = exp((an - c) / N0);
//			sum += probs[j];
//		}
//		for (unsigned j = 0; j < CONSTELLATION_SIZE; ++j)
//			probs[j] /= sum;
//
//		double mutualInf = 0.;
//
//		for (unsigned j = 0; j < Constellation.size(); ++j)
//			if (probs[j] != 0)
//				mutualInf += probs[j] * log2(probs[j]);
//		globalSum += mutualInf;
//	}
//
//	return -(globalSum / NUM_OF_ITERATIONS + log2(Constellation.size()));
//}
//

//vector<complex<double> > dlibMinSearch(unsigned ConstellationSize, unsigned ind, unsigned mode)
//{
//	auto tmp = GenerateNoisyConstellation(ConstellationSize, ind);
//		//GenerateRandomConstellation(ConstellationSize);
//	column_vector starting_point(ConstellationSize * 2);
//	for(unsigned i = 0; i < ConstellationSize; ++i)
//	{
//		starting_point(2 * i) = tmp[i].real();
//		starting_point(2 * i + 1) = tmp[i].imag();
//	}
//
//	column_vector zero_vector(ConstellationSize * 2);
//	for (unsigned i = 0; i < ConstellationSize * 2; ++i)
//		zero_vector(i) = -2;
//	column_vector mega_vector(ConstellationSize * 2);
//	for (unsigned i = 0; i < ConstellationSize * 2; ++i)
//		mega_vector(i) = 2.;
//	
//	switch(mode)
//	{
//	case 0:
//		dlib::find_min_using_approximate_derivatives(dlib::bfgs_search_strategy(),
//			dlib::objective_delta_stop_strategy(1e-7).be_verbose(),
//			CapacityApproxDlib,
//			starting_point, -10);
//		break;
//	case 1:
//		dlib::find_min_using_approximate_derivatives(dlib::lbfgs_search_strategy(100),
//			dlib::objective_delta_stop_strategy(1e-7).be_verbose(),
//			CapacityApproxDlib,
//			starting_point, -10);
//		break;
//	case 2:
//		dlib::find_min_using_approximate_derivatives(dlib::cg_search_strategy(),
//			dlib::objective_delta_stop_strategy(1e-7).be_verbose(),
//			CapacityApproxDlib,
//			starting_point, -10);
//		break;
//	default:
//		break;
//	}
//
//	for (unsigned i = 0; i < ConstellationSize; ++i)
//		tmp[i] = { starting_point(2 * i), starting_point(2 * i + 1) };
//
//	return tmp;
//}

