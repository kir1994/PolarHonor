#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <gsl/gsl_integration.h>
#include <numeric>
#include <algorithm>
#include <ctime>
#include <functional>

const double TARGET_CAPACITY = 3.6;
const unsigned NUM_OF_ITERATIONS = 10000;
const double EPS = 1e-3;
const unsigned CONSTELLATION_SIZE = 16;

const double TARGET_SNR = 12;
const double TARGET_SIGMA = sqrt(1. / (2 * pow(10., TARGET_SNR / 10)));

const double SCALING_FACTOR = 3;

using namespace std;


vector<complex<double> > qam16 =
{ { -3, -3 },{ -3, -1 },
{ -3, 1 },{ -3, 3 },
{ -1, -3 },{ -1, -1 },
{ -1, 1 },{ -1, 3 },
{ 1, -3 },{ 1, -1 },
{ 1, 1 },{ 1, 3 },
{ 3, -3 },{ 3, -1 },
{ 3, 1 },{ 3, 3 } };

//
//
//vector<complex<double> > qam16_4 =
//{ { 1, 1 },{ 1, 3 },
//{ 3, 1 },{ 3, 3 } };

vector<complex<double> > opt16 =
{ { 0.007, 0.767 },{ 0.126 ,0.106 },
{ 0.644, 0.545 },{ 1.279, 0.305 },
{ 0.906 ,-0.771 },{ -1.032, -0.103 },
{ -0.504, 0.332 },{ -0.611, 1.020 },
{ 0.758, -0.119 },{ -0.911, -0.772 },
{ -0.388, -0.329 },{ 0.245, -0.552 },
{ -0.272, -1.001 },{ 0.376, -1.215 },
{ -1.136, 0.571 },{ 0.512, 1.211 } };

double avgPower(const vector<complex<double> >& Constellation)
{
	double res = 0;

	for(auto elem : Constellation)
		res += pow(abs(elem), 2) / Constellation.size();
	
	return res;
}

uniform_real_distribution<double> distr;
auto gen = mt19937_64(time(nullptr));
uniform_real_distribution<double> mutation;
normal_distribution<double> constNoise(0., 0.02);
normal_distribution<double> mutationNoise(0., 0.5);

double CapacityApprox(const vector<complex<double> >& Constellation, double Sigma, unsigned NumOfIterations)
{
	uniform_int_distribution<int> distr(0., Constellation.size() - 1);
	normal_distribution<double> noiseD(0., Sigma);
	vector<double> probs(Constellation.size());

	double globalSum = 0;

	for(unsigned i = 0; i < NumOfIterations; ++i)
	{
		int elem = distr(gen);
		complex<double> noise = { noiseD(gen),  noiseD(gen) };
		complex<double> transmitted = Constellation[elem] + noise;


		double sum = 0;

		for(unsigned j = 0; j < Constellation.size(); ++j)
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

	while(fabs(rightSNR - leftSNR) > EPS)
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


void Normalize(vector<complex<double> >& Constellation)
{
	double P = avgPower(Constellation);
	for (auto&& elem : Constellation)
	{
		elem.real(elem.real() / sqrt(P));
		elem.imag(elem.imag() / sqrt(P));
	}
}


vector<complex<double> > GenerateRandomConstellation(unsigned ConstellationSize)
{
	vector<complex<double> > res(ConstellationSize);

	for(auto&& elem : res)
		elem = { (distr(gen) - 0.5) * SCALING_FACTOR,  (distr(gen) - 0.5) * SCALING_FACTOR };
	Normalize(res);
	return res;
}


vector<complex<double> > GenerateNoisyConstellation(unsigned ConstellationSize, unsigned Ind)
{
	if (Ind % 10 >= 2)
		return GenerateRandomConstellation(ConstellationSize);
	vector<complex<double> > res = (Ind % 2 == 0) ? qam16 : opt16;

	for (auto&& elem : res)
	{
		elem.real(elem.real() + constNoise(gen));
		elem.imag(elem.imag() + constNoise(gen));
	}
	Normalize(res);
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

#include <signal.h>



vector<complex<double> > ReconstructConstellation(const vector<complex<double> >& Constellation)
{
	vector<complex<double> > res;
	res.reserve(Constellation.size() * 4);

	for(auto&& elem : Constellation)
	{
		res.push_back({  elem.real(),  elem.imag() });
		res.push_back({  elem.real(), -elem.imag() });
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
		population[i].second = GenerateNoisyConstellation(ConstellationSize, i);//GenerateRandomConstellation(ConstellationSize);
		population[i].first = CapacityApprox(population[i].second, TARGET_SIGMA, NUM_OF_ITERATIONS);//FindSNR(population[i]);
	}

	sort(population.begin(), population.end(), std::greater<pair<double, vector<complex<double> > > >());
	

	vector<unsigned> permutation(PopulationSize);
	iota(permutation.begin(), permutation.end(), 0);

	vector<complex<double> > tempConstellation(ConstellationSize);
	double curBest = population.begin()->first;
	unsigned cnt = 0;
	bool flag = false;
	for(unsigned i = 0; F && i < NumOfGenerations; ++i)
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
					if(!UseBest)
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
			Normalize(tempConstellation);
			double tempMetric = CapacityApprox(tempConstellation, TARGET_SIGMA, NUM_OF_ITERATIONS);

			if (tempMetric > population[k].first)
			{
				population[k] = { tempMetric, tempConstellation };
				flag = true;
			}
		}
		sort(population.begin(), population.end(), std::greater<pair<double, vector<complex<double> > > >());
		if (i % 1 == 0)
			cout << "Generation " << i << ": best metric " << population.begin()->first << endl;//*max_element/*min_element*/(metrics.begin(), metrics.end()) << endl;
		/*if (curBest == population.begin()->first)
			cnt++;
		else
			cnt = 0;
		if(cnt == 20)*/		
		if (flag == false)
		{
			cout << "Generation " << i << ": convergence\n";
			break;
		}
	}
	//auto&& bestMetric = /*min_element*/max_element(metrics.begin(), metrics.end());
	return population.begin()->second;//population[distance(metrics.begin(), bestMetric)];
}

#include <fstream>

int main(int argc, char** argv)
{
	signal(SIGINT, handler);
	Normalize(qam16);
	//Normalize(qam16_4);
	Normalize(opt16);
	int population = stoi(argv[1]);
	float CR = stof(argv[2]);
	float F = stof(argv[3]);
	string mode = argv[4];
	//for (auto&& elem : opt16)
	//{
	//	cout << elem.real() << ' ' << elem.imag() << ' ';//endl;
	//}
	
	//cout << findSNR(qam16) << ' ' << findSNR(opt16);

	cout << "QAM16: " << CapacityApprox(qam16, TARGET_SIGMA, NUM_OF_ITERATIONS)
		<< ", OptimalConstellation: " << CapacityApprox(opt16, TARGET_SIGMA, NUM_OF_ITERATIONS) << endl;
	auto&& BestConstellation = DEConstellationSearch(CONSTELLATION_SIZE, mode == "best",
		population, CR, F, 10000);

	cout << CapacityApprox(BestConstellation, TARGET_SIGMA, NUM_OF_ITERATIONS);

	ofstream out("res" + to_string(CR) + "_" + to_string(F) + "_" + mode + ".txt");
	for (auto&& elem : BestConstellation)
		out << elem.real() << ' ' << elem.imag() << ' ';//endl;

	return 0;
}