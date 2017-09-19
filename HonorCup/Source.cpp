#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <numeric>
#include <algorithm>
#include <ctime>
#include <functional>
#include <signal.h>

const double TARGET_CAPACITY = 3.6;
const unsigned NUM_OF_ITERATIONS = 10000;
const double EPS = 1e-3;
const unsigned CONSTELLATION_SIZE = 16;

const double TARGET_SNR = 12;
const double TARGET_SIGMA = sqrt(1. / (2 * pow(10., TARGET_SNR / 10)));

const double SCALING_FACTOR = 3;
const double CAPACITY_SEED = 123456487;

using namespace std;
auto gen = mt19937_64(time(nullptr));

vector<double> noiseTable(NUM_OF_ITERATIONS * 2);
vector<unsigned> uniformTable(NUM_OF_ITERATIONS);
normal_distribution<double> noiseD(0., TARGET_SIGMA);
uniform_int_distribution<int> distr(0., CONSTELLATION_SIZE - 1);
uniform_real_distribution<double> distrR;
double probs[16];

vector<complex<double> > qam16 =
{ { -3, -3 },{ -3, -1 },
{ -3, 1 },{ -3, 3 },
{ -1, -3 },{ -1, -1 },
{ -1, 1 },{ -1, 3 },
{ 1, -3 },{ 1, -1 },
{ 1, 1 },{ 1, 3 },
{ 3, -3 },{ 3, -1 },
{ 3, 1 },{ 3, 3 } };

vector<complex<double> > apsk12_4 =
{ { 0.7071067812, 0.7071067812 },
{ -0.7071067812, 0.7071067812 },
{ -0.7071067812, -0.7071067812 },
{ 0.7071067812, -0.7071067812 },
{ 2.511407148, 0.6729295173 },
{ 1.838477631, 1.838477631 },
{ 0.6729295173, 2.511407148 },
{ -0.6729295173, 2.511407148 },
{ -1.838477631, 1.838477631 },
{ -2.511407148, 0.6729295173 },
{ -2.511407148, -0.6729295172 },
{ -1.838477631, -1.838477631 },
{ -0.6729295173, -2.511407148 },
{ 0.6729295173, -2.511407148 },
{ 1.838477631, -1.838477631 },
{ 2.511407148, -0.6729295173 } };

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

vector<complex<double> > hex16 =
{ {-1.183215957, 0},
{ -0.5070925528, 0 },
{ 0.1690308509, 0 },
{ 0.8451542547, 0 },
{ -0.8451542547, -0.5855400438 },
{ -0.1690308509, -0.5855400438 },
{ 0.5070925528, -0.5855400438 },
{ 1.183215957, -0.5855400438 },
{ -0.8451542547, 0.5855400438 },
{ -0.1690308509, 0.5855400438 },
{ 0.5070925528, 0.5855400438 },
{ 1.183215957, 0.5855400438 },
{ -0.5070925528, 1.171080088 },
{ 0.1690308509, 1.171080088 },
{ -0.5070925528, -1.171080088 },
{ 0.1690308509, -1.171080088 } };

vector<complex<double> > circle16 =
{ {0.000000000000000000000000000000, -0.783335257075577578989352063067},
{ -0.416424118540487841572566999821,  -0.663479523780077891110760966724 },
{ 0.416424118540487841572566999821,  -0.663479523780077891110760966724 },
{ -0.705416674059201028176189884850,  -0.340589842680189381976367195063 },
{ 0.705416674059201028176189884850,  -0.340589842680189381976367195063 },
{ -0.216664742924422421010647936933,  -0.278940013248970164633197783693 },
{ 0.216664742924422421010647936933,  -0.278940013248970164633197783693 },
{ -0.778541931420337049205105404241,   0.086525059941917500142787551612 },
{ 0.778541931420337049205105404241,   0.086525059941917500142787551612 },
{ -0.347827782219485848085061235425,   0.134062045376535804555624718135 },
{ 0.347827782219485848085061235425,   0.134062045376535804555624718135 },
{ 0.000000000000000000000000000000,   0.392500195080286289189873104776 },
{ -0.613422546404355735310803524023,   0.487162092675997572917936613924 },
{ 0.613422546404355735310803524023,   0.487162092675997572917936613924 },
{ -0.260587343786171975082805234959,   0.738720759987242156057605439640 },
{ 0.260587343786171975082805234959,   0.738720759987242156057605439640 }};

double avgPower(const vector<complex<double> >& Constellation)
{
	double res = 0;

	for(auto elem : Constellation)
		res += pow(abs(elem), 2) / Constellation.size();
	
	return res;
}

uniform_real_distribution<double> mutation;
normal_distribution<double> constNoise(0., 0.02);
normal_distribution<double> mutationNoise(0., 0.5);

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
	double sP = sqrt(P);
	for (auto&& elem : Constellation)
	{
		elem /= sP;
	}
}


vector<complex<double> > GenerateRandomConstellation(unsigned ConstellationSize)
{
	vector<complex<double> > res(ConstellationSize);

	for (auto&& elem : res)
		elem = { (distrR(gen) - 0.5) * SCALING_FACTOR,  (distrR(gen) - 0.5) * SCALING_FACTOR };
	Normalize(res);
	return res;
}

vector<complex<double> > GenerateNoisyConstellation(unsigned ConstellationSize, unsigned Ind)
{
	if (Ind % 10 >= 5)
		return GenerateRandomConstellation(ConstellationSize);
	vector<complex<double> > res;
	switch (Ind % 5) {
	case 0:
		res = opt16;
		break;
	case 1:
		res = hex16;
		break;
	case 2:
		res = circle16;
		break;
	case 3:
		res = apsk12_4;
		break;
	case 4:
		res = qam16;
		break;
	default:
		res = opt16;
	}

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
	Normalize(opt16);
	Normalize(hex16);
	Normalize(circle16);
	Normalize(apsk12_4);

	int population = stoi(argv[1]);
	float CR = stof(argv[2]);
	float F = stof(argv[3]);
	string mode = argv[4];

	gen.seed(CAPACITY_SEED);
	for (unsigned i = 0; i < 2 * NUM_OF_ITERATIONS; ++i) {
		noiseTable[i] = noiseD(gen);
	}
	for (unsigned i = 0; i < NUM_OF_ITERATIONS; ++i) {
		uniformTable[i] = distr(gen);
	}
	unsigned seed = time(nullptr);
	gen.seed(seed);

	cout << "QAM16: " << CapacityApprox(qam16, TARGET_SIGMA, NUM_OF_ITERATIONS)
		<< ", HEX16: " << CapacityApprox(hex16, TARGET_SIGMA, NUM_OF_ITERATIONS)
		<< ","<<endl<<" Circle16: " << CapacityApprox(circle16, TARGET_SIGMA, NUM_OF_ITERATIONS)
		<< ", APSK12_4: " << CapacityApprox(apsk12_4, TARGET_SIGMA, NUM_OF_ITERATIONS)
		<< ", Paper Opt: " << CapacityApprox(opt16, TARGET_SIGMA, NUM_OF_ITERATIONS) << endl;
	auto&& BestConstellation = DEConstellationSearch(CONSTELLATION_SIZE, mode == "best",
		population, CR, F, 10000);

	double bestCapacity = CapacityApprox(BestConstellation, TARGET_SIGMA, NUM_OF_ITERATIONS);
	cout << bestCapacity;

	ofstream out("res_" + to_string(seed) + "_" + to_string(population) + "_" + to_string(CR) + "_" + to_string(F) + "_" + mode + "_" + to_string(bestCapacity) + ".txt");
	for (auto&& elem : BestConstellation)
		out << elem.real() << ' ' << elem.imag() << ' ';//endl;

	return 0;
}