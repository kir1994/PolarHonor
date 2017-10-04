#include <iostream>
#include <vector>
#include <sstream>

#define  _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <algorithm>
#include <functional>
#include <numeric>
#include <random>
#include <ctime>
//#include <csignal>


using namespace std;

vector<double> a, b;
const unsigned K = 512;
vector<double> wk(K);
vector<complex<double>> Hk(K);
uniform_real_distribution<double> mutation;
mt19937_64 gen = mt19937_64(time(nullptr));

template<class T>
complex<double> H(double w, const vector<T>& a, const vector<T>& b)
{
	complex<double> rUp(0., 0.);
	for (unsigned i = 0; i < b.size(); ++i)
		rUp += complex<double>(b[i] * cos(-w * i), b[i] * sin(-w * i));

	complex<double> rDown(0., 0.);
	for (unsigned i = 0; i < a.size(); ++i)
		rDown += complex<double>(a[i] * cos(-w * i), a[i] * sin(-w * i));

	auto res = rUp / rDown;

	return res;
}

uniform_int_distribution<int> paramD(-255, 255);
pair<vector<double>, vector<double> > GenerateRandomParams()
{
	pair<vector<double>, vector<double> > tempFilterPair(vector<double>(b.size()), vector<double>(a.size()));
	for (unsigned i = 0; i < b.size(); ++i)
		tempFilterPair.first[i] = paramD(gen) / 16.; 
	for (unsigned i = 0; i < a.size(); ++i)
		tempFilterPair.second[i] = paramD(gen) / 16.;
	return tempFilterPair;
}

double MSE(const vector<double>& bq, const vector<double>& aq)
{
	double MSE = 0.;
	for (unsigned i = 0; i < K; ++i)
	{
		auto&& Hi = H(wk[i], aq, bq);
		double absH = abs(Hi);
		double offset = absH <= 1 ? 1. : abs(Hi);

		MSE += pow(abs(Hk[i] - Hi) * offset, 2.);
	}
	return MSE / K;
}

//volatile bool F = true;
//
//
//void handler(int sig)
//{
//	F = false;
//}

bool operator <(const pair<vector<double>, vector<double> >& a, const pair<vector<double>, vector<double> >& b)
{
	return a.first[0] < b.first[0];
}

pair<vector<double>, vector<double> > DEFilterSearch(
	bool UseBest,
	unsigned PopulationSize, double CrossoverProbability, double AmplificationFactor,
	unsigned NumOfGenerations)
{
	vector<pair<double, pair<vector<double>, vector<double> > > > population(PopulationSize);

	for (unsigned i = 0; i < PopulationSize; ++i)
	{
		population[i].second = GenerateRandomParams();
		population[i].first = MSE(population[i].second.first, population[i].second.second);
	}

	sort(population.begin(), population.end());


	vector<unsigned> permutation(PopulationSize);
	iota(permutation.begin(), permutation.end(), 0);

	pair<vector<double>, vector<double> > tempFilterPair(vector<double>(b.size()), vector<double>(a.size()));
	/*double curBest = population.begin()->first;
	unsigned cnt = 0;*/
	bool flag = false;
	for (unsigned i = 0; /*F &&*/ i < NumOfGenerations; ++i)
	{
		flag = false;
		for (unsigned k = 0; k < PopulationSize; ++k)
		{
			auto&& pos = find(permutation.begin(), permutation.end(), k);
			swap(*permutation.rbegin(), *pos);
			shuffle(permutation.begin(), prev(permutation.end()), gen);
			for (unsigned j = 0; j < b.size(); ++j)
			{
				if (mutation(gen) > CrossoverProbability)
					tempFilterPair.first[j] = population[k].second.first[j];
				else
				{
					if (!UseBest)
					{
						tempFilterPair.first[j] = population[permutation[0]].second.first[j] +
							AmplificationFactor * (population[permutation[1]].second.first[j] - population[permutation[2]].second.first[j]);
					}
					else
					{
						tempFilterPair.first[j] = population[0].second.first[j] +
							AmplificationFactor *
							(population[permutation[0]].second.first[j] + population[permutation[1]].second.first[j]
								- population[permutation[2]].second.first[j] - population[permutation[3]].second.first[j]);
					}
				}
				tempFilterPair.first[j] = round(tempFilterPair.first[j] * 16) / 16;		
				if (tempFilterPair.first[j] > 15.9375)
					tempFilterPair.first[j] = 15.9375;
				else if (tempFilterPair.first[j] < -15.9375)
					tempFilterPair.first[j] = 15.9375;
			}
			for (unsigned j = 0; j < a.size(); ++j)
			{
				if (mutation(gen) > CrossoverProbability)
					tempFilterPair.second[j] = population[k].second.second[j];
				else
				{
					if (!UseBest)
					{
						tempFilterPair.second[j] = population[permutation[0]].second.second[j] +
							AmplificationFactor * (population[permutation[1]].second.second[j] - population[permutation[2]].second.second[j]);
					}
					else
					{
						tempFilterPair.second[j] = population[0].second.second[j] +
							AmplificationFactor *
							(population[permutation[0]].second.second[j] + population[permutation[1]].second.second[j]
								- population[permutation[2]].second.second[j] - population[permutation[3]].second.second[j]);
					}
				}
				tempFilterPair.second[j] = round(tempFilterPair.second[j] * 16) / 16;
				if (tempFilterPair.second[j] > 15.9375)
					tempFilterPair.second[j] = 15.9375;
				else if (tempFilterPair.second[j] < -15.9375)
					tempFilterPair.second[j] = 15.9375;
			}

			double tempMetric = MSE(tempFilterPair.first, tempFilterPair.second);

			if (tempMetric < population[k].first)
			{
				population[k] = { tempMetric, tempFilterPair };
				flag = true;
			}
		}
		sort(population.begin(), population.end());

//		if (i % 1 == 0)
//			cout << "Generation " << i << ": best metric " << population.begin()->first << endl;
		if (flag == false)
		{
			//cout << "Generation " << i << ": convergence\n";
			break;
		}
	}
	return population.begin()->second;
}




int main()
{
	//signal(SIGINT, handler);
	string tmp;
	
	getline(cin, tmp);
	stringstream ss(tmp);
	double val;
	while (ss >> val)
	{
		b.push_back(val);
		ss.ignore(1);
	}

	getline(cin, tmp);
	stringstream ss2(tmp);
	while (ss2 >> val)
	{
		a.push_back(val);
		ss2.ignore(1);
	}
	for (unsigned k = 0; k < K; ++k)
	{
		wk[k] = k * M_PI / K;
		Hk[k] = H(wk[k], a, b);
	}


	auto && bestFilters = DEFilterSearch(true, (a.size() + b.size()) * 4, 0.71, 0.5, 50);


	/*for(auto&& w : wk)
	{
		auto res = H(w, bestFilters.second, bestFilters.first);
		cout << (abs(res) >= 1);
	}

	cout << endl << MSE(bestFilters.first, bestFilters.second) << endl;*/
	cout << int(bestFilters.first[0] * 16);
	for (unsigned i = 1; i < b.size(); ++i)
		cout << ',' << int(bestFilters.first[i] * 16);
	cout << endl;
	cout << int(bestFilters.second[0] * 16);
	for (unsigned i = 1; i < a.size(); ++i)
		cout << ',' << int(bestFilters.second[i] * 16) ;
}