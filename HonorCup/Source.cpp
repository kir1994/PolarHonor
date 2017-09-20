#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <ctime>
#include <functional>
#include <signal.h>

#include "Config.h"


#include <fstream>

normal_distribution<double> noiseD(0., TARGET_SIGMA);
uniform_int_distribution<int> distr(0., CONSTELLATION_SIZE - 1);

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