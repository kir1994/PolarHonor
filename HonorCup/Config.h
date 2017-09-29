#pragma once
#include <vector>
#include <complex>
#include <random>

using namespace std;

const double TARGET_CAPACITY = 3.6;
const unsigned NUM_OF_ITERATIONS = 10000;
const double EPS = 1e-3;
const unsigned CONSTELLATION_SIZE = 16;
const double SMALL_MODULE = 0.1;

const double TARGET_SNR = 12;
const double TARGET_SIGMA = sqrt(1. / (2 * pow(10., TARGET_SNR / 10)));

const double SCALING_FACTOR = 3;
const double CAPACITY_SEED = 123456487;


extern vector<complex<double> > circle16;
extern vector<complex<double> > hex16;
extern vector<complex<double> > opt16;
extern vector<complex<double> > apsk12_4;
extern vector<complex<double> > qam16;

extern mt19937_64 gen;
extern vector<double> noiseTable;
extern vector<unsigned> uniformTable;

enum handle_bad_points { DO_NOTHING, MOVE_POINT, MOVE_CONSTELLATION };
void Normalize(vector<complex<double> >& Constellation, handle_bad_points mode = DO_NOTHING, double eps = 0.0001);

double avgPower(const vector<complex<double> >& Constellation);
double CapacityApprox(const vector<complex<double> >& Constellation, double Sigma, unsigned NumOfIterations);
double FindSNR(const vector<complex<double> >& Constellation);

vector<complex<double> > DEConstellationSearch(
	unsigned ConstellationSize, bool UseBest,
	unsigned PopulationSize, double CrossoverProbability, double AmplificationFactor,
	unsigned NumOfGenerations);
vector<complex<double> > dlibMinSearch(unsigned ConstellationSize, unsigned ind, unsigned mode);


void PrintConstellation(const vector<complex<double> >& Constellation);

void handler(int sig);