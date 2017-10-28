#include <iostream>
#include <vector>
#include <sstream>

#define  _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <algorithm>
#include <functional>
#include <numeric>
//#include <ctime>
//#include <csignal>


using namespace std;

#define maxiter 500

const int maxdeg = 20;

int rootsMy(double *a, int n, double *wr, double *wi)
{
	double sq, b2, c, disc;
	int i, m, numroots;

	m = n;
	numroots = 0;
	while (m > 1) {
		b2 = -0.5*a[m - 2];
		c = a[m - 1];
		disc = b2*b2 - c;
		if (disc < 0.0) {                   // complex roots
			sq = sqrt(-disc);
			wr[m - 2] = b2;
			wi[m - 2] = sq;
			wr[m - 1] = b2;
			wi[m - 1] = -sq;
			numroots += 2;
		}
		else {                              // real roots
			sq = sqrt(disc);
			wr[m - 2] = fabs(b2) + sq;
			if (b2 < 0.0) wr[m - 2] = -wr[m - 2];
			if (wr[m - 2] == 0)
				wr[m - 1] = 0;
			else {
				wr[m - 1] = c / wr[m - 2];
				numroots += 2;
			}
			wi[m - 2] = 0.0;
			wi[m - 1] = 0.0;
		}
		m -= 2;
	}
	if (m == 1) {
		wr[0] = -a[0];
		wi[0] = 0.0;
		numroots++;
	}
	return numroots;
}

void deflate(double *a, int n, double *b, double *quad, double *err)
{
	double r, s;
	int i;

	r = quad[1];
	s = quad[0];

	b[1] = a[1] - r;

	for (i = 2; i <= n; i++) {
		b[i] = a[i] - r * b[i - 1] - s * b[i - 2];
	}
	*err = fabs(b[n]) + fabs(b[n - 1]);
}

void find_quad(double *a, int n, double *b, double *quad, double *err, int *iter)
{
	double *c, dn, dr, ds, drn, dsn, eps, r, s;
	int i;

	c = new double[n + 1];
	c[0] = 1.0;
	r = quad[1];
	s = quad[0];
	eps = 1e-6;
	*iter = 1;

	do {
		if (*iter > maxiter) break;
		if (((*iter) % 200) == 0) {
			eps *= 10.0;
		}
		b[1] = a[1] - r;
		c[1] = b[1] - r;

		for (i = 2; i <= n; i++) {
			b[i] = a[i] - r * b[i - 1] - s * b[i - 2];
			c[i] = b[i] - r * c[i - 1] - s * c[i - 2];
		}
		dn = c[n - 1] * c[n - 3] - c[n - 2] * c[n - 2];
		drn = b[n] * c[n - 3] - b[n - 1] * c[n - 2];
		dsn = b[n - 1] * c[n - 1] - b[n] * c[n - 2];

		if (fabs(dn) < 1e-6) {
			if (dn < 0.0) dn = -1e-8;
			else dn = 1e-8;
		}
		dr = drn / dn;
		ds = dsn / dn;
		r += dr;
		s += ds;
		(*iter)++;
	} while ((fabs(dr) + fabs(ds)) > eps);
	quad[0] = s;
	quad[1] = r;
	*err = fabs(ds) + fabs(dr);
	delete[] c;
}

void diff_poly(double *a, int n, double *b)
{
	double coef;
	int i;

	coef = (double)n;
	b[0] = 1.0;
	for (i = 1; i<n; i++) {
		b[i] = a[i] * ((double)(n - i)) / coef;
	}
}

void recurse(double *a, int n, double *b, int m, double *quad,
	double *err, int *iter)
{
	double *c, *x, rs[2], tst, e1, e2;

	if (fabs(b[m]) < 1e-7) m--;    // this bypasses roots at zero
	if (m == 2) {
		quad[0] = b[2];
		quad[1] = b[1];
		*err = 0;
		*iter = 0;
		return;
	}
	c = new double[m + 1];
	x = new double[n + 1];
	c[0] = x[0] = 1.0;
	rs[0] = quad[0];
	rs[1] = quad[1];
	*iter = 0;
	find_quad(b, m, c, rs, err, iter);
	tst = fabs(rs[0] - quad[0]) + fabs(rs[1] - quad[1]);
	if (*err < 1e-7) {
		quad[0] = rs[0];
		quad[1] = rs[1];
	}
	// tst will be 'large' if we converge to wrong root
	if (((*iter > 5) && (tst < 1e-4)) || ((*iter > 20) && (tst < 1e-1))) {
		diff_poly(b, m, c);
		recurse(a, n, c, m - 1, rs, err, iter);
		quad[0] = rs[0];
		quad[1] = rs[1];
	}
	delete[] x;
	delete[] c;
}


#include <random>

uniform_real_distribution<double> rnd(0, 1);
mt19937_64 gen = mt19937_64(5489u);//time(nullptr));

void get_quads(const double *a, int n, double *quad, double *x)
{
	double *b, *z, err, tmp;
	double xr, xs;
	int iter, i, m;


	if (n == 2) {
		x[0] = a[1];
		x[1] = a[2];
		return;
	}
	else if (n == 1) {
		x[0] = a[1];
		return;
	}
	m = n;
	b = new double[n + 1];
	z = new double[n + 1];
	b[0] = 1.0;
	for (i = 0; i <= n; i++) {
		z[i] = a[i];
		x[i] = 0.0;
	}
	do {
		if (n > m) {
			quad[0] = 3.14159e-1;
			quad[1] = 2.78127e-1;
		}
		do {                    // This loop tries to assure convergence
			for (i = 0; i<5; i++) {
				find_quad(z, m, b, quad, &err, &iter);
				if ((err > 1e-7) || (iter > maxiter)) {
					diff_poly(z, m, b);
					recurse(z, m, b, m - 1, quad, &err, &iter);
				}
				deflate(z, m, b, quad, &err);
				if (err < 0.001) break;

				quad[0] = rnd(gen);
				quad[1] = rnd(gen);
			}
			if (err > 0.01) {
				quad[0] = rnd(gen);
				quad[1] = rnd(gen);
			}
		} while (err > 0.01);
		x[m - 2] = quad[1];
		x[m - 1] = quad[0];
		m -= 2;
		for (i = 0; i <= m; i++) {
			z[i] = b[i];
		}
	} while (m > 2);
	if (m == 2) {
		x[0] = b[1];
		x[1] = b[2];
	}
	else x[0] = b[1];
	delete[] z;
	delete[] b;
}

void get_roots(int n, const double *a, vector<complex<double>>& vec)
{
	double x[maxdeg], quad[2], err, t, wr[maxdeg], wi[maxdeg];

	quad[0] = 2.71828e-1;//rnd(genR);//2.71828e-1;
	quad[1] = 3.14159e-1;//rnd(genR);//3.14159e-1;
	get_quads(a, n, quad, x);
	rootsMy(x, n, wr, wi);
	for (int i = 0; i < n; ++i)
	{
		vec[i] = { wr[i], wi[i] };
	}
}

vector<double> a, b;
vector<complex<double>> aRoots;
const unsigned K = 512;
vector<double> wk(K);
vector<complex<double> > Hk(K);
uniform_real_distribution<double> mutation;

const int MODULE = 2047;
const double DIV = 64.;
const double LIM = MODULE / DIV;

template<class T>
complex<double> H(double w, const vector<T>& a, const vector<T>& b, bool chk = false)
{
	complex<double> rUp(0., 0.);
	for (unsigned i = 0; i < b.size(); ++i)
		rUp += complex<double>(b[i] * cos(-w * i), b[i] * sin(-w * i));

	complex<double> rDown(0., 0.);
	for (unsigned i = 0; i < a.size(); ++i)
		rDown += complex<double>(a[i] * cos(-w * i), a[i] * sin(-w * i));

	if (chk && abs(rDown) <= 1e-6)
		return {10., 10.};

	auto res = rUp / rDown;

	return res;
}

void Normalize(double& param, double a)
{
	param /= a;
	if (param >= LIM)
		param = LIM;
	else if (param <= -LIM)
		param = LIM;
	else
		param = int(param * DIV) / DIV;
}

//normal_distribution<double> roots(0., 2);
uniform_int_distribution<int> paramD(-MODULE, MODULE);
uniform_real_distribution<double> roots(-0.8, 0.8);

uniform_real_distribution<double> addComplex(0, 1);

complex<double> poly[100];




pair<vector<double>, vector<double> > GenerateRandomParams(int ind)
{
	pair<vector<double>, vector<double> > tempFilterPair(vector<double>(b.size()), vector<double>(a.size()));

	for (unsigned i = 0; i < b.size(); ++i)
		tempFilterPair.first[i] = paramD(gen) / DIV;
	/*tempFilterPair.second[0] = 1;
	for (unsigned i = 1; i < a.size(); ++i)
		tempFilterPair.second[i] = paramD(gen) / DIV;
	return tempFilterPair;*/

	poly[0] = { 1., 0 };
	unsigned deg = 0;
	vector<complex<double >> rs;
	while(deg < a.size() - 1)
	{
		if(addComplex(gen) > 0.15 || a.size() - deg < 3)
		{
			double root = roots(gen);
			if (abs(root) >= 1.)
				root /= (abs(root) + 0.1);
			rs.push_back(root);
			deg++;
			poly[deg] = 1;
			for(int i = deg - 1; i > 0; --i)
				poly[i] = poly[i - 1] - poly[i] * root;
			poly[0] *= -root;
		}
		else
		{
			complex<double> root1 = { roots(gen), roots(gen) };
			if (abs(root1) >= 1.)
				root1 /= (abs(root1) + 0.1);
			complex<double> root2 = { root1.real(), -root1.imag() };
			rs.push_back(root1);
			rs.push_back(root2);
			deg++;
			poly[deg] = 1;
			for (int i = deg - 1; i > 0; --i)
				poly[i] = poly[i - 1] - poly[i] * root1;
			poly[0] *= -root1;
			deg++;
			poly[deg] = 1;
			for (int i = deg - 1; i > 0; --i)
				poly[i] = poly[i - 1] - poly[i] * root2;
			poly[0] *= -root2;
		}
	}
	vector<double> tmp(100);
	for (unsigned i = 0; i < b.size(); ++i)
	{
		tempFilterPair.second[i] = poly[b.size() - 1 - i].real();
		//tmp[i] = tempFilterPair.second[i];
		Normalize(tempFilterPair.second[i], 1.);		
	}
	/*vector<complex<double>> roots(tempFilterPair.second.size());
	get_roots(tempFilterPair.second.size() - 1, tmp.data(), roots);
	get_roots(tempFilterPair.second.size() - 1, tempFilterPair.second.data(), roots);*/

	double maxH = 0;
	for (unsigned i = 0; i < K; ++i)
		maxH = max(maxH, abs(H(wk[i], tempFilterPair.second, tempFilterPair.first)));

	for (unsigned i = 0; i < b.size(); ++i)
	{
		tempFilterPair.first[i] /= maxH;
	}
	return tempFilterPair;
}

double MSE(const vector<double>& bq, const vector<double>& aq)
{
	double BAD = 0;
	vector<complex<double>> roots(aq.size());
	get_roots(aq.size() - 1, aq.data(), roots);
	for(auto&& elem : roots)
	{
		if (abs(elem) >= 1.)
			BAD += 10;
	}
	
	double MSE = 0.;
	for (unsigned i = 0; i < K; ++i)
	{
		auto&& Hi = H(wk[i], aq, bq, true);
		double absH = abs(Hi);
		double W = (wk[i] <= 0.3 * M_PI) ? 2. : 1.;
		
		MSE += pow(abs((Hk[i] - Hi) * W), 2.);
	}
	return MSE / K + BAD;
}


bool operator <(const pair<vector<double>, vector<double> >& a, const pair<vector<double>, vector<double> >& b)
{
	return a.first[0] < b.first[0];
}

vector<double> tmpCoeff(100);

pair<vector<double>, vector<double> > DEFilterSearch(
	bool UseBest,
	unsigned PopulationSize, double CrossoverProbability, double AmplificationFactor,
	unsigned NumOfGenerations)
{
	vector<pair<double, pair<vector<double>, vector<double> > > > population(PopulationSize);

	for (unsigned i = 0; i < PopulationSize; ++i)
	{
		population[i].second = GenerateRandomParams(i);
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
			}
			tempFilterPair.second[0] = 1;
			for (unsigned j = 1; j < a.size(); ++j)
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
				Normalize(tempFilterPair.second[j], tempFilterPair.second[0]);
			}

			double tempMetric = MSE(tempFilterPair.first, tempFilterPair.second);

			if (tempMetric < population[k].first)
			{
				population[k] = { tempMetric, tempFilterPair };
				flag = true;
			}
		}
		sort(population.begin(), population.end());

		/*if (isnan(population.begin()->first))
			cout << "";*/

		//if (i % 1 == 0)
		//	cout << "Generation " << i << ": best metric " << population.begin()->first << endl;
		if (flag == false)
		{
			//cout << "Generation " << i << ": convergence\n";
			break;
		}
	}
	return population.begin()->second;
}


void printReducedPoint(double point)
{
	int discretePoint = int(point * DIV);
	if (discretePoint > MODULE)
		discretePoint = MODULE;
	else if (discretePoint < -MODULE)
		discretePoint = -MODULE;
	cout << discretePoint;
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

	/*aRoots.resize(a.size());
	get_roots(a.size() - 1, a.data(), aRoots); 
	for (auto&& elem : aRoots)
	{
		if(abs(elem) >= 1)
			cout << "FUCK" << endl;
	}*/

	auto && bestFilters = DEFilterSearch(false, /*(a.size() + b.size()) * 8*/50, 0.71, 0.7, 50);

	//for(auto&& w : wk)
	//{
	//	auto res = H(w, bestFilters.second, bestFilters.first);
	//	cout << (abs(res) /*>= 1*/) << endl;
	//}
	cout << endl << MSE(bestFilters.first, bestFilters.second) << endl;
	
	printReducedPoint(bestFilters.first[0]);
	for (unsigned i = 1; i < b.size(); ++i){
		cout << ", ";
		printReducedPoint(bestFilters.first[i]);
	}
	cout << endl;
	printReducedPoint(bestFilters.second[0]);
	for (unsigned i = 1; i < b.size(); ++i){
		cout << ", ";
		printReducedPoint(bestFilters.second[i]);
	}
	cout << endl;
}