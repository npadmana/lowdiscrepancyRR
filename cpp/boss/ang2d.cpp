#include "ang2d.h"
#include "src/misc/npvecmath.h"
#include <iostream>
#include "boost/format.hpp"

using namespace std;
using boost::format;



void Ang2D::vecshift(vector<double> &x, vector<double> &x0) {
	double tmp;
	transform(x0.begin(), x0.end(), x.begin(), x.begin(),
			[&tmp](double a, double b){return modf(a+b,&tmp);});
}

void Ang2D::OutputData::print() {

	double mean, stddev, err;
	int nrand1;
	dvector val1;

	for (const auto& v1 : stats) {
		tie(nrand1, ignore, mean, stddev, err) = v1;
		cout << format("Using %7i terms ... \n") % nrand1;
		cout << format("The mean is %13.10e +/- %13.10e with a scatter of %13.10e, a fractional error of %9.6f percent\n")
				% mean % err % stddev % (stddev/mean * 100);
	}


}

void Ang2D::OutputData::push_back(const int nrand, const dvector& vec) {
	results.push_back(make_tuple(nrand, vec));
}

void Ang2D::OutputData::finalize() {

	// MPI finalization will happen here


	// Compute the statistics
	{
		int nrand1, nsim;
		dvector val1;
		double mean, stddev, err;

		// Compute statistics
		stats.clear();
		for (const auto& v1 : results) {
			tie(nrand1, val1) = v1; // A copy, but this is small
			nsim = val1.size();
			mean = vecmean(val1);
			stddev = vecstddev(val1, mean);
			err = stddev/sqrt(nsim);
			stats.push_back(make_tuple(nrand1, nsim, mean, stddev, err));

		}
	}
}

void Ang2D::OutputData::save(ofstream& out) {
	int nrand1, nsim;
	dvector val1;

	int nelems = results.size();
	out.write((char*)&nelems, sizeof(int));

	for (const auto& v1 : results) {
		tie(nrand1, val1) = v1; // A copy, but this is small
		nsim = val1.size();
		out.write((char*)&nrand1, sizeof(int));
		out.write((char*)&nsim, sizeof(int));
		out.write((char*)&val1[0], sizeof(double)*nsim);
	}
}

double Ang2D::OutputData::mean() {
	return get<2>(stats.back());
}

double Ang2D::OutputData::error() {
	return get<3>(stats.back());
}
