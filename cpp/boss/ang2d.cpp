#include "ang2d.h"
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

	double mean, stddev, nsim;
	int nrand1;
	dvector val1;

	for (auto v1 : zip(nrands, vals)) {
		// This is an extra copy, but it makes the code cleaner
		boost::tie(nrand1, val1) = v1;
		nsim = val1.size();
		cout << format("Using %7i terms... \n") % nrand1;
		tie(mean, stddev) = stats(val1);
		cout << format("The mean is %13.10e +/- %13.10e with a scatter of %13.10e, a fractional error of %9.6f percent\n")
				% mean % (stddev/sqrt(nsim)) % stddev % (stddev/mean * 100);
	}


}

void Ang2D::OutputData::save(ofstream& out) {
	int nrand1, nsim;
	dvector val1;

	int nelems = vals.size();
	out.write((char*)&nelems, sizeof(int));

	for (auto v1 : zip(nrands, vals)) {
		boost::tie(nrand1, val1) = v1;
		nsim = val1.size();
		out.write((char*)&nrand1, sizeof(int));
		out.write((char*)&nsim, sizeof(int));
		out.write((char*)&val1[0], sizeof(double)*nsim);
	}
}
