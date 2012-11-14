/** RR calculation using quasi-Monte Carlo methods.
 *
 * Nikhil Padmanabhan, Yale
 * November 13, 2012
 */

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include "boost/format.hpp"
#include "boost/program_options.hpp"


#include "src/npgsl/npRandom.h"

using namespace std;
using boost::format;
using namespace chrono;

namespace po = boost::program_options;


const int dim = 6;
const double pi = 3.1415926535897932385;



vector<double> rreval(double rmin, double rmax, int nrand, int nsim) {

	// Set up the output vector
	vector<double> outvec(nsim);

	// Set up the random vector
	vector<double> x0(dim);
	npRandom rng(1234);

	// Temporary values
	double x1, y1, z1, r, cth, sth, phi, rlo3, dr3, out;
	rlo3 = pow(rmin,3.0);
	dr3 = pow(rmax, 3.0)-rlo3;

	for(int jj=0; jj<nsim; ++jj) {

		for_each(x0.begin(), x0.end(), [&rng](double &x){x=rng();});

		// Set up the quasi RNG
		npQuasiRandom qrng(dim);

		// Initialize the integrator
		out=0.0;

		// Loop over points
		for (int ii=0; ii<nrand; ++ii) {

			// Generate a random pseudo-random number and shift it mod 1
			vector<double> x = qrng();
			transform(x0.begin(), x0.end(), x.begin(), x.begin(),
					[&r](double a, double b){return modf(a+b,&r);});

			// Work out displacement
			phi = x[5]*2*pi;
			cth = x[4]*2 - 1;
			sth = sqrt(1-cth*cth);
			r = pow(x[3]*dr3 + rlo3, 1./3.);
			x1 = x[0] + r*sth*cos(phi);
			y1 = x[1] + r*sth*sin(phi);
			z1 = x[2] + r*cth;

			// Check to see if in box
			if ((x1 > 0.0) && (x1 < 1.0) && (y1 > 0.0) && (y1 < 1.0) && (z1 > 0.0) && (z1 < 1.0)) out++;

		}

		// Multiply back in the jacobian
		outvec[jj]  = out/(nrand) * 4 * pi * dr3/3.0;
	}
	return outvec;
}

int main(int argc, char **argv) {
	int nrand, nsim;
	double rmin, rmax;

	string fn;
	bool savefile=false;

	// Get the input parameters -- pull them into their own scope
	{
		try {
			po::options_description desc("Allowed options");
			desc.add_options()
	    				("help", "produce help message")
	    				("nrand", po::value<int>(&nrand)->default_value(10000), "Number of pseudo-random points")
	    				("rmin", po::value<double>(&rmin)->default_value(0.0), "Minimum r value for shell")
	    				("rmax", po::value<double>(&rmax)->default_value(0.5), "Maximum r value for shell")
	    				("nsim", po::value<int>(&nsim)->default_value(1), "Number of simulations")
	    				("save", po::value<string>(&fn), "save simulation results to file [optional]")
	    				;

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);
			po::notify(vm);

			if (vm.count("help")) {
				cout << desc << "\n";
				return 1;
			}
			if (vm.count("save")) {
				savefile=true;
			}

		}
		catch (exception &e) {
			cout << e.what() << "\n";
			return 1;
		}

	}

	// Print some informational messages
	cout << format("Running with bins from %1% to %2%...\n")%rmin%rmax;
	cout << format("Running with %1% pseudo-random numbers... \n")%nrand;
	cout << format("and %1% simulations\n")%nsim;
	if (savefile) {
		cout << format("Simulations will be saved in %1% \n")%fn;
	}

	steady_clock::time_point t1 = steady_clock::now();
	vector<double> val = rreval(rmin, rmax, nrand, nsim);
	steady_clock::time_point t2 = steady_clock::now();
	double mean=0.0, stddev=0.0;

	for (auto v1 : val) mean+=v1; mean /= val.size();
	for (auto v1 : val) stddev += (v1-mean)*(v1-mean); stddev = sqrt(stddev/val.size());

	cout << format("The mean is %13.10f +/- %13.10f with a scatter of %13.10f, a fractional error of %9.6f percent\n")
			% mean % (stddev/sqrt(nsim)) % stddev % (stddev/mean * 100);

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << format("Total evaluation time = %1% seconds \n")%(time_span.count());


	// Binary file format
	// int nsim, int nrand, double rmin, double rmax, double values[nrand]
	if (savefile) {
		ofstream ofs(fn.c_str(), ios::binary);
		if (!ofs) {
			cout << "ERROR! Unable to save file\n";
			return 1;
		}
		ofs.write((char*)&nsim, sizeof(int));
		ofs.write((char*)&nrand, sizeof(int));
		ofs.write((char*)&rmin, sizeof(double));
		ofs.write((char*)&rmax, sizeof(double));
		ofs.write((char*)&val[0], sizeof(double)*nsim);
		ofs.close();
	}

}
