/** Do a spherical cap integral, using the ang2d framework.
 *
 * Nikhil Padmanabhan, Yale.
 */

#include "ang2d.h"

#include <utility>
#include <iostream>
#include <fstream>
#include <chrono>
#include "boost/format.hpp"
#include "boost/program_options.hpp"

using namespace std;
using boost::format;
using namespace chrono;

namespace po = boost::program_options;


typedef pair<double, double> dpair;



/* Simple mask, centered on the N-pole at dec > dec0 */
class Cap {
public :
	double dec0;
	Cap(double dec0_) : dec0(dec0_) {};

	double operator()(double ra, double dec);

};

double Cap::operator()(double ra, double dec) {
	if (dec > dec0) {return 1.0;};
	return 0.0;
}


int main(int argc, char **argv) {
	int nrand, nsim;
	double decmin, thetamax;

	string fn;
	bool savefile=false;

	// Get the input parameters -- pull them into their own scope
	{
		try {
			po::options_description desc("Allowed options");
			desc.add_options()
	    				("help", "produce help message")
	    				("nrand", po::value<int>(&nrand)->default_value(10000), "Number of pseudo-random points")
	    				("decmin", po::value<double>(&decmin)->default_value(0.0), "Minimum dec value for cap")
	    				("thetamax", po::value<double>(&thetamax)->default_value(5.0), "Maximum theta value for RR")
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
	cout << format("Running with a cap centered on the N-pole down to dec = %1%...\n")%decmin;
	cout << format("Running with %1% pseudo-random numbers... \n")%nrand;
	cout << format("and %1% simulations\n")%nsim;
	if (savefile) {
		cout << format("Simulations will be saved in %1% \n")%fn;
	}

	// Define various bounds
	dpair RABounds(0, 360);
	dpair decBounds(decmin,90);
	dpair thetaBin(0, thetamax);
	Cap mask1(decmin);

	steady_clock::time_point t1 = steady_clock::now();
	// Actual call to code needs to go here.
	vector<double> val = Ang2D::rreval(RABounds, decBounds, thetaBin, nrand, nsim, mask1, mask1);
	steady_clock::time_point t2 = steady_clock::now();
	double mean=0.0, stddev=0.0;

	for (auto v1 : val) mean+=v1; mean /= val.size();
	for (auto v1 : val) stddev += (v1-mean)*(v1-mean); stddev = sqrt(stddev/val.size());

	cout << format("The mean is %13.10f +/- %13.10f with a scatter of %13.10f, a fractional error of %9.6f percent\n")
			% mean % (stddev/sqrt(nsim)) % stddev % (stddev/mean * 100);

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << format("Total evaluation time = %1% seconds \n")%(time_span.count());


	// Binary file format
	// int nsim, int nrand, double decmin, double thetamax, double values[nrand]
	if (savefile) {
		ofstream ofs(fn.c_str(), ios::binary);
		if (!ofs) {
			cout << "ERROR! Unable to save file\n";
			return 1;
		}
		ofs.write((char*)&nsim, sizeof(int));
		ofs.write((char*)&nrand, sizeof(int));
		ofs.write((char*)&decmin, sizeof(double));
		ofs.write((char*)&thetamax, sizeof(double));
		ofs.write((char*)&val[0], sizeof(double)*nsim);
		ofs.close();
	}

}

