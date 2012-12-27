/** Do the angular integral over the DEEP2 mask.
 *
 * Nikhil Padmanabhan, Yale
 */

#include "ang2d.h"
#include "BossMask.h"

#include <string>
#include <numeric>
#include <iostream>
#include <fstream>
#include <chrono>
#include "boost/format.hpp"
#include "boost/program_options.hpp"


using namespace std;
using boost::format;
using namespace chrono;
using namespace Boss;

namespace po = boost::program_options;


typedef pair<double, double> dpair;

int main(int argc, char **argv) {
	int nrand, nsim;
	double thetamin, thetamax;
	double cthresh; // Completeness threshold
	bool flatmask = false;
	string maskfn;

	string fn;
	bool savefile=false;
	bool use_prng=false;

	// Get the input parameters -- pull them into their own scope
	{
		try {
			po::options_description desc("Allowed options");
			desc.add_options()
	    				("help", "produce help message")
	    				("maskfn", po::value<string>(&maskfn), "mask file")
	    				("nrand", po::value<int>(&nrand)->default_value(10000), "Number of pseudo-random points")
	    				("thetamin", po::value<double>(&thetamin)->default_value(0.0), "Minimum theta value for RR")
	    				("thetamax", po::value<double>(&thetamax)->default_value(0.5), "Maximum theta value for RR")
	    				("cthresh", po::value<double>(&cthresh)->default_value(0.5), "Threshold completeness for mask [0.5]")
	    				("flatmask","Use completeness values for the mask in the RR calculation [false]")
	    				("nsim", po::value<int>(&nsim)->default_value(1), "Number of simulations")
	    				("save", po::value<string>(&fn), "save simulation results to file [optional]")
	    				("prng", "Use pseudo-random numbers instead of a low discrepancy sequence")
	    				;

			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);
			po::notify(vm);

			if (vm.count("help")) {
				cout << desc << "\n";
				return 1;
			}
			if (vm.count("maskfn")==0) {
				cout << "ERROR : A mask file is required\n";
				return 1;
			}
			if (vm.count("save")) {
				savefile=true;
			}
			if (vm.count("prng")) {
				use_prng=true;
			}
			if (vm.count("flatmask")) {
				flatmask=true;
			}


		}
		catch (exception &e) {
			cout << e.what() << "\n";
			return 1;
		}

	}

	// Print some informational messages
	cout << format("Running with a bin from %1% to %2%...\n")%thetamin%thetamax;
	cout << format("Running with %1% pseudo-random numbers... \n")%nrand;
	cout << format("and %1% simulations\n")%nsim;
	if (savefile) {
		cout << format("Simulations will be saved in %1% \n")%fn;
	}
	if (use_prng) {
		cout << "Using pseudo-random numbers instead of a low discrepancy sequence\n";
	}

	// Now define the mask
	BossMask mask1(maskfn,cthresh,flatmask);

	// Define various bounds
	dpair RABounds = mask1.RABounds;
	dpair decBounds = mask1.DecBounds;
	dpair thetaBin(thetamin, thetamax);

	// Test the cap area -- do just one simulation
	{
		vector<double> val = Ang2D::area(RABounds, decBounds, 1000000, 1, mask1, use_prng);
		cout << format("The mask integrates to %13.10f\n")%val[0];
	}


	steady_clock::time_point t1 = steady_clock::now();
	// Actual call to code needs to go here.
	vector<double> val = Ang2D::rreval(RABounds, decBounds, thetaBin, nrand, nsim, mask1, mask1, use_prng);
	steady_clock::time_point t2 = steady_clock::now();
	double mean=0.0, stddev=0.0;

	for (auto v1 : val) mean+=v1; mean /= val.size();
	for (auto v1 : val) stddev += (v1-mean)*(v1-mean); stddev = sqrt(stddev/val.size());

	cout << format("The mean is %13.10e +/- %13.10e with a scatter of %13.10e, a fractional error of %9.6f percent\n")
			% mean % (stddev/sqrt(nsim)) % stddev % (stddev/mean * 100);

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << format("Total evaluation time = %1% seconds \n")%(time_span.count());


	// Binary file format
	// int nsim, int nrand, double thetamin, double thetamax, double values[nrand]
	if (savefile) {
		ofstream ofs(fn.c_str(), ios::binary);
		if (!ofs) {
			cout << "ERROR! Unable to save file\n";
			return 1;
		}
		ofs.write((char*)&nsim, sizeof(int));
		ofs.write((char*)&nrand, sizeof(int));
		ofs.write((char*)&thetamin, sizeof(double));
		ofs.write((char*)&thetamax, sizeof(double));
		ofs.write((char*)&val[0], sizeof(double)*nsim);
		ofs.close();
	}



}
