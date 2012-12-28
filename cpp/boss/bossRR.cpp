/** Do the angular integral over a mangle mask
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
	Ang2D::InputParams p0;
	double thetamin, thetamax;
	string maskfn;

	string fn;
	bool savefile=false;

	// Get the input parameters -- pull them into their own scope
	{
		try {
			po::options_description desc("Allowed options");
			desc.add_options()
	    				("help", "produce help message")
	    				("maskfn", po::value<string>(&maskfn), "mask file")
	    				("nrand", po::value<int>(&p0.nrand)->default_value(10000), "Number of pseudo-random points")
	    				("thetamin", po::value<double>(&thetamin)->default_value(0.0), "Minimum theta value for RR")
	    				("thetamax", po::value<double>(&thetamax)->default_value(0.5), "Maximum theta value for RR")
	    				("nsim", po::value<int>(&p0.nsim)->default_value(1), "Number of simulations")
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
				p0.use_prng=true;
			}


		}
		catch (exception &e) {
			cout << e.what() << "\n";
			return 1;
		}

	}

	// Print some informational messages
	cout << format("Running with a bin from %1% to %2%...\n")%thetamin%thetamax;
	cout << format("Running with %1% pseudo-random numbers... \n")%p0.nrand;
	cout << format("and %1% simulations\n")%p0.nsim;
	if (savefile) {
		cout << format("Simulations will be saved in %1% \n")%fn;
	}
	if (p0.use_prng) {
		cout << "Using pseudo-random numbers instead of a low discrepancy sequence\n";
	}

	// Now define the mask
	BossMask mask1(maskfn);

	// Define various bounds
	p0.ramin = mask1.ramin; p0.ramax = mask1.ramax;
	p0.decmin = mask1.decmin; p0.decmax = mask1.decmax;

	// Define the save schedule
	{
		int n1=10000;
		while (n1 < p0.nrand) {
			p0.save_schedule.push_back(n1);
			n1 *= 2;
		}
	}


	steady_clock::time_point t1 = steady_clock::now();
	// Actual call to code needs to go here.
	Ang2D::OutputData out = Ang2D::rreval(mask1, mask1, p0, thetamin, thetamax);
	steady_clock::time_point t2 = steady_clock::now();
	out.finalize();

	out.print();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << format("Total evaluation time = %1% seconds \n")%(time_span.count());



	// Binary file format
	if (savefile) {
		ofstream ofs(fn.c_str(), ios::binary);
		if (!ofs) {
			cout << "ERROR! Unable to save file\n";
			return 1;
		}
		ofs.write((char*)&thetamin, sizeof(double));
		ofs.write((char*)&thetamax, sizeof(double));
		out.save(ofs);
		ofs.close();
	}



}
