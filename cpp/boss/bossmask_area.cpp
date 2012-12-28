/** Compute the area of a BOSS mask. This is a useful utility for appropriately normalizing masks
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
using namespace Padmanabhan;
namespace po = boost::program_options;



int main(int argc, char **argv) {
	if (argc < 2) throw runtime_error("Missing argument");

	Ang2D::InputParams p0(argv[1]);


	// Print some informational messages
	cout << format("Running with %1% pseudo-random numbers... \n")%p0.nrand;
	cout << format("and %1% simulations\n")%p0.nsim;
	if (p0._dump) {
		cout << format("Simulations will be saved in %1% \n")%p0.dumpfn;
	}
	if (p0.use_prng) {
		cout << "Using pseudo-random numbers instead of a low discrepancy sequence\n";
	}

	// Now define the mask
	BossMask mask1(p0.mask1fn);

	// Define various bounds
	p0.ramin = mask1.ramin; p0.ramax = mask1.ramax;
	p0.decmin = mask1.decmin; p0.decmax = mask1.decmax;


	steady_clock::time_point t1 = steady_clock::now();
	// Actual call to code needs to go here.
	Ang2D::OutputData out = Ang2D::area(mask1, p0);
	steady_clock::time_point t2 = steady_clock::now();
	out.finalize();

	out.print();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << format("Total evaluation time = %1% seconds \n")%(time_span.count());


	// Binary file format
	if (p0._dump) {
		ofstream ofs(p0.dumpfn, ios::binary);
		if (!ofs) {
			cout << "ERROR! Unable to save file\n";
			return 1;
		}
		out.save(ofs);
		ofs.close();
	}



}
