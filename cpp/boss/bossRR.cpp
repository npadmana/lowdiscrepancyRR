/** Do the angular integral over a mangle mask
 *
 * Nikhil Padmanabhan, Yale
 */

#include "ang2d.h"
#include "BossMask.h"

#include <string>
#include <list>
#include <numeric>
#include <iostream>
#include <fstream>
#include <chrono>
#include "boost/format.hpp"
#include "boost/program_options.hpp"

#include "mpi.h"

using namespace std;
using boost::format;
using namespace chrono;
using namespace Boss;

namespace po = boost::program_options;


typedef pair<double, double> dpair;

int main(int argc, char **argv) {
	// MPI initialization
	MPI_Init(&argc, &argv);
	int rank, nproc;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	if (argc < 2) throw runtime_error("Missing argument");

	Ang2D::InputParams p0(argv[1]);

	// Check to see that thetabins is defined
	int nbins = p0.thetabins.size()-1;
	if (nbins < 1) throw invalid_argument("thetabins needs to have length > 1");

	// Print some informational messages
	if (p0.verbose && (rank==0)) {
		cout << format("Running with bins from %1% to %2%...\n")%p0.thetabins[0]%p0.thetabins[nbins];
		cout << format("Running with %1% pseudo-random numbers... \n")%p0.nrand;
		cout << format("and %1% simulations\n")%p0.nsim;
		if (p0._dump) {
			cout << format("Simulations will be saved in %1% \n")%p0.dumpfn;
		}
		if (p0.use_prng) {
			cout << "Using pseudo-random numbers instead of a low discrepancy sequence\n";
		}
	}

	// Now define the mask
	BossMask mask1(p0.mask1fn);
	BossMask mask2(p0.mask2fn);

	// Define various bounds
	int nra, ndec;
	tie(ndec, nra) = Ang2D::partition(nproc);
	if (p0.verbose && (rank==0)) cout << "Partitioning into " << nra << " x " << ndec << endl;
	Ang2D::setBounds(nra,ndec,rank,mask1,p0);


	// Get ready to execute the loop over thetabins here
	list<Ang2D::OutputData> outlist;
	steady_clock::time_point t1 = steady_clock::now();

	for (int ibin=0; ibin<nbins; ++ibin) {
	// Actual call to code needs to go here.
		Ang2D::OutputData out = Ang2D::rreval(mask1, mask2, p0, p0.thetabins[ibin], p0.thetabins[ibin+1]);
		out.finalize();
		outlist.push_back(out);
	}
	steady_clock::time_point t2 = steady_clock::now();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);


	if (rank==0) {
		cout << format("Total evaluation time = %1% seconds \n")%(time_span.count());



		// Binary file format
		if (p0._dump) {
			ofstream ofs(p0.dumpfn, ios::binary);
			if (!ofs) {
				cout << "ERROR! Unable to save dump file\n";
				return 1;
			}

			// Save the number of bins
			ofs.write((char*)&nbins, sizeof(int));
			// Now save the data -- thetamin, thetamax, and then the dump
			auto ii = outlist.begin();
			for (int ibin=0;ibin < nbins; ii++, ++ibin) {
				ofs.write((char*)&p0.thetabins[ibin], sizeof(double));
				ofs.write((char*)&p0.thetabins[ibin+1], sizeof(double));
				ii->save(ofs);
			}
			ofs.close();
		}

		// Save file format
		if (p0._save) {
			ofstream ofs(p0.savefn);
			if (!ofs) {
				cout << "ERROR! Unable to save file\n";
				return 1;
			}

			auto ii = outlist.begin();
			ofs << format("# %15s %15s %15s %9s\n")%"theta_min"%"theta_max"%"RR"%"error";
			for (int ibin=0;ibin < nbins; ii++, ++ibin) {
				ofs << format("  %15.10e %15.10e %15.10e %9.6f\n")
									% p0.thetabins[ibin] % p0.thetabins[ibin+1] % ii->mean() % ii->error(true);
			}
			ofs.close();
		}
	}


	// Nothing goes below this
	MPI_Finalize();
}
