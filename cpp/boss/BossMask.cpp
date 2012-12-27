/*
 * BossMask.cpp
 *
 *  Created on: Dec 17, 2012
 *      Author: npadmana
 */

#include <iostream>
#include <vector>
#include <stdexcept>

#include "BossMask.h"
#include "ang2d.h"

// We will use boost::program_options to parse the file
#include "boost/program_options.hpp"

using namespace std;
namespace po = boost::program_options;

Boss::BossMask::BossMask(string configfn) :
		flat(false),
		thresh(0.0),
		area(1.0),
		RABounds(0.0, 360.0),
		DecBounds(-90.0,90.0)
{
	string acceptfn;
	double ramin, ramax, decmin, decmax;
	int flat_;

	// Get the mask configuration parameters
	{
		try {
			po::options_description desc("Allowed options");
			desc.add_options()
	    				("help", "produce help message")
	    				("acceptfn",po::value<string>(&acceptfn), "Acceptance mangle mask")
	    				("ra_min", po::value<double>(&ramin)->default_value(0), "Minimum RA value to consider [0]")
	    				("ra_max", po::value<double>(&ramax)->default_value(360), "Maximum RA value to consider [360.0]")
	    				("dec_min", po::value<double>(&decmin)->default_value(-90.0), "Minimum Dec value to consider [-90]")
	    				("dec_max", po::value<double>(&decmax)->default_value(90), "Maximum RA value to consider [90]")
	    				("thresh", po::value<double>(&thresh)->default_value(0.0), "Completeness threshold")
	    				("norm", po::value<double>(&area)->default_value(1.0), "Normalization constant [1.0]")
	    				("flatmask", po::value<int>(&flat_)->default_value(0), "Is the mask flat? [0]")
	    				;

			po::variables_map vm;
			ifstream ifs(configfn.c_str());
			if (!ifs) {
				cout << "Error opening configuration file " << configfn << endl;
				throw runtime_error("Error opening configuration file");
			}
			po::store(po::parse_config_file(ifs, desc), vm);
			po::notify(vm);
			ifs.close();

			if (vm.count("help")) {
				cout << desc << endl;
				throw runtime_error("Help requested, now aborting.");
			}
			if (vm.count("acceptfn")==0) {
				cout << "ERROR : An acceptance mask is required\n" << endl;
				throw invalid_argument("Missing acceptance mask");
			}

		}
		catch (exception &e) {
			cout << e.what() << endl;
			throw;
		}

	}


	acceptance.reset(new MaskClass(acceptfn));

	// The rejection mask code would go here.

	// Reset all values
	if (flat_ != 0) flat=true;
	RABounds = dpair(ramin, ramax);
	DecBounds = dpair(decmin, decmax);
}

double Boss::BossMask::operator ()(double ra, double dec) const {
	double retval = acceptance->completeness_radec(ra, dec);
	if (flat) {
		retval = retval > thresh ? 1.0 : 0.0;
	} else {
		retval = retval > thresh ? retval : 0.0;
	}

	// The rejection mask code would go here.
	// We would put in short-circuits to quickly return if any
	// of the masks return a true.

	return retval/area;
}

// Don't worry about pretty printing --- this is
void Boss::BossMask::print() {
	cout << "Normalization : " << area << endl;
	cout << "RA range : " << RABounds.first << " " << RABounds.second << endl;
}
