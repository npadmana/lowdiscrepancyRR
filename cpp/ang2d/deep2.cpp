/** Do the angular integral over the DEEP2 mask.
 *
 * Nikhil Padmanabhan, Yale
 */

#include "ang2d.h"

#include <string>
#include <numeric>
#include <iostream>
#include <fstream>
#include <chrono>
#include "boost/format.hpp"
#include "boost/program_options.hpp"


#include "fitsio.h"

using namespace std;
using boost::format;
using namespace chrono;

namespace po = boost::program_options;


typedef pair<double, double> dpair;




// Use the DEEP2 mask
class DEEP2Mask {
public:

	// RA, Dec limits and dra, ddec
	double ralim0, ralim1, declim0, declim1, dra, ddec;
	double weighted_area;
	int nra, ndec;

	// Mask
	vector<float> arr;

	// Constructor
	DEEP2Mask(const string &fn);

	// Call
	double operator()(double ra, double dec) const;

	// Weighted area calculation
	double calc_area();


};


DEEP2Mask::DEEP2Mask(const string &fn) {
	int status=0;
	fitsfile *fptr;

	// Open the DEEP2 mask file
	fits_open_file(&fptr, fn.c_str(), READONLY, &status);

	// Read in the header information
	fits_read_key(fptr, TDOUBLE, "RALIM0", &ralim0, NULL, &status);
	fits_read_key(fptr, TDOUBLE, "RALIM1", &ralim1, NULL, &status);
	fits_read_key(fptr, TDOUBLE, "DECLIM0", &declim0, NULL, &status);
	fits_read_key(fptr, TDOUBLE, "DECLIM1", &declim1, NULL, &status);
	fits_read_key(fptr, TINT, "NRA", &nra, NULL, &status);
	fits_read_key(fptr, TINT, "NDEC", &ndec, NULL, &status);
	dra = (ralim1-ralim0)/nra;
	ddec = (declim1-declim0)/ndec;

	// Allocate space
	arr.resize(nra*ndec,0.0);

	// Read in the image
	long fpixel[2] = {1,1};
	float nulval=0;
	int anynul;
	fits_read_pix(fptr, TFLOAT, fpixel, nra*ndec, &nulval, &arr[0], &anynul, &status);

	// Close the DEEP2 mask file
	fits_close_file(fptr, &status);

	if (status) throw("Error reading the FITS file");

	cout << format("RA limits from %10.6f to %10.6f in %i pixels\n") % ralim0 % ralim1 %nra;
	cout << format("Dec limits from %10.6f to %10.6f in %i pixels\n") % declim0 % declim1 %ndec;
	weighted_area = calc_area();
	cout << format("Weighted mask area = %10.6e\n")%weighted_area;

}

double DEEP2Mask::operator ()(double ra, double dec) const {
	int ira, idec;

	// ra, dec limits
	if ((ra < ralim0) || (ra >= ralim1)) return 0.0;
	if ((dec < declim0) || (dec >= declim1)) return 0.0;

	// positions
	ira = static_cast<int>((ra - ralim0)/dra);
	idec = static_cast<int>((dec-declim0)/ddec);

	// be paranoid
	if ((ira < 0) || (ira >= nra)) return 0.0;
	if ((idec < 0) || (idec >= ndec)) return 0.0;

	// Now return
	return static_cast<double>(arr[idec*nra + ira])/weighted_area;
}

double DEEP2Mask::calc_area() {
	double pixarea, dec0, dec1, area=0.0;

	// Loop over rows in dec
	for (int idec=0; idec < ndec; ++idec) {
		dec0 = idec*ddec; dec1 = dec0+ddec;
		pixarea =(sin(dec1*Ang2D::D2R) - sin(dec0*Ang2D::D2R)) * dra * Ang2D::D2R;
		area += accumulate(arr.begin()+idec*nra, arr.begin()+(idec+1)*nra, 0.0)*pixarea;
	}

	return area;

}


int main(int argc, char **argv) {
	int nrand, nsim;
	double thetamin, thetamax;
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
	DEEP2Mask mask1(maskfn);

	// Define various bounds
	dpair RABounds(mask1.ralim0, mask1.ralim1);
	dpair decBounds(mask1.declim0, mask1.declim1);
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
