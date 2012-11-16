/** Do the angular integral over the DEEP2 mask.
 *
 * Nikhil Padmanabhan, Yale
 */

#include "ang2d.h"

#include <string>
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
	int nra, ndec;

	// Mask
	vector<float> arr;

	// Constructor
	DEEP2Mask(const string &fn);

	// Call
	double operator()(double ra, double dec);


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
	ddec = (declim1-declim0)/ddec;

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

}

double DEEP2Mask::operator ()(double ra, double dec) {
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
	return static_cast<double>(arr[idec*nra + ira]);
}


int main(int argc, char **argv) {

	// Hard code for now.
	DEEP2Mask mask1("../../../data/deep2/windowf.31.fits");

}
