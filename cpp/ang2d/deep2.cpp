/** Do the angular integral over the DEEP2 mask.
 *
 * Nikhil Padmanabhan, Yale
 */

#include "ang2d.h"

#include <string>

#include "fitsio.h"

using namespace std;

// Use the DEEP2 mask
class DEEP2Mask {
public:

	// RA, Dec limits and dra, ddec
	double ralim0, ralim1, declim0, declim1, dra, ddec;


	//
	DEEP2Mask(const string &fn);


};


DEEP2Mask::DEEP2Mask(const string &fn) {
	// Open the DEEP2 mask file


	// Close the DEEP2 mask file

}

int main(int argc, char **argv) {

	// Hard code for now.
	DEEP2Mask mask1("../../../data/deep2/windowf.31.fits");

}
