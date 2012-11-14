/** Do a spherical cap integral, using the ang2d framework.
 *
 * Nikhil Padmanabhan, Yale.
 */

#include "ang2d.h"

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

