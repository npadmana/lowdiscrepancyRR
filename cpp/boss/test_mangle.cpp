/* A simple code that generates a grid in RA, Dec and runs them
 * through a simple mask.
 *
 * This program has no parameters, it is just designed to exercise the code.
 *
 * It assumes a mangle polygon file called boss_survey.ply is available in the
 * directory.
 */

#include <iostream>
#include <cmath>
#include "mangle.h"

using namespace std;
using namespace Mangle;

int main() {
	MaskClass boss("boss_survey.ply");

	// Start the loop
	double dec;
	for (double ra=0.0; ra<360.0; ++ra)
		for (double sindec=-0.95; sindec<0.95; sindec += 0.01) {
			dec = asin(sindec) * 180.0/M_PI;
			if (boss.completeness_radec(ra, dec)>0.9)
				cout << ra << " " << dec << endl;
		}

}
