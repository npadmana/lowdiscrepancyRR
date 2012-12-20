/*
 * BossMask.cpp
 *
 *  Created on: Dec 17, 2012
 *      Author: npadmana
 */

#include <iostream>

#include "BossMask.h"
#include "ang2d.h"

using namespace std;

Boss::BossMask::BossMask(string acceptfn, double thresh_, bool flat_) :
	thresh(thresh_),
	flat(flat_),
	RABounds(0.0, 360.0),
	DecBounds(-50.0, 90.0)  // Excessive, but what the heck
{
	acceptance.reset(new MaskClass(acceptfn));

	// The rejection mask code would go here.

	// Compute the area and cache it.
	area = 1.0;
	area = Ang2D::area(RABounds, DecBounds, 1000000, 1, *this, false)[0];
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

// Don't worry about pretty-printing, this is mostly a debugging function.
void Boss::BossMask::print() {
	cout << "The area of the acceptance mask is " << area << endl;
}
