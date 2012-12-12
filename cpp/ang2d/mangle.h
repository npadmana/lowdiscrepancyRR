/*
 * mangle.h
 *
 *   Extracted from the make_random.cpp file by Martin White.
 *
 *  Created on: Dec 12, 2012
 *      Author: npadmana
 */

#ifndef MANGLE_H_
#define MANGLE_H_

namespace Mangle {

/** Class definition for a cap.
 *
 * See mangle documentation for a detailed description of what all this means.
 */
class CapClass {
private:
	double	x,y,z,cm;
public:
	// Constructor
	CapClass(double x1, double y1, double z1, double cm1) : x(x1), y(y1), z(z1), cm(cm1) {};

	// Test if in cap
	bool incap(double x0, double y0, double z0);
	bool incap(double theta, double phi);
};	// CapClass





} // end namespace mangle




#endif /* MANGLE_H_ */
