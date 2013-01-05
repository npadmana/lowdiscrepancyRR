/*
 * mangle2.h
 *
 *  An attempt to rewrite the Mangle class to make it more efficient.
 *  Note, this is being optimized for the case where you need to do many, many lookups
 *
 *  All the code is in a Mangle2 namespace --- as far as possible, I will stick to the same
 *  function names for the MaskClass as what Martin did in his original code.
 *
 *  Created on: Jan 4, 2013
 *      Author: npadmana
 */

#ifndef MANGLE2_H_
#define MANGLE2_H_

#include <vector>
#include <string>

namespace Mangle2 {

using namespace std;

/* A simple POD to store cap data
 *
 */
class CapClass {
public :
	double x, y, z, cm;
};


/* A polygon class, again a simple POD. This does not have any of the
 * caps stored, and therefore only makes sense with an associated CapList
 *
 * Users should almost never need to use this.
 */
class PolygonClass {
public :
	long polyid, pixelid;
	long ncaps, capbegin; // Number of caps, and an index to the first cap.
	double wt;
};



class MaskClass {
	// A class to handle masks.  When invoked it reads an ascii mask file.
	// The main method returns a completeness, used for making random catalogs.
private:
	typedef vector<CapClass> CapList;
	typedef vector<PolygonClass> PolyList;
	typedef vector<PolygonClass>::iterator PolyIterator;

	// Store the information
	CapList caps;
	PolyList polys;

	// Pixel information
	vector<long> pixelIndex; // Lookup table for pixels
	vector<long> numPolyInPixel; // Second piece of the lookup table.
	int	pixelres;   // pixel resolution
	char pixeltype; // pixel type


	// Pixel number for simplepix
	long pixelnum(double theta, double phi);
	// Find polygon and returns its weight and polyid, else returns 0, and -1;
	double findpoly(PolyIterator start, PolyIterator end, double x0, double y0, double z0, long &polyid);
	// Parsing routine for reading ply files
	long parsepoly(std::string sbuf, long &ncap, double &weight, long &pixel);

	/* Check if in polygon or not
	 *
	 *  @param clist : CapList -- list of caps, associated with the polygon
	 *  @param poly : PolygonClass -- polygon to check
	 *  @param x0,y0,z0 : position
	 *
	 *  You need to make sure that CapList and PolygonClass are associated. The MaskClass does this automatically.
	 */
	bool inPolygon(const PolygonClass& poly, double x0, double y0, double z0);

	/* Check if in cap or not
	 *
	 *   @param cap : CapClass
	 *   @param x0, y0, z0 : position on sphere
	 *
	 *   This could really be a static method
	 */
	bool inCap(const CapClass &cap, double x0, double y0, double z0);


public:
	MaskClass(std::string fname);
	long npolygons();

	/* Completeness routines */
	double completeness(double theta, double phi, long &polyid);
	double completeness(double theta, double phi);

	/* Input in RA, Dec in degrees */
	double completeness_radec(double ra, double dec, long &polyid);
	double completeness_radec(double ra, double dec);

	// Mostly for debugging and tests
	int getPixelres() const;
	void setPixelres(int pixelres);
};


}




#endif /* MANGLE2_H_ */
