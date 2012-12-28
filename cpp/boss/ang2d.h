/** RR calculation using quasi-Monte Carlo methods.
 *
 * Nikhil Padmanabhan, Yale
 * November 13, 2012
 */

#ifndef ANG2D_H_
#define ANG2D_H_

#include "src/npmisc.h"
#include <fstream>
#include <list>
#include <algorithm>
#include <cmath>
#include <Eigen/Core>
#include <stdexcept>
#include "src/npgsl/npRandom.h"

namespace Ang2D {

using namespace std;
using namespace Eigen;
using namespace Padmanabhan;


/** An output container
 *
 */
class OutputData {
public :
	// Data
	list<int> nrands;
	list<dvector> vals;

	// Print out results to screen
	void print();

	// Write out results in a binary format
	void save(ofstream& out);
};


/** Define input parameters
 *
 *  Just a POD
 */
class InputParams {
public :
	int nrand, nsim; // Number of randoms, simulations
	double decmin, decmax, ramin, ramax; // Bounds for the simulation
	bool use_prng;

	vector<int> save_schedule;

	// Define the default constructor
	InputParams() : nrand(10000), nsim(1), decmin(-90.0), decmax(90.0), ramin(0.0), ramax(360.0), use_prng(false) {};
};



/** Convenience function -- shifts a random vector in a random direction mod 1
 *
 * @param x (vector<double>) : vector to be shifted (in place)
 * @param x0 (vector<double>) : shift vector
 *
 * The vectors should be the same size; no checking is done.
 */
void vecshift(vector<double> &x, vector<double> &x0);


/** Compute the weighted area of a mask
 *
 *	IMPORTANT NOTE : The code assumes RA runs from [0,360], and Dec from [-90,90]
 *
 *	The bounds (except for RAbounds, see below) are assumed to be increasing. If not, an exception is thrown.
 *  For the RABounds, if the second elt is lower, we assume we have wrapped around the sphere.
 *
 *
 *	 @param Mask1 (templated class) -- class where operator() takes RA, Dec (in degrees)
 *	                                    and returns a double (weight)
 *	 @param p [InputParams]
 *
 *	 Returns OutData element with the areas
 */
template <class Mask>
OutputData area(const Mask &Mask1, const InputParams& p) {
	const int DIM = 2;

	// Dec data validity
	if ((p.decmin < -90) || (p.decmin > 90)) throw invalid_argument("RABounds exception thrown");
	if ((p.decmax < -90) || (p.decmax > 90)) throw invalid_argument("RABounds exception thrown");
	if (p.decmax < p.decmin) throw invalid_argument("Dec Bounds exception thrown");

	// RA
	if ((p.ramin < 0) || (p.ramin > 360.0)) throw invalid_argument("RABounds exception thrown");
	if ((p.ramax < 0) || (p.ramax > 360.0)) throw invalid_argument("RABounds exception thrown");


	// Survey Bounds
	double cth1 = sin(p.decmin * D2R); // Cos theta_1
	double dcth = sin(p.decmax*D2R) - cth1; // Cos theta_2 - Cos theta_1
	double phi1 = p.ramin; // phi_1, leave in degrees
	double dphi = (p.ramax - p.ramin); // phi_2 - phi1
	if (p.ramax < p.ramin) dphi += 360; // Wrap around


	// Jacobian
	double jacobian = dcth * dphi * D2R;

	// Set up the output data
	OutputData out;

	// Set up the random number generators
	npRandom rng(1234);
	npQuasiRandom qrng(DIM);

	// Set up the random shift vectors
	// We don't need this if use_prng is set, but never mind
	vector<dvector> x0(p.nsim);
	for (dvector& x0tmp : x0) x0tmp = rng(DIM);

	// Temporary values
	dvector x, xsave;
	double ra1, dec1;
	dvector val(p.nsim, 0.0), val1(p.nsim);
	auto n1 = p.save_schedule.begin();

	// Start the loop over the points
	for (int ii=0; ii<p.nrand; ++ii) {
		if (!p.use_prng) {
			xsave = qrng(); // Save the low-discrepancy sequence
		}

		// Now loop over the simulations
		for (int jj=0; jj<p.nsim; ++jj) {

			// Generate a random pseudo-random number and shift it mod 1
			// use the regular random number generator if desired.
			if (p.use_prng) {
				// No reason to shift this
				x = rng(DIM);
			} else {
				x = xsave;
				vecshift(x,x0[jj]);
			}

			// Work out the RA, Dec of position 1
			ra1 = ((x[1] * dphi) + phi1 ); // in degrees
			if (ra1 > 360) ra1 = ra1-360;
			x[0] = x[0]*dcth + cth1;
 			dec1 = 90 - acos(x[0])/D2R;

			// Now do the actual calculation
			val[jj] += Mask1(ra1, dec1);
		}

		// Save if necessary
		if (n1 != p.save_schedule.end()) {
			if (*n1 == (ii+1)) {
				transform(val.begin(), val.end(), val1.begin(), [&ii,&jacobian](double x){return x/(ii+1) * jacobian;});
				out.nrands.push_back(ii+1);
				out.vals.push_back(val1);
				n1++;
			}
		}


	}

	// Multiply back in the jacobian
	transform(val.begin(), val.end(), val1.begin(), [&p,&jacobian](double x){return x/p.nrand * jacobian;});
	out.nrands.push_back(p.nrand);
	out.vals.push_back(val1);

	return out;
}


/** Compute the angular RR integral
 *
 *	IMPORTANT NOTE : The code assumes RA runs from [0,360], and Dec from [-90,90]
 *
 *	The bounds (except for RAbounds, see below) are assumed to be increasing. If not, an exception is thrown.
 *  For the RABounds, if the second elt is lower, we assume we have wrapped around the sphere.
 *
 *
 *	 @param Mask1 (templated class) -- class where operator() takes RA, Dec (in degrees)
 *	                                    and returns a double (weight)
 *	 @param Mask2 (templated class) -- class where operator() takes RA, Dec (in degrees)
 *	                                    and returns a double (weight)
 *
 *	 @param p [InputParams]
 *	 @param thetamin, thetamax -- range in angle to compute the integral in
 *
 *	 Returns OutData element with the areas
 */
template <class Mask>
OutputData rreval(const Mask &Mask1, const Mask &Mask2, const InputParams& p, double thetamin, double thetamax) {
	const int DIM = 4;

	// Dec data validity
	if ((p.decmin < -90) || (p.decmin > 90)) throw invalid_argument("RABounds exception thrown");
	if ((p.decmax < -90) || (p.decmax > 90)) throw invalid_argument("RABounds exception thrown");
	if (p.decmax < p.decmin) throw invalid_argument("Dec Bounds exception thrown");

	// RA
	if ((p.ramin < 0) || (p.ramin > 360.0)) throw invalid_argument("RABounds exception thrown");
	if ((p.ramax < 0) || (p.ramax > 360.0)) throw invalid_argument("RABounds exception thrown");


	// Survey Bounds
	double cth1 = sin(p.decmin * D2R); // Cos theta_1
	double dcth = sin(p.decmax*D2R) - cth1; // Cos theta_2 - Cos theta_1
	double phi1 = p.ramin; // phi_1, leave in degrees
	double dphi = (p.ramax - p.ramin); // phi_2 - phi1
	if (p.ramax < p.ramin) dphi += 360; // Wrap around
	dphi *= D2R; // Work in radians -- necessary for the rotations later

	// Scatter bounds
	double cDth1 = cos(thetamax * D2R); // Cos Delta theta_1 --- use the larger number first,
	// to avoid -ve's in the code
	double dcDth = cos(thetamin * D2R) - cDth1; // Cos Delta theta_2 - Cos Delta theta_1

	// Jacobian
	double jacobian = dcth * dphi * dcDth * 2 * PI;

	// Set up the output data
	OutputData out;

	// Set up the random number generators
	npRandom rng(1234);
	npQuasiRandom qrng(DIM);

	// Set up the random shift vectors
	// We don't need this if use_prng is set, but never mind
	vector<dvector> x0(p.nsim);
	for (dvector& x0tmp : x0) x0tmp = rng(DIM);

	// Temporary values
	dvector x, xsave;
	double ra1, dec1, ra2, dec2, tmp, cth_, sth_, phi_;
	Vector3d x1, x2; Matrix3d rotmat;;
	dvector val(p.nsim, 0.0), val1(p.nsim);
	auto n1 = p.save_schedule.begin();

	// Start the loop over the points
	for (int ii=0; ii<p.nrand; ++ii) {
		if (!p.use_prng) {
			xsave = qrng(); // Save the low-discrepancy sequence
		}

		// Now loop over the simulations
		for (int jj=0; jj<p.nsim; ++jj) {

			// Generate a random pseudo-random number and shift it mod 1
			// use the regular random number generator if desired.
			if (p.use_prng) {
				// No reason to shift this
				x = rng(DIM);
			} else {
				x = xsave;
				vecshift(x,x0[jj]);
			}

			// Work out the RA, Dec of position 1
			x[1] = ((x[1] * dphi) + phi1 );
			x[0] = x[0]*dcth + cth1;
			ra1 = x[1]/D2R; if (ra1 > 360) ra1 = ra1-360;
 			dec1 = 90 - acos(x[0])/D2R;

 			// Work out displacement vector
 			phi_ = x[3]*2 *PI;
 			cth_ = x[2]*dcDth + cDth1;
 			sth_ = sqrt(1-cth_*cth_);
 			x1 << sth_ * cos(phi_), sth_ * sin(phi_), cth_;

 			// Generate rotation matrix
 			tmp = sqrt(1-x[0]*x[0]); // sin
 			// x[0] = cos theta
 			// tmp = sin theta
 			// x[1] = phi
 			rotmat << x[0] * cos(x[1]), -sin(x[1]), tmp *cos(x[1]),
 					x[0] * sin(x[1]),  cos(x[1]), tmp *sin(x[1]),
 					tmp             ,          0, x[0];

 			x2.noalias() = rotmat * x1;

 			// Now convert back to ra, dec
 			dec2 = 90 - acos(x2.coeff(2))/D2R; // no range checking
 			ra2 = atan2(x2.coeff(1), x2.coeff(0))/D2R; // (-180, 180)
 			if (ra2 < 0) ra2 = 360+ra2;


			// Now do the actual calculation
			val[jj] += Mask1(ra1, dec1)*Mask2(ra2, dec2);
		}

		// Save if necessary
		if (n1 != p.save_schedule.end()) {
			if (*n1 == (ii+1)) {
				transform(val.begin(), val.end(), val1.begin(), [&ii,&jacobian](double x){return x/(ii+1) * jacobian;});
				out.nrands.push_back(ii+1);
				out.vals.push_back(val1);
				n1++;
			}
		}


	}

	// Multiply back in the jacobian
	transform(val.begin(), val.end(), val1.begin(), [&p,&jacobian](double x){return x/p.nrand * jacobian;});
	out.nrands.push_back(p.nrand);
	out.vals.push_back(val1);

	return out;
}




}


#endif
