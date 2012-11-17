/** RR calculation using quasi-Monte Carlo methods.
 *
 * Nikhil Padmanabhan, Yale
 * November 13, 2012
 */

#ifndef ANG2D_H_
#define ANG2D_H_

#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>
#include <Eigen/Core>
#include "src/npgsl/npRandom.h"

namespace Ang2D {

using namespace std;
using namespace Eigen;


const double PI = 3.1415926535897932385;
const double D2R = 0.017453292519943295769; // Degrees to radians

typedef pair<double, double> dpair;

/** Convenience function -- shifts a random vector in a random direction mod 1
 *
 * @param x (vector<double>) : vector to be shifted (in place)
 * @param x0 (vector<double>) : shift vector
 *
 * The vectors should be the same size; no checking is done.
 */
void vecshift(vector<double> &x, vector<double> &x0) {
	double tmp;
	transform(x0.begin(), x0.end(), x.begin(), x.begin(),
			[&tmp](double a, double b){return modf(a+b,&tmp);});
}

/** Compute the angular RR integral
 *
 *	IMPORTANT NOTE : The code assumes RA runs from [0,360], and Dec from [-90,90]
 *
 *	The bounds (except for RAbounds, see below) are assumed to be increasing. If not, an exception is thrown.
 *  For the RABounds, if the second elt is lower, we assume we have wrapped around the sphere.
 *
 *  The last element of thetaBin cannot be > 180.0, and the first element must be >= 0.
 *
 *  The mask() operator should return a quantity which has inverse area dimensions, since this
 *  allows the most easy transformation to the usual RR terms. The user is responsible for this,
 *  although we provide a QMC routine to calculate this (area).
 *
 *   @param RABounds (pair<double, double>)  -- RA bounds to consider, should correspond normally to
 *                                              mask 1 (degrees)
 *   @param DecBounds (pair<double, double>) -- DecBounds to consider (degrees)
 *   @param thetaBin (pair<double, double>) -- min, max of the bin (degrees)
 *	 @param nrand (int) -- number of randoms
 *	 @param nsim (int) -- number of simulations
 *	 @param const Mask1 (templated class) -- class where operator() takes RA, Dec (in degrees)
 *	                                    and returns a double (weight)
 *	 @param const Mask2 (templated class) -- class where operator() takes RA, Dec (in degrees)
 *	                                   and returns a double (weight)
 *
 *	 Returns vector<double> rr[nsim] --- nsim measurements of RR
 */
template <class Mask>
vector<double> rreval(dpair RABounds, dpair DecBounds, dpair thetaBin,
		int nrand, int nsim, const Mask &Mask1, const Mask &Mask2) {
	const int DIM = 4;

	// Dec data validity
	if ((DecBounds.first < -90) || (DecBounds.first > 90)) throw ("RABounds exception thrown");
	if ((DecBounds.second < -90) || (DecBounds.second > 90)) throw ("RABounds exception thrown");
	if (DecBounds.second < DecBounds.first) throw("Dec Bounds exception thrown");

	// theta data validity
	if (thetaBin.second < thetaBin.first) throw("thetaBin exception thrown");
	if ((thetaBin.first < 0) || (thetaBin.first > 180)) throw("thetaBin exception thrown");
	if ((thetaBin.second < 0) || (thetaBin.second > 180)) throw("thetaBin exception thrown");

	// RA
	if ((RABounds.first < 0) || (RABounds.first > 360.0)) throw ("RABounds exception thrown");
	if ((RABounds.second < 0) || (RABounds.second > 360.0)) throw ("RABounds exception thrown");
	if (RABounds.second < RABounds.first) RABounds.second += 360; // Wrap around


	// Survey Bounds
	double cth1 = sin(DecBounds.first * D2R); // Cos theta_1
	double dcth = sin(DecBounds.second*D2R) - cth1; // Cos theta_2 - Cos theta_1
	double phi1 = RABounds.first * D2R; // phi_1
	double dphi = (RABounds.second - RABounds.first) * D2R; // phi_2 - phi1

	// Scatter bounds
	double cDth1 = cos(thetaBin.second * D2R); // Cos Delta theta_1 --- use the larger number first,
	// to avoid -ve's in the code
	double dcDth = cos(thetaBin.first * D2R) - cDth1; // Cos Delta theta_2 - Cos Delta theta_1

	// Jacobian
	double jacobian = dcth * dphi * dcDth * 2 * PI;

	// Set up the output vector
	vector<double> outvec(nsim);

	// Set up the random vector
	vector<double> x0(DIM);
	npRandom rng(1234);

	// Temporary values
	double ra1, dec1, ra2, dec2, out, tmp, cth_, sth_, phi_;
	Vector3d x1, x2; Matrix3d rotmat;

	for(int jj=0; jj<nsim; ++jj) {

		for_each(x0.begin(), x0.end(), [&rng](double &x){x=rng();});

		// Set up the quasi RNG
		npQuasiRandom qrng(DIM);

		// Initialize the integrator
		out=0.0;

		// Loop over points
		for (int ii=0; ii<nrand; ++ii) {

			// Generate a random pseudo-random number and shift it mod 1
			vector<double> x = qrng();
			vecshift(x,x0);

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
			out += Mask1(ra1, dec1) * Mask2(ra2, dec2);

		}

		// Multiply back in the jacobian
		outvec[jj]  = out/(nrand) * jacobian;
	}
	return outvec;
}

/** Compute the weighted area of a mask
 *
 *	IMPORTANT NOTE : The code assumes RA runs from [0,360], and Dec from [-90,90]
 *
 *	The bounds (except for RAbounds, see below) are assumed to be increasing. If not, an exception is thrown.
 *  For the RABounds, if the second elt is lower, we assume we have wrapped around the sphere.
 *
 *   @param RABounds (pair<double, double>)  -- RA bounds to consider (degrees)
 *   @param DecBounds (pair<double, double>) -- DecBounds to consider (degrees)
 *	 @param nrand (int) -- number of randoms
 *	 @param nsim (int) -- number of simulations
 *	 @param const Mask1 (templated class) -- class where operator() takes RA, Dec (in degrees)
 *	                                    and returns a double (weight)
 *
 *	 Returns vector<double> rr[nsim] --- nsim measurements of RR
 */
template <class Mask>
vector<double> area(dpair RABounds, dpair DecBounds,
		int nrand, int nsim, const Mask &Mask1) {
	const int DIM = 2;

	// Dec data validity
	if ((DecBounds.first < -90) || (DecBounds.first > 90)) throw ("RABounds exception thrown");
	if ((DecBounds.second < -90) || (DecBounds.second > 90)) throw ("RABounds exception thrown");
	if (DecBounds.second < DecBounds.first) throw("Dec Bounds exception thrown");

	// RA
	if ((RABounds.first < 0) || (RABounds.first > 360.0)) throw ("RABounds exception thrown");
	if ((RABounds.second < 0) || (RABounds.second > 360.0)) throw ("RABounds exception thrown");
	if (RABounds.second < RABounds.first) RABounds.second += 360; // Wrap around


	// Survey Bounds
	double cth1 = sin(DecBounds.first * D2R); // Cos theta_1
	double dcth = sin(DecBounds.second*D2R) - cth1; // Cos theta_2 - Cos theta_1
	double phi1 = RABounds.first; // phi_1, leave in degrees
	double dphi = (RABounds.second - RABounds.first); // phi_2 - phi1

	// Jacobian
	double jacobian = dcth * dphi * D2R;

	// Set up the output vector
	vector<double> outvec(nsim);

	// Set up the random vector
	vector<double> x0(DIM);
	npRandom rng(1234);

	// Temporary values
	double ra1, dec1, out;

	for(int jj=0; jj<nsim; ++jj) {

		for_each(x0.begin(), x0.end(), [&rng](double &x){x=rng();});

		// Set up the quasi RNG
		npQuasiRandom qrng(DIM);

		// Initialize the integrator
		out=0.0;

		// Loop over points
		for (int ii=0; ii<nrand; ++ii) {

			// Generate a random pseudo-random number and shift it mod 1
			vector<double> x = qrng();
			vecshift(x,x0);

			// Work out the RA, Dec of position 1
			ra1 = ((x[1] * dphi) + phi1 ); // in degrees
			if (ra1 > 360) ra1 = ra1-360;
			x[0] = x[0]*dcth + cth1;
 			dec1 = 90 - acos(x[0])/D2R;

			// Now do the actual calculation
			out += Mask1(ra1, dec1);

		}

		// Multiply back in the jacobian
		outvec[jj]  = out/(nrand) * jacobian;
	}
	return outvec;
}




}


#endif
