/*
 * BossMask.h
 *
 *  Created on: Dec 17, 2012
 *      Author: npadmana
 */

#ifndef BOSSMASK_H_
#define BOSSMASK_H_

#include "mangle2.h"
#include <memory>
#include <string>
#include <utility>

namespace Boss {

using namespace std;
using namespace Mangle2;

/** BossMask --
 *    a generic class to handle BOSS window functions
 *    The window function is assumed to be specified by a single acceptance mask,
 *    and multiple rejection masks.
 *
 *    TODO : This is currently unimplemented.
 *
 *    The acceptance mask can also have completeness numbers, and we can
 *    specify a threshold below which the completeness is returned to be zero.
 *
 *
 */
class BossMask {
private :
	typedef std::unique_ptr<MaskClass> maskptr;

	// Define the threshold
	double thresh;
	bool flat;

	// Define the acceptance mask
	string acceptfn;
	maskptr acceptance;

	// Area
	double area;

public :
	// Define RA and Dec bounds for the survey
	// Make these publicly accessible.
	double ramin, ramax, decmin, decmax;

	/** Constructor
	 *
	 *  @param acceptfn [string] acceptance mask file name
	 *  @param threshold [double] threshold
	 */
	BossMask(string configfn);

	/** Return the completeness
	 *
	 */
	double operator()(double ra, double dec) const;

	/** Print information about the mask
	 *
	 */
	void print();
};

} /* namespace Boss */
#endif /* BOSSMASK_H_ */
