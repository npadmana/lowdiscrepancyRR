#include "gtest/gtest.h"
#include "../mangle.h"
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

#ifdef DEBUGMODE
int main() {
#else
TEST(Mangle, Test1) {
#endif
	// Open the file
	Mangle::MaskClass mask("lowz.ply");
	ifstream ff("lowz.random");
	string sbuf;
	double ra, dec, wt;
	long polyid_in {0}, polyid_out {0};

	// Stop as soon as something fails
	// This prevents the code from spewing all over the place.
	getline(ff, sbuf);
	while (getline(ff, sbuf) && (polyid_in==polyid_out)) {
		istringstream(sbuf) >> ra >> dec >> polyid_in;
		wt = mask.completeness_radec(ra, dec, polyid_out);
#ifndef DEBUGMODE
		EXPECT_EQ(polyid_in, polyid_out);
#endif
	}
	ff.close();
}

#ifdef DEBUGMODE
int main2() {
#else
TEST(Mangle, Test2) {
#endif
	// Open the file
	Mangle::MaskClass mask("lowz.ply");
	// Turn off pixelization
	mask.setPixelres(-1);
	ifstream ff("lowz.random");
	string sbuf;
	double ra, dec, wt;
	long polyid_in {0}, polyid_out {0};

	// Stop as soon as something fails
	// This prevents the code from spewing all over the place.
	getline(ff, sbuf);
	while (getline(ff, sbuf) && (polyid_in==polyid_out)) {
		istringstream(sbuf) >> ra >> dec >> polyid_in;
		wt = mask.completeness_radec(ra, dec, polyid_out);
#ifndef DEBUGMODE
		EXPECT_EQ(polyid_in, polyid_out);
#endif
	}
	ff.close();
}



#ifndef DEBUGMODE
int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	int retval = RUN_ALL_TESTS();
	return retval;
}
#endif
