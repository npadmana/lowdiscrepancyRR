#include "gtest/gtest.h"
#include "../ang2d.h"
#include <tuple>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;
using namespace Padmanabhan;

double computeCosTheta(double ra1, double dec1, double ra2, double dec2) {
	double cosTheta=0;
	// Do the dot product out by hand for simplicity.

	cosTheta += cos(ra1*D2R)*cos(ra2*D2R);
	cosTheta += sin(ra1*D2R)*sin(ra2*D2R);
	cosTheta *= cos(dec1*D2R)*cos(dec2*D2R);
	cosTheta += sin(dec1*D2R)*sin(dec2*D2R);

	return cosTheta;
}

TEST(ZeroDisplacement, Test1) {
	double ra1, dec1, ra2, dec2, phi;

	for (int ii=0; ii < 100000; ++ii) {
		ra1 = 360.0*drand48();
		dec1 = asin(drand48()*2.0-1.0)/D2R;
		phi = 360.0*drand48();
		tie(ra2,dec2) = Ang2D::rotatePoint(ra1, dec1, 1.0, phi);
		EXPECT_NEAR(ra1, ra2, 1.e-4);
		EXPECT_NEAR(dec1, dec2, 1.e-4);
	}
}

TEST(SpecificPoints, Test1) {
	double ra2, dec2;
	tie(ra2, dec2) = Ang2D::rotatePoint(0.0,0.0,cos(1.0*D2R),0.0);
	EXPECT_NEAR(0.0, ra2, 1.e-5);
	EXPECT_NEAR(-1.0, dec2, 1.e-5);
}

TEST(SpecificPoints, Test2) {
	double ra2, dec2;
	tie(ra2, dec2) = Ang2D::rotatePoint(0.0,0.0,cos(1.0*D2R),90.0);
	EXPECT_NEAR(1.0, ra2, 1.e-5);
	EXPECT_NEAR(0.0, dec2, 1.e-5);
}

TEST(SpecificPoints, Test3) {
	double ra2, dec2;
	tie(ra2, dec2) = Ang2D::rotatePoint(0.0,0.0,cos(1.0*D2R),180.0);
	EXPECT_NEAR(0.0, ra2, 1.e-5);
	EXPECT_NEAR(1.0, dec2, 1.e-5);
}


TEST(SpecificPoints, Test4) {
	double ra2, dec2;
	tie(ra2, dec2) = Ang2D::rotatePoint(0.0,0.0,cos(1.0*D2R),270.0);
	EXPECT_NEAR(359.0, ra2, 1.e-5);
	EXPECT_NEAR(0.0, dec2, 1.e-5);
}

TEST(SpecificPoints, Test5) {
	double ra2, dec2;
	tie(ra2, dec2) = Ang2D::rotatePoint(0.0,45.0,cos(1.0*D2R),0.0);
	EXPECT_NEAR(0.0, ra2, 1.e-5);
	EXPECT_NEAR(44.0, dec2, 1.e-5);
}



class Displacement : public ::testing::TestWithParam<double> {
public:
	const int nsim=10000;
	// Empty body
};

// More detailed specific tests...
TEST_P(Displacement, Dec0) {
	double ra2, dec2, phi, ctheta;

	for (int ii=0; ii < nsim; ++ii) {
		ctheta = drand48()*2.0-1.0;
		phi = 360.0*drand48();
		tie(ra2,dec2) = Ang2D::rotatePoint(GetParam(), 0.0, ctheta, phi);
		EXPECT_NEAR(ctheta, computeCosTheta(GetParam(),0.0,ra2,dec2), 1.e-4);
	}
}

TEST_P(Displacement, Dec45) {
	double ra2, dec2, phi, ctheta;

	for (int ii=0; ii < nsim; ++ii) {
		ctheta = drand48()*2.0-1.0;
		phi = 360.0*drand48();
		tie(ra2,dec2) = Ang2D::rotatePoint(GetParam(), 45.0, ctheta, phi);
		EXPECT_NEAR(ctheta, computeCosTheta(GetParam(), 45.0,ra2,dec2), 1.e-4);
	}
}

TEST_P(Displacement, DecNeg45) {
	double ra2, dec2, phi, ctheta;

	for (int ii=0; ii < nsim; ++ii) {
		ctheta = drand48()*2.0-1.0;
		phi = 360.0*drand48();
		tie(ra2,dec2) = Ang2D::rotatePoint(GetParam(), -45.0, ctheta, phi);
		EXPECT_NEAR(ctheta, computeCosTheta(GetParam(), -45.0,ra2,dec2), 1.e-4);
	}
}

// Stress test everything!
TEST(StressDisplacement,  Test1) {
	double ra1,dec1,ra2, dec2, phi, ctheta;
	int nsim = 100000;

	for (int ii=0; ii < nsim; ++ii) {
		ra1 = 360.0*drand48();
		dec1 = asin(drand48()*2.0-1.0)/D2R;
		ctheta = drand48()*2.0-1.0;
		phi = 360.0*drand48();
		tie(ra2,dec2) = Ang2D::rotatePoint(ra1, dec1, ctheta, phi);
		EXPECT_NEAR(ctheta, computeCosTheta(ra1,dec1,ra2,dec2), 1.e-4);
	}
}


INSTANTIATE_TEST_CASE_P(DisplacementTest1,
                        Displacement,
                        ::testing::Values(0.0));
                        //::testing::Values(0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0));


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	int retval = RUN_ALL_TESTS();
	return retval;
}



