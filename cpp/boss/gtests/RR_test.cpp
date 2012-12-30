#include "gtest/gtest.h"

#include <string>

#include "../ang2d.h"
#include "../BossMask.h"

#include "mpi.h"

using namespace std;
using namespace Boss;

class RRTest : public ::testing::TestWithParam<int> {
public:
	double calcRR(string fn, int nproc, double thetamin, double thetamax) {
		Ang2D::InputParams p0(fn);
		// Now define the mask
		BossMask mask1(p0.mask1fn);
		BossMask mask2(p0.mask2fn);
		// Define various bounds
		int nra, ndec;
		tie(ndec, nra) = Ang2D::partition(nproc);

		double RRsum = 0.0;
		for (int irank=0; irank < nproc; ++irank) {
			Ang2D::setBounds(nra,ndec,irank,mask1,p0);
			Ang2D::OutputData out = Ang2D::rreval(mask1, mask2, p0, thetamin, thetamax);
			out.finalize();
			RRsum += out.mean();
		}
		return RRsum;
	}
};

// This test checks to see that splitting gets the same results as the
// unsplit case. Note that it doesn't check to see if the first case is correct.
TEST_P(RRTest, RelativeTestSphereBin1) {
	string inputfn {"sphere_input.cfg"};
	double testbase = calcRR(inputfn, 1, 0.0, 1.0);
	EXPECT_NEAR(1.0, calcRR(inputfn, GetParam(), 0.0, 1.0)/testbase, 1.e-3);
}

TEST_P(RRTest, RelativeTestHemisphereSphereBin1) {
	string inputfn {"hemisphere_input.cfg"};
	double testbase = calcRR(inputfn, 1, 0.0, 1.0);
	EXPECT_NEAR(1.0, calcRR(inputfn, GetParam(), 0.0, 1.0)/testbase, 1.e-3);
}

TEST_P(RRTest, RelativeTestSphereBin2) {
	string inputfn {"sphere_input.cfg"};
	double testbase = calcRR(inputfn, 1, 3.0, 5.0);
	EXPECT_NEAR(1.0, calcRR(inputfn, GetParam(), 3.0, 5.0)/testbase, 1.e-3);
}

TEST_P(RRTest, RelativeTestHemisphereSphereBin2) {
	string inputfn {"hemisphere_input.cfg"};
	double testbase = calcRR(inputfn, 1, 3.0, 5.0);
	EXPECT_NEAR(1.0, calcRR(inputfn, GetParam(), 3.0, 5.0)/testbase, 1.e-3);
}

TEST_F(RRTest, AbsoluteTestSphereBin1) {
	string inputfn {"sphere_input.cfg"};
	double testbase = calcRR(inputfn, 1, 0.0, 1.0);
	double testabs = 0.000076152421804380421494; // Calculated from Mathematica
	EXPECT_NEAR(1.0, testabs/testbase, 1.e-3);
}


TEST_F(RRTest, AbsoluteTestSphereBin2) {
	string inputfn {"sphere_input.cfg"};
	double testbase = calcRR(inputfn, 1, 3.0, 5.0);
	double testabs = 0.0012174183314141707447; // Calculated from Mathematica
	EXPECT_NEAR(1.0, testabs/testbase, 1.e-3);
}


INSTANTIATE_TEST_CASE_P(MultiSectionRRTest,
                        RRTest,
                        ::testing::Values(2,4,8,16,32));


int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	::testing::InitGoogleTest(&argc, argv);
	int retval = RUN_ALL_TESTS();
	MPI_Finalize();
	return retval;
}
