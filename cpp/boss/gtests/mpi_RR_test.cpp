#include "gtest/gtest.h"

#include <string>

#include "../ang2d.h"
#include "../BossMask.h"

#include "mpi.h"

using namespace std;
using namespace Boss;

class RRTest : public ::testing::TestWithParam<int> {
public:
	double calcRR(string fn, double thetamin, double thetamax) {
		Ang2D::InputParams p0(fn);
		// Now define the mask
		BossMask mask1(p0.mask1fn);
		BossMask mask2(p0.mask2fn);
		// Define various bounds
		int nra, ndec, irank, nproc;
		MPI_Comm_rank(MPI_COMM_WORLD, &irank);
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);
		tie(ndec, nra) = Ang2D::partition(nproc);
		Ang2D::setBounds(nra,ndec,irank,mask1,p0);
		Ang2D::OutputData out = Ang2D::rreval(mask1, mask2, p0, thetamin, thetamax);
		out.finalize();
		return out.mean();
	}
};

TEST_F(RRTest, AbsoluteTestSphereBin1) {
	string inputfn {"sphere_input.cfg"};
	double testbase = calcRR(inputfn, 0.0, 1.0);
	double testabs = 0.000076152421804380421494; // Calculated from Mathematica
	EXPECT_NEAR(1.0, testabs/testbase, 1.e-3);
}


TEST_F(RRTest, AbsoluteTestSphereBin2) {
	string inputfn {"sphere_input.cfg"};
	double testbase = calcRR(inputfn, 3.0, 5.0);
	double testabs = 0.0012174183314141707447; // Calculated from Mathematica
	EXPECT_NEAR(1.0, testabs/testbase, 1.e-3);
}


TEST_F(RRTest, AbsoluteTestHemisphereBin1) {
	string inputfn {"hemisphere_input.cfg"};
	double testbase = calcRR(inputfn, 0.0, 1.0);
	double testabs = 0.000151736; // Calculated from Mathematica
	EXPECT_NEAR(1.0, testabs/testbase, 1.e-3);
}

TEST_F(RRTest, AbsoluteTestHemisphereBin2) {
	string inputfn {"hemisphere_input.cfg"};
	double testbase = calcRR(inputfn, 3.0, 5.0);
	double testabs = 0.0023797; // Calculated from Mathematica
	EXPECT_NEAR(1.0, testabs/testbase, 1.e-3);
}

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	::testing::InitGoogleTest(&argc, argv);
	int retval = RUN_ALL_TESTS();
	MPI_Finalize();
	return retval;
}

