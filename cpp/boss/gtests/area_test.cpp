#include "gtest/gtest.h"
#include "../ang2d.h"
#include "../BossMask.h"

#include "mpi.h"

using namespace std;
using namespace Boss;

class HemisphereAreaTest : public ::testing::TestWithParam<int> {
	// Empty test body
};

TEST_P(HemisphereAreaTest, MultiSection) {
	int nproc = GetParam();
	Ang2D::InputParams p0("hemisphere_input.cfg");
	// Now define the mask
	BossMask mask1(p0.mask1fn);
	// Define various bounds
	int nra, ndec;
	tie(ndec, nra) = Ang2D::partition(nproc);

	double areasum = 0.0;
	for (int irank=0; irank < nproc; ++irank) {
		Ang2D::setBounds(nra,ndec,irank,mask1,p0);
		Ang2D::OutputData out = Ang2D::area(mask1, p0);
		out.finalize();
		areasum += out.mean();
	}
	EXPECT_NEAR(1.0, areasum, 1.e-3);
}

INSTANTIATE_TEST_CASE_P(HemisphereAreas,
                        HemisphereAreaTest,
                        ::testing::Values(1,2,4,8,16,32));


int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	::testing::InitGoogleTest(&argc, argv);
	int retval = RUN_ALL_TESTS();
	MPI_Finalize();
	return retval;
}
