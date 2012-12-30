#include "gtest/gtest.h"

#include <string>

#include "../ang2d.h"
#include "../BossMask.h"

#include "mpi.h"

using namespace std;
using namespace Boss;

class AreaTest : public ::testing::TestWithParam<int> {
public:
	double calcArea(string fn, int nproc) {
		Ang2D::InputParams p0(fn);
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
		return areasum;
	}
};

TEST_P(AreaTest, Hemisphere) {
	EXPECT_NEAR(1.0, calcArea("hemisphere_input.cfg", GetParam()), 1.e-3);
}

TEST_P(AreaTest, Sphere) {
	EXPECT_NEAR(1.0, calcArea("sphere_input.cfg", GetParam()), 1.e-3);
}

INSTANTIATE_TEST_CASE_P(MultiSectionAreaTest,
                        AreaTest,
                        ::testing::Values(1,2,4,8,16,32));


int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	::testing::InitGoogleTest(&argc, argv);
	int retval = RUN_ALL_TESTS();
	MPI_Finalize();
	return retval;
}
