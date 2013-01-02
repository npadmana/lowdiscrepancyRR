#include "gtest/gtest.h"

#include <string>

#include "../ang2d.h"
#include "../BossMask.h"

#include "mpi.h"

using namespace std;
using namespace Boss;

class AreaTest : public ::testing::TestWithParam<int> {
public:
	double calcArea(string fn) {
		Ang2D::InputParams p0(fn);
		// Now define the mask
		BossMask mask1(p0.mask1fn);
		// Define various bounds
		int nra, ndec, irank, nproc;
		MPI_Comm_rank(MPI_COMM_WORLD, &irank);
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);
		tie(ndec, nra) = Ang2D::partition(nproc);
		Ang2D::setBounds(nra,ndec,irank,mask1,p0);
		Ang2D::OutputData out = Ang2D::area(mask1, p0);
		out.finalize();
		return out.mean();
	}
};

TEST_F(AreaTest, Hemisphere) {
	EXPECT_NEAR(1.0, calcArea("hemisphere_input.cfg"), 1.e-3);
}

TEST_F(AreaTest, Sphere) {
	EXPECT_NEAR(1.0, calcArea("sphere_input.cfg"), 1.e-3);
}


int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	::testing::InitGoogleTest(&argc, argv);
	int retval = RUN_ALL_TESTS();
	MPI_Finalize();
	return retval;
}
