#include "gtest/gtest.h"
#include "../ang2d.h"
#include <tuple>

using namespace std;

TEST(Partition, Test1) {
	int nra, ndec;
	tie(ndec, nra) = Ang2D::partition(1);
	EXPECT_EQ(ndec, 1);
	EXPECT_EQ(nra, 1);
}


TEST(Partition, Test2) {
	int nra, ndec;
	tie(ndec, nra) = Ang2D::partition(2);
	EXPECT_EQ(ndec, 1);
	EXPECT_EQ(nra, 2);
}


TEST(Partition, Test4) {
	int nra, ndec;
	tie(ndec, nra) = Ang2D::partition(4);
	EXPECT_EQ(ndec, 2);
	EXPECT_EQ(nra, 2);
}

TEST(Partition, Test8) {
	int nra, ndec;
	tie(ndec, nra) = Ang2D::partition(8);
	EXPECT_EQ(ndec, 2);
	EXPECT_EQ(nra, 4);
}

TEST(Partition, Test32) {
	int nra, ndec;
	tie(ndec, nra) = Ang2D::partition(32);
	EXPECT_EQ(ndec, 4);
	EXPECT_EQ(nra, 8);
}

TEST(Partition, Test64) {
	int nra, ndec;
	tie(ndec, nra) = Ang2D::partition(64);
	EXPECT_EQ(ndec, 8);
	EXPECT_EQ(nra, 8);
}

