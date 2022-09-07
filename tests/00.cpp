#include <iostream>
#include<gtest/gtest.h>

TEST(eq,eq1)
{
	EXPECT_EQ(1, 1);
}

int main()
{
	testing::InitGoogleTest();
	RUN_ALL_TESTS();
}