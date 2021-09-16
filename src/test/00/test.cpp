#include <iostream>
#include <type.hpp>

int main()
{
	Eigen::Vector3f vec(1, 2, 3);
	std::cout << vec.transpose();
}