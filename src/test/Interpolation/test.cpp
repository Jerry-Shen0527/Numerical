#include <iostream>
#include <type.hpp>
#include <Interpolation.h>

#include <Eigen/Eigen>
int main()
{
	std::vector<Point2> points;
	points.emplace_back(0, 0);
	points.emplace_back(1, 1);
	points.emplace_back(2, 2);
	points.emplace_back(33, 3);
	points.emplace_back(4, 48);

	NewtonPolynomial newton(points);
	LagrangianPolynomial lagrangian(points);
	newton.evaluate();
	lagrangian.evaluate();
	std::cout << newton(0.6);
}