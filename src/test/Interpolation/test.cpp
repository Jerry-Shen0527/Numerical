#include <iostream>
#include <type.hpp>
#include <Interpolation.h>

int main()
{

	std::vector<Point2> points;
	points.emplace_back(0, 0);
	points.emplace_back(1, 1);
	points.emplace_back(2, 0);

	LagrangianPolynomial lagrangian(points);
	lagrangian.evaluate();
	std::cout << lagrangian(0.5);
}