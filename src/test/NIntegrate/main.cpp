#include <NIntegrate/Integrate.hpp>
#include <NIntegrate/2DIntegrate.hpp>
#include <NIntegrate/MCIntegrate.hpp>

int main()
{
	//GaussIntegrate gauss;
	//auto func = [](Float x) {return sqrt(x); };

	//Interval interval(0, 1);
	//interval.SetPartitionCount(10);

	//std::cout << gauss(func, interval);

	//GaussIntegrate2D gauss;
	//TriangleDomain domain;
	//Float shit;
	//auto func = [](Eigen::Vector2d vector) { return vector.x()*vector.x()*vector.y(); };

	//std::cout << gauss(func, domain);

	auto func = [](Eigen::Vector2d vector)->Float {return sin(vector.x()) * cos(vector.y()) * cos(vector.y()); };

	NDifferential<2> differential;
	auto gradient = differential(func);

	std::cout << gradient(Eigen::Vector2d(0.5, 0.5));
}