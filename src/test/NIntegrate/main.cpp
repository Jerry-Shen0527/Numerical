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

	GaussIntegrate2D gauss;
	TriangleDomain domain;
	Float shit;
	std::cout << domain.RandomSample(shit);
	auto func = [](Eigen::Vector2f vector) { return static_cast<Float>(1.0); };

	std::cout << gauss(func, domain);
}
