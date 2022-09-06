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

	TriangleDomain triangle;

	auto func1 = [](Eigen::Vector3d vector)->Float
	{
		return vector.x();
	};
	auto func2 = [](Eigen::Vector3d vector)->Float { return vector.x(); };

	using Eigen::Vector3d;
	TriangleDomain triangle2(Vector3d(0, 0, 0), Vector3d(0.5, 0, 0), Vector3d(0, 0.5, 0));

	NDifferential<3> gradient;

	std::cout << L2InnerProduct2D<3>(gradient(func1), gradient(func2), triangle)<< std::endl;
	std::cout << L2InnerProduct2D<3>(gradient(triangle2.remap(triangle.scale(func1))), gradient(triangle2.remap(triangle.scale(func2))), triangle2);

	//std::cout << gradient(func1)(Eigen::Vector3d(0.5, 0.5,0));
}