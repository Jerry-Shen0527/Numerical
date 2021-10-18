#include <NIntegrate/Integrate.hpp>

int main()
{
	GaussIntegrate gauss;
	auto func = [](Float x) {return sqrt(x); };

	Interval interval(0, 1);
	interval.SetPartitionCount(10);

	std::cout << gauss(func, interval);
}