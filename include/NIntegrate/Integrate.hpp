#pragma once
#include "type.hpp"
#include "Domain.hpp"
#include <RNG.hpp>
#include <MathFunctions.h>

template<typename T>
class Integrate
{
public:
	Integrate();
	virtual ~Integrate();

	virtual Float operator() (const std::function<Float(T)>& func, const Domain<T>& domain) = 0;

protected:
	/**
	 * \brief Controls the upper bound of error
	 */
	Float precision;
};

template<typename T>
Integrate<T>::Integrate()
{
}

template<typename T>
Integrate<T>::~Integrate()
{
}

class SubdivIntegrate :public Integrate<Float>
{
public:
	Float operator()(const std::function<Float(Float)>& func, const Domain<Float>& domain) override
	{
		Interval interval = static_cast<const Interval&>(domain);
		Float ret = 0;

		for (int i = 0; i < interval.GetPartitionCount(); ++i)
		{
			ret += integrate(func, interval.SubInterval(i));
		}

		return ret;
	}

	virtual Float integrate(const std::function<Float(Float)>& func, const Interval& interval) const = 0;
};

class TrapeziumIntegrate :public SubdivIntegrate
{
private:
	Float integrate(const std::function<Float(Float)>& func, const Interval& interval) const override
	{
		return (func(interval.lerp(0.0)) + func(interval.lerp(1.0))) * interval.length() / 2.0;
	}
};

class SimpsonIntegrate :public SubdivIntegrate
{
	Float integrate(const std::function<Float(Float)>& func, const Interval& interval) const override
	{
		return (func(interval.lerp(0.0)) + 4.0 * func(interval.lerp(0.5)) + func(interval.lerp(1.0))) * interval.length() / 6.0;
	}
};

class GaussIntegrate :public SubdivIntegrate
{
public:
	GaussIntegrate() {}

	Float integrate(const std::function<Float(Float)>& func, const Interval& interval) const override
	{
		Float ret = 0;
		Float x1 = 0.538469310105683;
		Float x2 = 0.906179845938664;
		Float A0 = 0.568888888888889;
		Float A1 = 0.478628670499366;
		Float A2 = 0.236926885056189;

		ret += func(interval.lerp(0.5)) * A0;
		ret += A1 * (func(interval.lerp(0.5 + x1 / 2.0)) + func(interval.lerp(0.5 - x1 / 2.0)));
		ret += A2 * (func(interval.lerp(0.5 + x2 / 2.0)) + func(interval.lerp(0.5 - x2 / 2.0)));
		return  ret * interval.length() / 2.0;
	}
};

Float L2InnerProduct(const std::function<Float(Float)>& func1, const std::function<Float(Float)>& func2, const Interval& interval)
{
	auto func = [&](Float val) {return func1(val) * func2(val); };
	return GaussIntegrate()(func, interval);
}

Float WeightedL2InnerProduct(const std::function<Float(Float)>& func1, const std::function<Float(Float)>& func2, const std::function<Float(Float)>& weight, const Interval& interval)
{
	auto func = [&](Float val) {return func1(val) * func2(val) * weight(val); };
	return GaussIntegrate()(func, interval);
}