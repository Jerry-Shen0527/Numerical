#pragma once
#include "type.hpp"
#include <RNG.hpp>
#include <MathFunctions.h>

template<typename T>
class Domain
{
public:
	virtual bool Inside(const T& value) const = 0;

	//Random sample part
	virtual T RandomSample(Float& pdf) const = 0;
};

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

using Domain1D = Domain<Float>;

class Interval :public Domain1D
{
public:
	Interval(Float left, Float right) : left(left), right(right) { SetPartitionCount(1); }

	void SetPartitionCount(int partition)
	{
		n = partition;
		h = length() / Float(n);
	}

	int GetPartitionCount() const
	{
		return n;
	}

	bool Inside(const Float& value) const override
	{
		return value >= left && value < right;
	}

	Float RandomSample(Float& pdf) const override
	{
		pdf = 1.0 / (length());
		return  RandomFloat(left, right);
	}

	Float lerp(Float t) const
	{
		return Lerp(t, left, right);
	}

	Float length() const
	{
		return right - left;
	}

	Interval SubInterval(int idx) { return Interval(left + idx * h, left + (idx + 1) * h); }

	//Scale a f([left,right]) to f([0,1]), remember to multiply the factor when scaling derivatives
	std::function<Float(Float)> scale(const std::function<Float(Float)>& func, Float factor = 1.0)
	{
		return [factor, &func, this](Float val)
		{
			return factor * func((val - left) / length());
		};
	}

private:
	Float left;
	Float right;
	Float h = 0;
	int n = 1;

	Float GetMeasure() const
	{
		return length() / static_cast<Float>(n);
	}
};

template<typename T>
class Union :Domain<T>
{
public:
	Union(const Domain<T>& domain1, const Domain<T>& domain2) :domain_1(domain1), domain_2(domain2) {}
	bool Inside(const T& value) const override;
	T RandomSample(Float& pdf) const override;

private:
	Domain<T> domain_1;
	Domain<T> domain_2;
};

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