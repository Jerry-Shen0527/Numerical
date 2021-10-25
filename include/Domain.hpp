#pragma once
#include "MathFunctions.h"
#include "RNG.hpp"

template<typename T>
class Domain
{
public:
	virtual bool Inside(const T& value) const = 0;

	//Random sample part
	virtual T RandomSample(Float& pdf) const = 0;
};

using Domain1D = Domain<Float>;

class Interval :public Domain1D
{
public:
	Interval(Float left=0, Float right=1) : left(left), right(right) { SetPartitionCount(1); }

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

	Float scale(Float val)
	{
		return (val - left) / length();
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