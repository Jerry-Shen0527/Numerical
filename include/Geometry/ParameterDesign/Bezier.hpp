#pragma once
#include <vector>

#include "Approximation.h"
#include "Domain.hpp"
#include "type.hpp"
#include "MathFunctions.h"

class BezierApproximation :public Approximation
{
public:
	BezierApproximation(const std::vector<Point2d>& p) : Approximation(p)
	{
		evaluate();
	}

	Float operator()(Float t)
	{
		//Would rather not use the recursion algorithm
		if (!Points.empty())
		{
			auto a = Points;
			auto size = Points.size();
			for (int i = size - 1; i >= 0; --i)
			{
				for (int j = 0; j < i; ++j)
				{
					a[j] = Lerp( interval.scale(t), a[j], a[j + 1]);
				}
			}
			return a[0].y();
		}
		else return 0;
	}

	void evaluate() override
	{
		if (!Points.empty())
		{
			std::sort(Points.begin(), Points.end(), [](const Point2d& a, const Point2d& b) {	return a.x() < b.x(); });

			interval = Interval(Points[0].x(), Points.back().x());
		}
	}

	Interval interval;
private:
};

template<typename T>
class BezierBernstein
{
public:
	BezierBernstein(const std::vector<T>& p) :Points(p) {}

	T operator()(Float t)
	{
		//Would rather not use the recursion algorithm
		if (!Points.empty())
		{
			int n = Points.size() - 1;
			T ret = T::Zero();
			for (int i = 0; i <= n; ++i)
			{
				ret += berstein(n, i, t) * Points[i];
			}
			return ret;
		}

		return T();
	}

	static Float berstein(int n, int i, Float t)
	{
		return Binomial(n, i) * pow(t, n - i) * pow(1 - t, i);
	}

	static int Binomial(int n, int k)
	{
		if (k<0 || k>n) return 0;

		int ret = 1;

		for (int i = k + 1; i <= n; ++i)
		{
			ret *= i;
		}
		for (int i = 2; i <= n - k; ++i)
		{
			ret /= i;
		}
		return  ret;
	}

	std::vector<T> Points;

private:
};