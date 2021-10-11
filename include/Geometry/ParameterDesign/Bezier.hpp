#pragma once
#include <vector>

#include "type.hpp"
#include "MathFunctions.h"

template<typename T>
class Bezier
{
public:
	Bezier(const std::vector<T>& p) :points(p) {}

	T operator()(Float t)
	{
		//Would rather not use the recursion algorithm
		if (!points.empty())
		{
			auto a = points;
			auto size = points.size();
			for (int i = size - 1; i >= 0; --i)
			{
				for (int j = 0; j < i; ++j)
				{
					a[j] = Lerp((1 - t), a[j], a[j + 1]);
				}
			}
			return a[0];
		}

		return T();
	}

	std::vector<T> points;

private:
};

template<typename T>
class BezierBernstein
{
public:
	BezierBernstein(const std::vector<T>& p) :points(p) {}

	T operator()(Float t)
	{
		//Would rather not use the recursion algorithm
		if (!points.empty())
		{
			int n = points.size() - 1;
			T ret = T::Zero();
			for (int i = 0; i <= n; ++i)
			{
				ret += berstein(n, i, t) * points[i];
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

	std::vector<T> points;

private:
};