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