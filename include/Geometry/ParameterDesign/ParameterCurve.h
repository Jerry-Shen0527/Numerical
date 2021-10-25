#pragma once
#include "Bezier.hpp"
#include "MathFunctions.h"

template<typename T, int dim>
class ParameterCurve
{
protected:
	std::vector<T> control_points;

	virtual void evaluate() = 0;
	std::function<Float(Float)> function[dim];

public:
	virtual ~ParameterCurve() = default;

	ParameterCurve(const std::vector<T>& vec) :control_points(vec)
	{
	}

	virtual T operator()(Float t)
	{
		if (!control_points.empty())
		{
			t = Clamp(t);
			T ret;
			for (int i = 0; i < dim; ++i)
			{
				ret[i] = function[i](t);
			}
			return ret;
		}
		return T();
	}
};

class BezierCurve2D : public ParameterCurve<Point2, 2>
{
public:
	explicit BezierCurve2D(const std::vector<Point2>& points)
		: ParameterCurve<Point2, 2>(points)
	{
	}

	void evaluate() override
	{
		int count = control_points.size();
		for (int dim_i = 0; dim_i < 2; ++dim_i)
		{
			Float h = 1.0 / (count - 1);
			std::vector<Point2> points;
			for (int i = 0; i < count; ++i)
			{
				points.emplace_back(i * h, control_points[i][dim_i]);
			}

			BezierApproximation bezier(points);
			function[dim_i] = bezier;
		}
	}
};
