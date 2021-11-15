#pragma once
#include "Bezier.hpp"
#include "Interpolation.h"
#include "MathFunctions.h"

template<typename T, int dim>
class ParameterCurve
{
protected:
	std::vector<T> control_points;

	std::function<Float(Float)> function[dim];

public:
	virtual ~ParameterCurve() = default;
	virtual void evaluate() = 0;

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

class BezierCurve2D : public ParameterCurve<Point2f, 2>
{
public:
	explicit BezierCurve2D(const std::vector<Point2f>& points)
		: ParameterCurve<Point2f, 2>(points)
	{
	}

	void evaluate() override
	{
		int count = control_points.size();

		for (int dim_i = 0; dim_i < 2; ++dim_i)
		{
			Float h = 1.0 / (count - 1);
			std::vector<Point2f> points;
			for (int i = 0; i < count; ++i)
			{
				points.emplace_back(i * h, control_points[i][dim_i]);
			}

			BezierApproximation bezier(points);
			function[dim_i] = bezier;
		}
	}
};

class BezierSplineCurve : public ParameterCurve<Point2f, 2>
{
public:
	explicit BezierSplineCurve(const std::vector<Point2f>& points)
		: ParameterCurve<Point2f, 2>(points)
	{
	}

	void evaluate() override
	{
		int count = control_points.size();

		for (int dim_i = 0; dim_i < 2; ++dim_i)
		{
			Float h = 1.0 / (count - 1);
			if (count == 1)h = 0.0;
			std::vector<Point2f> points;
			for (int i = 0; i < count; ++i)
			{
				points.emplace_back(i * h, control_points[i][dim_i]);
			}

			BezierSplineInterpolation bezier(points);
			bezier.evaluate();
			function[dim_i] = bezier;
		}
	}
};

class BSplineCurve : public ParameterCurve<Point2f, 2>
{
public:
	explicit BSplineCurve(const std::vector<Eigen::Matrix<double, 2, 1, 0>>& matrices)
		: ParameterCurve<Eigen::Matrix<double, 2, 1, 0>, 2>(matrices)
	{
	}

	void evaluate() override
	{
		int count = control_points.size();

		for (int dim_i = 0; dim_i < 2; ++dim_i)
		{
			Float h = 1.0 / (count - 1);
			if (count == 1)h = 0.0;
			std::vector<Point2f> points;
			for (int i = 0; i < count; ++i)
			{
				points.emplace_back(i * h, control_points[i][dim_i]);
			}

			BSplineApproximation<3,false> bezier(points);
			bezier.evaluate();
			function[dim_i] = bezier;
		}
	}
};