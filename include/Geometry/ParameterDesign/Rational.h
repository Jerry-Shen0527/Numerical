#pragma once
#include "ParameterCurve.h"

template <int PointCount>
class RationalBSpline : public ParameterCurve<Eigen::Vector2d, PointCount>
{
public:
	explicit RationalBSpline(const std::vector<Point2d>& points, const std::vector<Float>& weights) :
		ParameterCurve<Eigen::Vector2d, PointCount>(points),
		weights(weights)
	{
		if (PointCount != points.size() || PointCount != weights.size())
		{
			throw std::runtime_error("Point size doesn't coincide!");
		}
	}

	void evaluate() override
	{
		bezier_points.clear();
		for (int i = 0; i < PointCount; ++i)
		{
			bezier_points.push_back(Eigen::Vector3d(control_points[i].x(), control_points[i].y(), weights[i]));
		}
		bezier = BezierCurveND<3>(bezier_points);
		bezier.evaluate();
	}

	Eigen::Vector2d operator()(Float t) override
	{
		auto p = bezier(t);
		return Eigen::Vector2d(p.x() / p.z(), p.y() / p.z());
	}

	BezierCurveND<3> bezier;

	std::vector<Float> weights;

	std::vector<Eigen::Vector3d> bezier_points;
};
