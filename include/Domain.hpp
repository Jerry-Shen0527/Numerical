#pragma once
#include "MathFunctions.h"
#include "RNG.hpp"
#include <set>

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
	Interval(Float left = 0, Float right = 1) : left(left), right(right) { SetPartitionCount(1); }

	void SetPartitionCount(int partition)
	{
		n = partition;
		h = length() / Float(n);
		knot_points.clear();
	}

	int GetPartitionCount() const
	{
		if (!knot_points.empty())
			return knot_points.size() + 1;
		return n;
	}

	bool Inside(const Float& value) const override
	{
		return value >= left && value < right;
	}

	void SetSubIntervalKnots(const std::vector<Float>& vector)
	{
		//Remove duplicated
		auto set = std::set<Float>(vector.begin(), vector.end());
		knot_points = std::vector<Float>(set.begin(), set.end());
		std::sort(knot_points.begin(), knot_points.end());
	}

	Float RandomSample(Float& pdf) const override
	{
		pdf = 1.0 / length();
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

	Interval SubInterval(int idx)
	{
		if (knot_points.empty())
		{
			return Interval(left + idx * h, left + (idx + 1) * h);
		}
		else
		{
			if (idx == 0)
			{
				return Interval(left, knot_points[0]);
			}
			if (idx == knot_points.size())
			{
				return Interval(knot_points.back(), right);
			}

			return Interval(knot_points[idx - 1], knot_points[idx]);
		}
	}

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

	std::vector<Float> knot_points;
};

template<typename T>
class Union :Domain<T>
{
public:
	Union(const Domain<T>& domain1, const Domain<T>& domain2) :domain_1(domain1), domain_2(domain2) {}

	bool Inside(const T& value) const override
	{
		return domain_1.Inside(value) || domain_2.Inside(value);
	}

	T RandomSample(Float& pdf) const override;

private:
	Domain<T> domain_1;
	Domain<T> domain_2;
};

class RectDomain :public Domain<Point2f>
{
public:
	RectDomain(Point2f p1, Point2f p2)
	{
		x_axis = Interval(p1.x(), p2.x());
		y_axis = Interval(p1.y(), p2.y());
	}
	bool Inside(const Point2f& value) const override
	{
		return x_axis.Inside(value.x()) && y_axis.Inside(value.y());
	}

	Point2f RandomSample(Float& pdf) const override
	{
		Float pdf_x;
		Float pdf_y;
		Float  x = x_axis.RandomSample(pdf_x);
		Float  y = x_axis.RandomSample(pdf_y);
		pdf = pdf_x * pdf_y;
		return Point2f(x, y);
	}
private:

	Interval x_axis;
	Interval y_axis;
};

class SphereDomain :public Domain<Point2f>
{
public:

	bool Inside(const Point2f& value) const override
	{
		return true;
	}

	/**
	 * \brief Uniform sampling in a sphere
	 * \param pdf
	 * \return (theta,phi)
	 */
	Point2f RandomSample(Float& pdf) const override
	{
		Eigen::Vector3f vec = UniformSampleSphere(rect_domain.RandomSample(pdf));
		pdf *= 1. / 4.0 / Pi;

		return Point2f(SphericalTheta(vec), SphericalPhi(vec));
	}
	static Eigen::Vector3f UniformSampleSphere(const Point2f& u) {
		Float z = 1 - 2 * u[0];
		Float r = std::sqrt(std::max((Float)0, (Float)1 - z * z));
		Float phi = 2 * Pi * u[1];
		return Eigen::Vector3f(r * std::cos(phi), r * std::sin(phi), z);
	}

	static Float SphericalTheta(const Eigen::Vector3f& v) {
		return std::acos(Clamp(v.z(), -1, 1));
	}

	static Float SphericalPhi(const Eigen::Vector3f& v) {
		Float p = std::atan2(v.y(), v.x());
		return (p < 0) ? (p + 2 * Pi) : p;
	}
	//For saving sample parameters for Rectdomain
	RectDomain rect_domain = RectDomain(Point2f(0, 0), Point2f(1, 1));
};