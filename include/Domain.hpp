#pragma once
#include "MathFunctions.h"
#include "RNG.hpp"
#include <set>

template <typename T>
class Domain
{
public:
	virtual bool Inside(const T& value) const = 0;

	//Random sample part
	virtual T RandomSample(Float& pdf) const = 0;
};

using Domain1D = Domain<Float>;

class Interval : public Domain1D
{
public:
	Interval(Float left = 0, Float right = 1) : left(left), right(right) { SetPartitionCount(1); }

	void SetPartitionCount(int partition)
	{
		n = partition;
		h = length() / static_cast<Float>(n);
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
		return RandomFloat(left, right);
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

	//Scale a f([left,right]) to f([0,1]), remember to multiply the factor when scaling derivatives
	std::function<Float(Float)> remap(const std::function<Float(Float)>& func, Float factor = 1.0)
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

template <typename T>
class Union : Domain<T>
{
public:
	Union(const Domain<T>& domain1, const Domain<T>& domain2) : domain_1(domain1), domain_2(domain2)
	{
	}

	bool Inside(const T& value) const override
	{
		return domain_1.Inside(value) || domain_2.Inside(value);
	}

	T RandomSample(Float& pdf) const override;

private:
	Domain<T> domain_1;
	Domain<T> domain_2;
};

class RectDomain : public Domain<Point2d>
{
public:
	RectDomain(Point2d p1, Point2d p2)
	{
		x_axis = Interval(p1.x(), p2.x());
		y_axis = Interval(p1.y(), p2.y());
	}

	bool Inside(const Point2d& value) const override
	{
		return x_axis.Inside(value.x()) && y_axis.Inside(value.y());
	}

	Point2d RandomSample(Float& pdf) const override
	{
		Float pdf_x;
		Float pdf_y;
		Float x = x_axis.RandomSample(pdf_x);
		Float y = x_axis.RandomSample(pdf_y);
		pdf = pdf_x * pdf_y;
		return Point2d(x, y);
	}

private:
	Interval x_axis;
	Interval y_axis;
};

class TriangleDomain : public Domain<Eigen::Vector3d>
{
public:
	TriangleDomain() : p0(0, 0, 0), p1(1, 0, 0), p2(0, 1, 0)
	{
	}

	TriangleDomain(Eigen::Vector3d p0, Eigen::Vector3d p1, Eigen::Vector3d p2) : p0(p0), p1(p1), p2(p2)
	{
	}

	bool Inside(const Eigen::Vector3d& value) const override
	{
		//TODO:Judge if a point is in the domain
		return true;
	}

	static Point2d UniformSampleTriangle(const Point2d& u)
	{
		Float su0 = std::sqrt(u[0]);
		return Point2d(1 - su0, u[1] * su0);
	}

	Float Area() const
	{
		Eigen::Vector3d l1 = p1 - p0;
		Eigen::Vector3d l2 = p2 - p0;
		return 0.5 * (l1.cross(l2)).norm();
	}

	Eigen::Vector3d RandomSample(Float& pdf) const override
	{
		Eigen::Vector2d u = rect_domain.RandomSample(pdf);

		Point2d b = UniformSampleTriangle(u);

		Eigen::Vector3d ret = b[0] * p0 + b[1] * p1 + (1 - b[0] - b[1]) * p2;
		pdf /= Area();

		return ret;
	}

	Eigen::Vector3d p0;
	Eigen::Vector3d p1;
	Eigen::Vector3d p2;

	/**
	 * \brief Remap barycentric space function to take world coordinates as inputs
	 * \param func A standard function defined with barycentric coordinates
	 * \param factor
	 * \return A new function that takes world space coordinates as inputs
	 */
	std::function<Float(Eigen::Vector3d)> remap(const std::function<Float(Eigen::Vector2d)>& func, Float factor = 1.0)
	{
		return[factor, func, *this](Eigen::Vector3d coord)
		{
			return factor * func(Barycentric(coord)) * factor;
		};
	}

	/**
	 * \brief Remap world space function to take barycentric as inputs
	 * \param func A standard function defined with barycentric coordinates
	 * \param factor
	 * \return A new function that takes barycentric coordinates as inputs
	 */
	std::function<Float(Eigen::Vector2d)> scale(const std::function<Float(Eigen::Vector3d)>& func, Float factor = 1.0)
	{
		return[factor, func, *this](Eigen::Vector2d bary)
		{
			return factor * func(Coordinates(bary)) * factor;
		};
	}

	//convert the Barycentric to coordinates.
	Eigen::Vector3d Coordinates(Eigen::Vector2d barycentric) const
	{
		Float alpha = barycentric.x();
		Float beta = barycentric.y();

		Eigen::Vector3d ret = alpha * p0 + beta * p1 + (1 - alpha - beta) * p2;
		return ret;
	}

	//convert the coordinates to Barycentric
	Eigen::Vector2d Barycentric(Eigen::Vector3d coord) const
	{
		Eigen::Vector3d normal = (p2 - p0).cross(p1 - p0).normalized();
		coord = coord - Eigen::Vector3d(coord - p0).dot(normal) * normal;
		//std::cerr << (coord - p0).dot((p2 - p0).cross(p1 - p0))<< std::endl;
		assert(abs((coord - p0).dot((p2 - p0).cross(p1 - p0))) < 1E-7);

		Float i = (-(coord.x() - p1.x()) * (p2.y() - p1.y()) + (coord.y() - p1.y()) * (p2.x() - p1.x())) /
			(-(p0.x() - p1.x()) * (p2.y() - p1.y()) + (p0.y() - p1.y()) * (p2.x() - p1.x()));
		Float j = (-(coord.x() - p2.x()) * (p0.y() - p2.y()) + (coord.y() - p2.y()) * (p0.x() - p2.x())) /
			(-(p1.x() - p2.x()) * (p0.y() - p2.y()) + (p1.y() - p2.y()) * (p0.x() - p2.x()));

		return Eigen::Vector2d(i, j);
	}

	RectDomain rect_domain = RectDomain(Point2d(0, 0), Point2d(1, 1));
};

class SphereDomain : public Domain<Point2d>
{
public:
	bool Inside(const Point2d& value) const override
	{
		return true;
	}

	/**
	 * \brief Uniform sampling in a sphere
	 * \param pdf
	 * \return (theta,phi)
	 */
	Point2d RandomSample(Float& pdf) const override
	{
		Eigen::Vector3f vec = UniformSampleSphere(rect_domain.RandomSample(pdf));
		pdf *= 1. / 4.0 / Pi;

		return Point2d(SphericalTheta(vec), SphericalPhi(vec));
	}

	static Eigen::Vector3f UniformSampleSphere(const Point2d& u)
	{
		Float z = 1 - 2 * u[0];
		Float r = std::sqrt(std::max(static_cast<Float>(0), static_cast<Float>(1) - z * z));
		Float phi = 2 * Pi * u[1];
		return Eigen::Vector3f(r * std::cos(phi), r * std::sin(phi), z);
	}

	static Float SphericalTheta(const Eigen::Vector3f& v)
	{
		return std::acos(Clamp(v.z(), -1, 1));
	}

	static Float SphericalPhi(const Eigen::Vector3f& v)
	{
		Float p = std::atan2(v.y(), v.x());
		return (p < 0) ? (p + 2 * Pi) : p;
	}

	//For saving sample parameters for Rectdomain
	RectDomain rect_domain = RectDomain(Point2d(0, 0), Point2d(1, 1));
};