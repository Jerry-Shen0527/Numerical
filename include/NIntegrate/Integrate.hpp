#pragma once
#include "type.hpp"
#include "Domain.hpp"
#include <RNG.hpp>
#include <MathFunctions.h>

template <typename T>
class Integrate
{
public:
	Integrate();
	virtual ~Integrate();

	virtual Float operator()(const std::function<Float(T)>& func, const Domain<T>& domain) = 0;

protected:
	/**
	 * \brief Controls the upper bound of error
	 */
	Float precision;
};

template <typename T>
Integrate<T>::Integrate()
{
}

template <typename T>
Integrate<T>::~Integrate()
{
}

class SubdivIntegrate : public Integrate<Float>
{
public:
	Float operator()(const std::function<Float(Float)>& func, const Domain<Float>& domain) override
	{
		auto interval = static_cast<const Interval&>(domain);
		Float ret = 0;

		for (int i = 0; i < interval.GetPartitionCount(); ++i)
		{
			ret += integrate(func, interval.SubInterval(i));
		}

		return ret;
	}

	virtual Float integrate(const std::function<Float(Float)>& func, const Interval& interval) const = 0;
};

class TrapeziumIntegrate : public SubdivIntegrate
{
private:
	Float integrate(const std::function<Float(Float)>& func, const Interval& interval) const override
	{
		return (func(interval.lerp(0.0)) + func(interval.lerp(1.0))) * interval.length() / 2.0;
	}
};

class SimpsonIntegrate : public SubdivIntegrate
{
	Float integrate(const std::function<Float(Float)>& func, const Interval& interval) const override
	{
		return (func(interval.lerp(0.0)) + 4.0 * func(interval.lerp(0.5)) + func(interval.lerp(1.0))) * interval.
			length() / 6.0;
	}
};

class GaussIntegrate : public SubdivIntegrate
{
public:
	GaussIntegrate()
	{
	}

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
		return ret * interval.length() / 2.0;
	}
};

inline Float L2InnerProduct(const std::function<Float(Float)>& func1, const std::function<Float(Float)>& func2,
	const Interval& interval)
{
	auto func = [&](Float val) { return func1(val) * func2(val); };
	return GaussIntegrate()(func, interval);
}

inline Float WeightedL2InnerProduct(const std::function<Float(Float)>& func1, const std::function<Float(Float)>& func2,
	const std::function<Float(Float)>& weight, const Interval& interval)
{
	auto func = [&](Float val) { return func1(val) * func2(val) * weight(val); };
	return GaussIntegrate()(func, interval);
}

class GaussIntegrate2D : public Integrate<Eigen::Vector3d>
{
public:
	Float operator()(const std::function<Float(Eigen::Vector3d)>& func,
		const Domain<Eigen::Vector3d>& domain) override
	{
		//TODO: Restrict this integration to triangles
		auto triangle = static_cast<const TriangleDomain&>(domain);

		auto bary_func = triangle.scale(func);

		Float ret = 0;

		for (int i = 0; i < 7; ++i)
		{
			ret += omegas[i] * bary_func(Eigen::Vector2d(gauss_x[i], gauss_y[i]));
		}
		ret *= triangle.Area();
		return ret;
	}

private:
	constexpr static Float omegas[7] = {
		0.2250000000000002, 0.1323941527885061, 0.1323941527885061, 0.1323941527885061, 0.1259391805448272,
		0.1259391805448272, 0.1259391805448272
	};
	constexpr static Float gauss_x[7] = {
		0.3333333333333333, 0.0597158717897698, 0.4701420641051151, 0.4701420641051151, 0.7974269853530873,
		0.1012865073234563, 0.1012865073234563
	};
	constexpr static Float gauss_y[7] = {
		0.3333333333333333, 0.4701420641051151, 0.4701420641051151, 0.0597158717897698, 0.1012865073234563,
		0.1012865073234563, 0.7974269853530873
	};
	constexpr static Float gauss_z[7] = {
		0.3333333333333333, 0.4701420641051151, 0.0597158717897698, 0.4701420641051151, 0.1012865073234563,
		0.7974269853530873, 0.1012865073234563
	};
};


inline Float L2InnerProduct2D(
	const std::function
	<Float(Eigen::Matrix<Float, 3, 1>)>& func1,
	const std::function
	<Float(Eigen::Matrix<Float, 3, 1>)>& func2,
	const TriangleDomain& triangle)
{
	auto func = [&](Eigen::Matrix<Float, 3, 1> pos) { return func1(pos) * func2(pos); };
	return GaussIntegrate2D()(func, triangle);
}

template<int dim>
inline Float L2InnerProduct2D(
	const std::function
	<Eigen::Matrix<Float, dim, 1>(Eigen::Matrix<Float, 3, 1>)>& func1,
	const std::function
	<Eigen::Matrix<Float, dim, 1>(Eigen::Matrix<Float, 3, 1>)>& func2,
	const TriangleDomain& triangle)
{
	auto func = [&](Eigen::Matrix<Float, 3, 1> pos) { return func1(pos).dot(func2(pos)); };
	return GaussIntegrate2D()(func, triangle);
}

inline Float WeightedL2InnerProduct2D(
	const std::function
	<Float(Eigen::Matrix<Float, 3, 1>)>& func1,
	const std::function
	<Float(Eigen::Matrix<Float, 3, 1>)>& func2,
	const std::function
	<Float(Eigen::Matrix<Float, 3, 1>)>& weight,
	const TriangleDomain& triangle)
{
	auto func = [&](Eigen::Matrix<Float, 3, 1> pos) { return func1(pos) * func2(pos) * weight(pos); };
	return GaussIntegrate2D()(func, triangle);
}

template<int dim>
inline Float WeightedL2InnerProduct2D(
	const std::function
	<Eigen::Matrix<Float, dim, 1>(Eigen::Matrix<Float, 3, 1>)>& func1,
	const std::function
	<Eigen::Matrix<Float, dim, 1>(Eigen::Matrix<Float, 3, 1>)>& func2,
	const std::function
	<Float(Eigen::Matrix<Float, 3, 1>)>& weight,
	const TriangleDomain& triangle)
{
	auto func = [&](Eigen::Matrix<Float, 3, 1> pos) { return func1(pos).dot(func2(pos)) * weight(pos); };
	return GaussIntegrate2D()(func, triangle);
}

template<int dim>
class NDifferential
{
	using Vector = Eigen::Matrix<Float, dim, 1>;
public:
	std::function<Vector(Vector)> operator()(std::function<Float(Vector)> function)
	{
		return [function, this](Vector pos)
		{
			Float h = 0.01;

			Vector temp_gradient;
			temp_gradient.setZero();
			for (int i = 0; i < dim; ++i)
			{
				Vector h_vec;
				h_vec.setZero();
				h_vec[i] = h;
				temp_gradient[i] = 1.0 / 12 / h * (function(pos - 2 * h_vec) - 8 * function(pos - h_vec) + 8 * function(pos + h_vec) - function(pos + 2 * h_vec));
			}

			Vector another_temp_gradient;
			another_temp_gradient.setZero();

			do
			{
				h *= 0.5;
				temp_gradient = another_temp_gradient;
				for (int i = 0; i < dim; ++i)
				{
					Vector h_vec;
					h_vec.setZero();
					h_vec[i] = h;
					another_temp_gradient[i] = 1.0 / 12 / h * (function(pos - 2 * h_vec) - 8 * function(pos - h_vec) + 8 * function(pos + h_vec) - function(pos + 2 * h_vec));
				}
			} while ((another_temp_gradient - temp_gradient).norm() > bound);

			return another_temp_gradient;
		};
	};

private:
	Float bound = 1E-8;
};