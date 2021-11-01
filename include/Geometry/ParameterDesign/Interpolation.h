#include <type.hpp>
#include <vector>

#include "Bezier.hpp"
#include "Domain.hpp"

enum class Result
{
	NotAvailable = 1,
	Success = 2
};

class Interpolation
{
public:
	Interpolation(const std::vector<Point2>& points) :Points(points) {  }

	virtual Result error_bound(Float& error_b) = 0;

	virtual void evaluate() = 0;
	virtual ~Interpolation() = default;

	virtual  Float operator() (Float in_val) = 0;
protected:
	std::vector<Point2> Points;
};

class RadialInterpolation : public Interpolation
{
public:
	using Interpolation::Interpolation;

	Result error_bound(Float& error_b) override { return Result::NotAvailable; }

	void evaluate() override
	{
		if (!Points.empty())
		{
			if (Points.size() > 1)
			{
				auto iter_max = std::max_element(Points.begin(), Points.end(), [](const Point2& a, const Point2& b) {return a.x() < b.x(); });
				auto iter_min = std::min_element(Points.begin(), Points.end(), [](const Point2& a, const Point2& b) {return a.x() < b.x(); });
				d = (iter_max->x() - iter_min->x()) / Points.size();
				d = d * d;
			}

			Float a = 0;
			for (int i = 0; i < Points.size(); ++i)
			{
				a += Points[i].y();
			}
			a /= Points.size();
			ave = a;
		}

		if (!Points.empty())
		{
			auto size = Points.size();
			matrix = Eigen::MatrixXd(size, size);

			Eigen::VectorXd b(size);
			for (int i = 0; i < size; ++i)
			{
				b[i] = Points[i].y() - ave;
			}

			for (int i = 0; i < size; ++i)
			{
				for (int j = 0; j < size; ++j)
				{
					matrix(i, j) = radial(Points[i].x(), Points[j].x(), d);
				}
			}

			rst = matrix.ldlt().solve(b);
		}
	}

	Float operator()(Float in_val) override
	{
		if (Points.empty())
		{
			return 0;
		}
		else
		{
			Float ret = ave;
			for (int i = 0; i < Points.size(); ++i)
			{
				ret += rst[i] * radial(in_val, Points[i].x(), d);
			}

			return ret;
		}
	}

	Float radial(Float x, Float center, Float d)
	{
		return 1.0 / ((x - center) * (x - center) + d);
	}

	Eigen::MatrixXd matrix;
	Eigen::VectorXd rst;
	Float ave;
	Float d = 1;
};

class LagrangianPolynomial : public Interpolation
{
public:
	using Interpolation::Interpolation;

	Float operator()(Float in_val) override;
	//TODO:: the coefficients can be precalculated.
	void evaluate() override {}
	Result error_bound(Float& error_b) override;
private:
};

inline Float LagrangianPolynomial::operator()(Float in_val)
{
	Float rst = 0;
	for (int i = 0; i < Points.size(); ++i)
	{
		Float temp = Points[i].y();
		if (temp == 0.0)
		{
			continue;
		}
		for (int j = 0; j < Points.size(); ++j)
		{
			if (j != i)
			{
				temp *= (in_val - Points[j].x()) / (Points[i].x() - Points[j].x());
			}
		}
		rst += temp;
	}
	return rst;
}

inline Result LagrangianPolynomial::error_bound(Float& error_b)
{
	//The maximum of the interpolated function's (n+1)th level derivative needs to be known.
	return  Result::NotAvailable;
}

class NewtonPolynomial :public Interpolation
{
public:
	using Interpolation::Interpolation;

	Result error_bound(Float& error_b) override;
	void evaluate() override;
	Float operator()(Float in_val) override;

	void addPoint(const Point2& point)
	{
		Points.push_back(point);
	}

private:
	std::vector<Float> coeff;
};

inline Result NewtonPolynomial::error_bound(Float& error_b)
{
	//The same with lagrangian.
	return Result::NotAvailable;
}

inline void NewtonPolynomial::evaluate()
{
	coeff.clear();
	coeff.resize(Points.size());
	for (int i = 0; i < coeff.size(); ++i)
	{
		coeff[i] = Points[i].y();
	}
	for (int i = 1; i < Points.size(); ++i)
	{
		for (int j = Points.size() - 1; j >= i; --j)
		{
			coeff[j] = (coeff[j] - coeff[j - 1]) / (Points[j].x() - Points[j - i].x());
		}
	}
}

inline Float NewtonPolynomial::operator()(Float in_val)
{
	Float rst = 0;
	for (int i = 0; i < Points.size(); ++i)
	{
		Float temp = coeff[i];
		for (int j = 0; j < i; ++j)
		{
			temp *= (in_val - Points[j].x());
		}
		rst += temp;
	}
	return rst;
}

//Hermite polynomial interpolation can be generalized to a lot of areas.
class HermitePolynomial :public Interpolation
{
public:

	HermitePolynomial(const std::vector<Point2>& points, const std::vector<Point2>& derivatives)
		:Interpolation(points), derivatives(derivatives) {	}
	Result error_bound(Float& error_b) override;
	void evaluate() override;
	Float operator()(Float in_val) override;

private:
	std::vector<Point2> derivatives;
};

class SplineInterpolation :public Interpolation
{
public:
	using Interpolation::Interpolation;

	Result error_bound(Float& error_b) override;
	void evaluate() override;
	Float operator()(Float in_val) override;
};

class BezierSplineInterpolation :public Interpolation
{
public:
	BezierSplineInterpolation(const std::vector<Point2>& points)
		: Interpolation(points)
	{
	}

	Result error_bound(Float& error_b) override
	{
		return Result::NotAvailable;
	}

	void evaluate() override
	{
		if (Points.empty())
		{
			return;
		}
		if (Points.size() == 1)
		{
			rst = Eigen::VectorXd(1);
			rst[0] = Points[0].y();
			return;
		}
		using namespace Eigen;
		int point_count = Points.size();
		int mat_size = point_count * 3 - 2;
		matrix = SparseMatrix<Float>(mat_size, mat_size);

		std::sort(Points.begin(), Points.end(), [](const Point2& a, const Point2& b) {return a.x() < b.x(); });

		Interval interval(Points[0].x(), Points.back().x());
		std::vector< Triplet<Float>> triplets;

		for (int i = 1; i < point_count - 1; ++i)
		{
			triplets.emplace_back(3 * i, 3 * i, 1);

			auto p = Points[i];
			auto x = p.x();
			auto diff_1 = (x - Points[i - 1].x()) / 3.0;
			auto diff_2 = (Points[i + 1].x() - x) / 3.0;

			//C1

			triplets.emplace_back(3 * i - 1, 3 * i, 1.0 / diff_1);
			triplets.emplace_back(3 * i - 1, 3 * i - 1, -1.0 / diff_1);
			//C1
			triplets.emplace_back(3 * i - 1, 3 * i, 1.0 / diff_2);
			triplets.emplace_back(3 * i - 1, 3 * i + 1, -1.0 / diff_2);
			//C2

			triplets.emplace_back(3 * i + 1, 3 * i + 0, 1.0 / diff_1 / diff_1);
			triplets.emplace_back(3 * i + 1, 3 * i - 1, -2.0 / diff_1 / diff_1);
			triplets.emplace_back(3 * i + 1, 3 * i - 2, 1.0 / diff_1 / diff_1);

			//C2

			triplets.emplace_back(3 * i + 1, 3 * i + 0, -1.0 / diff_2 / diff_2);
			triplets.emplace_back(3 * i + 1, 3 * i + 1, 2.0 / diff_2 / diff_2);
			triplets.emplace_back(3 * i + 1, 3 * i + 2, -1.0 / diff_2 / diff_2);
		}

		//C0
		int i = point_count - 1;

		auto p = Points[i];
		auto x = p.x();
		triplets.emplace_back(3 * i, 3 * i, 1);
		auto diff_1 = (x - Points[i - 1].x()) / 3.0;
		//C2 BC
		triplets.emplace_back(3 * i - 1, 3 * i - 0, 1.0 / diff_1 / diff_1);
		triplets.emplace_back(3 * i - 1, 3 * i - 1, -2.0 / diff_1 / diff_1);
		triplets.emplace_back(3 * i - 1, 3 * i - 2, 1.0 / diff_1 / diff_1);
		i = 0;

		p = Points[i];
		x = p.x();
		triplets.emplace_back(3 * i, 3 * i, 1);
		diff_1 = (Points[i + 1].x() - x) / 3.0;
		//C2 BC
		triplets.emplace_back(3 * i + 1, 3 * i + 0, 1.0 / diff_1 / diff_1);
		triplets.emplace_back(3 * i + 1, 3 * i + 1, -2.0 / diff_1 / diff_1);
		triplets.emplace_back(3 * i + 1, 3 * i + 2, 1.0 / diff_1 / diff_1);

		matrix.setFromTriplets(triplets.begin(), triplets.end());

		VectorXd rhs(mat_size);

		rhs.setZero();

		for (int i = 0; i < point_count; ++i)
		{
			rhs[3 * i] = Points[i].y();
		}

		SparseLU<Eigen::SparseMatrix<Float>> solver;
		solver.compute(matrix);

		rst = solver.solve(rhs);
	}

	Float operator()(Float in_val) override
	{
		if (Points.empty())
		{
			return 0;
		}
		for (int i = 0; i < Points.size(); ++i)
		{
			if (Points[i].x() < in_val)
				continue;
			{
				i--;
				if (i < 0)
				{
					return rst[0];
				}
				std::vector<Point2> bezier_points(4);

				auto diff = Points[i + 1].x() - Points[i].x();

				for (int j = 0; j < 4; ++j)
				{
					bezier_points[j] = Point2(Points[i].x() + j * (diff) / 3.0, rst[3 * i + j]);
				}

				BezierApproximation bezier(bezier_points);

				return bezier(in_val);
			}
		}

		return rst[3 * (Points.size() - 1)];
	}

	Eigen::SparseMatrix<Float> matrix;

	Eigen::VectorXd rst;

	std::vector<std::vector<Point2>> parts;
};
