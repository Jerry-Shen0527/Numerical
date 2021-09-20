#include <type.hpp>
#include <vector>

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
				auto iter_max = std::max_element(Points.begin(), Points.end(), [](const Point2& a, const Point2& b) {return a.x() > b.x(); });
				auto iter_min = std::min_element(Points.begin(), Points.end(), [](const Point2& a, const Point2& b) {return a.x() > b.x(); });
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
