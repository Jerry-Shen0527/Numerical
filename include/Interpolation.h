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
	Interpolation(const std::vector<Point2>& points) :Points(points) {}

	virtual Result error_bound(Float& error_b) = 0;

	virtual void evaluate() = 0;
	virtual ~Interpolation() = default;

	virtual  Float operator() (Float in_val) = 0;
protected:
	std::vector<Point2> Points;
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
		Float temp = Points[i].y;
		for (int j = 0; j < Points.size(); ++j)
		{
			if (j != i)
			{
				temp *= (in_val - Points[j].x) / (Points[i].x - Points[j].x);
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
		coeff[i] = Points[i].y;
	}
	for (int i = 1; i < Points.size(); ++i)
	{
		for (int j = Points.size() - 1; j >= i; --j)
		{
			coeff[j] = (coeff[j] - coeff[j - 1]) / (Points[j].x - Points[j - i].x);
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
			temp *= (in_val - Points[j].x);
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
