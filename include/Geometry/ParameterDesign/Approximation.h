#pragma once
#include "Domain.hpp"

class Approximation
{
public:
	Approximation(const std::vector<Point2>& points) :Points(points) {  }

	virtual void evaluate() = 0;
	virtual ~Approximation() = default;

	virtual  Float operator() (Float in_val) = 0;
protected:
	std::vector<Point2> Points;
};

template<int k>
class BSplineApproximation :public Approximation
{
public:
	explicit BSplineApproximation(const std::vector<Point2>& matrices)
		: Approximation(matrices)
	{
	}

	void evaluate() override
	{
		if (Points.size() <= 1)
		{
			return;
		}

		std::sort(Points.begin(), Points.end(), [](const Point2& a, const Point2& b) {return a.x() < b.x(); });

		Interval interval(Points.front().x(), Points.back().x());

		knot_vector.resize(Points.size());

		for (int i = 0; i < Points.size(); ++i)
		{
			knot_vector[i] = interval.lerp(Float(i) / (Points.size() - 1));
		}

		for (int i = 0; i < k - 1; ++i)
		{
			knot_vector.insert(knot_vector.begin(), knot_vector.front());
			knot_vector.push_back(knot_vector.back());
		}
	}

	Float operator()(Float t) override
	{
		if (Points.empty())
		{
			return 0;
		}

		if (t < Points.begin()->x())
		{
			return Points.begin()->y();
		}
		if (t > Points.back().x())
		{
			return Points.back().y();
		}

		auto iter = std::find_if(knot_vector.begin(), knot_vector.end(), [t](Float knot) {return  knot > t; });

		int r = std::distance(knot_vector.begin(), iter) - 1;

		std::vector<Float> tower(std::next(knot_vector.begin(), r - k + 1), std::next(knot_vector.begin(), r + 1));

		for (int j = 1; j < k - 1; ++j)
		{
			for (int i = r - k + 1 + j; i <= r; ++i)
			{
				int index = r - k + 1 + i;
				float alpha = (t - Points[i].x()) / (Points[i + k - j].x() - Points[i].x());
			}
		}

		return 0;
	}

	std::vector<Float> knot_vector;
};