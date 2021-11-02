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

//This is definitely a bad idea
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
		if (Points.size() < k)
		{
			return;
		}

		std::sort(Points.begin(), Points.end(), [](const Point2& a, const Point2& b) {return a.x() < b.x(); });

		Interval interval(Points.front().x(), Points.back().x());

		knot_vector.resize(Points.size() - k + 2);

		for (int i = 0; i < Points.size() - k + 2; ++i)
		{
			knot_vector[i] = interval.lerp(Float(i) / (Points.size() - k + 1));
		}

		for (int i = 0; i < k - 1; ++i)
		{
			knot_vector.insert(knot_vector.begin(), knot_vector.front());
			knot_vector.push_back(knot_vector.back());
		}
	}

	Float operator()(Float t) override
	{
		if (Points.size() < k)
		{
			return 0;
		}

		if (t <= Points.begin()->x())
		{
			return Points.front().y();
		}
		if (t >= Points.back().x())
		{
			return Points.back().y();
		}

		auto iter = std::find_if(knot_vector.begin(), knot_vector.end(), [t](Float knot) {return  knot > t; });

		int r = std::distance(knot_vector.begin(), iter) - 1;

		std::vector<Point2> tower(std::next(Points.begin(), r - k + 1), std::next(Points.begin(), r + 1));

		for (int j = 1; j <= k - 1; ++j)
		{
			for (int i = r; i >= r - k + 1 + j; --i)
			{
				int index = i - (r - k + 1);
				float alpha = (t - knot_vector[i]) / (knot_vector[i + k - j] - knot_vector[i]);
				tower[index] = Lerp(alpha, tower[index - 1], tower[index]);
			}
		}
		return tower.back().y();
	}

	std::vector<Float> knot_vector;
};