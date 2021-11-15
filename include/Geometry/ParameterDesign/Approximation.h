#pragma once
#include "Domain.hpp"

class Approximation
{
public:
	Approximation(const std::vector<Point2f>& points) :Points(points) {  }

	virtual void evaluate() = 0;
	virtual ~Approximation() = default;

	virtual  Float operator() (Float in_val) = 0;
protected:
	std::vector<Point2f> Points;
};

//This is definitely a bad idea
template<int k, bool Interpolate_EndPoints = true>
class BSplineApproximation :public Approximation
{
public:
	explicit BSplineApproximation(const std::vector<Point2f>& matrices)
		: Approximation(matrices)
	{
	}

	void evaluate() override
	{
		if (Points.size() <= k)
		{
			return;
		}

		std::sort(Points.begin(), Points.end(), [](const Point2f& a, const Point2f& b) {return a.x() < b.x(); });

		Interval interval(Points.front().x(), Points.back().x());

		if (Interpolate_EndPoints)
		{
			knot_vector.resize(Points.size() - k + 1);

			for (int i = 0; i < Points.size() - k + 1; ++i)
			{
				knot_vector[i] = interval.lerp(Float(i) / (Points.size() - k));
			}
			for (int i = 0; i < k - 1; ++i)
			{
				knot_vector.insert(knot_vector.begin(), knot_vector.front());
				knot_vector.push_back(knot_vector.back());
			}
		}
		else
		{
			knot_vector.resize(Points.size() + k - 1);
			for (int i = 0; i < Points.size() + k - 1; ++i)
			{
				knot_vector[i] = interval.lerp(Float(i) / (Points.size() + k - 2));
			}
		}
	}

	Float operator()(Float t) override
	{
		if (Points.size() <= k)
		{
			return 0;
		}

		if (t <= knot_vector[k - 1])
			t = knot_vector[k - 1] + std::numeric_limits<Float>::epsilon();
		if (t >= knot_vector[knot_vector.size() - k ])
			t = knot_vector[knot_vector.size() - k ];

		auto iter = std::find_if(knot_vector.begin(), knot_vector.end(), [t](Float knot) {return  knot >= t; });

		int r = std::distance(knot_vector.begin(), iter);

		std::vector<Point2f> tower(std::next(Points.begin(), r - k), std::next(Points.begin(), r + 1));

		for (int j = 1; j <= k; ++j)
		{
			for (int i = r; i >= r - k + j; --i)
			{
				int index = i - (r - k);
				float alpha = (t - knot_vector[i-1]) / (knot_vector[i + k - j] - knot_vector[i-1]);
				tower[index] = Lerp(alpha, tower[index - 1], tower[index]);
			}
		}
		return tower.back().y();
	}

	std::vector<Float> knot_vector;
};