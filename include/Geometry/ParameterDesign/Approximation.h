#pragma once

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
		if (Points.empty())
		{
			return;
		}

		std::sort(Points.begin(), Points.end(), [](const Point2& a, const Point2& b) {return a.x() < b.x(); });
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

		for (int i = 0; i < k - 1; ++i)
		{
			Points.insert(Points.begin(), Points.front());
			Points.push_back(Points.back());
		}

		auto iter = std::find_if(Points.begin(), Points.end(), [t](const Point2& point) {return  point.x() > t; });

		int r = std::distance(Points.begin(), iter) - 1;

		std::vector<Point2> tower(std::next(Points.begin(), r - k + 1), std::next(Points.begin(), r + 1));

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
};