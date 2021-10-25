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
