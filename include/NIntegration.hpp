#pragma once

template<typename T>
class Domain
{
	virtual bool inside(const T& val);

};


class Integrator
{
public:
	Integrator();
	~Integrator();



private:
};

Integrator::Integrator()
{
}

Integrator::~Integrator()
{
}