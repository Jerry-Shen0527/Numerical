#pragma once
#include "type.hpp"

template<typename T>
class Domain
{
	virtual bool Inside(const T& value) const = 0;
	virtual T RandomSample() { throw std::runtime_error("RandomSample Not available."); }

	void setSamples(int sample)
	{
		samples = sample;
	}
	virtual T UniformSample(int idx) = 0;

protected:

	int samples = 0;
};

template<typename T>
class Integrate
{
public:
	Integrate();
	~Integrate();

	virtual Float operator() (std::function<Float(T)> func, const Domain<T>& domain) = 0;

protected:
	/**
	 * \brief Controls the upper bound of error
	 */
	Float precision;
};

template<typename T>
Integrate<T>::Integrate()
{
}

template<typename T>
Integrate<T>::~Integrate()
{
}

class Domain1D :public Domain<Float>
{
public:
};

class TrapeziumIntegrate
{
public:
};