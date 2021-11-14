#pragma once

#include"Integrate.hpp"

template<typename T>
class MCIntegrate :public Integrate<T>
{
public:
	MCIntegrate(size_t sample_count, Float precision_bound = -1)
		: sample_count(sample_count),
		precision_bound(precision_bound)
	{
	}

	Float operator()(const std::function<Float(T)>& func, const Domain<T>& domain) override
	{
		Float ret = 0;
		for (size_t i = 0; i < sample_count; i++)
		{
			Float pdf;
			auto sample = domain.RandomSample(pdf);
			ret += sample / pdf;
		}
		ret /= (Float)sample_count;
		return ret;
	}
private:
	size_t sample_count;
	Float precision_bound;
};