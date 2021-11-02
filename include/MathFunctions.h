#pragma once
#include "type.hpp"

template<typename T>
inline T Lerp(Float t, T t1, T t2) { return (1 - t) * t1 + t * t2; }

inline Float Clamp(Float t, Float t1 = 0, Float t2 = 1) {
	if (t > t2)
	{
		return t2;
	}
	if (t < t1)
	{
		return t1;
	}
	return t;
}
