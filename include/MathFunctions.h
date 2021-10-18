#pragma once
#include "type.hpp"

template<typename T>
inline T Lerp(Float t, T t1, T t2) { return t * t1 + (1 - t) * t2; }