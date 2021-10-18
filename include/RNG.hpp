#pragma once
inline Float RandomFloat() { return rand() / (RAND_MAX + 1.0); }
inline Float RandomFloat(Float min, Float max) { return min + (max - min) * RandomFloat(); }