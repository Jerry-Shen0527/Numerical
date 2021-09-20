#pragma once
#include <type.hpp>
#include <functional>
template<typename T>
class FEM
{
public:
	FEM(std::function<std::vector<>> func);
	~FEM();

private:
};

template <typename T>
FEM<T>::~FEM()
{
}