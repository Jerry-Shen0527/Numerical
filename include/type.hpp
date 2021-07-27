#include <glm/glm.hpp>

#include <iostream>

template<typename T, int N>
std::ostream& operator<<(std::ostream& out, glm::vec<N, T>& vec)
{
	for (int i = 0; i < N; ++i)
	{
		out << vec[i] << ' ';
	}
	return out;
}

using Float = double;

using Point2 = glm::vec<2, Float>;