#include "FEM/FEM.hpp"

void StaticFEM1D::FillMatrix()
{
	std::vector<Eigen::Triplet<Float>> triplets;
	for (int i = 0; i < mat_size; ++i)
	{
		auto related_vec = RelatedFuncIdx(i);

		for (int j : related_vec)
		{
			triplets.emplace_back(i, j, GradientInnerProduct(i, j));
			triplets.emplace_back(i, j, SelfInnerProduct(i, j));
			triplets.emplace_back(i, j, GradientSelfInnerProduct(i, j));
		}
	}
	matrix_ = Eigen::SparseMatrix<Float>(mat_size, mat_size);
	matrix_.setFromTriplets(triplets.begin(), triplets.end());
}

void StaticFEM1D::FillRhs()
{
	rhs = Vector(mat_size);
	for (int i = 0; i < mat_size; ++i)
	{
		rhs(i) = RHSInnerProduct(i);
	}
}