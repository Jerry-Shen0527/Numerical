#include "FEM/FEM.hpp"

void StaticFEM1D::FillMatrix()
{
	std::vector<Eigen::Triplet<Float>> triplets;
	for (int i = 0; i < size; ++i)
	{
		auto related_vec = RelatedFuncIdx(i);

		for (int j : related_vec)
		{
			triplets.emplace_back(i, j, GradientInnerProduct(i, j));
		}
	}
	matrix_ = Eigen::SparseMatrix<Float>(size, size);
	matrix_.setFromTriplets(triplets.begin(), triplets.end());
}

void StaticFEM1D::FillRhs()
{
	rhs = Eigen::VectorXd(size);
	for (int i = 0; i < size; ++i)
	{
		rhs(i) = RHSInnerProduct(i);
	}
}