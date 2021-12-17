#include "FEM/FEM.hpp"

bool StaticFEMwithCG::DoConjugateGradient(bool full)
{
	bool solved = false;

	int size = rhs.rows();
	if (rst.rows() != size)
	{
		rst = Vector(size);
		rst.setZero();
	}

	Vector r = (rhs - matrix_ * rst);
	Vector p = r;

	for (int i = 0; i < (full ? size : 2); ++i)
	{
		Float r_sqr = r.dot(r);

		//std::cout << r_sqr << std::endl;
		//if (r_sqr>1E8)
		//{
		//	system("pause");
		//}
		if (sqrt(r_sqr) < 1E-9)
		{
			solved = true;
			break;
		}
		Float alpha = r_sqr / p.dot(matrix_ * p);
		rst += alpha * p;
		r -= alpha * matrix_ * p;

		Float beta = r.dot(r) / r_sqr;
		p = r + beta * p;
	}
	return solved;
}

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