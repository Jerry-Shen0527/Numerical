#pragma once
#include <type.hpp>
#include <functional>
#include "Eigen/Eigen"
class FEM
{
public:
	virtual ~FEM() {}

	virtual void evaluate()
	{
		FillMatrix();
		FillRhs();

		Eigen::SparseLU<Eigen::SparseMatrix<Float>> solver;

		solver.compute(matrix_);
		rst = solver.solve(rhs);
	}

	Eigen::VectorXd coeff_() { return rst; }

protected:
	//This is a premature design of interfaces, which might be refactored in large scale
	virtual void FillMatrix() = 0;
	virtual void FillRhs() = 0;

	Eigen::SparseMatrix<Float> matrix_;
	Eigen::VectorXd rhs;
	Eigen::VectorXd rst;
};

class StaticFEM1D	 :public FEM
{
protected:
	void FillMatrix() override;
	void FillRhs() override;
public:

};