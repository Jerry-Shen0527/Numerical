#pragma once
#include <type.hpp>
#include <functional>
#include "Eigen/Eigen"
class StaticFEM
{
public:
	virtual ~StaticFEM() {}

	virtual void evaluate()
	{
		FillMatrix();

		FillRhs();

		Eigen::SimplicialLLT<Eigen::SparseMatrix<Float>> solver;

		solver.compute(matrix_);
		rst = solver.solve(rhs);
	}

	Eigen::VectorXd& coeff_() { return rst; }

protected:
	//This is a premature design of interfaces, which might be refactored in large scale
	virtual void FillMatrix() = 0;
	virtual void FillRhs() = 0;

	virtual Float GradientInnerProduct(int i, int j) = 0;
	virtual Float RHSInnerProduct(int i) = 0;
	virtual std::vector<int> RelatedFuncIdx(int idx) = 0;

	size_t size;

	Eigen::SparseMatrix<Float> matrix_;
	Eigen::VectorXd rhs;
	Eigen::VectorXd rst;
};

class StaticFEM1D :public StaticFEM
{
protected:
	void FillMatrix() override;
	void FillRhs() override;
public:
	std::function<Float(Float)> Function;
	virtual Float Value(Float x) = 0;
};