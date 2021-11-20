#pragma once
#include <type.hpp>
#include <functional>
#include "Eigen/Eigen"
class StaticFEM
{
public:
	explicit StaticFEM()
	{
	}

	virtual ~StaticFEM() {}

	virtual void evaluate()
	{
		SetMatSize();
		FillMatrix();
		FillRhs();

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Float>> solver;

		solver.compute(matrix_);
		rst = solver.solve(rhs);
	}

	Vector& coeff_() { return rst; }

protected:
	//This is a premature design of interfaces, which might be refactored in large scale
	virtual void SetMatSize() = 0;
	virtual void FillMatrix() = 0;
	virtual void FillRhs() = 0;

	virtual Float GradientInnerProduct(int i, int j) = 0;
	virtual Float SelfInnerProduct(int i, int j) = 0;
	virtual Float GradientSelfInnerProduct(int i, int j) = 0;
	virtual Float RHSInnerProduct(int i) = 0;
	virtual std::vector<int> RelatedFuncIdx(int idx) = 0;

	size_t mat_size;

	Eigen::SparseMatrix<Float> matrix_;
	Vector rhs;
	Vector rst;
};

class StaticFEM1D :public StaticFEM
{
protected:
	StaticFEM1D()
	{
	}

	void FillMatrix() override;
	void FillRhs() override;
public:
	virtual Float Value(Float x) = 0;
};