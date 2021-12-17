#pragma once
#include <chrono>
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

	virtual void SolveMatrix()
	{
		Eigen::SparseLU<Eigen::SparseMatrix<Float>> solver;

		solver.compute(matrix_);
		rst = solver.solve(rhs);
	}

	virtual void evaluate()
	{
		SetMatSize();
		FillMatrix();
		FillRhs();

		SolveMatrix();
	}

	Vector& coeff_() { return rst; }

	void SetRhs(const Vector& rhs)
	{
		if (rhs.rows() == this->rhs.rows())
		{
			this->rhs = rhs;
			return;
		}
		throw std::runtime_error("RHS size wrong!");
	}

	Vector& GetRhs()
	{
		return  rhs;
	}

	void SetRst(const Vector& rst)
	{
		if (this->rst.rows() == 0 || rst.rows() == this->rst.rows())
		{
			this->rst = rst;
			return;
		}
		throw std::runtime_error("RST size wrong!");
	}

	Vector& GetRst()
	{
		return  rst;
	}

	Vector GetResidual()
	{
		if (rst.rows() == 0)
		{
			rst = Vector(matrix_.rows());
			rst.setZero();
		}
		Vector ret = rhs - matrix_ * rst;
		return ret;
	}

	//This is a premature design of interfaces, which might be refactored in large scale
	virtual void SetMatSize() = 0;
	virtual void FillMatrix() = 0;
	virtual void FillRhs() = 0;
protected:

	virtual Float GradientInnerProduct(int i, int j) = 0;
	virtual Float SelfInnerProduct(int i, int j) = 0;
	virtual Float GradientSelfInnerProduct(int i, int j) = 0;
	virtual Float RHSInnerProduct(int i) = 0;
	virtual std::vector<int> RelatedFuncIdx(int idx) = 0;

	size_t mat_size = 0;

	Eigen::SparseMatrix<Float> matrix_;
	Vector rhs;
	Vector rst;
};

class StaticFEMwithCG :public StaticFEM
{
public:
	bool DoConjugateGradient(bool full = true);

public:
	inline static Float time_last_cg = std::numeric_limits<Float>::infinity();

	void SolveMatrix() override
	{
		using namespace std::chrono;
		using std::cout;

		auto start = system_clock::now();

		int accu = 0;
		while (true)
		{
			if (DoConjugateGradient())
			{
				accu += rhs.rows() * rhs.rows();
				break;
			}
		}
		// do something...
		auto end = system_clock::now();
		auto duration = duration_cast<nanoseconds>(end - start).count();

		duration = accu;
		if (rhs.rows() > 10)
		{
			cout << double(duration) << "&" << log2(double(duration) / time_last_cg) << "\\\\" << std::endl;
			time_last_cg = double(duration);
		}
	}
};

class StaticFEM1D :public StaticFEMwithCG
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