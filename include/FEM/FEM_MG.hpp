#pragma once
#include "FEM.hpp"

#include <chrono>

class StaticFEM_Multigrid :public StaticFEM
{
protected:
	void SetMatSize() override
	{
		for (auto& solver : solvers)
		{
			solver->SetMatSize();
		}
	}

	void FillMatrix() override
	{
		for (auto& solver : solvers)
		{
			solver->FillMatrix();
		}
	}

	void FillRhs() override
	{
		for (auto& solver : solvers)
		{
			solver->FillRhs();
		}
	}

	Float GradientInnerProduct(int i, int j) override
	{
		throw std::runtime_error("Should not be called!");
	}
	Float SelfInnerProduct(int i, int j) override
	{
		throw std::runtime_error("Should not be called!");
	}
	Float GradientSelfInnerProduct(int i, int j) override
	{
		throw std::runtime_error("Should not be called!");
	}
	Float RHSInnerProduct(int i) override
	{
		throw std::runtime_error("Should not be called!");
	}
	std::vector<int> RelatedFuncIdx(int idx) override
	{
		throw std::runtime_error("Should not be called!");
	}

public:

	inline static Float time_last = std::numeric_limits<Float>::infinity();
	void evaluate() override
	{
		SetMatSize();
		FillMatrix();
		FillRhs();

		using namespace std::chrono;
		using std::cout;

		auto start = system_clock::now();
		MG(solvers);

		Vector	old_rst = solvers.back()->GetRst();

		int i = 0;

		while (true)
		{
			MG(solvers);
			rst = solvers.back()->GetRst();
			if ((old_rst - rst).norm() < 1E-8)
			{
				break;
			}
			old_rst = rst;
			i++;
		}

		//cout << "MG iteration times: " << i << std::endl;

		// do something...
		rst = solvers.back()->GetRst();

		auto end = system_clock::now();
		auto duration = duration_cast<nanoseconds>(end - start).count();

		duration = smooth_count;
		//cout << "MG takes " << double(duration) << " nanoseconds. " << std::endl;
		cout << rst.rows() + 1 << "&" << double(duration) << "&" << log2(double(duration) / time_last) << "&";
		time_last = double(duration);

		smooth_count = 0;
		//cout <<  double(duration) << ", ";
	}

	void MG(std::vector<std::shared_ptr<StaticFEMwithCG>> solvers)
	{
		//assert(solvers.size() > 1);
		if (solvers.size() == 1)
		{
			solvers[0]->SolveMatrix();
		}
		else
		{
			for (int i = 0; i < preSmoothPass; ++i)
			{
				Smooth(solvers.back());
			}

			for (int i = 0; i < MGPass; ++i)
			{
				std::vector sub_solvers(solvers.begin(), solvers.end() - 1);
				ProjectDown(solvers.back(), sub_solvers.back());
				MG(sub_solvers);
				ProjectUp(sub_solvers.back(), solvers.back());
			}

			for (int i = 0; i < postSmoothPass; ++i)
			{
				Smooth(solvers.back());
			}
		}
	}

	void Smooth(std::shared_ptr<StaticFEMwithCG> solver)
	{
		solver->DoConjugateGradient(false);
		smooth_count += solver->GetRst().rows();
	}
	inline static int smooth_count = 0;
protected:
	virtual void ProjectDown(std::shared_ptr<StaticFEMwithCG> up, std::shared_ptr<StaticFEMwithCG> down) = 0;
	virtual void ProjectUp(std::shared_ptr<StaticFEMwithCG> down, std::shared_ptr<StaticFEMwithCG> up) = 0;

private:
	const int preSmoothPass = 1;

	const int MGPass = 2;

	const int postSmoothPass = 1;

protected:
	std::vector<std::shared_ptr<StaticFEMwithCG>> solvers;
};

//int StaticFEM_Multigrid::smooth_count = 0;