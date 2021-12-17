#pragma once
#include "Domain.hpp"
#include "FEM.hpp"
#include "FEM_MG.hpp"
#include "Geometry/ParameterDesign/Interpolation.h"
#include "NIntegrate/Integrate.hpp"

class StaticFEM1DApp :public StaticFEM1D
{
public:
	StaticFEM1DApp(const std::function<Float(Float)>& rhs_func, const std::function<Float(Float)>& d_func, const std::function<Float(Float)>& b_func,
		const std::function<Float(Float)>& c_func, const Interval& interval)
		: RHS_func(rhs_func),
		d_func(d_func),
		b_func(b_func),
		c_func(c_func),
		interval(interval)
	{
	}

protected:
	Float GradientSelfInnerProduct(int i, int j) override;
	Float GradientInnerProduct(int i, int j) override;
	Float SelfInnerProduct(int i, int j) override;
	Float RHSInnerProduct(int i) override;

	std::vector<int> RelatedFuncIdx(int idx) override;

	virtual std::vector<int> IdxToMesh(int idx, std::vector<int>& shapeFuncId) = 0;
	virtual bool MeshToIdx(int mesh_idx, int shapefun_idx, int& idx) = 0;
public:
	Float Value(Float x) override;

public:
	std::vector<std::function<Float(Float)>> ShapeFunctions;
	std::vector<std::function<Float(Float)>> ShapeFunctionGradients;

	//- (d[x]u'[x])'+ b u'[x]+c u=RHSFunc
	std::function<Float(Float)> RHS_func;
	std::function<Float(Float)> d_func;
	std::function<Float(Float)> b_func;
	std::function<Float(Float)> c_func;

	Interval interval;
};

//*************************************************************//
//******************Polynomial FEM App*************************//
//*************************************************************//

std::function<Float(Float)> LagrangianBase(int N, int i);
std::function<Float(Float)> LagrangianBaseDerivative(int N, int i);

template<int N>
class PolynomialFEMApp :public StaticFEM1DApp
{
public:
	PolynomialFEMApp(const std::function<Float(Float)>& rhs_func,
		const std::function<Float(Float)>& d_func, const std::function<Float(Float)>& b_func, const std::function<Float(Float)>& c_func,
		const Interval& interval)
		: StaticFEM1DApp(rhs_func, d_func, b_func, c_func, interval)
	{
		static_assert(N > 0);

		for (int i = 0; i <= N; ++i)
		{
			ShapeFunctions.push_back(LagrangianBase(N, i));
			ShapeFunctionGradients.push_back(LagrangianBaseDerivative(N, i));
		}
	}

protected:
	std::vector<int> IdxToMesh(int idx, std::vector<int>& shapeFuncId) override
	{
		if ((1 + idx) % N != 0)
		{
			shapeFuncId = { (1 + idx) % N };
			return { idx / N };
		}
		else
		{
			shapeFuncId = { N,0 };
			return { idx / N,idx / N + 1 };
		}
	}

	bool MeshToIdx(int mesh_idx, int shapefunc_idx, int& idx) override
	{
		if ((mesh_idx == 0 && shapefunc_idx == 0) || (mesh_idx == interval.GetPartitionCount() - 1 && shapefunc_idx == N))
			return false;
		else
		{
			idx = N * mesh_idx + shapefunc_idx - 1;
		}

		return true;
	}

	void SetMatSize() override
	{
		mat_size = N * interval.GetPartitionCount() - 1;
	}

public:
};

class StaticFEM1DMG final :public StaticFEM_Multigrid
{
public:
	StaticFEM1DMG(const std::function<Float(Float)>& rhs_func,
		const std::function<Float(Float)>& d_func, const std::function<Float(Float)>& b_func, const std::function<Float(Float)>& c_func,
		const Interval& interval)
	{
		//Make sure the interval can be sub-partitioned.
		assert(interval.GetPartitionCount() % 2 == 0 && interval.GetPartitionCount() > 2);
		std::vector<Float> knots = interval.GetSubIntervalKnots();
		solvers.insert(solvers.begin(), std::make_shared<PolynomialFEMApp<1>>(rhs_func, d_func, b_func, c_func, interval));

		while (true)
		{
			if (knots.size() % 2 == 0 || knots.size() <= 1)
			{
				//When the interval has odd intervals, it's harder to map the points. However, it's not impossible.
				//TODO::Possibly fully subdivide the interval.
				break;
			}
			std::vector<Float> new_knots(knots.size() / 2);

			Interval new_interval(interval.lerp(0), interval.lerp(1));
			for (int i = 0; i < new_knots.size(); ++i)
			{
				new_knots[i] = knots[i * 2 + 1];
			}
			new_interval.SetSubIntervalKnots(new_knots);

			solvers.insert(solvers.begin(), std::make_shared<PolynomialFEMApp<1>>(rhs_func, d_func, b_func, c_func, new_interval));
			knots = new_knots;
		}
	}
	Float Value(Float x)
	{
		std::shared_ptr<PolynomialFEMApp<1>> ptr = std::static_pointer_cast<PolynomialFEMApp<1>>(solvers.back());
		return ptr->Value(x);
	}
private:

	Vector I_k_down(const Vector& vec)
	{
		Vector ret(vec.rows() / 2);

		auto size = ret.rows();

		for (int i = 0; i < size; ++i)
			ret[i] = vec[2 * i + 0] / 2 + vec[2 * i + 1] + vec[2 * i + 2] / 2;

		return ret;
	}

	void ProjectDown(std::shared_ptr<StaticFEMwithCG> up, std::shared_ptr<StaticFEMwithCG> down) override
	{
		auto down_rhs = I_k_down(up->GetResidual());

		//std::cout << "Project down:"<<down_rhs.transpose()<< std::endl;

		down->SetRhs(down_rhs);
		Vector zeros(down->GetRhs().rows());
		zeros.setZero();
		down->SetRst(zeros);
	}

	void ProjectUp(std::shared_ptr<StaticFEMwithCG> down, std::shared_ptr<StaticFEMwithCG> up) override
	{
		auto q = down->GetRst();

		Vector up_correction(2 * q.rows() + 1);
		std::shared_ptr<PolynomialFEMApp<1>> ptr = std::static_pointer_cast<PolynomialFEMApp<1>>(up);
		auto knots = ptr->interval.GetSubIntervalKnots();

		assert(knots.size() == 2 * q.rows() + 1);

		int lerp_count = 0;
		int direct_count = 0;
		for (int i = 0; i < 2 * q.rows() + 1; ++i)
		{
			if (i % 2 != 0)
			{
				up_correction[i] = q[i / 2];
				direct_count++;
			}
			else
			{
				lerp_count++;

				Float left = i - 1 < 0 ? ptr->interval.lerp(0) : knots[i - 1];
				Float left_val = i - 1 < 0 ? 0 : q[i / 2 - 1];
				Float right = i >= 2 * q.rows() ? ptr->interval.lerp(1) : knots[i + 1];
				Float right_val = i >= 2 * q.rows() ? 0 : q[i / 2];
				Float middle = knots[i];
				up_correction[i] = Lerp((middle - left) / (right - left), left_val, right_val);
				//up_correction[i]
			}
		}

		assert(lerp_count - direct_count == 1);
		up->GetRst() += up_correction;
	}
};