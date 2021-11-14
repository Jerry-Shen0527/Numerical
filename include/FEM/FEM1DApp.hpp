#pragma once
#include "Domain.hpp"
#include "FEM.hpp"
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
	Float GradientSelfInnerProduct(int i, int j) override
	{
		std::vector<int> i_id, j_id;
		auto i_mesh = IdxToMesh(i, i_id);
		auto j_mesh = IdxToMesh(j, j_id);

		Float ret = 0;

		for (int a = 0; a < i_mesh.size(); ++a)
		{
			for (int b = 0; b < j_mesh.size(); ++b)
			{
				if (i_mesh[a] == j_mesh[b])
				{
					auto sub_interval = interval.SubInterval(i_mesh[a]);
					ret += WeightedL2InnerProduct(sub_interval.scale(ShapeFunctions[i_id[a]]),
						sub_interval.scale(ShapeFunctionGradients[j_id[b]], 1.0 / sub_interval.length()), b_func, sub_interval);
				}
			}
		}

		return ret;
	}
	Float GradientInnerProduct(int i, int j) override
	{
		std::vector<int> i_id, j_id;
		auto i_mesh = IdxToMesh(i, i_id);
		auto j_mesh = IdxToMesh(j, j_id);

		Float ret = 0;

		for (int a = 0; a < i_mesh.size(); ++a)
		{
			for (int b = 0; b < j_mesh.size(); ++b)
			{
				if (i_mesh[a] == j_mesh[b])
				{
					auto sub_interval = interval.SubInterval(i_mesh[a]);
					ret += WeightedL2InnerProduct(sub_interval.scale(ShapeFunctionGradients[i_id[a]], 1.0 / sub_interval.length()),
						sub_interval.scale(ShapeFunctionGradients[j_id[b]], 1.0 / sub_interval.length()), d_func, sub_interval);
				}
			}
		}

		return ret;
	}

	Float SelfInnerProduct(int i, int j) override
	{
		std::vector<int> i_id, j_id;
		auto i_mesh = IdxToMesh(i, i_id);
		auto j_mesh = IdxToMesh(j, j_id);

		Float ret = 0;

		for (int a = 0; a < i_mesh.size(); ++a)
		{
			for (int b = 0; b < j_mesh.size(); ++b)
			{
				if (i_mesh[a] == j_mesh[b])
				{
					auto sub_interval = interval.SubInterval(i_mesh[a]);
					ret += WeightedL2InnerProduct(sub_interval.scale(ShapeFunctions[i_id[a]]), sub_interval.scale(ShapeFunctions[j_id[b]]), c_func, sub_interval);
				}
			}
		}

		return ret;
	}

	Float RHSInnerProduct(int i) override
	{
		std::vector<int> func_id;
		auto i_mesh = IdxToMesh(i, func_id);

		Float ret = 0;

		for (int a = 0; a < func_id.size(); ++a)
		{
			auto sub_interval = interval.SubInterval(i_mesh[a]);
			ret += L2InnerProduct(sub_interval.scale(ShapeFunctions[func_id[a]]), RHS_func, sub_interval);
		}

		return ret;
	}

	std::vector<int> RelatedFuncIdx(int idx) override
	{
		std::vector<int> ret;
		std::vector<int> foo_id;

		auto MeshIds = IdxToMesh(idx, foo_id);

		std::set<int> set_ret;

		for (auto mesh_id : MeshIds)
		{
			for (int i = 0; i < ShapeFunctions.size(); ++i)
			{
				int idx;
				if (MeshToIdx(mesh_id, i, idx))
				{
					set_ret.emplace(idx);
				}
			}
		}
		ret.assign(set_ret.begin(), set_ret.end());
		return ret;
	}

	virtual std::vector<int> IdxToMesh(int idx, std::vector<int>& shapeFuncId) = 0;
	virtual bool MeshToIdx(int mesh_idx, int shapefun_idx, int& idx) = 0;
public:
	Float Value(Float x) override
	{
		if (mat_size == 0)
		{
			return 0;
		}
		Float ret = 0;
		for (int i = 0; i < interval.GetPartitionCount(); ++i)
		{
			auto sub_interval = interval.SubInterval(i);
			if (sub_interval.Inside(x))
			{
				for (int j = 0; j < ShapeFunctions.size(); ++j)
				{
					int idx;
					if (MeshToIdx(i, j, idx))
					{
						ret += sub_interval.scale(ShapeFunctions[j])(x) * rst(idx);
					}
				}
			}
		}
		return ret;
	}

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

std::function<Float(Float)> LagrangianBase(int N, int i)
{
	std::vector<Point2f> points(N + 1);
	Float h = 1.0 / N;
	for (int i = 0; i <= N; ++i)
	{
		points[i] = Point2f(i * h, 0);
	}
	points[i] = Point2f(i * h, 1.0);
	return LagrangianPolynomial(points);
}

std::function<Float(Float)> LagrangianBaseDerivative(int N, int i)
{
	return [=](Float x)
	{
		Float ret = 0;
		for (int missing = 0; missing <= N; ++missing)
		{
			if (missing != i)
			{
				std::vector<Point2f> points;
				Float h = 1.0 / N;

				for (int j = 0; j <= N; ++j)
					if (j != missing)
						if (j == i)
							points.emplace_back(j * h, 1.0);
						else
							points.emplace_back(j * h, 0.0);

				ret += LagrangianPolynomial(points)(x) / (h * (i - missing));
			}
		}
		return ret;
	};
}

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