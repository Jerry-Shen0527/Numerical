#include <iostream>
#include <numeric>
#include <type.hpp>
#include <set>

#include <Eigen/Eigen>

#include "FEM/FEM.hpp"
#include "Geometry/ParameterDesign/Interpolation.h"
#include "imgui/implot.h"
#include "NIntegrate/Integrate.hpp"
#include "Visualization/Visualizer.h"

class StaticFEM1DApp :public StaticFEM1D
{
public:
	StaticFEM1DApp(int segment, const std::function<Float(Float)>& rhs_func, const std::function<Float(Float)>& d_func,
		const std::function<Float(Float)>& c_func, const Interval& interval)
		: StaticFEM1D(segment), RHS_func(rhs_func),
		d_func(d_func),
		c_func(c_func),
		interval(interval)
	{
		this->interval.SetPartitionCount(segment);
	}

protected:
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
						ret += sub_interval.scale(ShapeFunctions[j]
						)(x) * rst(idx);
					}
				}
			}
		}
		return ret;
	}

	std::vector<std::function<Float(Float)>> ShapeFunctions;
	std::vector<std::function<Float(Float)>> ShapeFunctionGradients;

	//-a u''+b u'+c u=RHSFunc
	std::function<Float(Float)> RHS_func;
	std::function<Float(Float)> d_func;
	std::function<Float(Float)> c_func;

	Interval interval;
};

std::function<Float(Float)> LagrangianBase(int N, int i)
{
	std::vector<Point2> points(N + 1);
	Float h = 1.0 / N;
	for (int i = 0; i <= N; ++i)
	{
		points[i] = Point2(i * h, 0);
	}
	points[i] = Point2(i * h, 1.0);
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
				std::vector<Point2> points;
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
	PolynomialFEMApp(int segment, const std::function<Float(Float)>& rhs_func,
		const std::function<Float(Float)>& d_func, const std::function<Float(Float)>& c_func,
		const Interval& interval)
		: StaticFEM1DApp(segment, rhs_func, d_func, c_func, interval)
	{
		static_assert(N > 0);
		mat_size = N * segment - 1;

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

public:
};

using Linear = PolynomialFEMApp<1>;
using Quadratic = PolynomialFEMApp<2>;

//For compatibility with Mathematica
#define Power pow
#define Sin sin
#define Cos cos

class FEM1DVisualizer :public Visualizer
{
protected:

	void evaluate()
	{
		int segement_ = segemnt;

		auto rhs = [](Float x) {return  -4 - x + Power(x, 2) - Power(x, 3) + Power(x, 4) + Cos(x) - 2 * x * Cos(x) - 2 * Sin(x); };
		auto a = [](Float x) {return sin(x) + 2; };
		auto c = [](Float x) {return x * x + 1; };

		Linear linear(segement_, rhs, a, c, Interval(0.0, 1.0));
		Quadratic quadratic(segement_, rhs, a, c, Interval(0.0, 1.0));
		linear.evaluate();
		quadratic.evaluate();
		for (int i = 0; i < Length; ++i)
		{
			quadratic_val[i] = quadratic.Value(1.0 / (Length - 1) * i);
			linear_val[i] = linear.Value(1.0 / (Length - 1) * i);

			quadratic_diff[i] = quadratic_val[i] - precise_val[i];
			linear_diff[i] = linear_val[i] - precise_val[i];
		}
	}

	void draw(bool* p_open) override;

	std::vector<Point2> points;
	bool updated = true;
	int segemnt = 1;

	void error(std::vector<Float>& ref, std::vector<Float>& eval, Float& L_1, Float& L_2, Float& L_inf)
	{
		assert(ref.size() == eval.size());

		std::vector<Float> minus(ref.size());

		for (int i = 0; i < ref.size(); ++i)
		{
			minus[i] = abs(ref[i] - eval[i]);
		}

		L_inf = *std::max_element(minus.begin(), minus.end(), [](Float a, Float b) {return a < b; });
		L_2 = sqrt(std::accumulate(minus.begin(), minus.end(), static_cast<Float>(0), [](Float r, Float a) {return r + a * a; }) / Float(ref.size()));
		L_1 = std::accumulate(minus.begin(), minus.end(), static_cast<Float>(0), [](Float r, Float a) {return r + a; }) / Float(ref.size());
	}

public:
	FEM1DVisualizer() {
		Float h = 1.0 / (Length - 1);
		for (int i = 0; i < Length; ++i)
		{
			xs[i] = i * h;
			rhs_f[i] = (xs[i] - 1) * sin(xs[i]);
			auto accurate_func = [](Float x) {return (-1 + x) * x; };
			precise_val[i] = accurate_func(xs[i]);
		}
		segemnt = 2;

		Float L1, L2, L_inf;

		do
		{
			evaluate();

			error(precise_val, linear_val, L1, L2, L_inf);

			pointcount.push_back(segemnt);
			linear_vec_L_1.push_back(L1);
			linear_vec_L_2.push_back(L2);
			linear_vec_L_inf.push_back(L_inf);
			error(precise_val, quadratic_val, L1, L2, L_inf);

			quadratic_vec_L_1.push_back(L1);
			quadratic_vec_L_2.push_back(L2);
			quadratic_vec_L_inf.push_back(L_inf);

			segemnt *= 2;
		} while (segemnt != 1024);
		segemnt = 2;
		evaluate();
	}
	const size_t Length = 2001;

	std::vector<Float> xs = std::vector<Float>(Length);
	std::vector<Float> rhs_f = std::vector<Float>(Length);
	std::vector<Float> precise_val = std::vector<Float>(Length);
	std::vector<Float> quadratic_val = std::vector<Float>(Length);
	std::vector<Float> linear_val = std::vector<Float>(Length);
	std::vector<Float> quadratic_diff = std::vector<Float>(Length);
	std::vector<Float> linear_diff = std::vector<Float>(Length);

	std::vector<Float> pointcount;
	std::vector<Float> linear_vec_L_1;
	std::vector<Float> linear_vec_L_2;
	std::vector<Float> linear_vec_L_inf;

	std::vector<Float> quadratic_vec_L_1;
	std::vector<Float> quadratic_vec_L_2;
	std::vector<Float> quadratic_vec_L_inf;
};

static inline ImVec2 operator-(const ImVec2& lhs, const ImVec2& rhs) { return ImVec2(lhs.x - rhs.x, lhs.y - rhs.y); }
void FEM1DVisualizer::draw(bool* p_open)
{
	if (ImGui::BeginTabBar("FEM 1D App")) {
		if (ImGui::BeginTabItem("FEM1D"))
		{
			if (ImPlot::BeginPlot("Line Plot", "x", "f(x)", ImGui::GetContentRegionAvail() - ImVec2(0, 100), ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMenus)) {
				//ImPlot::PlotLine("u", &xs[0], &ys1[0], Length);
				ImPlot::PlotLine("Precise solution", &xs[0], &precise_val[0], Length);
				//ImPlot::PlotLine("FEM Result", &xs[0], &ys3[0], Length);
				ImPlot::PlotLine("FEM Quadratic", &xs[0], &quadratic_val[0], Length);
				ImPlot::PlotLine("FEM Linear", &xs[0], &linear_val[0], Length);

				ImPlot::EndPlot();
			}
			if (ImGui::SliderInt("Number of segments", &segemnt, 1, 200))
			{
				evaluate();
			}
			ImGui::EndTabItem();
		}
		if (ImGui::BeginTabItem("FEM1D difference"))
		{
			if (ImPlot::BeginPlot("Line Plot", "x", "f(x)", ImGui::GetContentRegionAvail() - ImVec2(0, 100), ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMenus)) {
				//ImPlot::PlotLine("u", &xs[0], &ys1[0], Length);
				//ImPlot::PlotLine("FEM Result", &xs[0], &ys3[0], Length);
				ImPlot::PlotLine("FEM diff Linear", &xs[0], &linear_diff[0], Length);
				ImPlot::PlotLine("FEM diff Quadratic", &xs[0], &quadratic_diff[0], Length);

				ImPlot::EndPlot();
			}
			if (ImGui::SliderInt("Number of segments", &segemnt, 1, 200))
			{
				segemnt = segemnt < 1 ? 1 : segemnt;
				evaluate();
			}
			ImGui::EndTabItem();
		}

		if (ImGui::BeginTabItem("Error Plot"))
		{
			if (ImPlot::BeginPlot("Line Plot", "x", "f(x)", ImGui::GetContentRegionAvail() - ImVec2(0, 100), ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMenus, ImPlotAxisFlags_LogScale, ImPlotAxisFlags_LogScale)) {
				//ImPlot::PlotLine("u", &xs[0], &ys1[0], Length);
				//ImPlot::PlotLine("FEM Result", &xs[0], &ys3[0], Length);
				ImPlot::PlotLine("Linear L1    Error", &pointcount[0], &linear_vec_L_1[0], pointcount.size());
				ImPlot::PlotLine("Linear L2    Error", &pointcount[0], &linear_vec_L_2[0], pointcount.size());
				ImPlot::PlotLine("Linear L_inf Error", &pointcount[0], &linear_vec_L_inf[0], pointcount.size());

				ImPlot::PlotLine("Quadratic L1    Error", &pointcount[0], &quadratic_vec_L_1[0], pointcount.size());
				ImPlot::PlotLine("Quadratic L2    Error", &pointcount[0], &quadratic_vec_L_2[0], pointcount.size());
				ImPlot::PlotLine("Quadratic L_inf Error", &pointcount[0], &quadratic_vec_L_inf[0], pointcount.size());

				ImPlot::EndPlot();
			}

			ImGui::EndTabItem();
		}
		ImGui::EndTabBar();
	}

	ImGui::End();
}

int main()
{
	FEM1DVisualizer visualizer;
	visualizer.RenderLoop();
}