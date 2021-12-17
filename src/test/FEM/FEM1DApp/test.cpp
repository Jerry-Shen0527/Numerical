#include <iostream>
#include <numeric>
#include <type.hpp>

#include <Eigen/Eigen>

#include "FEM/FEM1DApp.hpp"
#include "imgui/implot.h"
#include "Visualization/Visualizer.h"

#include <numeric>

using Linear = PolynomialFEMApp<1>;
using Quadratic = PolynomialFEMApp<2>;

class FEM1DVisualizer :public Visualizer
{
protected:

	void evaluate()
	{
		int segement_ = segemnt;

		if (segement_ % 2 != 0)
		{
			segement_ += 1;
		}

		//Homework 2
		//auto rhs = [](Float x) {return  -4 - x + Power(x, 2) - Power(x, 3) + Power(x, 4) + Cos(x) - 2 * x * Cos(x) - 2 * Sin(x); };
		//auto a = [](Float x) {return sin(x) + 2; };
		//auto c = [](Float x) {return x * x + 1; };
		auto rhs = [](Float x) {return  2 * cos(1 - x) * cos(x) + 2 * sin(1 - x) * sin(x); };
		auto d = [this](Float x) {return 1; };
		auto b = [](Float x) {return 0;	};
		auto c = [](Float x) {return 0; };

		Interval interval(0.0, 1.0);

		interval.SetPartitionCount(segement_);

		StaticFEM1DMG mg_linear(rhs, d, b, c, interval);
		Linear cg_linear(rhs, d, b, c, interval);
		mg_linear.evaluate();
		cg_linear.evaluate();
		for (int i = 0; i < Length; ++i)
		{
			cg_val[i] = cg_linear.Value(1.0 / (Length - 1) * i);
			mg_val[i] = mg_linear.Value(1.0 / (Length - 1) * i);

			cg_diff[i] = cg_val[i] - precise_val[i];
			mg_diff[i] = mg_val[i] - precise_val[i];
		}
	}

	void Control_UI();
	void draw(bool* p_open) override;

	std::vector<Point2d> points;
	bool updated = true;
	int segemnt = 16;
	float epsilon = 1E-7;

	void error(std::vector<float>& ref, std::vector<float>& eval, Float& L_1, Float& L_2, Float& L_inf)
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
	void CalcAccurateRst()
	{
		Float h = 1.0 / (Length - 1);
		for (int i = 0; i < Length; ++i)
		{
			xs[i] = i * h;

			auto accurate_func = [this](Float x) {return sin(x) * sin(1 - x); };

			precise_val[i] = accurate_func(xs[i]);
		}
	}

	FEM1DVisualizer() {
		CalcAccurateRst();
		segemnt = 16;
		Float L1, L2, L_inf;
		do
		{
			//std::cout << "Segmentations: " << segemnt << std::endl;
			evaluate();

			error(precise_val, mg_val, L1, L2, L_inf);

			using std::cout;
			using std::endl;

			//if (segemnt == 16)
			//{
			//	cout << segemnt << '&' << L1 << '&' << '-' << '&' << L2 << '&' << '-' << '&' << L_inf << '&' << '-' << "\\\\" << endl;
			//}
			//else
			//	cout << segemnt << '&' << L1 << '&' <<- log2(L1 / linear_L1.back()) << '&' << L2 << '&' << -log2(L2 / mg_L2.back()) << '&' << L_inf << '&' << -log2(L_inf / mg_Linf.back()) << "\\\\" << endl;

			pointcount.push_back(segemnt);
			mg_L1.push_back(L1);
			mg_L2.push_back(L2);
			mg_Linf.push_back(L_inf);
			error(precise_val, cg_val, L1, L2, L_inf);

			//if (segemnt == 16)
			//{
			//	cout << segemnt << '&' << L1 << '&' << '-' << '&' << L2 << '&' << '-' << '&' << L_inf << '&' << '-' << "\\\\" << endl;
			//}
			//else
			//	cout << segemnt << '&' << L1 << '&' << -log2(L1 / cg_L1.back()) << '&' << L2 << '&' << -log2(L2 / cg_L2.back()) << '&' << L_inf << '&' << -log2(L_inf / cg_Linf.back()) << "\\\\" << endl;

			cg_L1.push_back(L1);
			cg_L2.push_back(L2);
			cg_Linf.push_back(L_inf);

			segemnt *= 2;
		} while (segemnt != 8192);
		segemnt = 16;

		evaluate();
	}
	const size_t Length = 10001;

	std::vector<float> xs = std::vector<float>(Length);
	std::vector<float> precise_val = std::vector<float>(Length);
	std::vector<float> cg_val = std::vector<float>(Length);
	std::vector<float> mg_val = std::vector<float>(Length);
	std::vector<float> cg_diff = std::vector<float>(Length);
	std::vector<float> mg_diff = std::vector<float>(Length);

	std::vector<float> pointcount;
	std::vector<float> mg_L1;
	std::vector<float> mg_L2;
	std::vector<float> mg_Linf;

	std::vector<float> cg_L1;
	std::vector<float> cg_L2;
	std::vector<float> cg_Linf;
};

static inline ImVec2 operator-(const ImVec2& lhs, const ImVec2& rhs) { return ImVec2(lhs.x - rhs.x, lhs.y - rhs.y); }

void FEM1DVisualizer::Control_UI()
{
	if (ImGui::SliderInt("Number of segments", &segemnt, 3, 200))
	{
		segemnt = segemnt < 2 ? 2 : segemnt;
		evaluate();
	}
	if (ImGui::SliderFloat("Epsilon", &epsilon, 1E-7, 1E-1, "%.8f", ImGuiSliderFlags_Logarithmic))
	{
		evaluate();
		CalcAccurateRst();
	}
}

void FEM1DVisualizer::draw(bool* p_open)
{
	if (ImGui::BeginTabBar("FEM 1D App")) {
		if (ImGui::BeginTabItem("FEM1D"))
		{
			if (ImPlot::BeginPlot("Line Plot", "x", "f(x)", ImGui::GetContentRegionAvail() - ImVec2(0, 100), ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMenus)) {
				ImPlot::PlotLine("Precise solution", &xs[0], &precise_val[0], Length);
				ImPlot::PlotLine("FEM Conjugate Gradient", &xs[0], &cg_val[0], Length);
				ImPlot::PlotLine("FEM Multigrid", &xs[0], &mg_val[0], Length);

				ImPlot::EndPlot();
			}
			Control_UI();
			ImGui::EndTabItem();
		}
		if (ImGui::BeginTabItem("FEM1D difference"))
		{
			if (ImPlot::BeginPlot("Line Plot", "x", "f(x)", ImGui::GetContentRegionAvail() - ImVec2(0, 100), ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMenus)) {
				ImPlot::PlotLine("FEM diff Multigrid", &xs[0], &mg_diff[0], Length);
				ImPlot::PlotLine("FEM diff Conjugate Gradient", &xs[0], &cg_diff[0], Length);

				ImPlot::EndPlot();
			}
			Control_UI();
			ImGui::EndTabItem();
		}

		if (ImGui::BeginTabItem("Error Plot"))
		{
			if (ImPlot::BeginPlot("Line Plot", "x", "f(x)", ImGui::GetContentRegionAvail() - ImVec2(0, 100), ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMenus, ImPlotAxisFlags_LogScale, ImPlotAxisFlags_LogScale)) {
				ImPlot::PlotLine("Multigrid L1    Error", &pointcount[0], &mg_L1[0], pointcount.size());
				ImPlot::PlotLine("Multigrid L2    Error", &pointcount[0], &mg_L2[0], pointcount.size());
				ImPlot::PlotLine("Multigrid L_inf Error", &pointcount[0], &mg_Linf[0], pointcount.size());

				ImPlot::PlotLine("Conjugate Gradient L1    Error", &pointcount[0], &cg_L1[0], pointcount.size());
				ImPlot::PlotLine("Conjugate Gradient L2    Error", &pointcount[0], &cg_L2[0], pointcount.size());
				ImPlot::PlotLine("Conjugate Gradient L_inf Error", &pointcount[0], &cg_Linf[0], pointcount.size());

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