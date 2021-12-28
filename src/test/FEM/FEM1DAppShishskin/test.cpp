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

using LinearDG = PolynomialFEMAppSD<1>;
using QuadraticDG = PolynomialFEMAppSD<2>;

class FEM1DVisualizer :public Visualizer
{
protected:

	void evaluate()
	{
		int segement_ = segemnt;

		//Homework 2
		//auto rhs = [](Float x) {return  -4 - x + Power(x, 2) - Power(x, 3) + Power(x, 4) + Cos(x) - 2 * x * Cos(x) - 2 * Sin(x); };
		//auto a = [](Float x) {return sin(x) + 2; };
		//auto c = [](Float x) {return x * x + 1; };
		auto rhs = [](Float x) {return  x; };
		auto d = [this](Float x) {return epsilon; };
		auto b = [](Float x) {return 1;	};
		auto c = [](Float x) {return 0; };

		Interval interval(0.0, 1.0);

		if (use_shishkin)
		{
			Interval interval1(0.0, 1 - 2 * epsilon * log(segement_));
			interval1.SetPartitionCount(segement_);
			Interval interval2(1 - 2 * epsilon * log(segement_), 1.0);
			interval2.SetPartitionCount(segement_);

			std::vector<Float> knot_vector(2 * segement_ - 1);

			for (int i = 0; i < segement_ - 1; ++i)
			{
				knot_vector[i] = interval1.SubInterval(i).lerp(1.0);
				knot_vector[segement_ + i] = interval2.SubInterval(i).lerp(1.0);
			}
			if (segement_ > 1)
				knot_vector[segement_ - 1] = interval1.SubInterval(segement_ - 1).lerp(1.0);

			interval.SetSubIntervalKnots(knot_vector);
		}
		else
			interval.SetPartitionCount(segement_);

		LinearDG linear(rhs, d, b, c, interval);
		QuadraticDG quadratic(rhs, d, b, c, interval);
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

	void Control_UI();
	void draw(bool* p_open) override;

	std::vector<Point2d> points;
	bool updated = true;
	int segemnt = 16;
	float epsilon = 1E-7;
	bool use_shishkin = false;

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

			auto accurate_func = [this](Float x) {return -(-1 + exp(x / epsilon) + pow(x, 2) - exp(1 / epsilon) * pow(x, 2) + 2 * epsilon * (-1 + exp(x / epsilon) + x - exp(1 / epsilon) * x)) / (2. * (-1 + exp(1 / epsilon))); };

			if (epsilon < 1E-2)
			{
				auto accurate_func = [this](Float x) {return -exp(1 / epsilon * (x - 1)) / 2. + x * epsilon + x * x / 2.; };
				precise_val[i] = accurate_func(xs[i]);
			}
			else
			{
				precise_val[i] = accurate_func(xs[i]);
			}
		}
	}

	FEM1DVisualizer() {
		CalcAccurateRst();
		segemnt = 16;
		Float L1, L2, L_inf;
		do
		{
			evaluate();

			error(precise_val, linear_val, L1, L2, L_inf);

			using std::cout;
			using std::endl;

			if (segemnt == 16)
			{
				cout << segemnt << '&' << L1 << '&' << '-' << '&' << L2 << '&' << '-' << '&' << L_inf << '&' << '-' << "\\\\" << endl;
			}
			else
				cout << segemnt << '&' << L1 << '&' <<- log2(L1 / linear_L1.back()) << '&' << L2 << '&' << -log2(L2 / linear_L2.back()) << '&' << L_inf << '&' << -log2(L_inf / linear_Linf.back()) << "\\\\" << endl;

			pointcount.push_back(segemnt);
			linear_L1.push_back(L1);
			linear_L2.push_back(L2);
			linear_Linf.push_back(L_inf);
			error(precise_val, quadratic_val, L1, L2, L_inf);

			//if (segemnt == 16)
			//{
			//	cout << segemnt << '&' << L1 << '&' << '-' << '&' << L2 << '&' << '-' << '&' << L_inf << '&' << '-' << "\\\\" << endl;
			//}
			//else
			//	cout << segemnt << '&' << L1 << '&' << -log2(L1 / quadratic_L1.back()) << '&' << L2 << '&' << -log2(L2 / quadratic_L2.back()) << '&' << L_inf << '&' << -log2(L_inf / quadratic_Linf.back()) << "\\\\" << endl;

			quadratic_L1.push_back(L1);
			quadratic_L2.push_back(L2);
			quadratic_Linf.push_back(L_inf);

			segemnt *= 2;
		} while (segemnt != 8192);
		segemnt = 16;

		evaluate();
	}
	const size_t Length = 20001;

	std::vector<float> xs = std::vector<float>(Length);
	std::vector<float> precise_val = std::vector<float>(Length);
	std::vector<float> quadratic_val = std::vector<float>(Length);
	std::vector<float> linear_val = std::vector<float>(Length);
	std::vector<float> quadratic_diff = std::vector<float>(Length);
	std::vector<float> linear_diff = std::vector<float>(Length);

	std::vector<float> pointcount;
	std::vector<float> linear_L1;
	std::vector<float> linear_L2;
	std::vector<float> linear_Linf;

	std::vector<float> quadratic_L1;
	std::vector<float> quadratic_L2;
	std::vector<float> quadratic_Linf;
};

static inline ImVec2 operator-(const ImVec2& lhs, const ImVec2& rhs) { return ImVec2(lhs.x - rhs.x, lhs.y - rhs.y); }

void FEM1DVisualizer::Control_UI()
{
	if (ImGui::SliderInt("Number of segments", &segemnt, 2, 200))
	{
		segemnt = segemnt < 2 ? 2 : segemnt;
		evaluate();
	}
	if (ImGui::SliderFloat("Epsilon", &epsilon, 1E-7, 1E-1, "%.8f", ImGuiSliderFlags_Logarithmic))
	{
		evaluate();
		CalcAccurateRst();
	}
	if (ImGui::Checkbox("Use Shishkin", &use_shishkin))
	{
		evaluate();
	}
}

void FEM1DVisualizer::draw(bool* p_open)
{
	if (ImGui::BeginTabBar("FEM 1D App")) {
		if (ImGui::BeginTabItem("FEM1D"))
		{
			if (ImPlot::BeginPlot("Line Plot", "x", "f(x)", ImGui::GetContentRegionAvail() - ImVec2(0, 100), ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMenus)) {
				ImPlot::PlotLine("Precise solution", &xs[0], &precise_val[0], Length);
				ImPlot::PlotLine("FEM Quadratic", &xs[0], &quadratic_val[0], Length);
				ImPlot::PlotLine("FEM Linear", &xs[0], &linear_val[0], Length);

				ImPlot::EndPlot();
			}
			Control_UI();
			ImGui::EndTabItem();
		}
		if (ImGui::BeginTabItem("FEM1D difference"))
		{
			if (ImPlot::BeginPlot("Line Plot", "x", "f(x)", ImGui::GetContentRegionAvail() - ImVec2(0, 100), ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMenus)) {
				ImPlot::PlotLine("FEM diff Linear", &xs[0], &linear_diff[0], Length);
				ImPlot::PlotLine("FEM diff Quadratic", &xs[0], &quadratic_diff[0], Length);

				ImPlot::EndPlot();
			}
			Control_UI();
			ImGui::EndTabItem();
		}

		if (ImGui::BeginTabItem("Error Plot"))
		{
			if (ImPlot::BeginPlot("Line Plot", "x", "f(x)", ImGui::GetContentRegionAvail() - ImVec2(0, 100), ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMenus, ImPlotAxisFlags_LogScale, ImPlotAxisFlags_LogScale)) {
				ImPlot::PlotLine("Linear L1    Error", &pointcount[0], &linear_L1[0], pointcount.size());
				ImPlot::PlotLine("Linear L2    Error", &pointcount[0], &linear_L2[0], pointcount.size());
				ImPlot::PlotLine("Linear L_inf Error", &pointcount[0], &linear_Linf[0], pointcount.size());

				ImPlot::PlotLine("Quadratic L1    Error", &pointcount[0], &quadratic_L1[0], pointcount.size());
				ImPlot::PlotLine("Quadratic L2    Error", &pointcount[0], &quadratic_L2[0], pointcount.size());
				ImPlot::PlotLine("Quadratic L_inf Error", &pointcount[0], &quadratic_Linf[0], pointcount.size());

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