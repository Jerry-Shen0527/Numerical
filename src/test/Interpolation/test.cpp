#include <iostream>
#include <type.hpp>
#include <Interpolation.h>

#include <Eigen/Eigen>

#include "imgui/implot.h"
#include "Visualization/Visualizer.h"

class InterpolationVisualizer :public Visualizer
{
protected:
	void draw(bool* p_open) override;
public:
};

void DrawPolynomial()
{
	static float xs1[1001], ys1[1001];
	for (int i = 0; i < 1001; ++i) {
		xs1[i] = i * 0.001f;
		ys1[i] = 0.5f + 0.5f * sinf(50 * (xs1[i] + (float)ImGui::GetTime() / 10));
	}
	static double xs2[11], ys2[11];
	for (int i = 0; i < 11; ++i) {
		xs2[i] = i * 0.1f;
		ys2[i] = xs2[i] * xs2[i];
	}
	ImGui::BulletText("Anti-aliasing can be enabled from the plot's context menu (see Help).");
	if (ImPlot::BeginPlot("Line Plot", "x", "f(x)")) {
		ImPlot::PlotLine("sin(x)", xs1, ys1, 1001);
		ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
		ImPlot::PlotLine("x^2", xs2, ys2, 11);
		ImPlot::EndPlot();
	}
}

void DrawRadial()
{
}

void InterpolationVisualizer::draw(bool* p_open)
{
	if (ImGui::BeginTabBar("ImPlotDemoTabs")) {
		if (ImGui::BeginTabItem("Homework1"))
		{
			if (ImGui::CollapsingHeader("DrawPolynomial"))
				DrawPolynomial();
			if (ImGui::CollapsingHeader("Filled Line Plots"))
				DrawRadial();

			ImGui::EndTabItem();
		}

		if (ImGui::BeginTabItem("DrawPolynomial"))
		{
			DrawPolynomial();
			ImGui::EndTabItem();
		}

		ImGui::EndTabBar();
	}
	ImGui::End();
}

int main()
{
	//std::vector<Point2> points;
	//points.emplace_back(0, 0);
	//points.emplace_back(1, 1);
	//points.emplace_back(2, 2);
	//points.emplace_back(33, 3);
	//points.emplace_back(4, 48);

	//NewtonPolynomial newton(points);
	//LagrangianPolynomial lagrangian(points);
	//newton.evaluate();
	//lagrangian.evaluate();
	//std::cout << newton(0.6);

	InterpolationVisualizer visualizer;
	visualizer.RenderLoop();
}