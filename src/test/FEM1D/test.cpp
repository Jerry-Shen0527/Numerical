#include <iostream>
#include <type.hpp>
#include <Interpolation.h>

#include <Eigen/Eigen>

#include "imgui/implot.h"
#include "Visualization/Visualizer.h"

class FEMV1Disualizer :public Visualizer
{
protected:
	void AddPoint();
	void DragPoint();
	void draw(bool* p_open) override;

	std::vector<Point2> points;
	bool updated = true;

	void evaluate_lagrangian()
	{
		LagrangianPolynomial lagrangian(points);
		lagrangian.evaluate();

		RadialInterpolation radial(points);
		radial.evaluate();
		for (int i = 0; i < Length; ++i)
		{
			ys1[i] = lagrangian(xs[i]);
			ys2[i] = radial(xs[i]);
		}
	}

public:
	FEMV1Disualizer() {
		for (int i = 0; i < Length; ++i)
		{
			xs[i] = i;
		}
	}
	const float Length = 1001;

	std::vector<float> xs = std::vector<float>(Length);
	std::vector<float> ys1 = std::vector<float>(Length);
	std::vector<float> ys2 = std::vector<float>(Length);
};

void FEMV1Disualizer::AddPoint()
{
	if (ImGui::IsMouseClicked(ImGuiMouseButton_Right))
	{
		auto pos = ImPlot::GetPlotMousePos();
		points.emplace_back(pos.x, pos.y);
		updated = true;
	}
}

void FEMV1Disualizer::DragPoint()
{
	for (int i = 0; i < points.size(); ++i)
	{
		auto& point = points[i];
		updated |= ImPlot::DragPoint(("Point " + std::to_string(i)).c_str(), &point.x(), &point.y());
	}
}

void FEMV1Disualizer::draw(bool* p_open)
{
	if (updated)
	{
		evaluate_lagrangian();
		updated = false;
	}
	if (ImGui::BeginTabBar("Homework 1")) {
		if (ImGui::BeginTabItem("Interpolation"))
		{
			if (ImPlot::BeginPlot("Line Plot", "x", "f(x)", ImGui::GetContentRegionAvail(), ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMenus)) {
				ImPlot::PlotLine("Polynomial", &xs[0], &ys1[0], Length);
				AddPoint();
				DragPoint();

				ImPlot::PlotLine("Radial", &xs[0], &ys2[0], Length);

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
	FEMV1Disualizer visualizer;
	visualizer.RenderLoop();
}