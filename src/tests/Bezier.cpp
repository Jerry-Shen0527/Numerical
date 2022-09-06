#include <iostream>
#include <type.hpp>

#include <Eigen/Eigen>

#include "Geometry/ParameterDesign/Bezier.hpp"
#include "imgui/implot.h"
#include "Visualization/Visualizer.h"

class BezierVisualizer :public Visualizer
{
protected:
	void AddPoint();
	void DragPoint();

	void draw(bool* p_open) override;

	std::vector<Point2d> points;
	bool updated = true;

	void evaluate_bezier()
	{
		BezierBernstein bezier(points);
		//Bezier bezier(points);
		for (int i = 0; i < Length; ++i)
		{
			xs[i] = bezier(ctr_points[i]).x();
			ys[i] = bezier(ctr_points[i]).y();
		}
		ctr_xs.resize(points.size());
		ctr_ys.resize(points.size());

		for (int i = 0; i < points.size(); ++i)
		{
			ctr_xs[i] = points[i].x();
			ctr_ys[i] = points[i].y();
		}
	}

public:
	BezierVisualizer() {
		for (int i = 0; i < Length; ++i)
		{
			ctr_points[i] = Float(i) / (Length - 1);
		}
	}
	const float Length = 1001;

	std::vector<float> ctr_points = std::vector<float>(Length);
	std::vector<float> xs = std::vector<float>(Length);
	std::vector<float> ys = std::vector<float>(Length);

	std::vector<float> ctr_xs = std::vector<float>(0);
	std::vector<float> ctr_ys = std::vector<float>(0);
};

void BezierVisualizer::AddPoint()
{
	if (ImGui::IsMouseClicked(ImGuiMouseButton_Right))
	{
		auto pos = ImPlot::GetPlotMousePos();
		points.emplace_back(pos.x, pos.y);
		updated = true;
	}
}

void BezierVisualizer::DragPoint()
{
	for (int i = 0; i < points.size(); ++i)
	{
		auto& point = points[i];
		updated |= ImPlot::DragPoint(("Point " + std::to_string(i)).c_str(), &point.x(), &point.y());
	}
}

void BezierVisualizer::draw(bool* p_open)
{
	if (updated)
	{
		evaluate_bezier();
		updated = false;
	}
	if (ImGui::BeginTabBar("Homework 2")) {
		if (ImGui::BeginTabItem("Interpolation"))
		{
			if (ImPlot::BeginPlot("Line Plot", "x", "f(x)", ImGui::GetContentRegionAvail(), ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMenus)) {
				ImPlot::PlotLine("Bezier", &xs[0], &ys[0], Length);
				if (!ctr_xs.empty())
				{
					ImPlot::PlotLine("ControlBox", &ctr_xs[0], &ctr_ys[0], ctr_xs.size());
				}
				AddPoint();
				DragPoint();

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
	BezierVisualizer visualizer;
	visualizer.RenderLoop();
}