#include <iostream>
#include <type.hpp>

#include <Eigen/Eigen>

#include "Geometry/ParameterDesign/Bezier.hpp"
#include "imgui/implot.h"
#include "Visualization/Visualizer.h"
#include "Geometry/ParameterDesign/ParameterCurve.h"
#include "Geometry/ParameterDesign/Rational.h"

class BezierVisualizer :public Visualizer
{
protected:
	void AddPoint();
	void DragPoint();

	void draw(bool* p_open) override;

	std::vector<Point2d> points;
	std::vector<Float> weights;
	bool updated = true;

	void evaluate_bezier()
	{
		using V = Eigen::Vector2d;
		auto temp_points = std::vector<Eigen::Vector2d>(3);

		points[0] = V(points[0].x(), 0);
		points[1] = V(0, points[1].y());
		temp_points[0] = points[0];
		temp_points[1] = V(points[0].x(), points[1].y());
		temp_points[2] = V(2 * points[0].x(), 2 * points[1].y());

		weights = { 1,1,0 };
		RationalBSpline<3> RBS(temp_points, weights);
		RBS.evaluate();

		for (int i = 0; i < Length; ++i)
		{
			xs[i] = RBS(ctr_points[i]).x();
			ys[i] = RBS(ctr_points[i]).y();
		}
		weights = { 1,-1,0 };
		temp_points[1] *= -1;
		RBS = RationalBSpline<3>(temp_points, weights);

		RBS.evaluate();
		for (int i = 0; i < Length; ++i)
		{
			xs[i + Length] = RBS(1.0 - ctr_points[i]).x();
			ys[i + Length] = RBS(1.0 - ctr_points[i]).y();
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
			//ctr_points[i + Length] = Float(i) / (Length - 1);
		}

		points.resize(2);
		points[0] = Eigen::Vector2d(2, 0);
		points[1] = Eigen::Vector2d(0, 1);
	}
	const float Length = 1001;

	std::vector<float> ctr_points = std::vector<float>(Length);
	std::vector<float> xs = std::vector<float>(2 * Length);
	std::vector<float> ys = std::vector<float>(2 * Length);

	std::vector<float> ctr_xs = std::vector<float>(0);
	std::vector<float> ctr_ys = std::vector<float>(0);
};

void BezierVisualizer::AddPoint()
{
	//if (ImGui::IsMouseClicked(ImGuiMouseButton_Right))
	//{
	//	auto pos = ImPlot::GetPlotMousePos();
	//	points.emplace_back(pos.x, pos.y);
	//	updated = true;
	//}
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
				int delta = 667;

				ImPlot::PlotLine("Bezier", &xs[0], &ys[0], Length + delta);

				ImPlot::PlotLine("Bezier", &xs[Length + delta], &ys[Length + delta], Length - delta);
				//if (!ctr_xs.empty())
				//{
				//	ImPlot::PlotLine("ControlBox", &ctr_xs[0], &ctr_ys[0], ctr_xs.size());
				//}
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