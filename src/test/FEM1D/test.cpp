#include <iostream>
#include <type.hpp>

#include <Eigen/Eigen>

#include "FEM/FEM.hpp"
#include "imgui/implot.h"
#include "Visualization/Visualizer.h"

static Float InnerProductLinear(Float x, Float xmin, Float xmax)
{
	return ((-2 + x * x + xmax - x * (1 + xmax)) * cos(x) + 2 * cos(xmax) + (1 - 2 * x + xmax) * sin(x) + (-1 + xmax) * sin(xmax))
		/ (x - xmax) + 1 / (x - xmin) * ((2 + x - x * x + (-1 + x) * xmin) * cos(x) - 2 * cos(xmin) - (1 - 2 * x + xmin) * sin(x) + sin(xmin) - xmin * sin(xmin));
}

static Float LinearFunc(Float x, Float x0, Float x_left, Float x_right)
{
	if (x<x_left || x>x_right)
	{
		return 0;
	}
	if (x < x0)
	{
		return (x - x_left) / (x0 - x_left);
	}
	else return (x - x_right) / (x0 - x_right);
}

class LinearBase1DFEM :public StaticFEM
{
public:
	LinearBase1DFEM(size_t size) :size(size)
	{
		mesh.resize(size);
		//A default devision is provided
		Float h = 1.0 / (size + 1);
		for (int i = 0; i < size; ++i)
		{
			mesh[i] = (i + 1) * h;
		}
	}
	std::vector<Float> mesh;

	Float value(Float x) {
		Float ret = 0;
		for (int i = 0; i < size; ++i)
		{
			ret += LinearFunc(x, Value(i), Value(i - 1), Value(i + 1)) * coeff_()(i);
		}
		return ret;
	}

protected:
	void FillRhs() override;
private:

	//return the indices of the functions related to the current function
	std::vector<int> RelatedFuncIdx(int idx)
	{
		std::vector<int> ret;
		ret.push_back(idx);
		int left = idx - 1, right = idx + 1;
		if (0 <= left && left < size)
		{
			ret.push_back(left);
		}
		if (0 <= right && right < size)
		{
			ret.push_back(right);
		}
		return ret;
	}

	Float Value(int index)
	{
		if (index < 0)
		{
			return 0;
		}
		if (index >= size)
		{
			return  1.0;
		}
		return mesh[index];
	}

	Float SelfGradientInnerProduct(int i, int j)
	{
		if (i == j)
		{
			return 1.0 / (Value(i) - Value(i - 1)) + 1.0 / (Value(i + 1) - Value(i));
		}
		else if (abs(i - j) == 1)
		{
			return -1.0 / abs(Value(i) - Value(j));
		}
		else return 0;//For robustness
	}
};


class FEM1DVisualizer :public Visualizer
{
protected:

	void evaluate()
	{
		LinearBase1DFEM linear(5);
		linear.evaluate();
		for (int i = 0; i < Length; ++i)
		{
			ys3[i] = linear.value(1.0 / Length * i);
		}
	}

	void AddPoint();
	void DragPoint();
	void draw(bool* p_open) override;

	std::vector<Point2> points;
	bool updated = true;

public:
	FEM1DVisualizer() {
		Float h = 1.0 / Length;
		for (int i = 0; i < Length; ++i)
		{
			xs[i] = i * h;
			ys1[i] = (xs[i] - 1) * sin(xs[i]);
			ys2[i] = 2 - 2 * xs[i] + 2 * xs[i] * cos(1) - 2 * cos(xs[i]) + sin(xs[i]) - xs[i] * sin(xs[i]);
		}
		evaluate();
	}
	const float Length = 1001;

	std::vector<float> xs = std::vector<float>(Length);
	std::vector<float> ys1 = std::vector<float>(Length);
	std::vector<float> ys2 = std::vector<float>(Length);
	std::vector<float> ys3 = std::vector<float>(Length);
};

void FEM1DVisualizer::AddPoint()
{
	if (ImGui::IsMouseClicked(ImGuiMouseButton_Right))
	{
		auto pos = ImPlot::GetPlotMousePos();
		points.emplace_back(pos.x, pos.y);
		updated = true;
	}
}

void FEM1DVisualizer::DragPoint()
{
	for (int i = 0; i < points.size(); ++i)
	{
		auto& point = points[i];
		updated |= ImPlot::DragPoint(("Point " + std::to_string(i)).c_str(), &point.x(), &point.y());
	}
}

void FEM1DVisualizer::draw(bool* p_open)
{
	if (ImGui::BeginTabBar("Homework 1")) {
		if (ImGui::BeginTabItem("FEM1D"))
		{
			if (ImPlot::BeginPlot("Line Plot", "x", "f(x)", ImGui::GetContentRegionAvail(), ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMenus)) {
				ImPlot::PlotLine("u", &xs[0], &ys1[0], Length);
				ImPlot::PlotLine("Precise solution", &xs[0], &ys2[0], Length);
				ImPlot::PlotLine("FEM Result", &xs[0], &ys3[0], Length);

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