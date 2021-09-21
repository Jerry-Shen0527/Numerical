#include <iostream>
#include <type.hpp>

#include <Eigen/Eigen>

#include "FEM/FEM.hpp"
#include "imgui/implot.h"
#include "Visualization/Visualizer.h"

static Float LinearInnerProduct(Float xmin, Float x_mid, Float xmax)
{
	return ((-2 + x_mid * x_mid + xmax - x_mid * (1 + xmax)) * cos(x_mid) + 2 * cos(xmax) + (1 - 2 * x_mid + xmax) * sin(x_mid) + (-1 + xmax) * sin(xmax))
		/ (x_mid - xmax) + 1 / (x_mid - xmin) * ((2 + x_mid - x_mid * x_mid + (-1 + x_mid) * xmin) * cos(x_mid) - 2 * cos(xmin) - (1 - 2 * x_mid + xmin) * sin(x_mid) + sin(xmin) - xmin * sin(xmin));
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

class LinearBaseFEM1D :public StaticFEM1D
{
public:
	LinearBaseFEM1D(size_t size)
	{
		this->size = size;
		mesh.resize(size);
		//A default devision is provided
		Float h = 1.0 / (size + 1);
		for (int i = 0; i < size; ++i)
		{
			mesh[i] = (i + 1) * h;
		}
	}
	std::vector<Float> mesh;

	Float Value(Float x) override
	{
		Float ret = 0;
		for (int i = 0; i < size; ++i)
		{
			ret += LinearFunc(x, MeshValue(i), MeshValue(i - 1), MeshValue(i + 1)) * coeff_()(i);
		}
		return ret;
	}

private:

	//return the indices of the functions related to the current function
	std::vector<int> RelatedFuncIdx(int idx) override
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

	Float MeshValue(int index)
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

	Float GradientInnerProduct(int i, int j) override
	{
		if (i == j)
		{
			return 1.0 / (MeshValue(i) - MeshValue(i - 1)) + 1.0 / (MeshValue(i + 1) - MeshValue(i));
		}
		else if (abs(i - j) == 1)
		{
			return -1.0 / abs(MeshValue(i) - MeshValue(j));
		}
		else return 0;//For robustness
	}

	Float RHSInnerProduct(int i) override
	{
		return LinearInnerProduct(MeshValue(i - 1), MeshValue(i), MeshValue(i + 1));
	}
};

static Float Quadratic(Float x, Float x_left, Float x_middle, Float x_right)
{
	if (x<x_left || x>x_right)
	{
		return 0;
	}
	return (x - x_left) * (x - x_right) / ((x_middle - x_left) * (x_middle - x_right));
}

static Float QuadraticInnerProduct(Float xmin, Float xmid, Float xmax)
{
	return (2 * (1 - 2 * xmax + xmin) * cos(xmax) - 2 * (1 + xmax - 2 * xmin) * cos(xmin) + 6 * (sin(xmax) - sin(xmin)) - (xmax - xmin) * ((-1 + xmax) * sin(xmax) + (-1 + xmin) * sin(xmin))) / ((xmax - xmid) * (xmid - xmin));
}

static Float QuadraticGradientInnerProduct(Float xmin, Float xmid, Float xmax, Float xmin2, Float xmid2, Float xmax2)
{
	auto left = std::max(xmin, xmin2);
	auto right = std::min(xmax, xmax2);

	if (left > right)
	{
		return  0;
	}
	auto PrimitiveFunction = [=](Float x) {return ((4 * x * x * x) / 3. + x * (xmax + xmin) * (xmax2 + xmin2) - x * x * (xmax + xmax2 + xmin + xmin2)) / ((xmax - xmid) * (xmax2 - xmid2) * (xmid - xmin) * (xmid2 - xmin2)); };
	return PrimitiveFunction(right) - PrimitiveFunction(left);
}

class QuadraticBaseFEM1D :public StaticFEM1D
{
public:
	QuadraticBaseFEM1D(size_t size) :mesh_size(size)
	{
		this->size = 2 * size + 1;
		mesh.resize(size);
		//A default division is provided
		Float h = 1.0 / (size + 1);
		for (int i = 0; i < size; ++i)
		{
			mesh[i] = (i + 1) * h;
		}
	}
	std::vector<Float> mesh;
	size_t mesh_size;

	Float Value(Float x) override
	{
		Float ret = 0;
		for (int i = 0; i < size; ++i)
		{
			Float xmin, xmid, xmax;
			idx_to_mesh(i, xmin, xmid, xmax);

			ret += Quadratic(x, xmin, xmid, xmax) * coeff_()(i);
		}
		return ret;
	}

private:

	//return the indices of the functions related to the current function
	std::vector<int> RelatedFuncIdx(int idx) override
	{
		std::vector<int> ret;
		if (idx < mesh_size)
		{
			ret.push_back(idx);
			int left = idx - 1, right = idx + 1;
			if (0 <= left && left < mesh_size)
			{
				ret.push_back(left);
			}
			if (0 <= right && right < mesh_size)
			{
				ret.push_back(right);
			}
			ret.push_back(idx + mesh_size);

			ret.push_back(idx + mesh_size + 1);
		}
		else
		{
			ret.push_back(idx);
			int left = idx - mesh_size - 1, right = idx - mesh_size;

			if (0 <= left && left < mesh_size)
			{
				ret.push_back(left);
			}
			if (0 <= right && right < mesh_size)
			{
				ret.push_back(right);
			}
		}
		return ret;
	}

	Float MeshValue(int index)
	{
		if (index < 0)
		{
			return 0;
		}
		if (index >= mesh_size)
		{
			return  1.0;
		}
		return mesh[index];
	}

	void idx_to_mesh(int idx, Float& xmin, Float& xmid, Float& xmax)
	{
		if (idx < mesh_size)
		{
			xmin = MeshValue(idx - 1);	xmid = MeshValue(idx);	xmax = MeshValue(idx + 1);
		}
		else
		{
			xmin = MeshValue(idx - 1 - mesh_size); xmid = (MeshValue(idx - mesh_size) + MeshValue(idx - 1 - mesh_size)) / 2.0; xmax = MeshValue(idx - mesh_size);
		}
	}

	Float GradientInnerProduct(int i, int j) override
	{
		Float xmin, xmid, xmax;
		Float xmin2, xmid2, xmax2;

		idx_to_mesh(i, xmin, xmid, xmax);
		idx_to_mesh(j, xmin2, xmid2, xmax2);

		return QuadraticGradientInnerProduct(xmin, xmid, xmax, xmin2, xmid2, xmax2);
	}

	Float RHSInnerProduct(int i) override
	{
		Float xmin, xmid, xmax;

		idx_to_mesh(i, xmin, xmid, xmax);
		return QuadraticInnerProduct(xmin, xmid, xmax);
	}
};

class FEM1DVisualizer :public Visualizer
{
protected:

	void evaluate()
	{
		int sizeeee = 64;
		LinearBaseFEM1D linear(sizeeee);
		QuadraticBaseFEM1D quadratic(sizeeee / 2);
		linear.evaluate();
		quadratic.evaluate();
		for (int i = 0; i < Length; ++i)
		{
			if (abs(ys2[i]) > 1e-5)
			{
				ys3[i] = (quadratic.Value(1.0 / (Length - 1) * i) - ys2[i])/ ys2[i];
				ys4[i] = (linear.Value(1.0 / (Length - 1) * i) - ys2[i])/ ys2[i];
			}
		}
	}

	void AddPoint();
	void DragPoint();
	void draw(bool* p_open) override;

	std::vector<Point2> points;
	bool updated = true;

public:
	FEM1DVisualizer() {
		Float h = 1.0 / (Length - 1);
		for (int i = 0; i < Length; ++i)
		{
			xs[i] = i * h;
			ys1[i] = (xs[i] - 1) * sin(xs[i]);
			ys2[i] = -(2 - 2 * xs[i] + 2 * xs[i] * cos(1) - 2 * cos(xs[i]) + sin(xs[i]) - xs[i] * sin(xs[i]));
		}
		evaluate();
	}
	const size_t Length = 2001;

	std::vector<float> xs = std::vector<float>(Length);
	std::vector<float> ys1 = std::vector<float>(Length);
	std::vector<float> ys2 = std::vector<float>(Length);
	std::vector<float> ys3 = std::vector<float>(Length);
	std::vector<float> ys4 = std::vector<float>(Length);
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
				//ImPlot::PlotLine("u", &xs[0], &ys1[0], Length);
				//ImPlot::PlotLine("Precise solution", &xs[0], &ys2[0], Length);
				//ImPlot::PlotLine("FEM Result", &xs[0], &ys3[0], Length);
				ImPlot::PlotLine("FEM Diff Quadratic", &xs[0], &ys3[0], Length);
				ImPlot::PlotLine("FEM Diff Linear", &xs[0], &ys4[0], Length);

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