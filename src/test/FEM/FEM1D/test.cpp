#include <iostream>
#include <numeric>
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
	LinearBaseFEM1D(size_t size) : StaticFEM1D(), size(size)
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

	Float Value(Float x) override
	{
		Float ret = 0;
		for (int i = 0; i < mat_size; ++i)
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
		if (0 <= left && left < mat_size)
		{
			ret.push_back(left);
		}
		if (0 <= right && right < mat_size)
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
		if (index >= mat_size)
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

protected:

	int size;

	Float SelfInnerProduct(int i, int j) override
	{
		return 0;
	}

	Float GradientSelfInnerProduct(int i, int j) override
	{
		return 0;
	}

	void SetMatSize() override
	{
		mat_size = size;
	}
};

static Float Quadratic(Float x, Float x_left, Float x_middle, Float x_right)
{
	if (x <= x_left || x >= x_right)
	{
		return 0.0;
	}

	return   (x - x_left) * (x - x_right) / (x_middle - x_left) / (x_middle - x_right);
}

static Float QuadraticInnerProduct(Float xmin, Float xmid, Float xmax)
{
	//return (2. * (1. - 2. * xmax + xmin) * cos(xmax) - 2. * (1. + xmax - 2. * xmin) * cos(xmin) + 6. * (sin(xmax) - sin(xmin)) - (xmax - xmin) * (xmax*sin(xmax) - sin(xmax) + xmin*sin(xmin) - sin(xmin))) / (xmax - xmid) / (xmid - xmin);
	return  (2 * cos(xmax) - 4 * xmax * cos(xmax) + 2 * xmin * cos(xmax) - 2 * cos(xmin) - 2 * xmax * cos(xmin) +
		4 * xmin * cos(xmin) + 6 * sin(xmax) + xmax * sin(xmax) - (xmax * xmax) * sin(xmax) -
		xmin * sin(xmax) + xmax * xmin * sin(xmax) - 6 * sin(xmin) + xmax * sin(xmin) -
		xmin * sin(xmin) - xmax * xmin * sin(xmin) + (xmin * xmin) * sin(xmin)) /
		((xmax - xmid) * (xmid - xmin));
}

static Float QuadraticGradientInnerProduct(Float xmin, Float xmid, Float xmax, Float xmin2, Float xmid2, Float xmax2)
{
	auto left = std::max(xmin, xmin2);
	auto right = std::min(xmax, xmax2);

	if (left > right)
	{
		return  0;
	}
	Float ret = -((left - right) * (4 * (left * left) + 4 * left * right + 4 * (right * right) - 3 * left * xmax -
		3 * right * xmax - 3 * left * xmax2 - 3 * right * xmax2 + 3 * xmax * xmax2 - 3 * left * xmin -
		3 * right * xmin + 3 * xmax2 * xmin - 3 * left * xmin2 - 3 * right * xmin2 + 3 * xmax * xmin2 +
		3 * xmin * xmin2)) / 3. / (xmax - xmid) / (xmax2 - xmid2) / (xmid - xmin) / (xmid2 - xmin2);

	return ret;
}

class QuadraticBaseFEM1D :public StaticFEM1D
{
public:
	QuadraticBaseFEM1D(size_t size) : StaticFEM1D(), mesh_size(size), mat_size_(2 * size + 1)
	{
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
	size_t mat_size_;

	Float Value(Float x) override
	{
		Float ret = 0;
		for (int i = 0; i < mat_size; ++i)
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
			xmin = MeshValue(idx - 1 - mesh_size);
			xmid = (MeshValue(idx - mesh_size) + MeshValue(idx - 1 - mesh_size)) / 2.0;
			xmax = MeshValue(idx - mesh_size);
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

protected:
	Float SelfInnerProduct(int i, int j) override
	{
		return 0;
	}

	Float GradientSelfInnerProduct(int i, int j) override
	{
		return 0;
	}

	void SetMatSize() override
	{
		mat_size = mat_size_;
	}
};

class FEM1DVisualizer :public Visualizer
{
protected:

	void evaluate()
	{
		int segement_ = segemnt;
		LinearBaseFEM1D linear(segement_);
		QuadraticBaseFEM1D quadratic(segement_ / 2);
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

	void AddPoint();
	void DragPoint();
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
			precise_val[i] = -(2 - 2 * xs[i] + 2 * xs[i] * cos(1) - 2 * cos(xs[i]) + sin(xs[i]) - xs[i] * sin(xs[i]));
		}
		evaluate();

		Float L1, L2, L_inf;

		segemnt = 1;
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
		segemnt = 1;
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
static inline ImVec2 operator-(const ImVec2& lhs, const ImVec2& rhs) { return ImVec2(lhs.x - rhs.x, lhs.y - rhs.y); }
void FEM1DVisualizer::draw(bool* p_open)
{
	if (ImGui::BeginTabBar("Homework 1")) {
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
			if (ImGui::SliderInt("Number of segments", &segemnt, 1, 512))
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
			if (ImGui::SliderInt("Number of segments", &segemnt, 1, 512))
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
			if (ImGui::SliderInt("Number of segments", &segemnt, 1, 100))
			{
				evaluate();
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