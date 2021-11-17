#include <iostream>
#include <numeric>
#include <type.hpp>

#include <Eigen/Eigen>

#include "FEM/FEM1DApp.hpp"
#include "FEM/FEM2DApp.hpp"
#include "imgui/implot.h"
#include "Visualization/Visualizer.h"

#include <numeric>

using Linear = PolynomialFEMApp<1>;
using Quadratic = PolynomialFEMApp<2>;

class FEM2DVisualizer :public Visualizer
{
protected:

	StaticFEM2DAppP1 app;
protected:

	void evaluate()
	{
		app.evaluate();
	}

	void Control_UI();
	void draw(bool* p_open) override;

	bool updated = true;

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
	}

	FEM2DVisualizer(const StaticFEM2DApp::FEM2DMesh& mesh) :app(mesh) {
		CalcAccurateRst();
		//segemnt = 16;
		//Float L1, L2, L_inf;
		//do
		//{
		//	evaluate();

		//	error(precise_val, linear_val, L1, L2, L_inf);

		//	using std::cout;
		//	using std::endl;

		//	//if (segemnt == 16)
		//	//{
		//	//	cout << segemnt << '&' << L1 << '&' << '-' << '&' << L2 << '&' << '-' << '&' << L_inf << '&' << '-' << "\\\\" << endl;
		//	//}
		//	//else
		//	//	cout << segemnt << '&' << L1 << '&' <<- log2(L1 / linear_L1.back()) << '&' << L2 << '&' << -log2(L2 / linear_L2.back()) << '&' << L_inf << '&' << -log2(L_inf / linear_Linf.back()) << "\\\\" << endl;

		//	pointcount.push_back(segemnt);
		//	linear_L1.push_back(L1);
		//	linear_L2.push_back(L2);
		//	linear_Linf.push_back(L_inf);
		//	error(precise_val, quadratic_val, L1, L2, L_inf);

		//	if (segemnt == 16)
		//	{
		//		cout << segemnt << '&' << L1 << '&' << '-' << '&' << L2 << '&' << '-' << '&' << L_inf << '&' << '-' << "\\\\" << endl;
		//	}
		//	else
		//		cout << segemnt << '&' << L1 << '&' << -log2(L1 / quadratic_L1.back()) << '&' << L2 << '&' << -log2(L2 / quadratic_L2.back()) << '&' << L_inf << '&' << -log2(L_inf / quadratic_Linf.back()) << "\\\\" << endl;

		//	quadratic_L1.push_back(L1);
		//	quadratic_L2.push_back(L2);
		//	quadratic_Linf.push_back(L_inf);

		//	segemnt *= 2;
		//} while (segemnt != 4096);
		//segemnt = 16;

		evaluate();
	}
	const size_t Length = 10001;

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

void FEM2DVisualizer::Control_UI()
{
	//if (ImGui::SliderInt("Number of segments", &segemnt, 2, 200))
	//{
	//	segemnt = segemnt < 2 ? 2 : segemnt;
	//	evaluate();
	//}
	evaluate();
	CalcAccurateRst();
}

void FEM2DVisualizer::draw(bool* p_open)
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
#include <OpenMesh/Core/IO/MeshIO.hh>
int main()
{
	using Mesh = StaticFEM2DApp::FEM2DMesh;
	StaticFEM2DApp::FEM2DMesh mesh;
	std::cout << OpenMesh::IO::read_mesh(mesh, "gd0.obj");

	FEM2DVisualizer visualizer(mesh);
	visualizer.RenderLoop();
}