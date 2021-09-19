#include"Visualization/Visualizer.h"

class p_Visualizer :public Visualizer
{
protected:
	void draw(bool* p_open) override;
};

void p_Visualizer::draw(bool* p_open)
{
	if (ImGui::BeginTabBar("ImPlotDemoTabs")) {
		if (ImGui::BeginTabItem("Homework1"))
		{
			ImGui::EndTabItem();
		}

		ImGui::EndTabBar();
	}
	ImGui::End();
}

int main()
{
	p_Visualizer visualizer;
	visualizer.RenderLoop();
}