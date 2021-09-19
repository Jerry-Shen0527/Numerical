#pragma once

#include "imgui/imgui_impl_vulkan.h"

class GLFWwindow;

class Visualizer
{
public:
	Visualizer();
	virtual ~Visualizer();

	virtual void Init() { initWindow(); };
	virtual void RenderLoop();

protected:
	void initWindow();
	void cleanUp();
	void drawDefault(bool* p_open);


	virtual void draw(bool* p_open) = 0;

private:
	bool show;
	GLFWwindow* window;
	ImGui_ImplVulkanH_Window* wd;
	ImVec4 clear_color;
};
