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

	virtual void draw()
	{
	}

private:
	bool show_demo_window;
	GLFWwindow* window;
	ImGui_ImplVulkanH_Window* wd;
	ImVec4 clear_color;
};
