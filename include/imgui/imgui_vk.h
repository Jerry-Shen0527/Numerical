#pragma once
#include "imgui/imgui_impl_vulkan.h"
#pragma once

#include <vulkan/vulkan.h>

struct VulkanContext
{
	static VkAllocationCallbacks* g_Allocator;
	static VkInstance               g_Instance;
	static VkPhysicalDevice         g_PhysicalDevice;
	static VkDevice                 g_Device;
	static uint32_t                 g_QueueFamily;
	static VkQueue                  g_Queue;
	static VkDebugReportCallbackEXT g_DebugReport;
	static VkPipelineCache          g_PipelineCache;
	static VkDescriptorPool         g_DescriptorPool;

	static ImGui_ImplVulkanH_Window g_MainWindowData;
	static int                      g_MinImageCount;
	static bool                     g_SwapChainRebuild;
};

void check_vk_result(VkResult err);
void SetupVulkan(const char** extensions, uint32_t extensions_count);
void glfw_error_callback(int error, const char* description);
void SetupVulkanWindow(ImGui_ImplVulkanH_Window* wd, VkSurfaceKHR surface, int width, int height);
void FrameRender(ImGui_ImplVulkanH_Window* wd, ImDrawData* draw_data);
void FramePresent(ImGui_ImplVulkanH_Window* wd);
void CleanupVulkan();
void CleanupVulkanWindow();


