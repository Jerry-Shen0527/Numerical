

message("vulkan include files at ${Vulkan_INCLUDE_DIRS}")

Ubpa_AddTarget(
  MODE STATIC
  LIB ${Vulkan_LIBRARIES} glfw
)