Ubpa_GetTargetName(imgui_lib "${PROJECT_SOURCE_DIR}/src/Imgui")

Ubpa_AddTarget(
  TEST
  MODE STATIC
  LIB ${imgui_lib} glfw
) 