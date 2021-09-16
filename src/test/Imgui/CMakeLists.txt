Ubpa_GetTargetName(imgui_lib "${PROJECT_SOURCE_DIR}/src/Imgui")

Ubpa_AddTarget(
  TEST  
  MODE EXE
  LIB ${imgui_lib} glfw
) 
