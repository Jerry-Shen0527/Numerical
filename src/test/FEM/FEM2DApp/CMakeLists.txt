Ubpa_GetTargetName(FEM "${PROJECT_SOURCE_DIR}/src/FEM")
Ubpa_GetTargetName(Tools "${PROJECT_SOURCE_DIR}/src/Tools")
Ubpa_GetTargetName(visualize "${PROJECT_SOURCE_DIR}/src/Visualization")

Ubpa_AddTarget(
  TEST
  MODE EXE
  LIB Eigen3::Eigen ${FEM} ${visualize} ${Tools}
) 