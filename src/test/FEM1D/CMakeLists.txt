Ubpa_GetTargetName(FEM "${PROJECT_SOURCE_DIR}/src/FEM")

Ubpa_AddTarget(
  TEST
  MODE EXE
  LIB Eigen3::Eigen ${FEM}
) 