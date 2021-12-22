message("Torch libraries: ${TORCH_LIBRARIES}")

Ubpa_AddTarget(
  TEST
  MODE EXE
  LIB ${TORCH_LIBRARIES}
) 