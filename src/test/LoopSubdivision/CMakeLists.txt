Ubpa_GetTargetName(FEM "${PROJECT_SOURCE_DIR}/src/FEM")
Ubpa_GetTargetName(Tools "${PROJECT_SOURCE_DIR}/src/Tools")
Ubpa_GetTargetName(visualize "${PROJECT_SOURCE_DIR}/src/Visualization")

Ubpa_GetTargetName(project_ubpa ${CMAKE_CURRENT_SOURCE_DIR})

set(project ${project_ubpa})

set(folder ${CMAKE_CURRENT_SOURCE_DIR})
message("target at ${project}_shaders")

#donut_compile_shaders(
#    TARGET ${project}_shaders
#    CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/shaders.cfg
#    SOURCES ${shaders}
#    FOLDER ${folder}
#    DXIL ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/shaders/${project}/dxil
#    SPIRV_DXC ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/shaders/${project}/spirv
#)

donut_compile_shaders_all_platforms(
    TARGET ${project}_shaders
    CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/shaders.cfg
    FOLDER ${folder}
    OUTPUT_BASE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/shaders/${project}
)

Ubpa_AddTarget(
  TEST
  MODE EXE
  LIB donut_app donut_render donut_engine Eigen3::Eigen ${FEM} ${visualize} ${Tools} 
) 

add_dependencies(${project_ubpa} ${project}_shaders donut_shaders)
set_target_properties(${project_ubpa} PROPERTIES FOLDER ${folder})