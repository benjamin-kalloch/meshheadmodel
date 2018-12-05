cmake_minimum_required(VERSION 2.6)

#
# NOTE: all paths are relative to the docker environment NOT the host system!
#		See: 'dockerrun.sh' for the entry point of the host system to the docker environment. 
#

project(meshHeadModel)
set (volmesh_VERSION_MAJOR 1)
set (volmesh_VERSION_MINOR 0)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/deploy)
set(GMSH_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gmsh-3.0.6-source/)
set (CMAKE_CXX_STANDARD 17)	# for std::remove to work properly

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

link_directories(${GMSH_DIRECTORY}/build_lib/)
include_directories(
 	${GMSH_DIRECTORY}/Common/
	${GMSH_DIRECTORY}/build_lib/Common/
	${GMSH_DIRECTORY}/Geo/
	${GMSH_DIRECTORY}/Numeric/
	${GMSH_DIRECTORY}/Mesh/
	${GMSH_DIRECTORY}/Post/
)

add_executable(meshHeadModel mesh_head_model.cpp)

target_link_libraries(meshHeadModel ${VTK_LIBRARIES} gmp mpfr CGAL CGAL_ImageIO Gmsh)
