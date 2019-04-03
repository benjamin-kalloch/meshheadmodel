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

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake/Modules/")

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


# Defaults to 'activated', make sure the path to the 'UseTBB.cmake' configuration file is correct for you system
# Strangely, the 'FindTBB.cmake' config in '/usr/share/cmake/Modules' does not set the full path of the
# UseTBB.cmake file in ${TBB_USE_FILE}...
option(CGAL_ACTIVATE_CONCURRENT_MESH_3 "Activate parallelism in Mesh_3" ON)
if( CGAL_ACTIVATE_CONCURRENT_MESH_3 OR ENV{CGAL_ACTIVATE_CONCURRENT_MESH_3} )
    add_definitions( -DCGAL_CONCURRENT_MESH_3 )
    find_package( TBB REQUIRED )
    if( TBB_FOUND )
        include(/usr/lib/x86_64-linux-gnu/cmake/CGAL/UseTBB.cmake)
    endif()
endif()

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

target_link_libraries(meshHeadModel ${VTK_LIBRARIES} ${TBB_LIBRARIES} gmp mpfr CGAL CGAL_ImageIO Gmsh)
