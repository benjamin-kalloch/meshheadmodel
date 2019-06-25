cmake_minimum_required(VERSION 2.6)

#
# NOTE: all paths are relative to the docker environment NOT the host system!
#		See: 'dockerrun.sh' for the entry point of the host system to the docker environment. 
#

project(meshHeadModel)
set (volmesh_VERSION_MAJOR 1)
set (volmesh_VERSION_MINOR 0)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/deploy)
set(GMSH_DIRECTORY /GMSH)
set (CMAKE_CXX_STANDARD 17)	# for std::remove to work properly

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/usr/share/cmake/Modules/")

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_library(GMSH_LIB gmsh ${GMSH_DIRECTORY}/build_lib/)
if(NOT GMSH_LIB)
    message(FATAL_ERROR "Could not find libgmsh!")
endif()

find_path(GMSH_INC gmsh.h ${GMSH_DIRECTORY}/api/)
if(NOT GMSH_INC)
    message(FATAL_ERROR "Could not find gmsh.h")
endif()

include_directories(${GMSH_INC})

# Defaults to 'activated', make sure the path to the 'UseTBB.cmake' configuration file is correct for your system
# Strangely, the 'FindTBB.cmake' config in '/usr/share/cmake/Modules' does not set the full path of the
# UseTBB.cmake file in ${TBB_USE_FILE}...
option(CGAL_ACTIVATE_CONCURRENT_MESH_3 "Activate parallelism in Mesh_3" ON)
if( CGAL_ACTIVATE_CONCURRENT_MESH_3 OR ENV{CGAL_ACTIVATE_CONCURRENT_MESH_3} )
    add_definitions( -DCGAL_CONCURRENT_MESH_3 )
    find_package( TBB REQUIRED )
    if( TBB_FOUND )
        include(/usr/local/lib/cmake/CGAL/UseTBB.cmake)
    endif()
endif()

add_executable(meshHeadModel mesh_head_model.cpp)

target_link_libraries(meshHeadModel ${VTK_LIBRARIES} ${TBB_LIBRARIES} gmp mpfr CGAL CGAL_ImageIO ${GMSH_LIB})
