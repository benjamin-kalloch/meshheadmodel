#!/bin/bash

BASE_DIR=/home/benny/Projects/tdcs-pipeline/meshheadmodel/
CONTAINER_NAME=cgal_and_gmsh_compiled_with_tbb_ubuntu2004
CGAL_LIBRARY_LOCATION=/usr/local/lib/
GMSH_IBRARY_LOCATION=/usr/local/lib/

docker run                  \
    --rm                    \
	--name=ubuntu           \
    -v ${BASE_DIR}:/shares  \
    $CONTAINER_NAME         \
    bash -c "cd /shares/deploy/ && export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CGAL_LIBRARY_LOCATION:$GMSH_IBRARY_LOCATION; ./runMeshingTool_lesions.sh"
