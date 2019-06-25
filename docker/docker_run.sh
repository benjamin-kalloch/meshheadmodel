#!/bin/bash

BASE_DIR=$HOME/Repositories/tdcs-pipeline/meshheadmodel/
CONTAINER_NAME=cgal_and_gmsh_compiled_with_tbb_ubuntu1904
CGAL_LIBRARY_LOCATION=/usr/local/lib/

docker run                  \
    --rm                    \
	--name=ubuntu           \
    -v ${BASE_DIR}:/shares  \
    $CONTAINER_NAME         \
    bash -c "cd /shares/deploy/ && export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CGAL_LIBRARY_LOCATION; ./runMeshingTool.sh"
