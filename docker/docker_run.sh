#!/bin/bash

BASE_DIR=$HOME/Repositories/tdcs-pipeline/meshheadmodel/

CONTAINER_NAME=cgal_tbb_ubuntu1904

docker run                  \
    --rm                    \
	--name=ubuntu           \
    -v ${BASE_DIR}:/shares  \
    $CONTAINER_NAME         \
    bash -c "apt purge libcgal-dev libcgal13 -y && cd /shares/deploy/CGAL-4.13.1/build/ && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j 12 && make install && cd /shares/deploy/ && export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/; ./runMeshingTool.sh"
