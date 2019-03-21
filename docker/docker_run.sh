#!/bin/bash

BASE_DIR=$HOME/Repositories/tdcs-pipeline/meshheadmodel
CONTAINER_NAME=cgal_ubuntu1904

docker run                  \
    --rm                    \
    --name=ubuntu           \
    -v ${BASE_DIR}:/shares  \
    $CONTAINER_NAME         \
    bash -c "cd /shares/deploy/ && ./runMeshingTool.sh"
