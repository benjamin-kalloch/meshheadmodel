#!/bin/bash

BASE_DIR=$HOME/Repositories/tdcs-pipeline/meshheadmodel/
CONTAINER_NAME=cgal_and_gmsh_compiled_with_tbb_ubuntu1904

docker run                      \
	--rm                        \
	--name=compile_meshingtool  \
	-v $BASE_DIR:/shares        \
	$CONTAINER_NAME             \
	/bin/bash -c "ls -l / && ls -l /usr/local/lib && cd /shares/ && cmake . && make VERBOSE=1"	
