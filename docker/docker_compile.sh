#!/bin/bash

BASE_DIR=/home/benny/Projects/tdcs-pipeline/meshheadmodel/
CONTAINER_NAME=cgal_and_gmsh_compiled_with_tbb_ubuntu2004

docker run                      \
	--rm                        \
	--name=compile_meshingtool  \
	-v $BASE_DIR:/shares        \
	$CONTAINER_NAME             \
	/bin/bash -c "ls -l /shares && cd /shares/ && cmake . && make VERBOSE=1"	
