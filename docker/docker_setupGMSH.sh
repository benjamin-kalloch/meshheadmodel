#!/bin/bash

BASE_DIR=$HOME/Repositories/tdcs-pipeline/meshheadmodel/

CONTAINER_NAME=cgal_ubuntu1904

docker run               \
	--rm                 \
	--name=ubuntu        \
	-v $BASE_DIR:/shares \
	$CONTAINER_NAME      \
	bash -c "ls -l /shares && cd /shares/deploy/gmsh-4.3.0-source/ && mkdir build_lib && cd build_lib && cmake -DDEFAULT=0 -DENABLE_PARSER=1 -DENABLE_BUILD_DYNAMIC=1 .. && make -j"	
