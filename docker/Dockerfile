FROM ubuntu:focal
ENV DEBIAN_FRONTEND noninteractive
RUN apt -y update 
RUN apt -y upgrade
RUN apt install -y libboost1.71-dev \
                   libboost-atomic1.71-dev \
                   libboost-atomic1.71.0 \
                   libboost-chrono1.71-dev \
                   libboost-chrono1.71.0 \
                   libboost-date-time1.71-dev \
                   libboost-date-time1.71.0 \
                   libboost-program-options-dev \
                   libboost-program-options1.71-dev \
                   libboost-program-options1.71.0 \
                   libboost-serialization1.71-dev \
                   libboost-serialization1.71.0 \
                   libboost-system-dev \
                   libboost-system1.71-dev \
                   libboost-thread-dev \
                   libboost-thread1.71-dev \
                   build-essential libgmp-dev libmpfr-dev libvtk6-dev libtbb-dev cmake
COPY gmsh /GMSH
COPY cgal /CGAL
RUN cd /GMSH && mkdir build_lib && cd build_lib && cmake -DDEFAULT=0 -DENABLE_PARSER=1 -DENABLE_BUILD_DYNAMIC=1 .. && make -j && make install
RUN cd /CGAL && mkdir build && cd build/ && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j 12 && make install 
