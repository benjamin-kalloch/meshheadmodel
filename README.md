# MeshHeadModel

This application uses the CGAL API to generate a tetrahedral volume mesh from a 3D label image (ANALYZE file format) and optionally surface files (OFF file format) of the outer shell (i.e. the skin compartment) and the electrodes attached to the skin. The application will be detailed in in the paper **A flexible open-source pipeline for simulating transcranial electric stimulation** by *Benjamin Kalloch, Pierre-Louis Bazin, Arno Villringer, Bernhard Sehm, and Mario Hlawitschka*.

## How to compile
You need the following prerequisite
- On Ubuntu 19.04 the following packages are required: *build-essential libgmp-dev libmpfr-dev libvtk6-dev cmake libboost-atomic1.67-dev libboost-atomic1.67.0 libboost-chrono1.67-dev libboost-chrono1.67.0 libboost-date-time1.67-dev libboost-date-time1.67.0 libboost-program-options-dev libboost-program-options1.67-dev libboost-program-options1.67.0 libboost-serialization1.67-dev libboost-serialization1.67.0 libboost-system-dev libboost-system1.67-dev libboost-thread-dev* (this will also install the dependencies for CGAL and GMSH)
- Download the GMSH 4.3.0 source code [3]
- Download the CGAL 4.13.1 source code [4]
*Note* While Ubuntu 19.04 contains CGAL v.4.13.1 as well, the version included in the Ubuntu repositories is lacking a fix concerning the mesh export causing arbitrary normal orientation for boundary faces in the exported mesh. Therefore, we must compile CGAL from the sources. A more convenient approach is explained in the second part of this README utilizing a Docker.


After ensuring that all prerequisites are met follow this procedure
1) Extract the gmsh source code and copy the entire directory into the deploy directory of this repository.
2) Compile GMSH as a dynamic link library:
    - After you have copied the source-code folder ('gmsh-4.3.0-source') into the 'deploy' directory, create a folder 'build_lib' in that source directory.
    - A minimal configuration can be crated using cmake `cmake -DDEFAULT=0 -DENABLE_PARSER=1 -DENABLE_BUILD_DYNAMIC=1 ..`
    - Run `make -j` to build the successfully configured GMSH
3) Extract the CGAL source code at any location.
4) Compile CGAL and install it globally:
    - Create the folder 'build' in the source directory ('CGAL-4.13.1').
    - Run `cmake -DCMAKE_BUILD_TYPE=Release ..` to create a minimal configuration.
    - Run `make -j 8 && make install` to build the successfully configured CGAL.
5) Issue the command `cmake .` in the base-directory of this repository (meshheadmodel).
6) Finally, run `make`to build the tool.
*Note* You must adapt the following paths in the 'CMakeLists.txt' file: `set(GMSH_DIRECTORY /GMSH)`

## How to run
We provide a bash script *runMeshingTool.sh* in the deploy directory with exemplary input parameters from the almi5 test case. You can adjust the parameters according to your input files.
The following parameters are supported by our meshing tool ([!]=mandatory parameter, [*]=optional parameter):
- **_imagefile_** = path to the segmented 3D label image of the subject's head, (in ANALYZE file format) [*]
- **_outfile_** = location where the generated volume mesh should be stored. [!]
- **_electrodefiles_** = comma-separated list surface descriptions files (in OFF file format) of the electrodes [*]
- **_tissuesurfaces_** = comma-separated list of paths to the surface descriptions of additional tissue sturctures (e.g. the smooth skin) of the head model. When providing tissue structures as a surface description make sure that the 3D label image does not contain the label for these structures anymore. [*]
- **_f_angbound_**, **_f_distbound_**, **_f_sizebound_**, **_c_reratio_**, **_c_sizebound_** = Metrics to control the mesh quality and element size  (refer to [5] for a description of these measures; the sample runMeshingTool file contains values which produced a good mesh quality) [*]
- **_lloyd_** = Enable the Lloyd global optimizer. [*]
- **_odt_** = Enable the global optimized Delaunay triangulation optimizer. [*]
- **_perturb_** = Enable local mesh optimization by mesh vertex perturbation. [*]
- **_exude_** = Enable local mesh optimization by exuding slivers. [*]
- **_simnibs_** = Export mesh in a SimNIBS compatible output file format. [*]
- **_help_** = display the help text [*]
*Note* If you did not install CGAL to the default location `/usr/local/lib/`, you must adapt the  variable 'CGAL_LIBRARY_LOCATION' in the runMeshingTool-script accordingly.

## Secondary information
To facilitate the setup of the proper environment to built the MeshHeadModel tool, we provide a Dockerfile to create a container based on Ubuntu 19.04 which contains all required libraries. 

Simply build the container with `docker build -t cgal_and_gmsh_compiled_with_tbb_ubuntu1904 .` in the docker directory of this repository.

Before compiling the meshing tool, you must copy the source code directories of GMSH and CGAL into the Docker directory and rename them to 'gmsh' and 'cgal' respectively.
They are compiled during the image-creation process and included in the final docker image. Therefore, you may remove these source-code directories after the image has been built.

When compiled in the container, the tool also must be launched from inside the container. You can do this by calling the *docker_run.sh* script which again invokes the *runMeshingTool.sh* script as mentioned before from inside the container. When using docker to run the tool consider that all paths within *runMeshingtTool.sh* (both input files and output files) must be specified according to the docker environment NOT according to the host machine. The best practice is to use the predefined input & output directories and to use paths relative to the MeshHeadModel executable. 

! Note that for both docker scripts you may have to adjust the variables $BASE_DIR and $CONTAINER_NAME if you use a path/name other then predefined in the file.

###### Useful links
[1] The CGAL project: https://doc.cgal.org/4.13/Manual/index.html

[2] THe GMSH project: http://gmsh.info/

[3] GMSH 4.3.0 source code: http://gmsh.info/bin/Linux/gmsh-4.3.0-Linux64.tgz

[4] CGAL 4.13.1 source code: https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.13.1/CGAL-4.13.1.tar.xz

[5] Explanation of the CGAL quality criteria: https://doc.cgal.org/4.13/Mesh_3/index.html#title11
