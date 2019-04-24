# MeshHeadModel

This application uses the CGAL API to generate a tetrahedral volume mesh from a 3D label image (ANALYZE file format) and optionally surface files (OFF file format) of the outer shell (i.e. the skin compartment) and the electrodes attached to the skin. The application is detailed in in the paper **A flexible open-source pipeline for simulating transcranial electric stimulation** by *Benjamin Kalloch, Pierre-Louis Bazin, Arno Villringer, Bernhard Sehm, and Mario Hlawitschka*.

## How to compile
You need the following prerequisite
- the CGAL development library v4.13
- On Ubuntu 19.04 the following packages are required: *build-essential libcgal-dev libgmp-dev libmpfr-dev libvtk6-dev cmake* (this will also install the CGAL development files of the required version)
- Download the GMSH 3.0.6 source code [3]

After ensuring that all prerequisites are met follow this procedure
1) Extract the gmsh source code and copy the entire directory into the deploy directory of this repository.
2) Compile GMSH as a dynamic link library:
    - After you have copied the source-code folder ('gmsh-3.0.6-source') into the 'deploy' directory, create a folder 'build_lib' in that source directory.
    - A minimal configuration can be crated using cmake `cmake -DDEFAULT=0 -DENABLE_PARSER=1 -DENABLE_BUILD_DYNAMIC=1 ..`
    - Run `make -j` to build the successfully configured GMSH
    - *Note* that any deviations in terms of location or naming scheme of the folder of the GMSH source code and the location of the built library necessitate an adaptation of the 'CMakeLists.txt' file. 
2) Issue the command `cmake .`.
3) Finally, run `make`to build the tool.

## How to run
We provide a bash script *runMeshingTool.sh* in the deploy directory with exemplary input parameters from the almi5 test case. You can adjust the parameters according to your input files.
The following parameters are supported by our meshing tool ([!]=mandatory parameter, [*]=optional parameter):
- **_imagefile_** = path to the segmented 3D label image of the subject's head, (in ANALYZE file format) [*]
- **_outfile_** = location where the generated volume mesh should be stored. [!]
- **_electrodefiles_** = comma-separated list surface descriptions files (in OFF file format) of the electrodes [*]
- **_tissuesurfaces_** = comma-separated list of paths to the surface descriptions of additional tissue sturctures (e.g. the smooth skin) of the head model. When providing tissue structures as a surface description make sure that the 3D label image does not contain the label for these structures anymore. [*]
- **_f_angbound_**, **_f_distbound_**, **_f_sizebound_**, **_c_reratio_**, **_c_sizebound_** = Metrics to control the mesh quality and element size  (refer to [4] for a description of these measures; the sample runMeshingTool file contains values which produced a good mesh quality) [*]
- **_lloyd_** = Enable the Lloyd global optimizer. [*]
- **_odt_** = Enable the global optimized Delaunay triangulation optimizer. [*]
- **_perturb_** = Enable local mesh optimization by mesh vertex perturbation. [*]
- **_exude_** = Enable local mesh optimization by exuding slivers. [*]
- **_help_** = display the help text [*]

## Secondary information
To facilitate the setup of the proper environment to built the MeshHeadModel tool, we provide a Dockerfile to create a container based on Ubuntu 18.10 which contains all required libraries. 

Simply build the container with `docker build -t cgal_ubuntu1904 .` in the docker directory of this repository.

Before compiling the meshing tool, GMSH must be compiled as a shared library. We prepared the bash script *docker_setupGMSH.sh* for this job.
Note that this script expects the GMSH source code *gmsh-3.0.6-source* to be in the *deploy* directory.
After the setup of GMSH, you may compile the meshing tool by running the script *docker_compile.sh*.

When compiled in the container, the tool also must be launched from inside the container. You can do this by calling the *docker_run.sh* script which again invokes the *runMeshingTool.sh* script as mentioned before from inside the container. When using docker to run the tool consider that all paths within *runMeshingtTool.sh* (both input files and output files) must be specified according to the docker environment NOT according to the host machine. The best practice is to use the predefined input & output directories and to use paths relative to the MeshHeadModel executable. 

! Note that for all three docker scripts you may have to adjust the variables $BASE_DIR and $CONTAINER_NAME if you use a path/name other then predefined in the file.

###### Useful links
[1] The CGAL project: https://doc.cgal.org/4.13/Manual/index.html

[2] THe GMSH project: http://gmsh.info/

[3] GMSH 3.0.6 source code: http://gmsh.info/src/gmsh-3.0.6-source.tgz

[4] Explanation of the CGAL quality criteria: https://doc.cgal.org/4.13/Mesh_3/index.html#title11
