# MeshHeadModel

This application uses the CGAL API to generate a tetrahedral volume mesh from a 3D label image (ANALYZE file format) and optionally surface files (OFF file format) of the outer shell (i.e. the skin compartment) and the electrodes attached to the skin. The application is detailed in in the paper **A flexible open-source pipeline for simulating transcranial electric stimulation** by *Benjamin Kalloch, Pierre-Louis Bazin, Arno Villringer, Bernhard Sehm, and Mario Hlawitschka*.

## How to compile
You need the following prerequisite
- the CGAL development library v4.12
- On Ubuntu 18.10 the following packages are required: *build-essential libcgal-dev libgmp-dev libmpfr-dev libvtk6-dev cmake* (this will also install the CGAL development files of the required version)
- Download the GMSH 3.0.6 source code [3]

After ensuring are prerequisits are met follow this procedure
1) Extract the gmsh source code and copy the entire directory into the deploy directory of this repository.
2) Issue the command `cmake .`.
3) Finally, run `make`to build the tool.

## How to run
We provide a bash script *runMeshingTool.sh* in the deploy directory with exemplary input parameters from the almi5 test case. You can adjust the parameters according to your input files.
The following parameters are supported by our meshing tool ([!]=mandatory parameter, [*]=optional parameter):
- **_imagefile_** = path to the segmented 3D label image of the subject's head, (in ANALYZE file format) [!]
- **_outfile_** = location where the generated volume mesh should be stored. [!]
- **_electrodefiles_** = comma-separated list surface descriptions files (in OFF file format) of the electrodes [*]
- **_smoothskin_** = location of the surface description of the outer skin surface of the head model, if providing the skin as an surface description make sure that the 3D label image does not contain the label for skin anymore. [*]
- **_f_angbound_**, **_f_distbound_**, **_f_sizebound_**, **_c_reratio_**, **_c_sizebound_** = Metrics to control the mesh quality and element size  (refer to [4] for a description of these measures; the sample runMeshingTool file contains values which produced a good mesh quality) [*]
- **_lloyd_** = Enable the Lloyd global optimizer. [*]
- **_odt_** = Enable the global optimized Delaunay triangulation optimizer. [*]
- **_perturb_** = Enable local mesh optimization by mesh vertex perturbation. [*]
- **_exude_** = Enable local mesh optimization by exuding slivers. [*]

## Secondary information
To facilitate the setup of the proper environment to built the MeshHeadModel tool, we provide a Dockerfile to create a container based on Ubuntu 18.10 which contains all required libraries. 

Simply build the container with `docker build -t cgal_ubuntu1810 .` in the docker directory of this repository.

To compile it is sufficient to run the bash script *docker_compile.sh*. Note that you may have to adjust the variables $BASE_DIR and $CONTAINER_NAME if you use a path/name other then predefined in the file.

When compiled in the container it must also be lauchned from inside the container. You can do this by calling the *docker_run.sh* script which again invokes the *runMeshingTool.sh* script as mentioned before from inside the container. Again you may have to adjust the name of the container and the location of the repository. When using docker to run the tool consider that all paths (both input files and output files) must be specified according to the docker environment NOT according to the host machine. The best pratice is to use the predefined input & output directories and to use paths releative to the MeshHeadModel executable. 

###### Useful links
[1] The CGAL project: https://doc.cgal.org/4.12/Manual/index.html

[2] THe GMSH project: http://gmsh.info/

[3] GMSH 3.0.6 source code: http://gmsh.info/src/gmsh-3.0.6-source.tgz

[4] Explanation of the CGAL quality criteria: https://doc.cgal.org/4.12/Mesh_3/index.html#title11
