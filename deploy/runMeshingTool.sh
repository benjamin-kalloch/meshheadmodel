#!/bin/bash

./meshHeadModel --imagefile 'input/almi_noskin.hdr' \
          --electrodefiles 'input/almi_dual_smooth_electrode0_origmesh_adapted.off,input/almi_dual_smooth_electrode1_origmesh_adapted.off' \
          --outfile 'output/almi_dual.msh' \
          --smoothskinfile 'input/almi_smooth_skin.off' \
		  --f_angbound=30 \
		  --f_distbound=0.4 \
		  --f_sizebound=2.0 \
		  --c_reratio=2.1 \
		  --c_sizebound=2.0  \
		  --lloyd \
          --exude
