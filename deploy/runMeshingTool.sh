#!/bin/bash

./meshHeadModel  \
          --imagefile 'input/almi5/almi5_allmaskscombined_noskin.hdr' \
          --electrodefiles 'input/almi5/almi_dual_smooth_electrode0.off,input/almi5/almi_dual_smooth_electrode1.off' \
          --tissuesurfaces 'input/almi5/almi_smooth_skin.off' \
          --outfile 'output/almi5_skull_improved_TEST.msh' \
          --f_angbound=30 \
          --f_distbound=0.4 \
          --f_sizebound=2.0 \
          --c_reratio=2.1 \
          --c_sizebound=2.0 \
          --simnibs
          #--lloyd \
          #--exude 
