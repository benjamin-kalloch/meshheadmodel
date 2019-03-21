#!/bin/bash

./meshHeadModel  \
          --imagefile 'input/lesioned_tissue.hdr' \
          --electrodefiles 'input/electrode0.off,input/electrode1.off' \
          --tissuesurfaces 'input/skin.off,input/skull.off,input/csf.off,input/gm.off,input/wm.off' \
          --outfile 'output/lesion_head_model.msh' \
          --f_angbound=30 \
          --f_distbound=0.8 \
          --f_sizebound=4.0 \
          --c_reratio=4.2 \
          --c_sizebound=4.0 \
          --lloyd \
          --exude 
