#!/bin/bash

echo "Before running this file, run model.sh, makeR.sh and zfp.sh in the folder above"

./farrmod_tilted_pos.sh

./farrmod_tilted_neg.sh

./fmute3D_tilted_pos.sh

./fmute3D_tilted_neg.sh

./marchenko3d_tilted_pos.sh

./marchenko3d_tilted_neg.sh

susum gplus_tilted_neg_zfpshot.su gmin_tilted_pos_zfpshot.su > green_total_zfpshot_neg.su
susum gplus_tilted_pos_zfpshot.su gmin_tilted_neg_zfpshot.su > green_total_zfpshot_pos.su
