#!/bin/bash

#Marchenko using time data
#../../../../marchenko3D file_shot=../reflx10y10.su verbose=2 \
#    file_tinv=farr_tilted_pos.su \
#    plane_wave=1 src_anglex=10 src_angley=5 src_velox=2170 src_veloy=2170 \
#    niter=10 shift=15 smooth=10 hw=5 \
#    file_green=green_tilted_pos_timeshot.su file_file_f2=f2_tilted_pos_timeshot.su \
#    file_gplus=gplus_tilted_pos_timeshot.su file_gmin=gmin_tilted_pos_timeshot.su

#Marchenko using frequency data
#../../../../marchenko3D file_shotw=../reflx10y10_W.bin verbose=2 \
#    file_tinv=farr_tilted_pos.su \
#    plane_wave=1 src_anglex=10 src_angley=5 src_velox=2170 src_veloy=2170 \
#    niter=10 shift=15 smooth=10 hw=5 \
#    file_green=green_tilted_pos_freqshot.su file_file_f2=f2_tilted_pos_freqshot.su
#    file_gplus=gplus_tilted_pos_freqshot.su file_gmin=gmin_tilted_pos_freqshot.su

#Marchenko using zfp data
../../../../marchenko3D file_shotzfp=../reflx10y10_zfp.bin verbose=2 \
    file_tinv=farr_tilted_pos.su \
    plane_wave=1 src_anglex=10 src_angley=5 src_velox=2170 src_veloy=2170 \
    niter=10 shift=15 smooth=10 hw=5 \
    file_green=green_tilted_pos_zfpshot.su file_file_f2=f2_tilted_pos_zfpshot.su \
    file_gplus=gplus_tilted_pos_zfpshot.su file_gmin=gmin_tilted_pos_zfpshot.su
