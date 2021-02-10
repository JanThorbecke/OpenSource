#!/bin/bash

#Marchenko using time data
#../../../../marchenko3D file_shot=../reflx10y10.su verbose=2 \
#    file_tinv=farr_plane.su \
#    plane_wave=1 src_anglex=0 src_angley=0 src_velox=2170 src_veloy=2170 \
#    niter=10 shift=15 smooth=10 hw=5 \
#    file_green=green_plane_timeshot.su file_file_f2=f2_plane_timeshot.su

#Marchenko using frequency data
#../../../../marchenko3D file_shotw=../reflx10y10_W.bin verbose=2 \
#    file_tinv=farr_plane.su \
#    plane_wave=1 src_anglex=0 src_angley=0 src_velox=2170 src_veloy=2170 \
#    niter=10 shift=15 smooth=10 hw=5 \
#    file_green=green_plane_freqshot.su file_file_f2=f2_plane_freqshot.su

#Marchenko using zfp data
../../../../marchenko3D file_shotzfp=../reflx10y10_zfp.bin verbose=2 \
    file_tinv=farr_plane.su \
    plane_wave=1 src_anglex=0 src_angley=0 src_velox=2170 src_veloy=2170 \
    niter=10 shift=15 smooth=10 hw=5 \
    file_green=green_plane_zfpshot.su file_file_f2=f2_plane_zfpshot.su
