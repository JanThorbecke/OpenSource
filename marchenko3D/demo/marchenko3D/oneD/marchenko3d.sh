#!/bin/bash

#Marchenko using time data
#marchenko3D file_shot=reflx10y10.su verbose=2 \
#    file_tinv=farr.su \
#    niter=10 shift=15 smooth=10 hw=5 \
#    file_green=green_timeshot.su file_file_f2=f2_timeshot.su

#Marchenko using frequency data
#marchenko3D file_shotw=reflx10y10_W.bin verbose=2 \
#    file_tinv=farr.su \
#    niter=10 shift=15 smooth=10 hw=5 \
#    file_green=green_freqshot.su file_file_f2=f2_freqshot.su

#Marchenko using zfp data
marchenko3D file_shotzfp=reflx10y10_zfp.bin verbose=2 \
    file_tinv=farr.su \
    niter=10 shift=15 smooth=10 hw=5 \
    file_green=green_zfpshot.su file_file_f2=f2_zfpshot.su
