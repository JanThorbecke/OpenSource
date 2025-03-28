#!/bin/bash

export PATH=$HOME/src/OpenSource/bin:$PATH

dx=2.5
dt=0.0005

#the model upto 900 m depth, deeper reflections are not needed to model the direct transmission response

makemod sizex=10000 sizez=1400 dx=$dx dz=$dx cp0=1800 ro0=1000 \
        orig=-5000,0 file_base=modelup.su verbose=2 \
        intt=def x=-5000,5000 z=400,400 poly=0 cp=2300 ro=3000 \
        intt=def x=-5000,5000 z=700,700 poly=0 cp=2000 ro=1100 

makewave fp=25 dt=$dt file_out=wave.su nt=4096 t0=0.1 scale=1

export OMP_NUM_THREADS=4
angle=0

fdelmodc \
    file_cp=modelup_cp.su ischeme=1 iorder=4 \
    file_den=modelup_ro.su \
    file_src=wave.su \
    file_rcv=iniFocusPlane_a${angle}.su \
    src_type=1 \
    src_orient=1 \
    src_injectionrate=1 \
    rec_type_vz=0 \
    rec_type_p=1 \
    rec_int_vz=2 \
    dtrcv=0.004 \
    rec_delay=0.1 \
    verbose=2 \
    tmod=2.144 \
    dxrcv=5 \
    plane_wave=1 nsrc=4001 src_angle=${angle} src_velo=1500 \
    xrcv1=-2250 xrcv2=2250 \
    zrcv1=0 zrcv2=0 \
    xsrc=0 zsrc=900 \
    ntaper=101 \
    left=2 right=2 top=2 bottom=2 

exit;

fdelmodc \
    file_cp=model10_cp.su ischeme=1 iorder=4 \
    file_den=model10_ro.su \
    file_src=wave.su \
    file_rcv=RefFocusPlane.su \
    src_type=1 \
    src_orient=1 \
    src_injectionrate=1 \
    rec_type_vz=0 \
    rec_type_p=1 \
    rec_int_vz=2 \
    dtrcv=0.004 \
    rec_delay=0.1 \
    verbose=2 \
    tmod=2.144 \
    dxrcv=5 \
plane_wave=1 nsrc=4001 src_angle=3 src_velo=1500 \
    xrcv1=-2250 xrcv2=2250 \
    zrcv1=0 zrcv2=0 \
    xsrc=0 zsrc=540 \
    ntaper=101 \
    left=2 right=2 top=2 bottom=2 

