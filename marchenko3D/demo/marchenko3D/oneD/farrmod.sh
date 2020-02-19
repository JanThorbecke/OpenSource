#!/bin/bash

dt=0.001

makewave fp=15 fmax=22 dt=${dt} file_out=wavefpmod.su nt=8192 t0=0.2 scale=1 scfft=0

export OMP_NUM_THREADS=8

#Model shot record in middle of model
fdelmodc3D \
    file_cp=cp3d.su ischeme=1 iorder=4 \
    file_den=ro3d.su \
    file_src=wavefpmod.su \
    file_rcv=farrmod.su \
    src_type=1 \
    src_orient=1 \
    src_injectionrate=1 \
    rec_type_vz=0 \
    rec_type_p=1 \
    rec_int_vz=2 \
    dtrcv=0.004 \
    rec_delay=0.2 \
    verbose=2 \
    tmod=1.2 \
    dxrcv=10.0 \
    xrcv1=-1000 xrcv2=1000 \
    dyrcv=10.0 \
    yrcv1=-300 yrcv2=300 \
    zrcv1=0 zrcv2=0 \
    xsrc=0 ysrc=0 zsrc=900 \
    ntaper=61 \
    left=4 right=4 top=4 bottom=4 front=4 back=4 \
