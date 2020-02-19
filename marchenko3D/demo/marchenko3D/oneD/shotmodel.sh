#!/bin/bash

dt=0.001

#define wavelet for modeling R
makewave w=fw fmin=0 flef=5 frig=25 fmax=30 dt=$dt file_out=wavefw.su nt=4096 t0=0.3 scale=0 scfft=1

#Set the number of threads for parallel jobs. More is not always better
export OMP_NUM_THREADS=8

#Model shot record in middle of model
fdelmodc3D \
    file_cp=cp3d.su ischeme=1 iorder=4 \
    file_den=ro3d.su \
    file_src=wavefw.su \
    file_rcv=shotx10y10.su \
    src_type=7 \
    src_orient=1 \
    src_injectionrate=1 \
    rec_type_vz=0 \
    rec_type_p=1 \
    rec_int_vz=2 \
    dtrcv=0.004 \
    rec_delay=0.3 \
    verbose=2 \
    tmod=4.392 \
    dxrcv=10.0 \
    xrcv1=-2000 xrcv2=2000 \
    dyrcv=10.0 \
    yrcv1=-600 yrcv2=600 \
    zrcv1=0 zrcv2=0 \
    xsrc=0 ysrc=0 zsrc=0 \
    ntaper=41 \
    left=4 right=4 top=4 bottom=4 front=4 back=4

#Model direct wave only in middle of model
fdelmodc3D \
    file_cp=cphom3d.su ischeme=1 iorder=4 \
    file_den=rohom3d.su \
    file_src=wavefw.su \
    file_rcv=shotx10y10hom.su \
    src_type=7 \
    src_orient=1 \
    src_injectionrate=1 \
    rec_type_vz=0 \
    rec_type_p=1 \
    rec_int_vz=2 \
    dtrcv=0.004 \
    rec_delay=0.3 \
    verbose=2 \
    tmod=4.392 \
    dxrcv=10.0 \
    xrcv1=-2000 xrcv2=2000 \
    dyrcv=10.0 \
    yrcv1=-600 yrcv2=600 \
    zrcv1=0 zrcv2=0 \
    xsrc=0 ysrc=0 zsrc=0 \
    ntaper=41 \
    left=4 right=4 top=4 bottom=4 front=4 back=4

#subtract direct wave from full model shot record: this defines R
sudiff shotx10y10_rp.su shotx10y10hom_rp.su > shotx10y10.su
