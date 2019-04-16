#!/bin/bash


dt=0.0004
../../../../utils/makewave fp=25 dt=$dt file_out=wave.su nt=4096 t0=0.1 scale=1

export OMP_NUM_THREADS=80
../../../fdelmodc3D par=model10_data.dat


