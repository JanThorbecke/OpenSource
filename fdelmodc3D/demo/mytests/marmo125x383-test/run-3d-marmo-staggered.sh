#!/bin/bash


dt=0.0005
../../../utils/makewave fp=3 dt=$dt file_out=wave.su nt=4096 t0=0.0 scale=1 shift=1

../../fdelmodc3D par=marmo3d-data.dat


