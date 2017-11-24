#!/bin/bash

filename_in=/vardim/home/boulleng/Canada_Data/Lalor_shots_clean_4000ms_agc_sortgdelsdel_sdel136out.su
filename_out=virshots_fmax160_FROM_all.su


ulimit -c unlimited
../corrvir \
	file_shots=${filename_in} \
	file_out=${filename_out} \
	nsources=908 \
	src_sel=0 \
	fmax=160 \
	nbm=1 \
	nreceivers=2685 \
	normsrc=0 \
	normalize=0 \
	cohr=0 \
	causal=1 \
	verbose=4
