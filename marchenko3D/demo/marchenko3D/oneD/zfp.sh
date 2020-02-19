#!/bin/bash

# Pre-transform the data to the frequency domain
#TWtransform file_in=reflx10y10.su file_out=reflx10y10_W.bin verbose=1 fmin=0 fmax=30 zfp=0 tolerance=1e-7

# Pre-transform the data to the frequency domain and apply zfp compression
TWtransform file_in=reflx10y10.su file_out=reflx10y10_zfp.bin verbose=1 fmin=0 fmax=30 zfp=1 tolerance=1e-7

rm reflx10y10.su
