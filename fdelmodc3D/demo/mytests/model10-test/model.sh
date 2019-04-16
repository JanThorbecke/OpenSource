#!/bin/bash

#adjust this PATH to where the code is installed
export PATH=$HOME/src/OpenSource/bin:$PATH:

dx=2.5

	#define gridded model for FD computations
makemod sizex=2000 sizez=1400 dx=$dx dz=$dx cp0=1800 ro0=1000 \
	orig=-1000,0 file_base=model10.su verbose=2 \
        intt=def x=-1000,1000 z=400,400 poly=0 cp=2300 ro=3000 \
        intt=def x=-1000,1000 z=700,700 poly=0 cp=2000 ro=1100 \
        intt=def x=-1000,1000 z=1100,1100 poly=0 cp=2500 ro=4000

