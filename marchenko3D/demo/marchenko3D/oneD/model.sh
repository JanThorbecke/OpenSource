#!/bin/bash

#adjust this PATH to where the code is installed
export PATH=$HOME/src/OpenSource/bin:$PATH:

dx=10
ny=181
dy=10000

#define gridded model for FD computations
makemod sizex=5000 sizez=1400 dx=$dx dz=$dx cp0=1800 ro0=1000 \
        orig=-2500,0 file_base=model5.su verbose=2 \
        intt=def x=-2500,2500 z=400,400 poly=0 cp=2300 ro=3000 \
        intt=def x=-2500,2500 z=700,700 poly=0 cp=2000 ro=1100 \
        intt=def x=-2500,2500 z=1100,1100 poly=0 cp=2500 ro=4000

#define homogenoeus model to compute direct wave only
makemod sizex=5000 sizez=300 dx=$dx dz=$dx cp0=1800 ro0=1000 \
        orig=-2500,0 file_base=hom.su verbose=2

rm cp3d.su
rm ro3d.su
rm cphom3d.su
rm rohom3d.su

for (( i=0; i<$ny; i++ ));do

    (( ypos = -1*($ny-1)/2 + $i ))
    (( yloc = $dy*$ypos ))
    echo $yloc

    sushw < model5_cp.su key=gy,scalel a=$yloc,-1000 >> cp3d.su
    sushw < model5_ro.su key=gy,scalel a=$yloc,-1000 >> ro3d.su
    sushw < hom_cp.su key=gy,scalel a=$yloc,-1000 >> cphom3d.su
    sushw < hom_ro.su key=gy,scalel a=$yloc,-1000 >> rohom3d.su

done



