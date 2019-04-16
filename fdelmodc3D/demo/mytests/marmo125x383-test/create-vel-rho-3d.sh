#!/bin/bash

# Create 2.5-d binary model
echo "CREATING 3D BINARY VELmodel"
velName=vel125x383
nz=125
nx=383
ny=104
twoD_to_twoHalfD.exe nz=$nz nx=$nx ny=$ny input2d=$velName.bin output3d=$velName-3d.bin
echo ""

# Create rho from vel3d
echo "CREATING 3D BINARY RHOmodel"
rhoName=rho125x383
nsamp=$(($nz*$nx*$ny))
vel2rho.exe nsamp=$nsamp velfile=$velName-3d.bin rhofile=$rhoName-3d.bin
echo ""

# Put SU headers on vel3d and rho3d
echo "PUTTING SU HEADERS"
dz=24.0
dx=24.0
dy=24.0
bin2su.exe binmodel=$velName-3d.bin sumodel=$velName-3d.su nz=$nz nx=$nx ny=$ny dz=$dz dx=$dx dy=$dy opendt=1
bin2su.exe binmodel=$rhoName-3d.bin sumodel=$rhoName-3d.su nz=$nz nx=$nx ny=$ny dz=$dz dx=$dx dy=$dy opendt=1
echo ""






