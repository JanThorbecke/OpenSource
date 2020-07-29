#!/bin/bash

echo "WARNING! The data size of the 3D reflection data is very large!!!"
echo "When you are positive your machine has enough memory,"
echo "delete the exit; statement from makeR.sh"

exit;

#create the reflection data that was modeled in 1 1D medium
makeR1D file_in=shotx10y10.su file_out=reflx10y10.su verbose=2 nxrcv=201 nyrcv=61 nxsrc=201 nysrc=61

