#!/bin/bash

# Check inputs
if [ "$#" -ne 3 ]; then
	echo "2dSU_to_3dSU_v2: wrong number of input parameters. Exiting."
	echo "Creates a 3D SU model by replicating 2D SU model in the y-direction."
	echo "Usage: "
	echo "./2dSU_to_3dSU_v2.sh file.su ny dy"
	echo "Ex: "
	echo "./2dSU_to_3dSU_v2.sh file.su 120 2.5" 
	echo "(outputs file-3d.su)"
	#return #for function
	exit #for script
fi

# Functions
function gethdrval(){
	val=$(cat $1 | grep $2) # get parameter
	val=${val#$2} # remove prefix
	val="${val#"${val%%[![:space:]]*}"}" #remove whitespaces
	
	echo "$val"
	#return "$val";
}

# Main 
name="${1%%.*}" #strip the extension
ny=$2
dy=$3

# Get header values
surange <$1 >hdr.txt
nz=$(gethdrval hdr.txt ns)
nx=$(gethdrval hdr.txt trwf)
dz=$(gethdrval hdr.txt d1)
dx=$(gethdrval hdr.txt d2)


echo "nz is $nz"
echo "nx is $nx"
echo "ny is $ny"
echo "dz is $dz"
echo "dx is $dx"
echo "dy is $dy"

echo -e "\n Running sustrip"
sustrip <$1 >$name.bin outpar=hdr.txt

echo -e "\n Running twoD_to_twoHalfD "
twoD_to_twoHalfD.exe nz=$nz nx=$nx ny=$ny input2d=$name.bin output3d=$name-3d.bin

echo -e "\n Running bin2su "
bin2su.exe binmodel=$name-3d.bin sumodel=$name-3d.su nz=$nz nx=$nx ny=$ny dz=$dz dx=$dx dy=$dy odtkey=1


# Clean up
rm hdr.txt
rm $name.bin
rm $name-3d.bin






