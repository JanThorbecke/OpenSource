#!/bin/bash

# Check inputs
if [ "$#" -ne 4 ]; then
	echo "2dSU_to_3dSU_v1.sh: wrong number of input parameters. Exiting."
	echo "Creates a 3D SU model by replicating 2D SU model in the y-direction."
	echo "Usage: "
	echo "./2dSU_to_3dSU_v1.sh dx dy ny file.su"
	echo "Ex: "
	echo "./2dSU_to_3dSU_v1.sh 2.5 2.5 120 file.su" 
	echo "(outputs file-3d.su)"
	#return #for function
	exit #for script
fi

dx=$1
dy=$2
ny=$3
indata=$4

name="${4%%.*}" #strips extension from name
outdata=$name-3d.su

tmp=tmp1.su
tmp2=tmp2.su

rm $outdata
rm $tmp

t1=`date +%s`
for (( iy=1; iy<=$ny; iy++ ))
do  
	echo "Adding y-slice $iy "

	# Set y-slice with correct fldr
	cp $indata $tmp
	sushw <$tmp >$tmp2 key=fldr a=$iy b=0 c=0

	# Cat it to our 2.5D
	cat $tmp2 >>$outdata
done
t2=`date +%s`
echo "Extending data with cat took $((t2-t1)) seconds"

rm $tmp2

# Add cdp (x coordinate for odtect)
echo "Setting CDP key"
#suchw <$outdata >$tmp key1=cdp key2=tracl a=1
time suchw <$outdata >$tmp key1=cdp key2=tracl a=1

# Set sx
echo "Setting sx key"
factor=$( echo "$dx*1000" | bc)
#suchw <$tmp >$outdata key1=sx key2=cdp a=0 b=$factor 
time suchw <$tmp >$outdata key1=sx key2=cdp a=0 b=$factor 

# Set sy
echo "Setting sy key"
factor=$( echo "$dy*1000" | bc)
#suchw <$outdata >$tmp key1=sy key2=fldr a=0 b=$factor 
time suchw <$outdata >$tmp key1=sy key2=fldr a=0 b=$factor 

# Set dt
echo "Setting dt key"
#suchw <$tmp >$outdata key1=dt key2=d1 b=1000
time suchw <$tmp >$outdata key1=dt key2=d1 b=1000

rm $tmp
#mv $tmp $outdata







