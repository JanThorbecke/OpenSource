#!/bin/bash
#PBS -N fdelmod
#PBS -q verylong
#PBS -l nodes=1
#PBS -k eo
#PBS -j eo

set -x
#-----------------------------------------------------------------------------
# Modeling of acoustic response of marine-type acquisition
# where a horizontal and a slanted cable is modeled simultaneously
# 
# Author: Eric Verschuur, Delft University of Technology
# Date  : July 2013
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# define the grid size of the modeling
# based on the grid, the other parameters can be define
#-----------------------------------------------------------------------------

grid=5
grid=10

#-----------------------------------------------------------------------------
# define the source wavelet:
# for gridsize 10 m use fp=9.7 (such that fmax<30)
# for gridsize 5 m use fp=22 (such that fmax<60) and
#-----------------------------------------------------------------------------

if [ $grid -eq 10 ]; then
   dt=0.00020
   fp=9.7
   vmax=3000
fi
if [ $grid -eq 5 ]; then
   dt=0.00010
   fp=22
   vmax=3000
fi

makewave file_out=wavelet.su dt=$dt nt=1024 fp=$fp shift=1 w=g2 verbose=1

#-----------------------------------------------------------------------------
# define the velocity model; in this case a flat water bottom and
# a curved reflector below the water bottom
#-----------------------------------------------------------------------------

makemod file_base=model.su \
   cp0=1500 ro0=1400 sizex=5000 sizez=1200 \
   dx=$grid dz=$grid orig=0,0 \
   \
   x=0,5000 z=500,500 \
   intt=def poly=2 cp=2000 ro=1800 \
   \
   x=0,2000,5000 z=1000,800,1000 \
   intt=def poly=2 cp=$vmax ro=2300 \
   \
   verbose=1

# display the models (velocity and density)

filecp=model_cp.su
filero=model_ro.su

suximage < $filecp wbox=800 hbox=300 title="Vp model" xbox=0 legend=1 &
suximage < $filero wbox=800 hbox=300 title="Rho model" xbox=800 legend=1 &

sleep 1

#-----------------------------------------------------------------------------
# define the source location and make receivers dependent on the source
# define length of cable and the vertical slant of the cable, which
# is defined as the relative depth at the end of the cable
# the receiver spacing is taken as twice the FD grid size
#-----------------------------------------------------------------------------

xsrc=500
cable=4000
slant=200

xrcv1=`dadd $xsrc 103`
xrcv2=`dadd $xrcv1 $cable`
zrcv1=`dadd $grid 3`
zrcv2=`dadd $zrcv1 $slant`
zrcv2=`dnint $zrcv2`
dxrcv=`dmul $grid 2`
xrcv1=`dnint $xrcv1`
xrcv2=`dnint $xrcv2`
dxrcv=`dnint $dxrcv`

#-----------------------------------------------------------------------------
# generate a list of coorindates for the slanted cable
# finally, these coordinates are stored in the file rcv.par
#-----------------------------------------------------------------------------

xrcv=$xrcv1
zrcv=$zrcv1
dzrcv=`echo $xrcv1 $xrcv2 $dxrcv $slant | awk '{nx=int(1.5+($2-$1)/$3);print ($4/(nx-1))}'`
/bin/rm rcv.list
while [ $xrcv -le $xrcv2 ]
do
   echo "$xrcv $zrcv" >> rcv.list
   xrcv=`expr $xrcv + $dxrcv`
   zrcv=`dadd $zrcv $dzrcv`
done 
mkparfile < rcv.list string1=xrcva string2=zrcva > rcv.par
/bin/rm rcv.list

#-----------------------------------------------------------------------------
# now do the actual modeling with free surface and monopole src/rcv
# we use both the regular rcv defintions for a horizontal cable and the
# extended option with xrcva,zrcva for the slanted cable locations
#-----------------------------------------------------------------------------

../fdelmodc \
	file_cp=$filecp file_den=$filero \
	ischeme=1 \
	file_src=wavelet.su verbose=4 \
	file_rcv=rec.su \
	file_snap=snap.su \
	xrcv1=$xrcv1 xrcv2=$xrcv2 dxrcv=$dxrcv \
	zrcv1=$zrcv1 zrcv2=$zrcv1 \
	par=rcv.par \
	sinkdepth=0 \
	rec_type_p=1 rec_type_vz=1 rec_int_vz=3 \
	dtrcv=0.004 \
        src_type=1 src_orient=1 xsrc=$xsrc zsrc=$grid nshot=1 \
	tsnap1=0.1 tsnap2=3.0 dtsnap=0.1 \
	sna_type_p=1 sna_type_vz=0 \
	top=1 bottom=4 left=4 right=4 ntaper=100 tapfact=0.3 \
	tmod=4.0 nzmax=300 nxmax=300

#-----------------------------------------------------------------------------
# show a movie of the wavefield snapshots
#-----------------------------------------------------------------------------

suxmovie < snap_sp.su clip=10 width=800 height=300 loop=1 \
   title="Shot at x=$xsrc - P comp. - snapshot at t=%f s." sleep=200000 fframe=0.1 dframe=0.1 &

#-----------------------------------------------------------------------------
# the 1st part of the receivers is the slanted cable, the 2nd part the flat one
#-----------------------------------------------------------------------------

nxtot=`surange < rec_rp.su | head -1 | awk '{print $1}'`
nx=`ddiv $nxtot 2`
nx=`dnint $nx`

# split into two files; on-the-fly repair dt from 3999 to 4000 us

sushw key=fldr < rec_rp.su a=1 c=1 j=$nx | \
suchw key1=dt key2=dt a=1 d=10 | \
suchw key1=dt key2=dt b=10 | \
sushw key=ntr a=0 | \
file_distribute file_base=shot verbose=1

suximage < shot1.su wbox=800 hbox=600 xbox=0 ybox=500 perc=99 \
   title="Shot record xsrc=$xsrc - slanted (x1,z1)=$xrcv1,$zrcv1 (x2,z2)=$xrcv2,$zrcv2" &

suximage < shot2.su wbox=800 hbox=600 xbox=800 ybox=500 perc=99 \
   title="Shot record xsrc=$xsrc - flat cable  (x1,z1)=$xrcv1,$zrcv1 (x2,z2)=$xrcv2,$zrcv1" &

#-----------------------------------------------------------------------------
# end of demo, remove some tmp files
#-----------------------------------------------------------------------------

sleep 1
#/bin/rm rcv.par
