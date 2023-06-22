#!/bin/bash
#PBS -N fdelmod
#PBS -q long
#PBS -l nodes=1
#PBS -k eo
#PBS -j eo

export PATH=$HOME/src/OpenSource/bin:$PATH:

dt=0.0005
makewave w=fw fmin=0 flef=5 frig=80 fmax=100  dt=$dt file_out=wavefw.su nt=4096 t0=0.3 scale=0 scfft=1

#./model.scr

mkdir -p shots
mkdir -p jobs

dxshot=10
ishot=0
nshots=601
zsrc=0

cat << EOF > jobs/slurm.job 
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --array=0-600
#SBATCH -J mod2D
#SBATCH -V
#SBATCH -p max30m

(( xsrc = -3000 + \${SLURM_ARRAY_TASK_ID}*${dxshot} ))
export PATH=\$HOME/src/OpenSource/bin:\$PATH:

echo ishot=\$SLURM_ARRAY_TASK_ID xsrc=\$xsrc zsrc=$zsrc
export OMP_NUM_THREADS=4
file_rcv=shots/shots_\${xsrc}.su

fdelmodc \
  	file_cp=ge_cp.su ischeme=1 iorder=4 \
   	file_den=ge_ro.su \
   	file_src=wavefw.su \
   	file_rcv=\$file_rcv \
	src_type=7 \
	src_orient=1 \
	src_injectionrate=1 \
   	rec_type_vz=0 \
   	rec_type_p=1 \
   	rec_int_vz=2 \
	rec_delay=0.3 \
   	dtrcv=0.004 \
   	verbose=2 \
   	tmod=4.394 \
   	dxrcv=10.0 \
   	xrcv1=-3000 xrcv2=3000 \
   	zrcv1=0 zrcv2=0 \
   	xsrc=\$xsrc zsrc=$zsrc \
   	ntaper=200 \
   	left=2 right=2 top=2 bottom=2
EOF

#sbatch jobs/slurm.job


