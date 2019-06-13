
export CWPROOT=/home/users/jan/SeisUnix/
export PATH=.:$CWPROOT/bin:/lus/scratch/$USER/OpenSource/bin:$PATH:
alias lt='ls -lart'

module swap PrgEnv-cray PrgEnv-intel
module list

mkdir -p /lus/scratch/$USER

cd /lus/scratch/$USER
rsync -av /lus/scratch/jan/OpenSource .
#git clone https://github.com/JanThorbecke/OpenSource.git

