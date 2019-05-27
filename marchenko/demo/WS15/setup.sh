
export CWPROOT=/home/users/jan/SeisUnix/
export PATH=.:$CWPROOT/bin:/lus/scratch/$USER/OpenSource/bin:$PATH:

module swap PrgEnv-cray PrgEnv-intel
module list

mkdir -p /lus/scratch/$USER

cd /lus/scratch/$USER
#git clone https://github.com/JanThorbecke/OpenSource.git

