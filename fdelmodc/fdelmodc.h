#include<stdlib.h>
#include<stdio.h>
#include<math.h>

typedef struct _compType { /* Receiver Type */
	int vz;
	int vx;
	int p;
	int txx;
	int tzz;
	int txz;
	int pp;
	int ss;
	int ud;
} compType;

typedef struct _receiverPar { /* Receiver Parameters */
	char *file_rcv;
	compType type;
	int n;
	int nt;
	int delay;
	int skipdt;
	int *z;
	int *x;
	float *zr;
	float *xr;
	int int_vx;
	int int_vz;
	int scale;
	int sinkdepth;
	int sinkvel;
	float cp;
	float rho;
} recPar;

typedef struct _snapshotPar { /* Snapshot Parameters */
	char *file_snap;
	char *file_beam;
	compType type;
	int nsnap;
	int delay;
	int skipdt;
	int skipdz;
	int skipdx;
	int nz;
	int nx;
	int z1;
	int z2;
	int x1;
	int x2;
	int vxvztime;
	int beam;
} snaPar;

typedef struct _modelPar { /* Model Parameters */
	int iorder;
	int ischeme;
	int grid_dir;
	int sh;
	char *file_cp;
	char *file_ro;
	char *file_cs;
	char *file_qp;
	char *file_qs;
	float dz;
	float dx;
	float dt;
	float tmod;
	int nt;
	float z0;
	float x0;
	int nz;
	int nx;
	int naz;
	int nax;
	/* Vx: rox */
	int ioXx;
	int ioXz;
	int ieXx;
	int ieXz;
	/* Vz: roz */
	int ioZx;
	int ioZz;
	int ieZx;
	int ieZz;
	/* P, Txx, Tzz: lam, l2m */
	int ioPx;
	int ioPz;
	int iePx;
	int iePz;
	/* Txz: muu */
	int ioTx;
	int ioTz;
	int ieTx;
	int ieTz;
	float Qp;
	float Qs;
	float fw;
} modPar;

typedef struct _waveletPar { /* Wavelet Parameters */
	char *file_src;
	int nt;
	int nx;
	float dt;
	float fmax;
	int random;
	int seed;
	int nst;
	size_t *nsamp;
} wavPar;

typedef struct _sourcePar { /* Source Array Parameters */
	int n;
	int type;
	int orient;
	int *z;
	int *x;
	int single;	
	int plane;
	int circle;
	int array;
	int random;
	float *tbeg;
	float *tend;
	int multiwav;
	float angle;
	float velo;
	float amplitude;
	int distribution;
	int window;
    int injectionrate;
	int sinkdepth;
} srcPar;

typedef struct _shotPar { /* Shot Parameters */
	int n;
	int *z;
	int *x;
} shotPar;

typedef struct _boundPar { /* Boundary Parameters */
	int top;
	int bot;
	int lef;
	int rig;
	float *tapz;
	float *tapx;
	float *tapxz;
	int cfree;
	int ntap;
	int *surface;
    int pml;
    float *pml_Vx;
    float *pml_nzVx;
    float *pml_nxVz;
    float *pml_nzVz;
    float *pml_nxP;
    float *pml_nzP;

} bndPar;


#if __STDC_VERSION__ >= 199901L
  /* "restrict" is a keyword */
#else
#define restrict 
#endif

