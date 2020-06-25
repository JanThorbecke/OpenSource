#include<stdlib.h>
#include<stdio.h>
#include<math.h>

typedef struct _compType { /* Receiver Type */
	long vz;
	long vx;
	long vy;
	long p;
	long txx;
	long tyy;
	long tzz;
	long txz;
	long tyz;
	long txy;
	long pp;
	long ss;
	long ud;
} compType;

typedef struct _receiverPar { /* Receiver Parameters */
	char *file_rcv;
	compType type;
	long n;
	long nt;
	long delay;
	long skipdt;
	long max_nrec;
	long *z;
	long *y;
	long *x;
	float *zr;
	float *yr;
	float *xr;
	long int_p;
	long int_vx;
	long int_vy;
	long int_vz;
	long scale;
	long sinkdepth;
	long sinkvel;
	float cp;
	float rho;
} recPar;

typedef struct _snapshotPar { /* Snapshot Parameters */
	char *file_snap;
	char *file_beam;
	compType type;
	long nsnap;
	long delay;
	long skipdt;
	long skipdz;
	long skipdy;
	long skipdx;
	long nz;
	long ny;
	long nx;
	long z1;
	long z2;
	long x1;
	long x2;
	long y1;
	long y2;
	long vxvztime;
	long beam;
	long withbnd;
} snaPar;

typedef struct _modelPar { /* Model Parameters */
	long iorder;
	long ischeme;
	long grid_dir;
	long sh;
	char *file_cp;
	char *file_ro;
	char *file_cs;
	char *file_qp;
	char *file_qs;
	float dz;
	float dy;
	float dx;
	float dt;
	float tmod;
	long nt;
	float z0;
	float y0;
	float x0;
	/* medium max/min values */
	float cp_min;
	float cp_max;
	float cs_min;
	float cs_max;
	float ro_min;
	float ro_max;
	long nz;
	long ny;
	long nx;
	long naz;
	long nay;
	long nax;
	long nfz;
	long nfy;
	long nfx;
	/* Vx: rox */
	long ioXx;
	long ioXy;
	long ioXz;
	long ieXx;
	long ieXy;
	long ieXz;
	/* Vy: roy */
	long ioYx;
	long ioYy;
	long ioYz;
	long ieYx;
	long ieYy;
	long ieYz;
	/* Vz: roz */
	long ioZx;
	long ioZy;
	long ioZz;
	long ieZx;
	long ieZy;
	long ieZz;
	/* P, Txx, Tyy, Tzz: lam, l2m */
	long ioPx;
	long ioPy;
	long ioPz;
	long iePx;
	long iePy;
	long iePz;
	/* Txz, Txy, Tyz: muu */
	long ioTx;
	long ioTy;
	long ioTz;
	long ieTx;
	long ieTy;
	long ieTz;
	/* attenuation / dissipative medium */
	float Qp;
	float Qs;
	float fw;
	float qr;
} modPar;

typedef struct _waveletPar { /* Wavelet Parameters */
	char *file_src; /* general source */
	long nsrcf;
	long nt;
	long ns;
	long nx;
	float dt;
	float ds;
	float fmax;
	long random;
	long seed;
	long nst;
	size_t *nsamp;
} wavPar;

typedef struct _sourcePar { /* Source Array Parameters */
	long n;
	long nx;
	long ny;
	long type;
	long orient;
	long *z;
	long *y;
	long *x;
	long single;	
	long plane;
	long circle;
	long array;
	long random;
	float Mxx;
	float Mxy;
	float Mxz;
	float Myy;
	float Myz;
	float Mzz;
	float *tbeg;
	float *tend;
	long multiwav;
	float angle;
	float velo;
	float amplitude;
	float dip;
	float strike;
	float rake;
	long distribution;
	long nxwindow;
	long nywindow;
	long window;
    long injectionrate;
	long sinkdepth;
	long src_at_rcv; /* Indicates that wavefield should be injected at receivers */
} srcPar;

typedef struct _shotPar { /* Shot Parameters */
	long n;
	long *z;
	long *y;
	long *x;
} shotPar;

typedef struct _boundPar { /* Boundary Parameters */
	long top;
	long bot;
	long lef;
	long rig;
	long fro;
	long bac;
	float *tapz;
	float *tapy;
	float *tapx;
	float *tapxz;
	float *tapxyz;
	long cfree;
	long ntap;
	long *surface;
    long npml;
    float R; /* reflection at side of model */
    float m; /* scaling order */
    float *pml_Vx;
		float *pml_Vy;
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

