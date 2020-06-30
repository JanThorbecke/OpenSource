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
	int max_nrec;
	int *z;
	int *x;
	float *zr;
	float *xr;
	int int_p;
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
	int withbnd;
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
	/* medium max/min values */
	float cp_min;
	float cp_max;
	float cs_min;
	float cs_max;
	float ro_min;
	float ro_max;
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
	/* attenuation / dissipative medium */
	float Qp;
	float Qs;
	float fw;
	float qr;
} modPar;

typedef struct _waveletPar { /* Wavelet Parameters */
	char *file_src; /* general source */
	int nsrcf;
	int nt;
	int ns;
	int nx;
	float dt;
	float ds;
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
	float Mxx;
	float Mzz;
	float Mxz;
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
	float dip;
	float strike;
	int distribution;
	int window;
    int injectionrate;
	int sinkdepth;
	int src_at_rcv; /* Indicates that wavefield should be injected at receivers */
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
    int npml;
    float R; /* reflection at side of model */
    float m; /* scaling order */
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

