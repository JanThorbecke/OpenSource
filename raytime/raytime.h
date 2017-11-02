#include<stdlib.h>
#include<stdio.h>
#include<math.h>

typedef struct _receiverPar { /* Receiver Parameters */
	char *file_rcv;
	int n;
	int nt;
	int max_nrec;
	int *z;
	int *x;
	float *zr;
	float *xr;
	int scale;
	int sinkdepth;
	int sinkvel;
	float cp;
	float rho;
} recPar;

typedef struct _modelPar { /* Model Parameters */
	int sh;
	char *file_cp;
	float dz;
	float dx;
	float dt;
	float z0;
	float x0;
	/* medium max/min values */
	float cp_min;
	float cp_max;
	int nz;
	int nx;
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
	int single;	
	int plane;
	int circle;
	int array;
	int random;
	int multiwav;
	float angle;
	float velo;
	float amplitude;
	int distribution;
	int window;
    int injectionrate;
	int sinkdepth;
	int src_at_rcv; /* Indicates that wavefield should be injected at receivers */
} srcPar;

typedef struct _shotPar { /* Shot Parameters */
	int n;
	int nx;
	int nz;
	int *z;
	int *x;
} shotPar;

typedef struct _raypar { /* ray-tracing parameters */
    int smoothwindow;
    int useT2;
    int geomspread;
    int nray;
} rayPar;

#if __STDC_VERSION__ >= 199901L
  /* "restrict" is a keyword */
#else
#define restrict 
#endif

