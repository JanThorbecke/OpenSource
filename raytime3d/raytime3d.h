#include<stdlib.h>
#include<stdio.h>
#include<limits.h>
#include<float.h>
#include<math.h>

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
typedef struct _dcomplexStruct { /* complex number */
    double r,i;
} dcomplex;
#endif/* complex */


typedef struct _icoord { /* 3D coordinate integer */
    int z;
    int x;
    int y;
} icoord;

typedef struct _fcoord { /* 3D coordinate float */
    float z;
    float x;
    float y;
} fcoord;

struct s_ecount {
  int       corner,corner_min,side;
};

typedef struct _receiverPar { /* Receiver Parameters */
	char *file_rcv;
	int n;
	int nt;
	int max_nrec;
	int *z;
	int *x;
	int *y;
	float *zr;
	float *xr;
	float *yr;
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
	float dy;
	float dt;
	float z0;
	float x0;
	float y0;
	/* medium max/min values */
	float cp_min;
	float cp_max;
	int nz;
	int nx;
	int ny;
} modPar;

typedef struct _waveletPar { /* Wavelet Parameters */
	char *file_src; /* general source */
	int nsrcf;
	int nt;
	int ns;
	int nx;
	int ny;
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
	int *y;
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
	int ny;
	int nx;
	int nz;
	int *z;
	int *x;
	int *y;
} shotPar;

typedef struct _raypar { /* ray-tracing parameters */
    int smoothwindow;
    int useT2;
    int geomspread;
    int nray;
} rayPar;

#ifndef TRUE
#  define TRUE 1
#endif

#ifndef FALSE
#  define FALSE 0
#endif

#define equal(x,y) !strcmp(x,y)
#define min2(a,b) (((a) < (b)) ? (a) : (b))
#define max2(a,b) (((a) > (b)) ? (a) : (b))

#define Infinity FLT_MAX

#if __STDC_VERSION__ >= 199901L
  /* "restrict" is a keyword */
#else
#define restrict 
#endif

