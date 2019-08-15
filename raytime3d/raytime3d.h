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
    long z;
    long x;
    long y;
} icoord;

typedef struct _fcoord { /* 3D coordinate float */
    float z;
    float x;
    float y;
} fcoord;

struct s_ecount {
  long       corner,corner_min,side;
};

typedef struct _receiverPar { /* Receiver Parameters */
	char *file_rcv;
	long n;
	long nx;
	long ny;
	long nt;
	long max_nrec;
	long *z;
	long *x;
	long *y;
	float *zr;
	float *xr;
	float *yr;
	long scale;
	long sinkdepth;
	long sinkvel;
	float cp;
	float rho;
} recPar;

typedef struct _modelPar { /* Model Parameters */
	long sh;
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
	long nz;
	long nx;
	long ny;
} modPar;

typedef struct _waveletPar { /* Wavelet Parameters */
	char *file_src; /* general source */
	long nsrcf;
	long nt;
	long ns;
	long nx;
	long ny;
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
	long type;
	long orient;
	long *z;
	long *x;
	long *y;
	long single;	
	long plane;
	long circle;
	long array;
	long random;
	long multiwav;
	float angle;
	float velo;
	float amplitude;
	long distribution;
	long window;
    long injectionrate;
	long sinkdepth;
	long src_at_rcv; /* Indicates that wavefield should be injected at receivers */
} srcPar;

typedef struct _shotPar { /* Shot Parameters */
	long n;
	long ny;
	long nx;
	long nz;
	long *z;
	long *x;
	long *y;
	float *zs;
	float *xs;
	float *ys;
} shotPar;

typedef struct _raypar { /* ray-tracing parameters */
    long smoothwindow;
    long useT2;
    long geomspread;
    long nray;
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

