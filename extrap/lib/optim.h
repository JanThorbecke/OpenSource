#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
typedef struct _dcomplexStruct { /* complex number */
    double r,i;
} dcomplex;
#define COMPLEX
#endif

#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define ISODD(n) ((n) & 01)
#define ABS(x) ((x) < 0 ? -(x) : (x))

/* #define REM  defines Remez exchange option */
/* #define CFSQP_ON  defines CFSQP option */

/*********************** Double complex definitions *********************/

/*********************** 1D Extrapolation Operators **********************/

void forwExtr(complex *oper, float k, float dx, float dz, int nkx);
void forwExtr_smooth(complex *oper, float k, float dx, float dz, int nkx, float k1, float k2);
void invExtr(complex *oper, float k, float dx, float dz, int nkx);
void forwExtr_ph(complex *oper, float k, float dx, float dz, float alfa, 
	int nkx);
void invExtr_ph(complex *oper, float k, float dx, float dz, float alfa, 
	int nkx);
void onewayextr(complex *oper, float dx, int n, float dz, float k);


/*********************** 2D Extrapolation Operators *********************/

void circ_op(complex **hopx, int hoplx, int hoply, float dx, float dy, 
	float k, float dz);
void Woper3d(float k, float dx, float dy, float dz, int nkx, 
	int nky, complex **oper);
void Woper3d_ph(float k, float dx, float dy, float dz, int nkx2, 
	int nky2, float alpha, complex **oper);
void Kzoper3d(float om, float c, float dx, float dy, float dz, int nkx, 
	int nky, complex **oper);
void W3d(float kp, float dx, float dy, float dz, int nkx, int nky, 
	float alpha, complex **oper, float *wfacto);


/*********************** 1D optimization methods *********************/
/*--------- 1) Truncation and Windowing */

void KaiserWindow(complex *kxwoper, int nkx, complex *xwop, int oplength, 
	float beta);
void GaussWindow(complex *kxwoper, float dx, int nkx, complex *xwop, 
	int oplength, float end);
void phase_k0(complex *xwoper, int oplength, float *phase);
void trunc1D(complex *kxwop, complex *xwop, int oplength, int nkx);

/*--------- 2) L_2 norm optimization */

void weightfunct(float *weight, float k, float dx, int nkx, float alfa1, 
	float alfa2, float scale);
void kxwfilter(complex *data, float k, float dx, int nkx, float alfa1, 
	float alfa2, float perc);
void shoperror(complex *opkx, complex *opx, float *weight, int nkx, 
	int oplength, float dx, float *err2);
void shortoper(complex *kxwop, int nkx, complex *xwop, int opl, float dx,
	float kf, float alfa1_f, float alfa2_f, float perc, 
	float kw, float alfa1_w, float alfa2_w, float scale, int filter);
void shortoptime(complex *wavelet, int nt, float *optwave, int oplength, 
	float dt, float fmin, float fmax, float scale);
void wlsq1d(complex *kxop, int nkx, complex *xop, float *weight, int hopl);
void fsqp_oper(complex *oper, int nkx, float dx, float dz, float alpha, 
	int oplength, float k);

/*--------- 3) L_infinity norm optimization */

void remez(complex *opx, int oplx, float k, float alpha, float dx, float dz);
void Re_remez(float *oper, int nkx, int nkc, float mrp, int ntrans, 
	float *opx, int oplx);
void remez_diff(float *opx, int hoplx, float dx, float kmax);

/*********************** 2D optimization methods *********************/
/*--------- 1) Truncation and Windowing */

void trunc2D(complex **hopkx, int nkx, int nky, complex **opx, 
	int hoplx, int hoply);

/*--------- 2) L_2 norm optimization */

void wlsq2dr(float **opkx, int nkx, int nky, float **opx, int hoplx, 
	int hoply, float dx, float dy, float k, float alpha, float wfact);
void wlsq2dc(complex **opkx, int nkx, int nky, complex **opx, int hoplx, 
	int hoply, float dx, float dy, float k, float alpha, float wfact);
void wlsq2d8c(complex **opkx, int nkx, complex **opx, int hoplx, 
	float dx, float k, float alpha, float wfact);
void wlsq2d_kkx(float *opkx, int nkx, int nk, float *opx, int oplx, int oplt,
	float dx, float dt, float alpha, float c, float perc, float wfact);
void operror(complex **kxop, complex **hopkx, float k, float alpha, int nkx,
	int nky, float dx, float *ErrL2, float *EaLin, float *EpLin);

/*--------- 3) Circular symmetrical operators */

void circ2D(complex *opx1d, int nkx, complex **hopx, int oplx, float dx);
void circ2Dr(float *hopkxi, int nkx, int nky, float **hopx, int hoplx, 
	int hoply, float dx, float dy, float k, float alpha, float wfact);

/*--------- 4) McClellan Transformation */

void McCoeff(int order, int ncnt, float maxk, float *coeff);
void HazRedCoeff(float w1c, float *coeff, float *wc);
void InvHazRed(float kc, float *kxc);

/*--------- 5) Series Expansion */

void coeffKz(float **oper, int nkx, int nky, float dx, float dy, 
	float dz, float k, float alpha, complex *x, int order, float wfact);
void spectrumKz(float **opx, int oplx, int oply, int nkx, int nky, 
	float **opkx);
void Kzseries(float k, float dx, int nkx, int oplx, float dy, int nky, 
	int oply, float dz, float alpha, int order, complex *a_m, 
	float **opxkz, float wfacto, float wfacts);

/*********************** Solution of Matrix Equations *********************/

void cSolveAxb(complex **a, int nrow, int ncol, complex *b, complex *x);
void SolveAxb(float **a, int nrow, int ncol, float *b, float *x);


/*********************** miscelaneous functions *********************/
void spline3(float x1, float x2, float z1, float z2, float dzdx1,
	     float dzdx2, float *a, float *b, float *c, float *d);
complex froot(float x);
double wallclock_time(void);
