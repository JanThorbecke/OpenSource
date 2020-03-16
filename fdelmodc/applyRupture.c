#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc.h"

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

float stresscalc(float Tu, float Ts, float D, float dD, float D0, float V0);
void slipcalc(float *D, float *dD, float *vx, float dx, int izsrc, int ixsrc1, int ixsrc2, int n1);

void applyRupture(modPar mod, int itime, float *vx, float *txz, int verbose)
{

    float           *D, *dD, Tout;
    int             irup;
    static float    zrup, xrup1, xrup2, dT, T0, Tu, Ts, D0, V0, dx, trup;
    static int      izrup, ixrup1, ixrup2, itrup, nrup, n1, first=1;

    if (first==1) {
        if(!getparfloat("zrup", &zrup)) zrup = 0.0;
        if(!getparfloat("xrup1", &xrup1)) xrup1 = 0.0;
        if(!getparfloat("xrup2", &xrup2)) xrup2 = 0.0;
        if(!getparfloat("trup", &trup)) trup = 0.0;
        if(!getparfloat("V0", &V0)) V0 = 0.0;
        if(!getparfloat("D0", &D0)) D0 = 0.0;
        if(!getparfloat("Tu", &Tu)) Tu = 0.0;
        if(!getparfloat("Ts", &Ts)) Ts = 0.0;
        if(!getparfloat("dT", &dT)) dT = 0.0;

        dx = mod.dx;
        n1 = mod.naz;

        ixrup1 = NINT(xrup1/dx);
        ixrup2 = NINT(xrup2/dx);
        itrup = NINT(trup/mod.dt);
        nrup = ixrup2 - ixrup1;
        first = 0;
    }

    if (itime>itrup) {
        D   = (float *)calloc(nrup,sizeof(float));
        dD  = (float *)calloc(nrup,sizeof(float));

        slipcalc(D, dD, vx, dx, izrup, ixrup1, ixrup2, n1);

        for (irup=ixrup1; irup<ixrup2; irup++){
            vmess("D[%d]=%.3e",D[irup-ixrup1]);
            if (fabs(dD[irup-ixrup1])>0.001) {
                Tout = stresscalc(Tu, Ts, D[irup-ixrup1], dD[irup-ixrup1], D0, V0);
                txz[irup*n1+izrup-1] += Tout - T0;
                txz[irup*n1+izrup] += Tout - T0;
            }
            else {
                vx[irup*n1+izrup] = 0.0;
            }
        }

        free(D); free(dD);
    }
    else if (itime==itrup) {
        for (irup=ixrup1; irup<ixrup2; irup++){
            txz[irup*n1+izrup-1] += dT;
            txz[irup*n1+izrup] += dT;
        }
    }

}

float stresscalc(float Tu, float Ts, float D, float dD, float D0, float V0)
{
	float   T1, T2, Tout;

    T1 = Tu*(1-D/D0);
    T2 = Ts*(V0/(V0+dD));

    if (fabs(T1) > fabs(T2)) Tout=T1;
    else Tout = T2;
	
    return Tout;
}

void slipcalc(float *D, float *dD, float *vx, float dx, int izsrc, int ixsrc1, int ixsrc2, int n1)
{
    int     ix;

    for (ix=ixsrc1; ix<ixsrc2; ix++) {
        dD[ix-ixsrc1] = vx[ix*n1+izsrc+1] - vx[ix*n1+izsrc-1];
    }

    for (ix=ixsrc1+1; ix<ixsrc2; ix++) {
        D[ix-ixsrc1] = (dx/2.0)*(dD[ix-ixsrc1-1]+dD[ix-ixsrc1]);
    }
    D[0] = (dx/2.0)*dD[0];

}