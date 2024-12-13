#include <stdio.h>
#include <complex.h>
#include <math.h>

//SUBROUTINE ZBESH (ZR, ZI, FNU, KODE, M, N, CYR, CYI, NZ, IERR)
void cbesh_(float complex *k, float *FNU, int *KODE, int *M, int *N, float complex *CY, int *NZ, int *IERR);
void zbesh_(double *kr, double *ki, double *FNU, int *KODE, int *M, int *N, double *CYR, double *CYI, int *NZ, int *IERR);


void H12(float complex *cdata, float complex *cwave, int nx, float *xi, float *zi, int iomin, int iomax, int nfreq, float deltom, float rho, float c, float maxdip, float xsrc, float zsrc, float te, float ts)
{
    int ix, iom, KODE, M, N, NZ, IERR;
    double x, z, r, phi, scl, om, FNU, kr, ki, CYR, CYI;
    double complex Mc, k, CY, tmp;

    //if (verbose) vmess("near and far P field of monopole");
    FNU=1.0;
    KODE=1;
    M=2;
    N=1;
    for (ix = 0; ix < nx; ix++) {
        x = xi[ix] - xsrc;
        z = fabs(zi[ix] - zsrc);
        r = sqrt(x*x + z*z);
        if (r != 0) phi = acos(z/r);
        else phi = M_PI/2;
        scl = 0.25*rho;
        if (fabs(phi) < maxdip*M_PI/180.0) {
            for (iom = iomin; iom <= iomax; iom++) {
                om = iom*deltom;
                Mc = csqrt(c*c*(CMPLX(1.0,om*te)/CMPLX(1.0,om*ts)));
                k  = r*om/Mc;
                kr = creal(k);
                ki = cimag(k);
                zbesh_(&kr, &ki, &FNU, &KODE, &M, &N, &CYR, &CYI, &NZ, &IERR);
                //fprintf(stderr,"om=%f ierr=%d CY= %e %e\n", om, IERR, (CYR), (CYI));
                tmp = -scl*CMPLX(CYR, CYI);
                cdata[ix*nfreq+iom] = tmp*cwave[iom];
            }
        }
        else {
            for (iom = iomin; iom <= iomax; iom++) {
                cdata[ix*nfreq+iom] = 0.0;
            }
        }
    }

    return;
}
