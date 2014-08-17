#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "genfft.h"
/*
To compile
gcc  -O3 -ffast-math -DDOUBLE -I../include -I. -o pmod PMOD.c -L/Users/jan/src/OpenSource/lib -ldgenfft -lm
gcc  -O3 -ffast-math -DDOUBLE -I../include -I. -o pmod PMOD.c -L/Users/jan/src/OpenSource/lib -lgenfft -lm
*/

REAL gauss1freq(REAL f, REAL freq)
{
    REAL   value;

    value = f*f/(2.0*freq*freq);
    value = sqrt(2.0*exp(1))*f*exp(-value)/(sqrt(2.0)*freq);

    return value;
}


void Woper3d(REAL k, REAL dx, REAL dy, REAL dz, int nkx, int nky, complex *oper)
{
       int       ikx, iky;
       REAL      kx, kx2, kz2, ky, ky2, k2;
       REAL      dkx, dky, kz;

       k2      = k*k;
       dkx = 2.0*M_PI/(nkx*dx);
       dky = 2.0*M_PI/(nky*dy);

       for (iky = 0; iky < nky/2; iky++) {
               ky   = (iky)*dky;
               ky2  = ky*ky;
               for (ikx = 0; ikx < nkx/2; ikx++) {
                       kx   = (ikx)*dkx;
                       kx2  = kx*kx;
                       kz2 = k2 - (kx2 + ky2);
                       if (kz2 >= 0.0) {
                               kz = dz*sqrt(kz2);
                               oper[iky*nkx+ikx].r = cos(kz);
                               oper[iky*nkx+ikx].i = 1.0*sin(kz);

                       }
                       else {
                               kz = dz*sqrt(-kz2);
                               oper[iky*nkx+ikx].r = exp(-kz);
                               oper[iky*nkx+ikx].i = 0.0;
                       }
					oper[iky*nkx+nkx-1-ikx] = oper[iky*nkx+ikx];
					oper[(nky-iky-1)*nkx+ikx] = oper[iky*nkx+ikx];
					oper[(nky-iky-1)*nkx+nkx-1-ikx] = oper[iky*nkx+ikx];
               }
       }

       return;
}

int main(int argc, char **argv)
{
        int      i, iz, np, nk, nlayers, npi, ip, ik, ikx, iky, ipz, iom, nom, iom2, ix, iy, nkx, nky, nt;
        REAL     dp, dk, dx1, dx2, dkx, dky, p, dom, omI, kappa, dt, dx, df, f;
        REAL     *mu, *lambda;
        complex  om, *qP, *qSV, *qSH, *cdum;
        complex  *L1SH, *L2SH, *L1TSH, *L2TSH, *LSH;
        complex  *L1PSV, *L2PSV, *L1TPSV, *L2TPSV, *LPSV;
        complex  *WpluslaySH, *WminlaySH, *WplussrcSH, *WminsrcSH, *WplusrcvSH, *WminrcvSH;
        complex  *WpluslayPSV, *WminlayPSV, *WplussrcPSV, *WminsrcPSV, *WplusrcvPSV, *WminrcvPSV;
        complex  *RpSH, *RmSH, *Pup, *Pdown, *PupplusSH, *PupminSH, *PdownplusSH, *PdownminSH;
        complex  *RpPSV, *RmPSV, *PupplusPSV, *PupminPSV, *PdownplusPSV, *PdownminPSV, *PPSV, *PSH, *P;
        complex  value;
        REAL     fp, t0;
        complex  *f1, *v1f1;
        REAL     kint, pint, a, b, kx, ky;
        int      icomp, sign;
        char     *filename1, *filename2;
        complex  kinter, *field, *ctrace;
        REAL     *trace;
        FILE     *fp_in, *fp_out, *fp_out1, *fp_out2;
        size_t   nwrite;
        REAL     k, dy, dz;
        REAL     y, scale, period, coi;
        int      jtot, nfreq;
        complex  *wavelet, *tmp, ctmp;
        REAL     fmin, fleft, fright, fmax, power;

    nt = 256; // total number of timesamples
    nom = nt/2+1;
    dt = 4e-3; // time sampling step
    dom = 2.0*M_PI/(nt*dt);
    df = 1.0/(nt*dt);
    omI = 0.5*dom;
    omI = 0.0;
    nkx = 512; 
    nky = nkx;


//Memory allocation
filename1 = (char *)calloc(100,sizeof(char));
filename2 = (char *)calloc(100,sizeof(char));
//kinter = (complex *)calloc(4*4*np*nom,sizeof(complex));
field = (complex *)calloc(nkx*nky*nom,sizeof(complex));
ctrace = (complex *)calloc(nom,sizeof(complex));
trace = (REAL *)calloc(nt,sizeof(REAL));
wavelet = (complex *)calloc(nom,sizeof(complex));
tmp = (complex *)calloc(nkx*nkx,sizeof(complex));

fprintf(stderr,"writing to disk\n");
for (icomp=10; icomp<11; icomp++) {  //ncomp=16!!!
    sprintf(filename1,"PPSV_FK_%d.bin",icomp+1);
    sprintf(filename2,"PPSV_time_%d.bin",icomp+1);
    // open output file for component icomp
    FILE *fp_out1 = fopen(filename1, "w");
    FILE *fp_out2 = fopen(filename2, "w");
    for (iom=0; iom<nom; iom++) {
        om.r = (iom)*dom;
        om.i = -omI;

        //inverse FFT

        // Phaseshift-Function
	k = om.r/2000.0;
	dx = 10.0; //dx = 5.0;
	dy = 10.0; //dy = 5.0;
	dz = 1000.0;
    dkx = 2.0*M_PI/(nkx*dx);
    dky = 2.0*M_PI/(nky*dy);
        
        
	Woper3d(k, dx, dy, dz, nkx, nky, &field[iom*nkx*nky]);
	//field[iom*nkx*nky].r = 1.0;

        sign = -1;
        cc2dfft(&field[iom*nkx*nky], nkx, nkx, nky, sign);
        for (iky = 0; iky<nky; iky++) {
            for (ikx = 0; ikx<nkx; ikx++) {
                tmp[ikx+iky*nkx].r = field[iom*nkx*nky+ikx+nkx*iky].r*dkx*dkx/(4*M_PI*M_PI);
                tmp[ikx+iky*nkx].i = field[iom*nkx*nky+ikx+nkx*iky].i*dkx*dkx/(4*M_PI*M_PI);
	    }
        }
/*
        for (iky = 0; iky<nky/2; iky++) {
            for (ikx = 0; ikx<nkx/2; ikx++) {

                field[iom*nkx*nky+nkx/2+ikx+(nky/2+iky)*nkx].r = tmp[ikx+nkx*iky].r;
                field[iom*nkx*nky+nkx/2+ikx+(nky/2+iky)*nkx].i = tmp[ikx+nkx*iky].i;
                field[iom*nkx*nky+ikx+(nky/2+iky)*nkx].r = tmp[ikx+nkx/2+nkx*iky].r;
                field[iom*nkx*nky+ikx+(nky/2+iky)*nkx].i = tmp[ikx+nkx/2+nkx*iky].i;
                field[iom*nkx*nky+ikx+nkx/2+nkx*iky].r = tmp[ikx+(nky/2+iky)*nkx].r;
                field[iom*nkx*nky+ikx+nkx/2+nkx*iky].i = tmp[ikx+(nky/2+iky)*nkx].i;
                field[iom*nkx*nky+ikx+nkx*iky].r = tmp[nkx/2+ikx+(nky/2+iky)*nkx].r;
                field[iom*nkx*nky+ikx+nkx*iky].i = tmp[nkx/2+ikx+(nky/2+iky)*nkx].i;
            }
        }
*/

        /* write to disk */
        nwrite = fwrite(&field[iom*nkx*nky].r, sizeof(complex), nkx*nky, fp_out1);
        assert(nwrite == nkx*nky);
    }
    fclose(fp_out1);
        fp=10; //peak value of frequency [Hz]
        t0=0.0;
        sign = -1;
    for (iy = 0; iy<nky; iy++) {
        for (ix = 0; ix<nkx; ix++) {
            for (iom=0; iom<nom; iom++) {
               om.r = (iom)*dom;
	       f = iom*df;
               om.i = 0.0;
               wavelet[iom].r = gauss1freq(f, fp);
//               wavelet[iom].r = 1.0;
//               wavelet[iom].i = 0.0;

               ctrace[iom] = field[ix+nkx*iy+iom*nkx*nky];
               // add wavelet to data
               ctmp.r = ctrace[iom].r*wavelet[iom].r -  ctrace[iom].i*wavelet[iom].i; //Multiplication defined as (a, b) * (c, d) = (ac - bd, ad + bc)
               ctmp.i = ctrace[iom].r*wavelet[iom].i +  ctrace[iom].i*wavelet[iom].r;
               ctrace[iom] = ctmp;
//               ctrace[iom] = wavelet[iom];
            }
	    cr1fft(&ctrace[0], trace, nt, sign);
	    /* write to disk */
//		for (i=0; i<nt; i++) trace[i] = 0.0;
//		trace[ix]=1.0;
            nwrite = fwrite(trace, sizeof(REAL), nt, fp_out2);
            assert(nwrite == nt);
	}
    }
    /* close output file of icomp */
    fclose(fp_out2);
}    
   return;
}

