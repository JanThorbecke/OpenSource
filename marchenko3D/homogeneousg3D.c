#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <genfft.h>

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

double wallclock_time(void);

void scl_data(float *data, long nsam, long nrec, float scl, float *datout, long nsamout);
void pad_data(float *data, long nsam, long nrec, long nsamout, float *datout);
void pad2d_data(float *data, long nsam, long nrec, long nsamout, long nrecout, float *datout);
void conjugate(float *data, long nsam, long nrec, float dt);

void corr(float *data1, float *data2, float *cov, long nrec, long nsam, float dt, long shift);
void convol(float *data1, float *data2, float *con, long nrec, long nsam, float dt, long shift);

long readSnapData3D(char *filename, float *data, segy *hdrs, long nsnaps, long nx, long ny, long nz, long sx, long ex, long sy, long ey, long sz, long ez);
long getFileInfo3D(char *filename, long *n1, long *n2, long *n3, long *ngath, float *d1, float *d2, float *d3, float *f1, float *f2, float *f3, float *sclsxgxsygy, long *nxm);

//void kxwfilter(float *data, long nt, long nx, float dt, float dx, float fmin, float fmax, float angle, float cp, float perc);
void timeDiff(float *data, long nsam, long nrec, float dt, float fmin, float fmax, long opt);
void depthDiff(float *data, long nsam, long nrec, float dt, float dx, float fmin, float fmax, float c, long opt);

void homogeneousg3D(float *HomG, float *green, float *f2p, float *zsyn, long nx, long ny, long nt, float dx, float dy,
    float dt, long Nfoc, long *sx, long *sy, long *sz, long verbose)
{
    char    *file_inp;
    long    i, j, k, l, ret, scheme, count=0, n_source;
    long    is, kxwfilt, n1, n2, n3, ngath, ntraces, zsrc;
    float   d1, d2, d3, f1, f2, f3, scl;
	float   *conv, fmin, fmax, alpha, cp, rho, perc;
	float   *greenjkz, *input, *inputjkz, *tmp1, *tmp2;
    double  t0, t2, tfft;
    segy    *hdr_inp;

    if (!getparstring("file_inp", &file_inp)) file_inp = NULL;
	if (!getparlong("scheme", &scheme)) scheme = 0;
	if (!getparlong("kxwfilt", &kxwfilt)) kxwfilt = 0;
	if (!getparfloat("fmin", &fmin)) fmin = 0.0;
	if (!getparfloat("fmax", &fmax)) fmax = 100.0;
	if (!getparfloat("alpha", &alpha)) alpha = 65.0;
	if (!getparfloat("cp", &cp)) cp = 1000.0;
	if (!getparfloat("rho", &rho)) rho = 1000.0;
	if (!getparfloat("perc", &perc)) perc = 0.15;

	tfft = 0.0;
	ret = 0;
    t0   = wallclock_time();

    /* Read in the source function input and check the size */
    n_source=0;
    ret = getFileInfo3D(file_inp, &n1, &n2, &n3, &n_source, &d1, &d2, &d3, &f1, &f2, &f3, &scl, &ntraces);

    if (n1!=nt) verr("number of t-samples between input (%li) and retrieved (%li) is not equal",n1,nt);
    if (n2!=nx) verr("number of x-samples between input (%li) and retrieved (%li) is not equal",n2,nx);
    if (n3!=ny) verr("number of y-samples between input (%li) and retrieved (%li) is not equal",n3,ny);
    if (NINT(d1*1000.0)!=NINT(dt*1000.0)) verr("dt sampling between input (%.3e) and retrieved (%.3e) is not equal",d1,dt);
    if (NINT(d2*1000.0)!=NINT(dx*1000.0)) verr("dx sampling between input (%.3e) and retrieved (%.3e) is not equal",d2,dx);
    if (NINT(d3*1000.0)!=NINT(dy*1000.0)) verr("dy sampling between input (%.3e) and retrieved (%.3e) is not equal",d3,dy);

    input   = (float *)calloc(n_source*nx*ny*nt,sizeof(float));
    hdr_inp = (segy *)calloc(n_source*nx*ny,sizeof(segy));
    readSnapData3D(file_inp, input, hdr_inp, n_source, nx, ny, nt, 0, nx, 0, ny, 0, nt);

    zsrc = labs(hdr_inp[0].selev);
    for (l = 0; l < nx*ny; l++) {
        sx[l] = hdr_inp[0].sx;
        sy[l] = hdr_inp[0].sy;
        sz[l] = hdr_inp[0].sdepth;
    }
    for (l = 0; l < n_source; l++) {
        sx[l] = hdr_inp[l*nx*ny].sx;
        sy[l] = hdr_inp[l*nx*ny].sy;
        sz[l] = hdr_inp[l*nx*ny].sdepth;
    }

	if (scheme==1) { // Scale the Green's functions if the classical scheme is used
		if (verbose) vmess("Classical Homogeneous Green's function retrieval");
		greenjkz	= (float *)calloc(Nfoc*nx*ny*nt,sizeof(float));
		inputjkz	= (float *)calloc(n_source*nx*ny*nt,sizeof(float));

		for (l = 0; l < Nfoc; l++) {
            for (k = 0; k < ny; k++) {
			    depthDiff(&greenjkz[l*ny*nx*nt+k*nx*nt], nt, nx, dt, dx, fmin, fmax, cp, 1);
            }
        }

        for (l = 0; l < n_source; l++) {
            for (i = 0; i < ny*nx*nt; i++) {
                inputjkz[l*ny*nx*nt+i] = input[l*ny*nx*nt+i];
            }
            conjugate(&inputjkz[l*ny*nx*nt], nt, nx*ny, dt);
            conjugate(&input[l*ny*nx*nt], nt, nx*ny, dt);
            for (k = 0; k < ny; k++) {
                depthDiff(&inputjkz[l*ny*nx*nt+k*nx*nt], nt, nx, dt, dx, fmin, fmax, cp, 1);
            }
        }
	}
	else if (scheme==2) {
		if (verbose) vmess("Marchenko Green's function retrieval with G source");
	}
    else if (scheme==6) {
		if (verbose) vmess("Marchenko Green's function retrieval with f2 source");
	}
    else if (scheme==7) {
		if (verbose) vmess("Marchenko Green's function retrieval with source depending on position");
	}
	else if (scheme==3) {
		if (verbose) vmess("Marchenko Homogeneous Green's function retrieval with multiple sources");
        if (verbose) vmess("Looping over %li source positions",n_source);
	}
    else if (scheme==4) {
		if (verbose) vmess("Back propagation with multiple sources");
        if (verbose) vmess("Looping over %li source positions",n_source);
	}
	else if (scheme==5) {
		if (verbose) vmess("Marchenko Homogeneous Green's function retrieval with f2 source");
	}
	else {
		if (verbose) vmess("Marchenko Homogeneous Green's function retrieval with G source");
	}

#pragma omp parallel default(shared) \
  private(i,j,k,is,conv,tmp1,tmp2,zsrc) 
{	
	conv	= (float *)calloc(nx*nt,sizeof(float));
	if (scheme==1) {
		tmp1	= (float *)calloc(nx*nt,sizeof(float));
		tmp2	= (float *)calloc(nx*nt,sizeof(float));
	}
	if (scheme==3) tmp1 = (float *)calloc(nx*nt,sizeof(float));

#pragma omp for schedule(guided,1)
	for (l = 0; l < Nfoc; l++) {

		count+=1;
        zsrc = hdr_inp[0].selev;

        if (verbose>2) vmess("Creating Homogeneous G at location %li out of %li",l,Nfoc);

		if (scheme==0) { //Marchenko representation with G source
            for (k = 0; k < ny; k++) {
                depthDiff(&f2p[l*ny*nx*nt+k*nx*nt], nt, nx, dt, dx, fmin, fmax, cp, 1);
                convol(&input[k*nx*nt], &f2p[l*ny*nx*nt+k*nx*nt], conv, nx, nt, dt, 0);
                timeDiff(conv, nt, nx, dt, fmin, fmax, -3);
                if (kxwfilt) {
                    //kxwfilter(conv, nt, nx, dt, dx, fmin, fmax, alpha, cp, perc);
                }
                for (i=0; i<nx; i++) {
                    for (j=0; j<nt/2; j++) {
                        HomG[(j+nt/2)*Nfoc+l] += conv[i*nt+j]/rho;
                        HomG[j*Nfoc+l] += conv[i*nt+(j+nt/2)]/rho;
                    }
                }
            }
		}
        else if (scheme==5) { //Marchenko representation with f2 source
            for (k = 0; k < ny; k++) {
                depthDiff(&green[l*ny*nx*nt+k*nx*nt], nt, nx, dt, dx, fmin, fmax, cp, 1);
                convol(&input[k*nx*nt], &green[l*ny*nx*nt+k*nx*nt], conv, nx, nt, dt, 0);
                timeDiff(conv, nt, nx, dt, fmin, fmax, -3);
                if (kxwfilt) {
                    //kxwfilter(conv, nt, nx, dt, dx, fmin, fmax, alpha, cp, perc);
                }
                for (i=0; i<nx; i++) {
                    for (j=0; j<nt/2; j++) {
                        HomG[(j+nt/2)*Nfoc+l] += conv[i*nt+j]/rho;
                        HomG[j*Nfoc+l] += conv[i*nt+(j+nt/2)]/rho;
                    }
                }
            }
		}
		else if (scheme==1) { //classical representation
            for (k = 0; k < ny; k++) {
                convol(&greenjkz[l*ny*nx*nt+k*nx*nt], &input[k*nx*nt], tmp1, nx, nt, dt, 0);
                convol(&green[l*ny*nx*nt+k*nx*nt], &inputjkz[k*nx*nt], tmp2, nx, nt, dt, 0);
                for (i = 0; i < nx; i++) {
                    for (j = 0; j < nt; j++) {
                        conv[i*nt+j] = tmp1[i*nt+j]+tmp2[i*nt+j];
                    }
                }
                timeDiff(conv, nt, nx, dt, fmin, fmax, -1);
                for (i=0; i<nx; i++) {
                    for (j=0; j<nt/2; j++) {
                        HomG[(j+nt/2)*Nfoc+l] += conv[i*nt+j]/rho;
                        HomG[j*Nfoc+l] += conv[i*nt+(j+nt/2)]/rho;
                    }
                }
            }
		}
		else if (scheme==2) { //Marchenko representation without time-reversal G source
            for (k = 0; k < ny; k++) {
                depthDiff(&f2p[l*ny*nx*nt+k*nx*nt], nt, nx, dt, dx, fmin, fmax, cp, 1);
                convol(&input[k*nx*nt], &f2p[l*ny*nx*nt+k*nx*nt], conv, nx, nt, dt, 0);
                timeDiff(conv, nt, nx, dt, fmin, fmax, -1);
                for (i=0; i<nx; i++) {
                    for (j=0; j<nt/2; j++) {
                        HomG[(j+nt/2)*Nfoc+l] += 2*conv[i*nt+j]/rho;
                        HomG[j*Nfoc+l] += 2*conv[i*nt+(j+nt/2)]/rho;
                    }
                }
            }
		}
		else if (scheme==6) { //Marchenko representation without time-reversal f2 source
            for (k = 0; k < ny; k++) {
                depthDiff(&green[l*ny*nx*nt+k*nx*nt], nt, nx, dt, dx, fmin, fmax, cp, 1);
                convol(&input[k*nx*nt], &green[l*ny*nx*nt+k*nx*nt], conv, nx, nt, dt, 0);
                timeDiff(conv, nt, nx, dt, fmin, fmax, -1);
                for (i=0; i<nx; i++) {
                    for (j=0; j<nt/2; j++) {
                        HomG[(j+nt/2)*Nfoc+l] += 2*conv[i*nt+j]/rho;
                        HomG[j*Nfoc+l] += 2*conv[i*nt+(j+nt/2)]/rho;
                    }
                }
            }
		}
        else if (scheme==7) { //Marchenko representation without time-reversal using varying sources
            for (k = 0; k < ny; k++) {
                depthDiff(&f2p[l*ny*nx*nt+k*nx*nt], nt, nx, dt, dx, fmin, fmax, cp, 1);
                depthDiff(&green[l*ny*nx*nt+k*nx*nt], nt, nx, dt, dx, fmin, fmax, cp, 1);
                for (is=0; is<(n_source/2); is+=2) {
                    zsrc = labs(hdr_inp[is].selev);
                    if (zsrc > NINT(1000.0*zsyn[l])) {
                        if (verbose > 1) vmess("Homogeneous Green's function at %li uses G source (zsrc=%li)",NINT(1000.0*zsyn[l]));
                        convol(&input[is*ny*nx*nt+k*nx*nt], &f2p[l*ny*nx*nt+k*nx*nt], conv, nx, nt, dt, 0);
                    }
                    else {
                        if (verbose > 1) vmess("Homogeneous Green's function at %li uses f2 source (zsrc=%li)",NINT(1000.0*zsyn[l]));
                        convol(&input[(is+1)*ny*nx*nt+k*nx*nt], &green[l*ny*nx*nt+k*nx*nt], conv, nx, nt, dt, 0);
                    }
                    timeDiff(conv, nt, nx, dt, fmin, fmax, -1);
                    for (i=0; i<nx; i++) {
                        for (j=0; j<nt/2; j++) {
                            HomG[(j+nt/2)*Nfoc+l] += 2*conv[i*nt+j]/rho;
                            HomG[j*Nfoc+l] += 2*conv[i*nt+(j+nt/2)]/rho;
                        }
                    }
                }
            }
		}
		else if (scheme==3) { //Marchenko representation with multiple shot gathers
            for (k = 0; k < ny; k++) {
                depthDiff(&f2p[l*ny*nx*nt+k*nx*nt], nt, nx, dt, dx, fmin, fmax, cp, 1); 
                for (is=0; is<n_source; is++) {
                    convol(&input[is*ny*nx*nt+k*nx*nt], &f2p[l*ny*nx*nt+k*nx*nt], conv, nx, nt, dt, 0);
                    timeDiff(conv, nt, nx, dt, fmin, fmax, -3);
                    if (kxwfilt) {
                        //kxwfilter(conv, nt, nx, dt, dx, fmin, fmax, alpha, cp, perc);
                    }
                    for (i=0; i<nx; i++) {
                        for (j=0; j<nt/2; j++) {
                            HomG[is*nt*Nfoc+(j+nt/2)*Nfoc+l] += conv[i*nt+j]/rho;
                            HomG[is*nt*Nfoc+j*Nfoc+l] += conv[i*nt+(j+nt/2)]/rho;
                        }
                    }
                }
            }
        }
		else if (scheme==4) { //Marchenko representation with multiple shot gathers without time-reversal
            for (k = 0; k < ny; k++) {
                depthDiff(&f2p[l*ny*nx*nt+k*nx*nt], nt, nx, dt, dx, fmin, fmax, cp, 1);
                for (is=0; is<n_source; is++) {
                    convol(&input[is*ny*nx*nt+k*nx*nt], &f2p[l*ny*nx*nt+k*nx*nt], conv, nx, nt, dt, 0);
                    timeDiff(conv, nt, nx, dt, fmin, fmax, -1);
                    if (kxwfilt) {
                        //kxwfilter(conv, nt, nx, dt, dx, fmin, fmax, alpha, cp, perc);
                    }
                    for (i=0; i<nx; i++) {
                        for (j=0; j<nt/2; j++) {
                            HomG[is*nt*Nfoc+(j+nt/2)*Nfoc+l] += conv[i*nt+j]/rho;
                            HomG[is*nt*Nfoc+j*Nfoc+l] += conv[i*nt+(j+nt/2)]/rho;
                        }
                    }
                }
            }
        }
        
	}
    free(conv);
	if (scheme==1) { 
		free(tmp1);
		free(tmp2);
	}
	if (scheme==3) free(tmp1);
    if (scheme==4) free(tmp1);
}
	if (scheme==1) {
		free(input);
		free(inputjkz);
        free(greenjkz);
	}
    free(hdr_inp);
		
    t2 = wallclock_time();
    if (verbose) {
        vmess("Total Homogeneous G time = %.3f", t2-t0);
    }

    return;
}

void timeDiff(float *data, long nsam, long nrec, float dt, float fmin, float fmax, long opt)
{
    long     optn, iom, iomin, iomax, nfreq, ix, sign;
    float   omin, omax, deltom, om, df, *rdata, scl;
    complex *cdata, *cdatascl;

    optn = optncr(nsam);
    nfreq = optn/2+1;
    df    = 1.0/(optn*dt);

    cdata = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (cdata == NULL) verr("memory allocation error for cdata");

    rdata = (float *)malloc(optn*nrec*sizeof(float));
    if (rdata == NULL) verr("memory allocation error for rdata");

    /* pad zeroes until Fourier length is reached */
    pad_data(data,nsam,nrec,optn,rdata);

    /* Forward time-frequency FFT */
    sign = -1;
    rcmfft(&rdata[0], &cdata[0], optn, nrec, optn, nfreq, sign);

    deltom = 2.*PI*df;
    omin   = 2.*PI*fmin;
    omax   = 2.*PI*fmax;
    iomin  = (long)MIN((omin/deltom), (nfreq));
    iomin  = MAX(iomin, 1);
    iomax  = MIN((long)(omax/deltom), (nfreq));

    cdatascl = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (cdatascl == NULL) verr("memory allocation error for cdatascl");

    for (ix = 0; ix < nrec; ix++) {
        for (iom = 0; iom < iomin; iom++) {
            cdatascl[ix*nfreq+iom].r = 0.0;
            cdatascl[ix*nfreq+iom].i = 0.0;
        }
        for (iom = iomax; iom < nfreq; iom++) {
            cdatascl[ix*nfreq+iom].r = 0.0;
            cdatascl[ix*nfreq+iom].i = 0.0;
        }
        if (opt == 1) {
            for (iom = iomin ; iom < iomax ; iom++) {
                om = deltom*iom;
                cdatascl[ix*nfreq+iom].r = -om*cdata[ix*nfreq+iom].i;
                cdatascl[ix*nfreq+iom].i = om*cdata[ix*nfreq+iom].r;
            }
        }
        else if (opt == -1) {
            for (iom = iomin ; iom < iomax ; iom++) {
                om = 1.0/(deltom*iom);
                cdatascl[ix*nfreq+iom].r = om*cdata[ix*nfreq+iom].i;
                cdatascl[ix*nfreq+iom].i = -om*cdata[ix*nfreq+iom].r;
            }
        }
		else if (opt == -2) {
            for (iom = iomin ; iom < iomax ; iom++) {
                om = 4.0/(deltom*iom);
                cdatascl[ix*nfreq+iom].r = om*cdata[ix*nfreq+iom].r;
                cdatascl[ix*nfreq+iom].i = om*cdata[ix*nfreq+iom].i;
            }
        }
		else if (opt == -3) {
            for (iom = iomin ; iom < iomax ; iom++) {
                om = 1.0/(deltom*iom);
                cdatascl[ix*nfreq+iom].r = 2*om*cdata[ix*nfreq+iom].i;
                cdatascl[ix*nfreq+iom].i = 0.0;
            }
        }
    }
    free(cdata);

    /* Inverse frequency-time FFT and scale result */
    sign = 1;
    scl = 1.0/(float)optn;
    crmfft(&cdatascl[0], &rdata[0], optn, nrec, nfreq, optn, sign);
    scl_data(rdata,optn,nrec,scl,data,nsam);

    free(cdatascl);
    free(rdata);

    return;
}

void depthDiff(float *data, long nsam, long nrec, float dt, float dx, float fmin, float fmax, float c, long opt)
{
    long    optn, iom, iomin, iomax, nfreq, ix, ikx, nkx, ikxmax;
    float   omin, omax, deltom, df, dkx, *rdata, kx, scl;
    float   kx2, kz2, kp2, kp;
    complex *cdata, *cdatascl, kz, kzinv;

    optn  = optncr(nsam);
    nfreq = optncr(nsam)/2+1;
    df    = 1.0/(optn*dt);
    nkx   = optncc(nrec);
    dkx   = 2.0*PI/(nkx*dx);
    cdata = (complex *)malloc(nfreq*nkx*sizeof(complex));
    if (cdata == NULL) verr("memory allocation error for cdata");

    rdata = (float *)malloc(optn*nkx*sizeof(float));
    if (rdata == NULL) verr("memory allocation error for rdata");

    /* pad zeroes in 2 directions to reach FFT lengths */
    pad2d_data(data,nsam,nrec,optn,nkx,rdata);

    /* double forward FFT */
    xt2wkx(&rdata[0], &cdata[0], optn, nkx, optn, nkx, 0);

    deltom = 2.*PI*df;
    omin   = 2.*PI*fmin;
    omax   = 2.*PI*fmax;

    iomin  = (long)MIN((omin/deltom), nfreq);
    iomin  = MAX(iomin, 0);
    iomax  = MIN((long)(omax/deltom), nfreq);

    cdatascl = (complex *)malloc(nfreq*nkx*sizeof(complex));
    if (cdatascl == NULL) verr("memory allocation error for cdatascl");

    for (iom = 0; iom < iomin; iom++) {
        for (ix = 0; ix < nkx; ix++) {
            cdatascl[iom*nkx+ix].r = 0.0;
            cdatascl[iom*nkx+ix].i = 0.0;
        }
    }
    for (iom = iomax; iom < nfreq; iom++) {
        for (ix = 0; ix < nkx; ix++) {
            cdatascl[iom*nkx+ix].r = 0.0;
            cdatascl[iom*nkx+ix].i = 0.0;
        }
    }
    if (opt > 0) {
        for (iom = iomin ; iom <= iomax ; iom++) {
            kp = (iom*deltom)/c;
            kp2 = kp*kp;

            ikxmax = MIN((long)(kp/dkx), nkx/2);

            for (ikx = 0; ikx < ikxmax; ikx++) {
                kx  = ikx*dkx;
                kx2 = kx*kx;
                kz2 = kp2 - kx2;
                kz.r  = 0.0;
                kz.i  = sqrt(kz2);
                cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kz.r-cdata[iom*nkx+ikx].i*kz.i;
                cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kz.r+cdata[iom*nkx+ikx].r*kz.i;

            }
            for (ikx = ikxmax; ikx <= nkx-ikxmax+1; ikx++) {
                cdatascl[iom*nkx+ikx].r = 0.0;
                cdatascl[iom*nkx+ikx].i = 0.0;
            }
            for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
                kx  = (ikx-nkx)*dkx;
                kx2 = kx*kx;
                kz2 = kp2 - kx2;
                kz.r  = 0.0;
                kz.i  = sqrt(kz2);
                cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kz.r-cdata[iom*nkx+ikx].i*kz.i;
                cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kz.r+cdata[iom*nkx+ikx].r*kz.i;
            }
        }
    }
    else if (opt < 0) {
        for (iom = iomin ; iom < iomax ; iom++) {
            kp = iom*deltom/c;
            kp2 = kp*kp;
            ikxmax = MIN((long)(kp/dkx), nkx/2);
            for (ikx = 0; ikx < ikxmax; ikx++) {
                kx = ikx*dkx;
                kx2  = kx*kx;
                kz2 = kp2 - kx2;
                kzinv.r  = 0.0;
                kzinv.i  = -sqrt(kz2)/kz2;
                cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kzinv.r-cdata[iom*nkx+ikx].i*kzinv.i;
                cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kzinv.r+cdata[iom*nkx+ikx].r*kzinv.i;
            }
            for (ikx = ikxmax; ikx <= nkx-ikxmax+1; ikx++) {
                cdatascl[iom*nkx+ikx].r = 0.0;
                cdatascl[iom*nkx+ikx].i = 0.0;
            }
            for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
                kx = (ikx-nkx)*dkx;
                kx2  = kx*kx;
                kz2 = kp2 - kx2;
                kzinv.r  = 0.0;
                kzinv.i  = -sqrt(kz2)/kz2;
                cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kzinv.r-cdata[iom*nkx+ikx].i*kzinv.i;
                cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kzinv.r+cdata[iom*nkx+ikx].r*kzinv.i;
            }
        }
    }
    free(cdata);

    /* inverse double FFT */
    wkx2xt(&cdatascl[0], &rdata[0], optn, nkx, nkx, optn, 0);
    /* select original samples and traces */
    scl = 1.0;
    scl_data(rdata,optn,nrec,scl,data,nsam);

    free(cdatascl);
    free(rdata);

    return;
}

void pad2d_data(float *data, long nsam, long nrec, long nsamout, long nrecout, float *datout)
{
    long it,ix;
    for (ix=0;ix<nrec;ix++) {
        for (it=0;it<nsam;it++)
            datout[ix*nsam+it]=data[ix*nsam+it];
        for (it=nsam;it<nsamout;it++)
            datout[ix*nsam+it]=0.0;
    }
    for (ix=nrec;ix<nrecout;ix++) {
        for (it=0;it<nsamout;it++)
            datout[ix*nsam+it]=0.0;
    }
}

void conjugate(float *data, long nsam, long nrec, float dt)
{
    long     optn,  nfreq, j, ix, it, sign, ntdiff;
    float   *rdata, scl;
    complex *cdata;

    optn  = optncr(nsam);
    ntdiff = optn-nsam;
    nfreq = optn/2+1;

    cdata = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (cdata == NULL) verr("memory allocation error for cdata");

    rdata = (float *)malloc(optn*nrec*sizeof(float));
    if (rdata == NULL) verr("memory allocation error for rdata");

    /* pad zeroes until Fourier length is reached */
    pad_data(data,nsam,nrec,optn,rdata);

    /* Forward time-frequency FFT */
    sign = -1;
    rcmfft(&rdata[0], &cdata[0], optn, nrec, optn, nfreq, sign);

    /* take complex conjugate */
    for(ix = 0; ix < nrec; ix++) {
        for(j = 0; j < nfreq; j++) cdata[ix*nfreq+j].i = -cdata[ix*nfreq+j].i;
    }

    /* Inverse frequency-time FFT and scale result */
    sign = 1;
    scl = 1.0/(float)optn;
    crmfft(&cdata[0], &rdata[0], optn, nrec, nfreq, optn, sign);
    for (ix = 0; ix < nrec; ix++) {
        for (it = 0 ; it < nsam ; it++)
            data[ix*nsam+it] = scl*rdata[ix*optn+it+ntdiff];
    }
    //scl_data(rdata,optn,nrec,scl,data,nsam);
        
	free(cdata);
    free(rdata);

    return;
}
