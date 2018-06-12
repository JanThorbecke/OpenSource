#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <genfft.h>
#include "marchenko.h"
#include "raytime.h"

int omp_get_max_threads(void);
int omp_get_num_threads(void);
void omp_set_num_threads(int num_threads);

void applyMute( float *data, int *mute, int smooth, int above, int Nsyn, int nxs, int nt, int *xrcvsyn, int npossyn, int shift, int pad, int nt0);
int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);
double wallclock_time(void);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);

void scl_data(float *data, int nsam, int nrec, float scl, float *datout, int nsamout);
void pad_data(float *data, int nsam, int nrec, int nsamout, float *datout);
void convol(float *data1, float *data2, float *con, int nrec, int nsam, float dt, int shift);
void corr(float *data1, float *data2, float *cov, int nrec, int nsam, float dt, int shift);

void synthesis(complex *Refl, complex *Fop, float *Top, float *iRN, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, double *tfft, int first, int verbose);

void kxwfilter(float *data, int nt, int nx, float dt, float dx, float fmin, float fmax, float angle, float cp, float perc);
void timeDiff(float *data, int nsam, int nrec, float dt, float fmin, float fmax, int opt);
void depthDiff(float *data, int nsam, int nrec, float dt, float dx, float fmin, float fmax, float c, int opt);
void pad2d_data(float *data, int nsam, int nrec, int nsamout, int nrecout, float *datout);

void homogeneousg(float *HomG, float *green, complex *Refl, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, float *pmin, float *f1min, float *f1plus, float *f2p, float *G_d, int *muteW, int smooth, int shift, int above, int pad, int nt0, int *synpos, int verbose)
{
    int     i, j, l, ret, scheme,count=0;
    int     iter, niter, ix, kxwfilt;
	float   scl, *conv, *greenf2, fmin, fmax, alpha, cp, perc;
	float   *greenjkz, *greenf2jkz, *tmp1, *tmp2;
    double  t0, t2, tfft;
	FILE	*fp;

	if (!getparint("scheme", &scheme)) scheme = 0;
	if (!getparint("kxwfilt", &kxwfilt)) kxwfilt = 0;
	if (!getparint("niter",&niter)) niter=10;
	if (!getparfloat("fmin", &fmin)) fmin = 0.0;
	if (!getparfloat("fmax", &fmax)) fmax = 100.0;
	if (!getparfloat("alpha", &alpha)) alpha = 65.0;
	if (!getparfloat("cp", &cp)) cp = 1500.0;
	if (!getparfloat("perc", &perc)) perc = 0.15;

	tfft = 0.0;
	ret = 0;
    t0   = wallclock_time();
	scl	= 1.0/((float)npossyn);

	if (scheme==1) {
		if (verbose) vmess("Classical Homogeneous Green's function retrieval");
		greenf2   	= (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));
		greenf2jkz	= (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));
		greenjkz	= (float *)calloc(nxs*ntfft,sizeof(float));

		for (l = 0; l < Nsyn; l++) {
			if (niter<1) {
				for (i = 0; i < npossyn; i++) {
                    j = 0;
					ix = ixpossyn[i];
                    greenf2[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j];
                    greenf2jkz[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j];
                    for (j = 1; j < nts; j++) {
                        greenf2[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+nts-j];
                        greenf2jkz[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+nts-j];
                    }
                }
			}
			else {
            	for (i = 0; i < npossyn; i++) {
                	j = 0;
                	/* set green to zero if mute-window exceeds nt/2 */
                	if (muteW[l*nxs+ixpossyn[i]] >= nts/2) {
                    	memset(&greenf2[l*nxs*nts+i*nts],0, sizeof(float)*nt);
						memset(&greenf2jkz[l*nxs*nts+i*nts],0, sizeof(float)*nt);
                    	continue;
                	}
                	greenf2[l*nxs*nts+i*nts+j] = f2p[l*nxs*nts+i*nts+j] + pmin[l*nxs*nts+i*nts+j];
					greenf2jkz[l*nxs*nts+i*nts+j] = f2p[l*nxs*nts+i*nts+j] + pmin[l*nxs*nts+i*nts+j];
                	for (j = 1; j < nts; j++) {
                    	greenf2[l*nxs*nts+i*nts+j] = f2p[l*nxs*nts+i*nts+nts-j] + pmin[l*nxs*nts+i*nts+j];
						greenf2jkz[l*nxs*nts+i*nts+j] = f2p[l*nxs*nts+i*nts+nts-j] + pmin[l*nxs*nts+i*nts+j];
                	}
				}
            }
			depthDiff(&greenf2jkz[l*nxs*nts], ntfft, nshots, dt, dx, fmin, fmax, cp, 1);
        }
		applyMute(greenf2, muteW, smooth, 4, Nsyn, nxs, nts, ixpossyn, npossyn, shift, pad, nt0);
		applyMute(greenf2jkz, muteW, smooth, 4, Nsyn, nxs, nts, ixpossyn, npossyn, shift, pad, nt0);

		for (i = 0; i < npossyn; i++) {
			for (j = 0; j < nts; j++) {
				greenjkz[i*nts+j] = green[i*nts+j];
			}
		}
		depthDiff(greenjkz, ntfft, nshots, dt, dx, fmin, fmax, cp, 1);
	}
	else {
		if (verbose) vmess("Marchenko Homogeneous Green's function retrieval");
	}

#pragma omp parallel default(shared) \
  private(i,j,conv,tmp1,tmp2) 
{	
	conv	= (float *)calloc(nxs*ntfft,sizeof(float));
	if (scheme==1) {
		tmp1	= (float *)calloc(nxs*ntfft,sizeof(float));
		tmp2	= (float *)calloc(nxs*ntfft,sizeof(float));
	}

#pragma omp for
	for (l = 0; l < Nsyn; l++) {

		count+=1;

		if (verbose > 2) vmess("Creating Homogeneous G at location %d out of %d",count,Nsyn);
		if (scheme==0) { //Marchenko representation
			convol(green, &f2p[l*nxs*nts], conv, nxs, nts, dt, 0);
			if (kxwfilt) {
				kxwfilter(conv, ntfft, nshots, dt, dx, fmin, fmax, alpha, cp, perc);
			}
			for (i=0; i<npossyn; i++) {
            	j=0;
            	HomG[(j+nts/2)*Nsyn+synpos[l]] += scl*(conv[i*nts+j] + conv[i*nts+j]);
            	for (j=1; j<nts/2; j++) {
                	HomG[(j+nts/2)*Nsyn+synpos[l]] += scl*(conv[i*nts+j] + conv[i*nts+nts-j]);
					HomG[j*Nsyn+synpos[l]] += scl*(conv[i*nts+(j+nts/2)] + conv[i*nts+nts-(j+nts/2)]);
            	}
        	}
		}
		else if (scheme==1) { //classical representation
			corr(&greenf2[l*nxs*nts], greenjkz, tmp1, nxs, nts, dt, 0);
			corr(&greenf2jkz[l*nxs*nts], green, tmp2, nxs, nts, dt, 0);
			for (i = 0; i < npossyn; i++) {
				for (j = 0; j < nts; j++) {
					conv[i*nts+j] = -(tmp2[i*nts+j]-tmp1[i*nts+j]);
				}
			}
			timeDiff(conv, ntfft, nshots, dt, fmin, fmax, -1);
			for (i=0; i<npossyn; i++) {
                j=0;
                HomG[(j+nts/2)*Nsyn+synpos[l]] += scl*(conv[i*nts+j]);
                for (j=1; j<nts/2; j++) {
                    HomG[(j+nts/2)*Nsyn+synpos[l]] += scl*conv[i*nts+j];
					HomG[j*Nsyn+synpos[l]] += scl*conv[i*nts+(j+nts/2)];
                }
            }
		}
	}
    free(conv);
	if (scheme==1) free(tmp1);free(tmp2);
}
	if (scheme==1) free(greenf2);free(greenf2jkz);
		
    t2 = wallclock_time();
    if (verbose) {
        vmess("Total Homogeneous G time = %.3f", t2-t0);
    }

    return;
}

void corr(float *data1, float *data2, float *cov, int nrec, int nsam, float dt, int shift)
{
    int     i, j, n, optn, nfreq, sign;
    float   df, dw, om, tau, scl;
    float   *qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2;
    complex *cdata1, *cdata2, *ccov, tmp;

    optn = optncr(nsam);
    nfreq = optn/2+1;

    cdata1 = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (cdata1 == NULL) verr("memory allocation error for cdata1");
    cdata2 = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (cdata2 == NULL) verr("memory allocation error for cdata2");
    ccov = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (ccov == NULL) verr("memory allocation error for ccov");

    rdata1 = (float *)malloc(optn*nrec*sizeof(float));
    if (rdata1 == NULL) verr("memory allocation error for rdata1");
    rdata2 = (float *)malloc(optn*nrec*sizeof(float));
    if (rdata2 == NULL) verr("memory allocation error for rdata2");

    /* pad zeroes until Fourier length is reached */
    pad_data(data1, nsam, nrec, optn, rdata1);
    pad_data(data2, nsam, nrec, optn, rdata2);

    /* forward time-frequency FFT */
    sign = -1;
    rcmfft(&rdata1[0], &cdata1[0], optn, nrec, optn, nfreq, sign);
    rcmfft(&rdata2[0], &cdata2[0], optn, nrec, optn, nfreq, sign);

    /* apply correlation */
    p1r = (float *) &cdata1[0];
    p2r = (float *) &cdata2[0];
    qr  = (float *) &ccov[0].r;
    p1i = p1r + 1;
    p2i = p2r + 1;
    qi = qr + 1;
    n = nrec*nfreq;
    for (j = 0; j < n; j++) {
        *qr = (*p1r * *p2r + *p1i * *p2i);
        *qi = (*p1i * *p2r - *p1r * *p2i);
        qr += 2;
        qi += 2;
        p1r += 2;
        p1i += 2;
        p2r += 2;
        p2i += 2;
    }
    free(cdata1);
    free(cdata2);

    /* shift t=0 to middle of time window (nsam/2)*/
    if (shift) {
        df = 1.0/(dt*optn);
        dw = 2*PI*df;
        tau = dt*(nsam/2);

        for (j = 0; j < nrec; j++) {
            om = 0.0;
            for (i = 0; i < nfreq; i++) {
                tmp.r = ccov[j*nfreq+i].r*cos(om*tau) + ccov[j*nfreq+i].i*sin(om*tau);
                tmp.i = ccov[j*nfreq+i].i*cos(om*tau) - ccov[j*nfreq+i].r*sin(om*tau);
                ccov[j*nfreq+i] = tmp;
                om += dw;
            }
        }
    }

    /* inverse frequency-time FFT and scale result */
    sign = 1;
    scl = 1.0/(float)optn;
    crmfft(&ccov[0], &rdata1[0], optn, nrec, nfreq, optn, sign);
    scl_data(rdata1,optn,nrec,scl,cov,nsam);

    free(ccov);
    free(rdata1);
    free(rdata2);
    return;
}

void timeDiff(float *data, int nsam, int nrec, float dt, float fmin, float fmax, int opt)
{
    int     optn, iom, iomin, iomax, nfreq, ix, sign;
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
    iomin  = (int)MIN((omin/deltom), (nfreq));
    iomin  = MAX(iomin, 1);
    iomax  = MIN((int)(omax/deltom), (nfreq));

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
        if (opt > 0) {
            for (iom = iomin ; iom < iomax ; iom++) {
                om = deltom*iom;
                cdatascl[ix*nfreq+iom].r = -om*cdata[ix*nfreq+iom].i;
                cdatascl[ix*nfreq+iom].i = om*cdata[ix*nfreq+iom].r;
            }
        }
        else if (opt < 0) {
            for (iom = iomin ; iom < iomax ; iom++) {
                om = 1.0/(deltom*iom);
                cdatascl[ix*nfreq+iom].r = om*cdata[ix*nfreq+iom].i;
                cdatascl[ix*nfreq+iom].i = -om*cdata[ix*nfreq+iom].r;
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

void depthDiff(float *data, int nsam, int nrec, float dt, float dx, float fmin, float fmax, float c, int opt)
{
    int     optn, iom, iomin, iomax, nfreq, ix, ikx, nkx, ikxmax;
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

    iomin  = (int)MIN((omin/deltom), nfreq);
    iomin  = MAX(iomin, 0);
    iomax  = MIN((int)(omax/deltom), nfreq);

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

            ikxmax = MIN((int)(kp/dkx), nkx/2);

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
            ikxmax = MIN((int)(kp/dkx), nkx/2);
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

void pad2d_data(float *data, int nsam, int nrec, int nsamout, int nrecout, float *datout)
{
    int it,ix;
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
