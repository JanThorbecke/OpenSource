#include "optim.h"
#include "genfft.h"
#include "par.h"

void kwExtr(float *data, int nx, int nt, float dt, float *velmod, 
	int ndepth, float fmin, float fmax, int *xi, int *zi, int nrec, 
	int ntap, int conjg, float *extr, int nxm, 
	float ox, float dxm, float dz, float *xrcv, int verbose)
{
	int    	iomin, iomax, iom, ix, d, hopl, hopl2, i, j, ikx, nkx;
	int	index1, i1, i2, nfreq, optn, lenx, ixrcv,hoplen;
	int endkx;
	float  	dom, om, c, cprev, df, scl, dkx, scl2, k, k2, kx, kx2;
	float	*rdata, *taper, *pdata, kz2;
	complex	*cdata, *cextr, *tmp1, tmp, ez;

	optn  = optncr(nt);
	nfreq = optn/2 + 1;
	tmp1 = (complex *)malloc(nx*nfreq*sizeof(complex));

	/* pad zero's to get fourier length */
	if( nt != optn ) {
		if (verbose >1) vmess("kwExtr: padding zeros to data from %d to %d",nt,optn);
			pdata = (float *)calloc(optn*nx,sizeof(float));
			for (i=0; i<nx; i++) {
				for (j=0; j<nt; j++) pdata[i*optn+j] = data[i*nt+j];
				for ( ; j<optn; j++) pdata[i*optn+j] = 0.0;
			}
	}
	else {
		pdata = &data[0];
	}

	xt2wx(pdata, tmp1, optn, nx, optn, nx);

	if( nt != optn ) free(pdata);

/* positioning of shot record into velocity model */

	if (conjg) scl = -1.0;
	else scl = 1.0;

	cdata = (complex *)calloc(nxm*nfreq, sizeof(complex));
	for (iom = 0; iom < nfreq; iom++) {
		for (ix = 0; ix < nx; ix++) {
			ixrcv = NINT((xrcv[ix]-ox)/dxm);
			cdata[iom*nxm+ixrcv].r += tmp1[iom*nx+ix].r;
			cdata[iom*nxm+ixrcv].i += tmp1[iom*nx+ix].i*scl;
		}
	}
	free(tmp1);

/* define some constants */

	nx    = nxm;
	df    = 1.0/(optn*dt);
	dom   = 2.*PI*df;
	iomin = (int)MIN((fmin*dt*optn), (nfreq-1));
	iomin = MAX(iomin, 1);
	iomax = MIN((int)(fmax*dt*optn), (nfreq-1));
	nkx   = optncc(2*ntap+nxm);
	lenx  = nkx;
	ntap  = (nkx-nxm)/2;
	scl2  = 1.0/nkx;
	dkx   = 2.0*M_PI/(nkx*dxm);


	if(verbose) {
		fprintf(stderr,"    kwExtr: number of depth steps = %d\n", ndepth);
		fprintf(stderr,"    kwExtr: number of frequencies = %d\n", iomax-iomin+1);
	}

	cextr = (complex *)calloc(nrec*nfreq, sizeof(complex));
	taper = (float *)malloc(lenx*sizeof(float));
	for (ix = 0; ix < lenx; ix++) taper[ix] = 1.0;
	for (ix = 0; ix < ntap; ix++) {
		taper[ix] = exp(-1.0*(pow((0.4*(ntap-ix)/ntap), 2)));
		taper[lenx-1-ix] = taper[ix];
	}

	tmp1  = (complex *)calloc(lenx, sizeof(complex));

	for (iom = iomin; iom <= iomax; iom++) {
		om = iom*dom;
		cprev = 0;
		for (j = 0; j < nx; j++) {
			tmp1[ntap+j] = cdata[iom*nx+j];
		}
		for (j=0; j<nrec; j++) 
			if (zi[j] == 0) cextr[iom*nrec+j] = tmp1[ntap+xi[j]];


		for (d = 0; d < ndepth; d++) {

			/* transform to wavenumber domain */

			cc1fft(tmp1, nkx, 1);

			c = 0.0;
			for (ix = 0; ix < nx; ix++) c += velmod[d*nxm+ix];
			k = nx*om/c;
			k2 = k*k;

			/* kx = 0 */
			ez.r = cos(k*dz);
			ez.i = -sin(k*dz);

			tmp.r  = ez.r*tmp1[0].r;
			tmp.r += ez.i*tmp1[0].i;
			tmp.i  = ez.r*tmp1[0].i;
			tmp.i -= ez.i*tmp1[0].r;
			tmp1[0] = tmp;

			/* kx != 0 */
			endkx = MIN((int)(k/dkx),nkx/2);
        	for (ikx = 1; ikx <= endkx; ikx++) {
				kx  = ikx*dkx;
				kx2 = kx*kx;
				kz2 = k2 - kx2;

				ez.r = cos(sqrt(kz2)*dz);
				ez.i = -sin(sqrt(kz2)*dz);

				tmp.r  = ez.r*tmp1[ikx].r;
				tmp.r += ez.i*tmp1[ikx].i;
				tmp.i  = ez.r*tmp1[ikx].i;
				tmp.i -= ez.i*tmp1[ikx].r;
				tmp1[ikx] = tmp;

				tmp.r  = ez.r*tmp1[nkx-ikx].r;
				tmp.r += ez.i*tmp1[nkx-ikx].i;
				tmp.i  = ez.r*tmp1[nkx-ikx].i;
				tmp.i -= ez.i*tmp1[nkx-ikx].r;
				tmp1[nkx-ikx] = tmp;
	       	}

			cc1fft(tmp1, nkx, -1);

			for (j = 0; j < nkx; j++) {
				tmp1[j].r *= scl2*taper[j];
				tmp1[j].i *= scl2*taper[j];
			}
			for (j=0; j<nrec; j++) {
				if (zi[j] == (d+1)) {
					cextr[iom*nrec+j] = tmp1[ntap+xi[j]];
				}
			}

		}

	}

	free(tmp1);

	rdata = (float *)malloc(nrec*optn*sizeof(float));
	wx2xt(cextr, rdata, optn, nrec, nrec, optn);

	for (ix = 0; ix < nrec; ix++) {
		for (j = 0; j < nt; j++) extr[ix*nt+j] = rdata[ix*optn+j];
	}

	free(cdata);
	free(cextr);
	free(taper);
	free(rdata);

	return;
}
