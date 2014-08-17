#include "optim.h"
#include "genfft.h"
#include "par.h"

void readtable_opt(complex *oper, float k, int *hopl);

void xwVSP(float *data, int nx, int nt, float dt, float *xrcv, float *velmod, int nxm, int ndepth, float ox, float dxm, float fmin, float fmax, int opl, int ntap, int *xi, int *zi, int nvsp, int ispr, int nrec, float *vsp, int verbose)
{
	int    	iomin, iomax, iom, ix, d, hopl, hopl2, hoplen, i, j;
	int	index1, i1, i2, nfreq, optn, lenx, ixp, izp, ixrcv;
	float  	dom, om, c, cprev, df;
	float	*rdata, *taper, *pdata;
	complex	*opx, *cdata, *tmp1, *tmp2, wa, *cvsp;

	optn  = optncr(nt);
	nfreq = optn/2 + 1;
	tmp1  = (complex *)malloc(nx*nfreq*sizeof(complex));

	/* pad zero's to get fourier length */
	if( nt != optn ) {
		if (verbose >1) vmess("xwVSP: padding zeros to data from %d to %d",nt,optn);
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

	cdata = (complex *)calloc(nxm*nfreq, sizeof(complex));
	for (iom = 0; iom < nfreq; iom++) {\
		for (ix = 0; ix < nx; ix++) {
			ixrcv = NINT((xrcv[ix]-ox)/dxm);
			cdata[iom*nxm+ixrcv].r += tmp1[iom*nx+ix].r;
			cdata[iom*nxm+ixrcv].i += tmp1[iom*nx+ix].i;
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
	hopl  = (opl+1)/2;
	hopl2 = hopl-1;
	lenx  = 2*hopl2+nx;

	if(verbose) {
		fprintf(stderr,"    xwVSP: number of depth steps = %d\n",ndepth);
		fprintf(stderr,"    xwVSP: number of frequencies = %d\n",iomax-iomin+1);
	}
	cvsp  = calloc(nrec*nfreq*nvsp, sizeof(complex));

	taper  = (float *)malloc(nx*sizeof(float));
	for (ix = 0; ix < nx; ix++) taper[ix] = 1.0;
	for (ix = 0; ix < ntap; ix++) {
		taper[ix] = exp(-1.0*(pow((0.4*(ntap-ix)/ntap), 2)));
		taper[nx-1-ix] = taper[ix];
	}

	for (iom = 0; iom < iomin; iom++) {
		for (ix = 0; ix < nx; ix++) {
			cdata[iom*nx+ix].r = 0.0;
			cdata[iom*nx+ix].i = 0.0;
		}
	}
	for (iom = iomax; iom < nfreq; iom++) {
		for (ix = 0; ix < nx; ix++) {
			cdata[iom*nx+ix].r = 0.0;
			cdata[iom*nx+ix].i = 0.0;
		}
	}

#if defined (SGI)
#pragma parallel
#pragma shared(cdata, cvsp)
#pragma byvalue(hopl, hopl2, velmod, lenx, iomin, iomax, taper, dom)
#pragma byvalue(nx, ndepth, xi, zi, nvsp, ispr, nrec, nfreq, verbose)
#pragma local(iom, tmp1, tmp2, wa, cprev, d, ix, j, c, om, index1)
#pragma local(i1, i2, opx, ixp, izp)
    { /* start of parallel region */
#endif

	tmp1  = (complex *)calloc(lenx, sizeof(complex));
	tmp2  = (complex *)calloc(lenx, sizeof(complex));
	opx   = (complex *)calloc(hopl, sizeof(complex));

#if defined(SGI)
#pragma pfor iterate(iom=iomin;iomax;1) schedtype (simple)
#endif
	for (iom = iomin; iom <= iomax; iom++) {
		om = iom*dom;
		cprev = 0;
		if (verbose>=2) fprintf(stderr,"    xwVSP: processing frequency %.3f\n", om/(2*PI));
		for (j = 0; j < nx; j++) {
			tmp1[hopl2+j].r = cdata[iom*nx+j].r*taper[j];
			tmp1[hopl2+j].i = cdata[iom*nx+j].i*taper[j];
		}

		for (d = 0; d < ndepth; d++) {

			for (j = 0; j < nvsp; j++) {
				for (ix = 0; ix < nrec; ix++) {
					ixp = xi[ix] + j*ispr;
					izp = zi[ix];
					if (izp == d) {
						cvsp[j*nfreq*nrec+iom*nrec+ix] = tmp1[hopl2+ixp];
					}
				}
			}


			for (ix = 0; ix < nx; ix++) {
				c = velmod[d*nxm+ix];
				if (c != cprev) {
					readtable_opt(opx, om/c, &hoplen);
					cprev = c;
				}

				wa.r = wa.i = 0.0;
				index1 = ix + hopl2;
                for (j = 0; j < hoplen; j++) {
					i1 = index1+j;
					i2 = index1-j;
					wa.r += (tmp1[i1].r + tmp1[i2].r)*opx[j].r;
					wa.r -= (tmp1[i1].i + tmp1[i2].i)*opx[j].i;
					wa.i += (tmp1[i1].i + tmp1[i2].i)*opx[j].r;
					wa.i += (tmp1[i1].r + tmp1[i2].r)*opx[j].i;
                }
				tmp2[index1] = wa;
			}
			for (ix = 0; ix < lenx; ix++) tmp1[ix] = tmp2[ix];
		}
		for (ix = 0; ix < nx; ix++) cdata[iom*nx+ix] = tmp2[hopl2+ix];
	}

	free(opx);
	free(tmp1);
	free(tmp2);

#if defined(SGI)
} /* end of parallel region */
#endif

	rdata  = (float *)malloc(nx*optn*sizeof(float));
	wx2xt(cdata, rdata, optn, nx, nx, optn);

	for (ix = 0; ix < nx; ix++) {
		for (j = 0; j < nt; j++) data[ix*nt+j] = rdata[ix*optn+j];
	}

	for (j = 0; j < nvsp; j++) {
		wx2xt(&cvsp[j*nfreq*nrec], (float *)&vsp[j*optn*nrec], optn, nrec, nrec, optn);
	}

	free(cdata);
	free(cvsp);
	free(taper);
	free(rdata);

	return;
}
