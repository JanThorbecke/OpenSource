#include "optim.h"
#include "genfft.h"
#include "par.h"

void readtable_opt(complex *oper, float k, int *hopl);

int extrapEdge(complex *data, complex *opx, complex *tmp, int ntap, int hoplen, int mode, int i0);   

void xwCFP(float *data, int nx, int nt, float dt, float *velmod, 
	float fmin, float fmax, float *wavelet, complex *source, 
	int opl_max, int ntap, 
	int id0, int id1, int ds, int *xi, int *zi, int nrec, int izmax, int izmin,
	int plane_waves, int ixs, float dt_int, 
	int beam, float *beams, int verbose)
{
	int    	iomin, iomax, iom, ix, d, hopl, hopl2, j, k, i1, i2;
	int		index1, nfreq, optn, lenx, hoplen;
	float  	dom, om, c, cprev, df, dxfact, tom, scale;
	float	*rdata, *taper;
	complex	*opx, *cdata, *cextr, *tmp1, *tmp2, *cwave, tmp, wa;

	optn  = optncr(nt);
	nfreq = optn/2 + 1;
	df    = 1.0/(optn*dt);
	dom   = 2.*PI*df;
	iomin = MIN((int)(fmin*dt*optn), (nfreq-1));
	iomin = MAX(iomin, 1);
	iomax = MIN((int)(fmax*dt*optn), (nfreq-1));
	hopl  = (opl_max+1)/2;
	hopl2 = hopl-1;
	lenx  = 2*hopl2+nx+2*ntap;
	scale = 1.0/(float)(iomax-iomin+1);


	cwave = (complex *)malloc(nfreq*sizeof(complex));
	cdata = (complex *)calloc(nx*nfreq, sizeof(complex));
	cextr = (complex *)calloc(nrec*nfreq, sizeof(complex));
	tmp1  = (complex *)calloc(lenx, sizeof(complex));
	tmp2  = (complex *)calloc(lenx, sizeof(complex));
	opx   = (complex *)calloc(hopl, sizeof(complex));

	if(verbose>=2) {
		fprintf(stderr,"    xwCFP: number of depth steps = %d\n", abs(id1-id0));
		fprintf(stderr,"    xwCFP: number of frequencies = %d\n", iomax-iomin+1);
	}

	if (ntap) {
		taper = (float *)malloc(ntap*sizeof(float));
		for (ix = 0; ix < ntap; ix++) {
			taper[ix] = exp(-1.0*(pow((2.0*(ntap-ix)/ntap), 2)));
		}
	}

	if (plane_waves) {
		for (iom = iomin; iom <= iomax; iom++) {
			tom = dt_int*iom*dom;
			for (ix = 0; ix < nx; ix++) {
				k = ixs - ix;
				cdata[iom*nx+ix].r = cos(k*tom);
				cdata[iom*nx+ix].i = -sin(k*tom);
			}
			for (j = 0; j < ntap; j++) {
				cdata[iom*nx+j].r *= taper[j];
				cdata[iom*nx+j].i *= taper[j];
				cdata[iom*nx+nx-1-j].r *= taper[j];
				cdata[iom*nx+nx-1-j].i *= taper[j];
			}
		}
	}


	for (iom = iomin; iom <= iomax; iom++) {
		om = iom*dom;
		cprev = 0;

		for (j = 0; j < nx; j++) {
			tmp1[hopl2+ntap+j] = cdata[iom*nx+j];
			tmp2[hopl2+ntap+j] = cdata[iom*nx+j];
		}

/* first: extrapolate from deepest source up to lowest receiver level */

		for (d = izmax; d >= id0; d--) {

			if (!plane_waves) {
				for (j = 0; j < nx; j++) {
					tmp1[hopl2+ntap+j].r += source[d*nx+j].r;
					tmp1[hopl2+ntap+j].i += source[d*nx+j].i;
				}
			}

			for (j=0; j<nrec; j++) {
				if (zi[j] == izmax) {
					cextr[iom*nrec+j].r += tmp1[hopl2+ntap+xi[j]].r;
					cextr[iom*nrec+j].i += tmp1[hopl2+ntap+xi[j]].i;
				}
			}

			if (beam) {
            	for (ix = 0; ix < nx; ix++) {
                	beams[d*nx+ix] += sqrt(tmp1[hopl2+ntap+ix].r*tmp1[hopl2+ntap+ix].r+
                       	tmp1[hopl2+ntap+ix].i*tmp1[hopl2+ntap+ix].i)*scale;
            	}
			}

/* Extrapolation of data at the ntap edges of the model */

			c = velmod[d*nx+0];
			readtable_opt(opx, om/c, &hoplen);
			if (ntap) {
				extrapEdge(tmp1, opx, tmp2, ntap, hoplen, 1, hopl2);
			}

			for (ix = 0; ix < nx; ix++) {

				c = velmod[d*nx+ix];
				if (c != cprev) {
					readtable_opt(opx, om/c, &hoplen);
					cprev = c;
				}

				wa.r = wa.i = 0.0;
				index1 = ix + hopl2 + ntap;
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

/* Extrapolation of data at the ntap edges of the model */

			if (ntap) {
				extrapEdge(tmp1, opx, tmp2, ntap, hoplen, 1, nx+ntap+hopl2);
			}
			for (j = 0; j < lenx; j++) tmp1[j] = tmp2[j];

/* select receivers out of the extrapolated data */

			for (j=0; j<nrec; j++) {
				if (zi[j] == (d-1)) {
					cextr[iom*nrec+j].r += tmp2[hopl2+ntap+xi[j]].r;
					cextr[iom*nrec+j].i += tmp2[hopl2+ntap+xi[j]].i;
				}
			}

			for (j = 0; j < ntap; j++) {
				tmp1[hopl2+j].r *= taper[j];
				tmp1[hopl2+j].i *= taper[j];
				tmp1[lenx-j-hopl2-1].r *= taper[j];
				tmp1[lenx-j-hopl2-1].i *= taper[j];
			}
		} /* end of first depth loop */

/* second: extrapolate from lowest source to deepest receiver */ 

		memset(tmp1,0,lenx*sizeof(float));
		memset(tmp2,0,lenx*sizeof(float));

		for (j = 0; j < nx; j++) {
			tmp1[hopl2+ntap+j] = cdata[iom*nx+j];
			tmp2[hopl2+ntap+j] = cdata[iom*nx+j];
		}

		for (d = izmin; d <= id1; d++) {

			if (!plane_waves) {
				for (j = 0; j < nx; j++) {
					tmp1[hopl2+ntap+j].r += source[d*nx+j].r;
					tmp1[hopl2+ntap+j].i += source[d*nx+j].i;
				}
			}

			if (beam) {
            	for (ix = 0; ix < nx; ix++) {
                	beams[d*nx+ix] += sqrt(tmp1[hopl2+ntap+ix].r*tmp1[hopl2+ntap+ix].r+
                       	tmp1[hopl2+ntap+ix].i*tmp1[hopl2+ntap+ix].i)*scale;
            	}
			}

/* Extrapolation of data at the ntap edges of the model */

			c = velmod[d*nx+0];
			readtable_opt(opx, om/c, &hoplen);
			if (ntap) {
				extrapEdge(tmp1, opx, tmp2, ntap, hoplen, 1, hopl2);
			}

			for (ix = 0; ix < nx; ix++) {

				c = velmod[d*nx+ix];
				if (c != cprev) {
					readtable_opt(opx, om/c, &hoplen);
					cprev = c;
				}

				wa.r = wa.i = 0.0;
				index1 = ix + hopl2 + ntap;
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

/* Extrapolation of data at the ntap edges of the model */

			if (ntap) {
				extrapEdge(tmp1, opx, tmp2, ntap, hoplen, 1, nx+ntap+hopl2);
			}
			for (j = 0; j < lenx; j++) tmp1[j] = tmp2[j];

/* select receivers out of the extrapolated data */

			for (j=0; j<nrec; j++)
				if (zi[j] == (d+1)) {
					cextr[iom*nrec+j].r += tmp2[hopl2+ntap+xi[j]].r;
					cextr[iom*nrec+j].i += tmp2[hopl2+ntap+xi[j]].i;
				}

			for (j = 0; j < ntap; j++) {
				tmp1[hopl2+j].r *= taper[j];
				tmp1[hopl2+j].i *= taper[j];
				tmp1[lenx-j-hopl2-1].r *= taper[j];
				tmp1[lenx-j-hopl2-1].i *= taper[j];
			}
		} /* end of second depth loop */

	} /* end of iom loop */

	free(opx);
	free(tmp1);
	free(tmp2);

/* transform wavelet and convolve data with it*/ 

	rc1_fft(wavelet, cwave, optn, -1);

	for (iom = iomin; iom <= iomax; iom++) {
		for (ix = 0; ix < nrec; ix++) {
			tmp = cextr[iom*nrec+ix];
			cextr[iom*nrec+ix].r = tmp.r*cwave[iom].r - tmp.i*cwave[iom].i;
			cextr[iom*nrec+ix].i = tmp.r*cwave[iom].i + tmp.i*cwave[iom].r;
		}
	}

	rdata = (float *)malloc(optn*nrec*sizeof(float));
	wx2xt(cextr, rdata, optn, nrec, nrec, optn);

	dxfact = sqrt((zi[1]-zi[0])*(zi[1]-zi[0])+(xi[1]-xi[0])*(xi[1]-xi[0]));
	for (ix = 0; ix < nrec; ix++) {
		for (j = 0; j < nt; j++) data[ix*nt+j] = rdata[ix*optn+j]*dxfact;
	}

	free(cwave);
	if (ntap) free(taper);
	free(rdata);
	free(cdata);
	free(cextr);

	return;
}
