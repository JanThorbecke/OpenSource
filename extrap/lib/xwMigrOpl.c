#include "optim.h"
#include "genfft.h"
#include "par.h"
#include <float.h>

void readtable_opt(complex *oper, float k, int *hopl);

int extrapEdge(complex *data, complex *opx, complex *tmp, int ntap, int hoplen, int mode, int i0);

void xwZoMigr(float *data, int nx, int nt, float dt, float *velmod, int nxm, int nzm, int ixa, int ixb, float fmin, float fmax, float *xrcv, int izrcv, float ox, float dxm, int opl_max, int ntap, int conjg, int ndepth, float *image, int verbose, float *exrcv, int ndepthex);

void xwMigrOpl(float *data, int nx, int nt, float dt, float *velmod1, float *velmod2, int nxm, int nzm, int ixa, int ixb, float fmin, float fmax, float *wavelet, int ntw, int nxw, float *xareal, int izsrc, float *xrcv, int izrcv, float ox, float dxm, int opl_max, int ntap, int conjg, int conjgs, int ndepth, float eps_r, float eps_a, float *image, int imc, int verbose, float *exsrc, float *exrcv, int ndepthex, int zomigr)
{
	int     iomin, iomax, iom, niom, ix, jx, d, hopl, hopl2, i, j;
	int     index1, i1, i2, nfreq, optn, lenx, hoplen, sign, slen;
	int     ixrcv, ixsrc, ixmin, ixmax, ixo, ixn;
	float   dom, om, c, cprev, df, sr, max_e, eps, smooth;
	float   *taper, scl, *locima, *locima2, *tmpim, *tmpim2, *pdata;
	float 	*trace, *smoother;
    float   c1prev, c2prev;
    int     hoplen1, hoplen2;
	complex *ctrace;
	float   t0, t1;
	complex *opx, *opxs, *cdata, *csrc, *tmp1, *tmp2;
	complex wa, *ctmp, da, *locdat, *locsrc;
	complex *cexsrc=(complex *) exsrc;
	complex *cexrcv=(complex *) exrcv;

	if (zomigr) {
		xwZoMigr(data, nx, nt, dt, velmod1, nxm, nzm, ixa, ixb, fmin, fmax, 
		        xrcv, izrcv, ox, dxm, opl_max, ntap, conjg, ndepth, image, 
                verbose, exrcv, ndepthex);
		return;
	}

/* define some constants */

	optn  = optncr(MAX(nt, ntw));
	nfreq = optn/2 + 1;
	df    = 1.0/(optn*dt);
	dom   = 2.*M_PI*df;
	iomin = (int)MIN((fmin*dt*optn), (nfreq-1));
	iomin = MAX(iomin, 1);
	iomax = MIN((int)(fmax*dt*optn), (nfreq-1));
	niom  = iomax-iomin+1;
	hopl  = (opl_max+1)/2;
	hopl2 = hopl-1;

/* transformation of shot record to frequency domain  */

    trace  = (float *)calloc(optn,sizeof(float));
    ctrace = (complex *)malloc(optn*sizeof(complex));
	cdata  = (complex *)calloc(nxm*nfreq, sizeof(complex));
	csrc   = (complex *)calloc(nxm*nfreq, sizeof(complex));

	if (conjg) scl = -1.0;
	else scl = 1.0;

    sign  = -1;
	max_e = 0.0;
	for (ix = 0; ix < nx; ix++) {

		memcpy(trace,&data[ix*nt],nt*sizeof(float));
        if (optn > nt) 
			memset( &trace[nt], 0, sizeof(float)*(optn-nt) );
            
        rc1fft(trace,ctrace,optn,sign);

   		ixrcv = NINT((xrcv[ix]-ox)/dxm);
		if ((ixrcv < 0 || ixrcv > nxm-1)) {
			if (verbose>1) fprintf(stderr,"xwMigrOpl: ixrcv %f (%d) outside model\n", xrcv[ix], ixrcv);
			continue;
		}
        for (iom=0; iom<nfreq; iom++) {
			/* Determine total energy in shot record */
			max_e += ctrace[iom].r*ctrace[iom].r;
			max_e += ctrace[iom].i*ctrace[iom].i;

			/* positioning of shot record into velocity model */
			cdata[iom*nxm+ixrcv].r = ctrace[iom].r;
			cdata[iom*nxm+ixrcv].i = ctrace[iom].i*scl;
        }
	}

	if (max_e < 100*FLT_MIN) {
		vwarn("xwMigrOpl: shot record does not contain any energy and is skipped !");
		return;
	}

	eps = max_e*eps_r + eps_a;

	if (verbose >=2) {
		vmess("xwMigrOpl: energy in shot = %.3e", max_e);
		vmess("xwMigrOpl: eps_r = %.2e eps_a = %.2e eps = %.2e", eps_r, eps_a, eps);
		vmess("xwMigrOpl: source depth izsrc = %d receiver depth izrcv = %d", izsrc, izrcv);
	}

/* transform wavelet and place at position within src array*/ 

	if (conjgs) scl = -1.0;
	else scl = 1.0;

	ixo = nxm;
	ixn = 0;
	max_e = 0.0;
	for (ix = 0; ix < nxw; ix++) {
		memcpy(trace,&wavelet[ix*ntw],ntw*sizeof(float));
        if (optn > ntw) 
            memset( &trace[ntw], 0, sizeof(float)*(optn-ntw) );
            
        rc1fft(trace,ctrace,optn,sign);

		ixsrc = NINT((xareal[ix]-ox)/dxm);
		if (ixsrc < ixo) ixo = MAX(0,ixsrc);
		if (ixsrc > ixn) ixn = MIN(nxm-1,ixsrc);
		if (ixsrc < 0 || ixsrc > nxm-1) {
			fprintf(stderr,"xwMigrOpl: ixsrc %f (%d) outside model\n", xareal[ix], ixsrc);
			continue;
		}
        for (iom=0; iom<nfreq; iom++) {
			/* Determine total energy in source record */
			max_e += ctrace[iom].r*ctrace[iom].r;
			max_e += ctrace[iom].i*ctrace[iom].i;

			csrc[iom*nxm+ixsrc].r = ctrace[iom].r;
			csrc[iom*nxm+ixsrc].i = ctrace[iom].i*scl;
        }
	}
    free(ctrace);
    free(trace);

	if (max_e < 100*FLT_MIN) {
		vwarn("xwMigrOpl: source record does not contain any energy and is skipped !");
		return;
	}

/* determine aperture to be calculated */

	for (ix = 0; ix < nx; ix++) {
		ixrcv = NINT((xrcv[ix]-ox)/dxm);
		if (ixrcv < ixo) ixo = ixrcv;
		if (ixrcv > ixn) ixn = ixrcv;
	}
	ixmin = MAX(0, ixo-ixb-1);
	ixmax = MIN(ixn+ixa+1, nxm-1);
	nx    = (ixmax-ixmin)+1;
	lenx  = 2*hopl2+nx;
	scl   = 2.0/nfreq;
	scl   = 1.0/(dt*dxm*dxm);

	if (verbose) {
		vmess("xwMigrOpl: calculation aperture: %.2f (%d) <--> %.2f (%d) (%d positions)", ixmin*dxm+ox, ixmin, ixmax*dxm+ox, ixmax, nx);
	}

/* allocate taper and local image arrays */

	if (ntap) {
		taper = (float *)malloc(ntap*sizeof(float));
		for (ix = 0; ix < ntap; ix++) {
			taper[ix] = exp(-1.0*(pow((2.0*(ntap-ix)/ntap), 2)));
//        	taper[ntap-ix-1] =(cos(PI*(ix)/(ntap))+1)/2.0;
//			fprintf(stderr,"taper %d = %f\n",ix, taper[ix]);
		}
		lenx  += 2*ntap;
	}
	if (imc == 2 || imc == 5) {
		tmpim  = (float *)calloc(nxm*ndepth, sizeof(float));
		tmpim2 = (float *)calloc(nxm*ndepth, sizeof(float));
	}

/* calculate image at depth = 0 */

	if(izsrc==izrcv && izsrc==0 ) {

	for (ix = ixmin; ix <= ixmax; ix++) {
		sr = 0.0; df = 0.0;
		for (iom = iomin; iom <= iomax; iom++) {
			wa = cdata[iom*nxm+ix];
			da = csrc[iom*nxm+ix];
			if (imc == 0) 
				image[ix*nzm+0] += scl*(da.r*wa.r+da.i*wa.i);
			else if (imc == 1) 
				image[ix*nzm+0] += scl*(da.r*wa.r+da.i*wa.i)/(da.r*da.r+da.i*da.i+eps);
			else if (imc == 2){ 
				df += (da.r*wa.r+da.i*wa.i);
				sr += da.r*da.r+da.i*da.i;
			}
		}
		if (imc == 2) 
			image[ix*nzm+0] += scl*df/(sr+eps);
	}
	}

	t0 = wallclock_time();

	tmp1  = (complex *)calloc(lenx, sizeof(complex));
	tmp2  = (complex *)calloc(lenx, sizeof(complex));
	opx   = (complex *)calloc(hopl, sizeof(complex));
	opxs  = (complex *)calloc(hopl, sizeof(complex));
	locdat  = (complex *)calloc(lenx, sizeof(complex));
	locsrc  = (complex *)calloc(lenx, sizeof(complex));
	locima  = (float *)calloc(nzm*nx, sizeof(float));
	if (imc == 2 || imc == 5) locima2 = (float *)calloc(nzm*nx, sizeof(float));

/* start extrapolation for all frequencies, depths and x-positions */

	for (iom = iomin; iom <= iomax; iom++) {
		memset(locdat, 0, sizeof(complex)*(lenx));
		memset(locsrc, 0, sizeof(complex)*(lenx));
		memset(tmp1, 0, sizeof(complex)*(lenx));
		memset(tmp2, 0, sizeof(complex)*(lenx));
		for (j = 0, ix = ixmin; j < nx; j++, ix++) {
			locdat[ntap+hopl2+j] = cdata[iom*nxm+ix];
			locsrc[ntap+hopl2+j] = csrc[iom*nxm+ix];
		}
		om    = iom*dom;
		d     = MIN(izsrc,izrcv);
		cprev = 0;
		c1prev = 0;
		c2prev = 0;

		/* if source at depth, first process only receivers */
		if (izsrc>izrcv) {
		if (verbose && iom==iomin) {
			vmess("xwMigrOpl: extrapolate receivers (%d) to src level (%d)", izrcv, izsrc);
		}

		for (d = izrcv; d < izsrc; d++) {

/* Extrapolation of data at the ntap edges of the model */

			c = velmod1[d*nxm+0];
			if (c!=0.0) readtable_opt(opx, om/c, &hoplen);
			if (ntap && c!=0.0) {
				 extrapEdge(locdat, opx, tmp1, ntap, hoplen, -1, hopl2);
			}

/* Extrapolation of data */

			for (ix = ixmin, jx = 0; ix <= ixmax; ix++, jx++) {

/* read operator from operator table (static array)*/

				c = velmod1[d*nxm+ix];
				if (c != cprev && c!=0.0) {
					readtable_opt(opx, om/c, &hoplen);
					cprev = c;
				}
				if (c==0.0) hoplen=0;

				wa.r = wa.i = 0.0;
				index1 = jx + hopl2 + ntap;
				for (j = 0; j < hoplen; j++) {
					i1 = index1+j; 
					i2 = index1-j;
					wa.r += (locdat[i1].r+locdat[i2].r)*opx[j].r;
					wa.r += (locdat[i1].i+locdat[i2].i)*opx[j].i;
					wa.i += (locdat[i1].i+locdat[i2].i)*opx[j].r;
					wa.i -= (locdat[i1].r+locdat[i2].r)*opx[j].i;
				}
				if (hoplen != 0) tmp1[index1] = wa;
				else tmp1[index1]=locdat[index1];
			}

/* Extrapolation of data at the ntap edges of the model */

//			c = velmod1[d*nxm+nxm-1];
//			readtable_opt(opx, om/c, &hoplen);
			if (ntap) {
				extrapEdge(locdat, opx, tmp1, ntap, hoplen, -1, hopl2+ntap+nx);
			}
			for (j = 0; j < lenx; j++) locdat[j] = tmp1[j];

			/* taper edges of extrapolated data */

			for (j = 0; j < ntap; j++) {
				locdat[hopl2+j].r *= taper[j];
				locdat[hopl2+j].i *= taper[j];

				locdat[lenx-j-hopl2-1].r *= taper[j];
				locdat[lenx-j-hopl2-1].i *= taper[j];
			}

			if(d==ndepthex) {
				if ( cexrcv ) {
					for (ix = 0; ix < nx; ix++) {
						cexrcv[iom*nxm+ix+ixmin] = locdat[hopl2+ntap+ix];
					}
				}
			}
		}
		}

		/* if receivers deeper than source, first extrapolate source */
		else if(izsrc<izrcv) {
			if (verbose && iom==iomin) {
				vmess("xwMigrOpl: extrapolate src (%d) to receiver level (%d)", izsrc, izrcv);
			}

		for (d = izsrc; d < izrcv; d++) {

/* Extrapolation of data at the ntap edges of the model */

			c = velmod2[d*nxm+0];
			if (c!=0.0) readtable_opt(opx, om/c, &hoplen);
			if (ntap && c!=0.0) {
				extrapEdge(locsrc, opx, tmp2, ntap, hoplen, 1, hopl2);
			}

/* Extrapolation of source */

			for (ix = ixmin, jx = 0; ix <= ixmax; ix++, jx++) {

/* read operator from operator table (static array)*/

				c = velmod2[d*nxm+ix];
				if (c != cprev && c!=0.0) {
					readtable_opt(opx, om/c, &hoplen);
					cprev = c;
				}
				if (c==0.0) hoplen=0;

				wa.r = wa.i = 0.0;
				index1 = jx + hopl2 + ntap;
				for (j = 0; j < hoplen; j++) {
					i1 = index1+j; 
					i2 = index1-j;
					wa.r += (locsrc[i1].r+locsrc[i2].r)*opx[j].r;
					wa.r -= (locsrc[i1].i+locsrc[i2].i)*opx[j].i;
					wa.i += (locsrc[i1].r+locsrc[i2].r)*opx[j].i;
					wa.i += (locsrc[i1].i+locsrc[i2].i)*opx[j].r;
				}
				if (hoplen != 0) tmp2[index1] = wa;
				else tmp2[index1]=locsrc[index1];
			}

/* Extrapolation of data at the ntap edges of the model */

//			c = velmod2[d*nxm+nxm-1];
//			readtable_opt(opx, om/c, &hoplen);
			if (ntap) {
				extrapEdge(locsrc, opx, tmp2, ntap, hoplen, 1, hopl2+ntap+nx);
			}
			for (j = 0; j < lenx; j++) locsrc[j] = tmp2[j];

			/* taper edges of extrapolated data */

			for (j = 0; j < ntap; j++) {
				locsrc[hopl2+j].r *= taper[j];
				locsrc[hopl2+j].i *= taper[j];

				locsrc[lenx-j-hopl2-1].r *= taper[j];
				locsrc[lenx-j-hopl2-1].i *= taper[j];
			}

			if(d==ndepthex) {
				if ( cexsrc ) {
					for (ix = 0; ix < nx; ix++) {
						cexsrc[iom*nxm+ix+ixmin] = locsrc[hopl2+ntap+ix];
					}
				}
			}

		}
		}

		if(ndepthex==0) {
			if ( cexsrc ) {
			    for (ix = 0; ix < nx; ix++) {
					cexsrc[iom*nxm+ix+ixmin] = locsrc[hopl2+ntap+ix];
				}
			}
			if ( cexrcv ) {
				for (ix = 0; ix < nx; ix++) {
					cexrcv[iom*nxm+ix+ixmin] = locdat[hopl2+ntap+ix];
				}
			}
		}

/* start extrapolation of both src and receiver arrays */


		for (; d < ndepth; d++) {
			 
	/* Extrapolation of data at the ntap edges of the model */

			c = velmod1[d*nxm+0];
			if (c!=0.0) readtable_opt(opx, om/c, &hoplen1);
			if (ntap && c!=0.0) {
				 extrapEdge(locdat, opx, tmp1, ntap, hoplen1, -1, hopl2);
			}
			c = velmod2[d*nxm+0];
			if (c!=0.0) readtable_opt(opxs, om/c, &hoplen2);
			if (ntap && c!=0.0) {
				 extrapEdge(locsrc, opxs, tmp2, ntap, hoplen2, 1, hopl2);
			}

			for (ix = ixmin, jx = 0; ix <= ixmax; ix++, jx++) {

				/* read operator from operator table (static array)*/

				c = velmod1[d*nxm+ix];
				if (c != c1prev && c != 0.0) {
					readtable_opt(opx, om/c, &hoplen1);
					c1prev = c;
				}
				if (c==0.0) hoplen1=0;

				c = velmod2[d*nxm+ix];
				if (c != c2prev && c != 0.0) {
					readtable_opt(opxs, om/c, &hoplen2);
					c2prev = c;
				}
				if (c==0.0) hoplen2=0;

				/* Extrapolation of data and source */

				wa.r = wa.i = 0.0;
				da.r = da.i = 0.0;
				index1 = jx + hopl2 + ntap;

    			for (j = 0; j < hoplen1; j++) {
        			i1 = index1+j;
        			i2 = index1-j;

        			wa.r += (locdat[i1].r+locdat[i2].r)*opx[j].r;
        			wa.r += (locdat[i1].i+locdat[i2].i)*opx[j].i;
        			wa.i += (locdat[i1].i+locdat[i2].i)*opx[j].r;
        			wa.i -= (locdat[i1].r+locdat[i2].r)*opx[j].i;
    			}

    			for (j = 0; j < hoplen2; j++) {
        			i1 = index1+j;
        			i2 = index1-j;

        			da.r += (locsrc[i1].r+locsrc[i2].r)*opxs[j].r;
        			da.r -= (locsrc[i1].i+locsrc[i2].i)*opxs[j].i;
        			da.i += (locsrc[i1].r+locsrc[i2].r)*opxs[j].i;
        			da.i += (locsrc[i1].i+locsrc[i2].i)*opxs[j].r;
    			}


				if (hoplen1 != 0) tmp1[index1] = wa;
				else tmp1[index1]=locdat[index1];
				if (hoplen2 != 0) tmp2[index1] = da;
				else tmp2[index1]=locsrc[index1];

				/* imaging condition */
	
				sr = (da.r*wa.r+da.i*wa.i);
				if (imc == 0) {
					locima[d*nx+jx] += sr;
				}
				else if (imc == 1) {
					locima[d*nx+jx] += sr/(da.r*da.r+da.i*da.i+eps);
				}
				else if (imc == 2 || imc == 5) {
					locima[d*nx+jx] += sr;
					locima2[d*nx+jx] += da.r*da.r+da.i*da.i;
				}
				else if (imc == 3) {
					locima[d*nx+jx] += da.r+wa.r;
				}
			}
			for (j = 0; j < lenx; j++) {
				locdat[j] = tmp1[j];
				locsrc[j] = tmp2[j];
			}
			if (ntap) {
				 extrapEdge(locdat, opx, tmp1, ntap, hoplen1, -1, hopl2+ntap+nx);
				 extrapEdge(locsrc, opxs, tmp2, ntap, hoplen2, 1, hopl2+ntap+nx);
			}

			for (j = 0; j < lenx; j++) {
				locdat[j] = tmp1[j];
				locsrc[j] = tmp2[j];
			}

			/* taper edges of extrapolated data */

			for (j = 0; j < ntap; j++) {
				locsrc[hopl2+j].r *= taper[j];
				locsrc[hopl2+j].i *= taper[j];
				locdat[hopl2+j].r *= taper[j];
				locdat[hopl2+j].i *= taper[j];

				locsrc[lenx-j-hopl2-1].r *= taper[j];
				locsrc[lenx-j-hopl2-1].i *= taper[j];
				locdat[lenx-j-hopl2-1].r *= taper[j];
				locdat[lenx-j-hopl2-1].i *= taper[j];
			}


			if(d==ndepthex-1) {
				if ( cexsrc ) {
					for (ix = 0; ix < nx; ix++) 
						cexsrc[iom*nxm+ix+ixmin] = locsrc[hopl2+ntap+ix];
				}
				if ( cexrcv ) {
					for (ix = 0; ix < nx; ix++) 
						cexrcv[iom*nxm+ix+ixmin] = locdat[hopl2+ntap+ix];
				}
			}

		}  /* end of depth loop */
	} /* end of iom loop */

	if (imc < 2) {
		for (d = 0; d < ndepth; d++) {
			for (ix = ixmin, jx = 0; ix <= ixmax; ix++, jx++) {
				image[ix*nzm+d+1] += scl*locima[d*nx+jx];
			}
		}
	}
	else if (imc == 2) {
		for (d = 0; d < ndepth; d++) {
			for (ix = ixmin, jx = 0; ix <= ixmax; ix++, jx++) {
				image[ix*nzm+d+1] += scl*locima[d*nx+jx]/(locima2[d*nx+jx]+eps);
			}
		}
	}
	else if (imc == 3) {
		for (d = 0; d < ndepth; d++) {
			for (ix = ixmin, jx = 0; ix <= ixmax; ix++, jx++) {
				image[ix*nzm+d+1] += scl*locima[d*nx+jx];
			}
		}
	}
	else if (imc == 5) {
		slen = hopl;
		smoother  = (float *)calloc((2*slen+1), sizeof(float));

		for (ix = 0; ix <=slen; ix++) {
        	smoother[slen+ix] =(cos(PI*(ix)/(slen))+1)/2.0;
		}
		for (ix = 0; ix <slen; ix++) {
        	smoother[ix] =smoother[2*slen-ix];
		}

		for (d = 0; d < ndepth; d++) {
			for (ix = ixmin, jx = 0; ix <= ixmax; ix++, jx++) {
				smooth = 0.0;
				i1 = MAX(jx-slen,0);
				i2 = MIN(jx+slen,nx);
				for (j=i1; j<i2; j++) {
					smooth += smoother[j-i1]*locima2[d*nx+j];
				}
				if (smooth > 1e3*FLT_MIN) {
					sr = locima[d*nx+jx];
					image[ix*nzm+d+1] += scl*locima[d*nx+jx]/(smooth+eps);
				}
			}
		}
		free(smoother);
	}

	free(opx);
	free(opxs);
	free(tmp1);
	free(tmp2);
	free(locdat);
	free(locsrc);
	free(locima);
	if (imc == 2 || imc == 5) {
		free(tmpim);
		free(tmpim2);
		free(locima2);
	}

	if(exsrc) wx2xt(cexsrc, exsrc, optn, nxm, nxm, optn);
	if(exrcv) wx2xt(cexrcv, exrcv, optn, nxm, nxm, optn);

	free(cdata);
	free(csrc);
	if (ntap) {
		free(taper);
	}

	t1 = wallclock_time();
	vmess("xwMigrOpl took: %f seconds", t1-t0);

	return;
}


