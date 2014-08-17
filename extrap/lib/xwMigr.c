#include "optim.h"
#include "../FFT/genfft.h"
#include "par.h"

void xwMigr(float **data, int nx, int nt, float dt, float **velmod, int nxm, int ixa, int ixb, float fmin, float fmax, float **wavelet, int ntw, int nxw, float *xareal, int izsrc, float *xrcv, int izrcv, float ox, float dxm, int opl, int ntap, int conjg, int conjgs, int ndepth, float eps_r, float eps_a, float *image, int imc, int verbose, float *exsrc, float *exrcv, int ndepthex)
{
	int     iomin, iomax, iom, ix, jx, d, hopl, hopl2, i, j;
	int     index1, i1, i2, nfreq, optn, lenx;
	int     ixrcv, ixsrc, ixmin, ixmax, ixo, ixn;
	float   dom, om, c, cprev, df, sr, max_e, eps;
	float   *taper, scl, *locima, *locima2, *tmpim, *tmpim2, *pdata;
	float   t0, t1;
	complex *opx, *cdata, *csrc, *tmp1, *tmp2;
	complex wa, *ctmp, da, *locdat, *locsrc;

	complex *cexsrc=(complex *) exsrc;
	complex *cexrcv=(complex *) exrcv;

#if defined(SGI)
	int     np = mp_suggested_numthreads(0);

	if (verbose >=2) vmess("xwMigr: number of CPU's = %d", np);
#endif

/* transformation of shot record to frequency domain  */

	optn  = optncr(MAX(nt, ntw));
	nfreq = optn/2 + 1;
	ctmp  = (complex *)malloc(nx*nfreq*sizeof(complex));

        /* pad zero's if not a fourier length */
        if( nt != optn ) {
		if (verbose >1) vmess("xwMigr: padding zeros to data from %d to %d\n",nt,optn);
        		pdata = (float *)calloc(optn*nx,sizeof(float));
                for (i=0; i<nx; i++) {
                        for (j=0; j<nt; j++) pdata[i*optn+j] = data[i][j];
                        for (; j<optn; j++) pdata[i*optn+j] = 0.0;
                        }
                }
        else {
                pdata = &data[0][0];
                }

        xt2wx(pdata, ctmp, optn, nx, optn, nx);

        if( nt != optn ) free(pdata);

/* Determine maximum energy in wavefield */

	max_e = 0.0;
	for (iom = 0; iom < nfreq; iom++) {
		for (ix = 0; ix < nx; ix++) {
			max_e += ctmp[iom*nx+ix].r*ctmp[iom*nx+ix].r;
			max_e += ctmp[iom*nx+ix].i*ctmp[iom*nx+ix].i;
		}
	}
	eps = max_e*eps_r + eps_a;

	if (verbose >=2) {
		vmess("xwMigr: energy in shot = %.3e", max_e);
		vmess("xwMigr: eps_r = %.2e eps_a = %.2e eps = %.2e", eps_r, eps_a, eps);
	}

/* positioning of shot record into velocity model */

	if (conjg) scl = -1.0;
	else scl = 1.0;

	cdata = (complex *)calloc(nxm*nfreq, sizeof(complex));

	for (iom = 0; iom < nfreq; iom++) {
		for (ix = 0; ix < nx; ix++) {
   			ixrcv = NINT((xrcv[ix]-ox)/dxm);
			cdata[iom*nxm+ixrcv].r = ctmp[iom*nx+ix].r;
			cdata[iom*nxm+ixrcv].i = ctmp[iom*nx+ix].i*scl;
		}
	}
	free(ctmp);

/* transform wavelet and place at position within src array*/ 

	csrc = (complex *)calloc(nxm*nfreq, sizeof(complex));

	if (conjgs) scl = -1.0;
	else scl = 1.0;

	ctmp  = (complex *)malloc(nxw*nfreq*sizeof(complex));

        /* pad zero's if not a fourier length */
        if( ntw != optn ) {
		if (verbose >1) vmess("xwMigr: padding zeros to source from %d to %d\n",ntw,optn);
        		pdata = (float *)calloc(optn*nxw,sizeof(float));
                for (i=0; i<nxw; i++) {
                        for (j=0; j<ntw; j++) pdata[i*optn+j] = wavelet[i][j];
                        for (  ; j<optn; j++) pdata[i*optn+j] = 0.0;
                        }
                }
        else {
                pdata = &wavelet[0][0];
                }

        xt2wx(pdata, ctmp, optn, nxw, optn, nxw);

        if( ntw != optn ) free(pdata);


	ixo = nxm;
	ixn = 0;
	for (ix = 0; ix < nxw; ix++) {
		ixsrc = NINT((xareal[ix]-ox)/dxm);
		if (ixsrc < ixo) ixo = MAX(0,ixsrc);
		if (ixsrc > ixn) ixn = MIN(nxm-1,ixsrc);
		if (ixsrc < 0 || ixsrc > nxm-1) continue;
		for (iom = 0; iom < nfreq; iom++) {
			csrc[iom*nxm+ixsrc].r = ctmp[iom*nxw+ix].r;
			csrc[iom*nxm+ixsrc].i = ctmp[iom*nxw+ix].i*scl;
		}
	}
	free(ctmp);

/* determine aperture to be calculated */

	for (ix = 0; ix < nx; ix++) {
		ixrcv = NINT((xrcv[ix]-ox)/dxm);
		if (ixrcv < ixo) ixo = ixrcv;
		if (ixrcv > ixn) ixn = ixrcv;
	}
	ixmin = MAX(0, ixo-ixb-1);
	ixmax = MIN(ixn+ixa+1, nxm-1);
	nx    = (ixmax-ixmin)+1;

	if (verbose>=2) {
		vmess("xwMigr: calculation aperture: %.2f (%d) <--> %.2f (%d) (%d positions)\n", ixmin*dxm+ox, ixmin, ixmax*dxm+ox, ixmax, nx);
	}

/* define some constants */

	df    = 1.0/(optn*dt);
	dom   = 2.*PI*df;
	iomin = (int)MIN((fmin*dt*optn), (nfreq-1));
	iomin = MAX(iomin, 1);
	iomax = MIN((int)(fmax*dt*optn), (nfreq-1));
	hopl  = (opl+1)/2;
	hopl2 = hopl-1;
	lenx  = 2*hopl2+nx;
	scl   = 2.0/nfreq;
	scl   = 1.0/(dt*dxm*dxm);

	taper = (float *)malloc(ntap*sizeof(float));

	for (ix = 0; ix < ntap; ix++) {
		taper[ix] = exp(-1.0*(pow((2.0*(ntap-ix)/ntap), 2)));
	}

	if (imc == 2) {
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
				image[ix*ndepth+0] += scl*(da.r*wa.r+da.i*wa.i);
			else if (imc == 1) 
				image[ix*ndepth+0] += scl*(da.r*wa.r+da.i*wa.i)/(da.r*da.r+da.i*da.i+eps);
			else if (imc == 2){ 
				df += (da.r*wa.r+da.i*wa.i);
				sr += da.r*da.r+da.i*da.i;
			}
		}
		if (imc == 2) 
			image[ix*ndepth+0] += scl*df/(sr+eps);
	}
	}

	t0 = wallclock_time();
#if defined (SGI)
#pragma parallel
#pragma shared(image, tmpim, tmpim2, cexsrc, cexrcv)
#pragma byvalue(hopl, hopl2, velmod, lenx, iomin, iomax, taper, dom)
#pragma byvalue(nx, ndepth, ixmin, ixmax, eps, imc, csrc, cdata, scl)
#pragma byvalue(ndepthex,izsrc,izrcv)
#pragma local(iom, tmp1, tmp2, cprev, d, ix, j, c, om, index1)
#pragma local(sr, opx, jx, wa, da)
#pragma local(locdat, locsrc, locima, locima2, i1, i2)
	{ /* start of parallel region */
#endif

	tmp1  = (complex *)calloc(lenx, sizeof(complex));
	tmp2  = (complex *)calloc(lenx, sizeof(complex));
	opx   = (complex *)calloc(hopl, sizeof(complex));
	locdat  = (complex *)calloc(lenx, sizeof(complex));
	locsrc  = (complex *)calloc(lenx, sizeof(complex));
	locima  = (float *)calloc(ndepth*nx, sizeof(float));
	if (imc == 2) locima2 = (float *)calloc(ndepth*nx, sizeof(float));

/* start extrapolation for all frequencies, depths and x-positions */

#if defined(SGI)
#pragma pfor iterate(iom=iomin;iomax;1) schedtype (simple)
#endif
	for (iom = iomin; iom <= iomax; iom++) {
		for (j = 0, ix = ixmin; j < nx; j++, ix++) {
			locdat[hopl2+j] = cdata[iom*nxm+ix];
			locsrc[hopl2+j] = csrc[iom*nxm+ix];
		}
		om    = iom*dom;
		d     = MIN(izsrc,izrcv);
		cprev = 0;

	/* if source at depth, first process only receivers */
		if (izsrc>izrcv) {
		for (d = izrcv; d < izsrc; d++) {
			for (ix = ixmin, jx = 0; ix <= ixmax; ix++, jx++) {

/* read operator from operator table (static array)*/

				c = velmod[d][ix];
				if (c != cprev) {
					readtable(opx, om/c, hopl);
					cprev = c;
				}

/* Extrapolation of data */
				wa.r = wa.i = 0.0;
				index1 = jx + hopl2;

				for (j = 0; j < hopl; j++) {
					i1 = index1+j; 
					i2 = index1-j;
					wa.r += (locdat[i1].r+locdat[i2].r)*opx[j].r;
					wa.r += (locdat[i1].i+locdat[i2].i)*opx[j].i;
					wa.i += (locdat[i1].i+locdat[i2].i)*opx[j].r;
					wa.i -= (locdat[i1].r+locdat[i2].r)*opx[j].i;
				}
				tmp1[index1] = wa;
			}

			for (j = 0; j < lenx; j++) locdat[j] = tmp1[j];

			/* taper edges of extrapolated data */

			for (j = 0; j < ntap; j++) {
				locdat[j].r *= taper[j];
				locdat[j].i *= taper[j];

				locdat[lenx-j-1].r *= taper[j];
				locdat[lenx-j-1].i *= taper[j];
			}

			if(d==ndepthex) {
			  if ( cexrcv ) {
			    for (ix = 0; ix < nx; ix++) {
				cexrcv[iom*nxm+ix+ixmin] = locdat[ix];
				}
			    }
			  }

			}
		}

		/* if receivers deeper than source, first extrapolate source */
		else if(izsrc<izrcv) {

		for (d = izsrc; d < izrcv; d++) {
			for (ix = ixmin, jx = 0; ix <= ixmax; ix++, jx++) {

/* read operator from operator table (static array)*/

				c = velmod[d][ix];
				if (c != cprev) {
					readtable(opx, om/c, hopl);
					cprev = c;
				}

/* Extrapolation of source */
				da.r = da.i = 0.0;
				index1 = jx + hopl2;

				for (j = 0; j < hopl; j++) {
					i1 = index1+j; 
					i2 = index1-j;
					da.r += (locsrc[i1].r+locsrc[i2].r)*opx[j].r;
					da.r -= (locsrc[i1].i+locsrc[i2].i)*opx[j].i;
					da.i += (locsrc[i1].r+locsrc[i2].r)*opx[j].i;
					da.i += (locsrc[i1].i+locsrc[i2].i)*opx[j].r;
				}
				tmp2[index1] = da;
			}

			for (j = 0; j < lenx; j++) locsrc[j] = tmp2[j];

			/* taper edges of extrapolated data */

			for (j = 0; j < ntap; j++) {
				locsrc[j].r *= taper[j];
				locsrc[j].i *= taper[j];

				locsrc[lenx-j-1].r *= taper[j];
				locsrc[lenx-j-1].i *= taper[j];
			}

			if(d==ndepthex) {
			  if ( cexsrc ) {
			    for (ix = 0; ix < nx; ix++) {
				cexsrc[iom*nxm+ix+ixmin] = locsrc[ix];
				}
			    }
			  }

			}

		}

		for (; d < ndepth; d++) {

			for (ix = ixmin, jx = 0; ix <= ixmax; ix++, jx++) {

/* read operator from operator table (static array)*/

				c = velmod[d][ix];
				if (c != cprev) {
					readtable(opx, om/c, hopl);
					cprev = c;
				}

/* Extrapolation of data and source */

				wa.r = wa.i = 0.0;
				da.r = da.i = 0.0;
				index1 = jx + hopl2;

    			for (j = 0; j < hopl; j++) {
        			i1 = index1+j;
        			wa.r += locdat[i1].r*opx[j].r;
        			wa.r += locdat[i1].i*opx[j].i;
        			wa.i += locdat[i1].i*opx[j].r;
        			wa.i -= locdat[i1].r*opx[j].i;
        			da.r += locsrc[i1].r*opx[j].r;
        			da.r -= locsrc[i1].i*opx[j].i;
        			da.i += locsrc[i1].r*opx[j].i;
        			da.i += locsrc[i1].i*opx[j].r;
			
        			i2 = index1-j;
        			wa.r += locdat[i2].r*opx[j].r;
        			wa.r += locdat[i2].i*opx[j].i;
        			wa.i += locdat[i2].i*opx[j].r;
        			wa.i -= locdat[i2].r*opx[j].i;
        			da.r += locsrc[i2].r*opx[j].r;
        			da.r -= locsrc[i2].i*opx[j].i;
        			da.i += locsrc[i2].r*opx[j].i;
        			da.i += locsrc[i2].i*opx[j].r;
    			}

				tmp1[index1] = wa;
				tmp2[index1] = da;
		
/* imaging condition */
	
				sr = (da.r*wa.r+da.i*wa.i);
				if (imc == 0) {
					locima[d*nx+jx] += sr;
				}
				else if (imc == 1) {
					locima[d*nx+jx] += sr/(da.r*da.r+da.i*da.i+eps);
				}
				else if (imc == 2) {
					locima[d*nx+jx] += sr;
					locima2[d*nx+jx] += da.r*da.r+da.i*da.i;
				}
			}

			for (j = 0; j < lenx; j++) {
				locdat[j] = tmp1[j];
				locsrc[j] = tmp2[j];
			}
			/* taper edges of extrapolated data */

			for (j = 0; j < ntap; j++) {
				locdat[j].r *= taper[j];
				locdat[j].i *= taper[j];
				locsrc[j].r *= taper[j];
				locsrc[j].i *= taper[j];

				locdat[lenx-j-1].r *= taper[j];
				locdat[lenx-j-1].i *= taper[j];
				locsrc[lenx-j-1].r *= taper[j];
				locsrc[lenx-j-1].i *= taper[j];
			}

			if(d==ndepthex) {
			  if ( cexsrc ) {
			    for (ix = 0; ix < nx; ix++) {
				cexsrc[iom*nxm+ix+ixmin] = locsrc[ix];
				}
			    }
			  if ( cexrcv ) {
			    for (ix = 0; ix < nx; ix++) {
				cexrcv[iom*nxm+ix+ixmin] = locdat[ix];
				}
			    }
			  }

		}  /* end of depth loop */
	} /* end of iom loop */
#if defined(SGI)
#pragma critical
{
#endif
	if (imc < 2) {
		for (d = 0; d < ndepth; d++) {
			for (ix = ixmin, jx = 0; ix <= ixmax; ix++, jx++) {
				image[ix*ndepth+d+1] += scl*locima[d*nx+jx];
			}
		}
	}
	else if (imc == 2) {
		for (d = 0; d < ndepth; d++) {
			for (ix = ixmin, jx = 0; ix <= ixmax; ix++, jx++) {
				tmpim[d*nx+jx] += locima[d*nx+jx];
				tmpim2[d*nx+jx] += locima2[d*nx+jx];
			}
		}
	}
#if defined(SGI)
}
#endif

	free(opx);
	free(tmp1);
	free(tmp2);
	free(locdat);
	free(locsrc);
	free(locima);
	if (imc == 2) free(locima2);
#if defined(SGI)
} /* end of parallel region */
#endif

	if (imc == 2) {
		for (d = 0; d < ndepth; d++) {
			for (ix = ixmin, jx = 0; ix <= ixmax; ix++, jx++) {
				image[ix*ndepth+d+1] += scl*tmpim[d*nx+jx]/(tmpim2[d*nx+jx]+eps);
			}
		}
		free(tmpim);
		free(tmpim2);
	}
	t1 = wallclock_time();
	vmess("parallel region: %f seconds", t1-t0);

	if(exsrc) wx2xt(cexsrc, exsrc, optn, nxm, nxm, optn);
	if(exrcv) wx2xt(cexrcv, exrcv, optn, nxm, nxm, optn);

	free(cdata);
	free(csrc);
	free(taper);

	return;
}


