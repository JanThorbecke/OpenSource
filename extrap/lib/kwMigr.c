#include "optim.h"
#include "genfft.h"
#include "par.h"

void kwZoMigr(float *data, int nx, int nt, float dt, float *velmod, int nxm, int nzm, int ixa, int ixb, float fmin, float fmax, float *xrcv, int izrcv, float ox, float dxm, float dz, int ntap, int conjg, int ndepth, float *image, int verbose, float *exrcv, int ndepthex);

void kwMigr(float *data, int nx, int nt, float dt, float *velmod, int nxm, int nzm, int ixa, int ixb, float fmin, float fmax, float *wavelet, int ntw, int nxw, float *xareal, int izsrc, float *xrcv, int izrcv, float ox, float dxm, float dz, int ntap, int conjg, int conjgs, int ndepth, float eps_r, float eps_a, float *image, int imc, int verbose, float *exsrc, float *exrcv, int ndepthex, int zomigr)
{
	int     iomin, iomax, iom, ix, d, i, j;
	int     nfreq, optn, nkx, sign, endkx;
	int     ixrcv, ixsrc, ixmin, ixmax, ixo, ixn, ikx;
	float   k, k2, kz2, kx, kx2;
	float   dom, om, c, dkx, df, sr, max_e, eps;
	float   *taper, scl, scl2, *locima, *locima2, *tmpim, *tmpim2, *pdata;
	float 	*trace;
	complex *ctrace;
	float   t0, t1;
	complex *cdata, *csrc, tmp, ez;
	complex wa, *ctmp, da, *locdat, *locsrc;
	complex *cexsrc=(complex *) exsrc;
	complex *cexrcv=(complex *) exrcv;

	if (zomigr) {
		kwZoMigr(data, nx, nt, dt, velmod, nxm, nzm, ixa, ixb, fmin, fmax, 
			xrcv, izrcv, ox, dxm, dz, ntap, conjg, ndepth, image, 
			verbose, exrcv, ndepthex);
		return;
	}

/* define some constants */

	optn  = optncr(MAX(nt, ntw));
	nfreq = optn/2 + 1;
	df    = 1.0/(optn*dt);
	dom   = 2.0*M_PI*df;
	iomin = (int)MIN((fmin*dt*optn), (nfreq-1));
	iomin = MAX(iomin, 1);
	iomax = MIN((int)(fmax*dt*optn), (nfreq-1));

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
		if (ixrcv < 0 || ixrcv > nxm-1) {
			fprintf(stderr,"kwMigr: ixrcv %f (%d) outside model\n", xrcv[ix], ixrcv);
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

	eps = max_e*eps_r + eps_a;

	if (verbose >=2) {
		vmess("kwMigr: energy in shot = %.3e", max_e);
		vmess("kwMigr: eps_r = %.2e eps_a = %.2e eps = %.2e", eps_r, eps_a, eps);
	}

/* transform wavelet and place at position within src array*/ 

	if (conjgs) scl = -1.0;
	else scl = 1.0;

	ixo = nxm;
	ixn = 0;
	for (ix = 0; ix < nxw; ix++) {
		memcpy(trace,&wavelet[ix*ntw],ntw*sizeof(float));
        if (optn > ntw) 
            memset( &trace[ntw], 0, sizeof(float)*(optn-ntw) );
            
        rc1fft(trace,ctrace,optn,sign);

		ixsrc = NINT((xareal[ix]-ox)/dxm);
		if (ixsrc < ixo) ixo = MAX(0,ixsrc);
		if (ixsrc > ixn) ixn = MIN(nxm-1,ixsrc);
		if (ixsrc < 0 || ixsrc > nxm-1) {
			fprintf(stderr,"kwMigr: ixsrc %f (%d) outside model\n", xareal[ix], ixsrc);
			continue;
		}
        for (iom=0; iom<nfreq; iom++) {
			csrc[iom*nxm+ixsrc].r = ctrace[iom].r;
			csrc[iom*nxm+ixsrc].i = ctrace[iom].i*scl;
        }
	}
    free(ctrace);
    free(trace);

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
		vmess("kwMigr: calculation aperture: %.2f (%d) <--> %.2f (%d) (%d positions)", ixmin*dxm+ox, ixmin, ixmax*dxm+ox, ixmax, nx);
	}

	/* allocate taper and local image arrays */

	scl   = 2.0/nfreq;
    scl   = 1.0/(dt*dxm*dxm);
	nkx   = optncc(2*ntap+nxm);
	ntap  = (nkx-nxm)/2;
	scl2  = 1.0/nkx;
	dkx   = 2.0*M_PI/(nkx*dxm);

	taper = (float *)malloc(ntap*sizeof(float));
	for (ix = 0; ix < ntap; ix++) {
		taper[ix] = exp(-1.0*(pow((0.4*(ntap-ix)/ntap), 2)));
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

	locdat  = (complex *)calloc(nkx, sizeof(complex));
	locsrc  = (complex *)calloc(nkx, sizeof(complex));
	locima  = (float *)calloc(ndepth*nxm, sizeof(float));
	if (imc == 2) locima2 = (float *)calloc(ndepth*nxm, sizeof(float));

/* start extrapolation for all frequencies, depths and x-positions */

	for (iom = iomin; iom <= iomax; iom++) {
		for (ix = 0; ix < nkx; ix++) {
			locdat[ix].r = locdat[ix].i = 0.0;
			locsrc[ix].r = locsrc[ix].i = 0.0;
		}

		for (ix = ixmin; ix <= ixmax; ix++) {
			locdat[ntap+ix] = cdata[iom*nxm+ix];
			locsrc[ntap+ix] = csrc[iom*nxm+ix];
		}
		om = iom*dom;
		d  = MIN(izsrc,izrcv);

		/* if source at depth, first process only receivers */
		if(izsrc>izrcv) {
		if (verbose) {
			vmess("kwMigr: extrapolate receivers (%d) to src level (%d)\n", izrcv, izsrc);
		}

		cc1fft(locdat, nkx, 1);

		for (d = izrcv; d < izsrc; d++) {
			if(verbose)vmess("recv processing depth %d",d);
			c = 0.0;
			for (ix = ixmin; ix <= ixmax; ix++) c += velmod[d*nxm+ix];
			k = nx*om/c;
			k2 = k*k;

			/* kx = 0 */
			ez.r = cos(k*dz);
			ez.i = -sin(k*dz);

			tmp.r  = ez.r*locdat[0].r;
			tmp.r += ez.i*locdat[0].i;
			tmp.i  = ez.r*locdat[0].i;
			tmp.i -= ez.i*locdat[0].r;
			locdat[0] = tmp;

			/* kx != 0 */
			endkx = MIN((int)(k/dkx),nkx/2);
        	for (ikx = 1; ikx <= endkx; ikx++) {
				kx  = ikx*dkx;
				kx2 = kx*kx;
				kz2 = k2 - kx2;

				ez.r = cos(sqrt(kz2)*dz);
				ez.i = -sin(sqrt(kz2)*dz);

				tmp.r  = ez.r*locdat[ikx].r;
				tmp.r += ez.i*locdat[ikx].i;
				tmp.i  = ez.r*locdat[ikx].i;
				tmp.i -= ez.i*locdat[ikx].r;
				locdat[ikx] = tmp;

				tmp.r  = ez.r*locdat[nkx-ikx].r;
				tmp.r += ez.i*locdat[nkx-ikx].i;
				tmp.i  = ez.r*locdat[nkx-ikx].i;
				tmp.i -= ez.i*locdat[nkx-ikx].r;
				locdat[nkx-ikx] = tmp;
	       		}
				/* save extrapolated field at requested depth */
				if(d==ndepthex) {
					if ( cexrcv ) {
			  			for (ix = 0; ix < nxm; ix++) {
							cexrcv[iom*nxm+ix] = locdat[ntap+ix];
						}
			  		}
				}

			}
			cc1fft(locdat, nkx, -1);
		}

		/* if receivers deeper than source, first extrapolate source */
		else if(izsrc<izrcv) {
			vmess("kwMigr: extrapolate src (%d) to receiver level (%d)\n", izsrc, izrcv);

              cc1fft(locsrc, nkx, 1);


		for (d = izsrc; d < izrcv; d++) {
			c = 0.0;
			for (ix = ixmin; ix <= ixmax; ix++) c += velmod[d*nxm+ix];
			k = nx*om/c;
			k2 = k*k;

			/* kx = 0 */
			ez.r = cos(k*dz);
			ez.i = -sin(k*dz);

			tmp.r  = ez.r*locsrc[0].r;
			tmp.r -= ez.i*locsrc[0].i;
			tmp.i  = ez.r*locsrc[0].i;
			tmp.i += ez.i*locsrc[0].r;
			locsrc[0] = tmp;

			/* kx != 0 */
			endkx = MIN((int)(k/dkx),nkx/2);
        	for (ikx = 1; ikx <= endkx; ikx++) {
				kx  = ikx*dkx;
				kx2 = kx*kx;
				kz2 = k2 - kx2;

				ez.r = cos(sqrt(kz2)*dz);
				ez.i = -sin(sqrt(kz2)*dz);

				tmp.r  = ez.r*locsrc[ikx].r;
				tmp.r -= ez.i*locsrc[ikx].i;
				tmp.i  = ez.r*locsrc[ikx].i;
				tmp.i += ez.i*locsrc[ikx].r;
				locsrc[ikx] = tmp;

				tmp.r  = ez.r*locsrc[nkx-ikx].r;
				tmp.r -= ez.i*locsrc[nkx-ikx].i;
				tmp.i  = ez.r*locsrc[nkx-ikx].i;
				tmp.i += ez.i*locsrc[nkx-ikx].r;
				locsrc[nkx-ikx] = tmp;
	       		}
				/* save extrapolated field at requested depth */
				if(d==ndepthex) {
					if ( cexsrc ) {
						for (ix = 0; ix < nxm; ix++) {
							cexsrc[iom*nxm+ix] = locsrc[ntap+ix];
						}
					}
				}

			}
			cc1fft(locsrc, nkx, -1);

		}

		/* start extrapolation of both src and receiver arrays */

		for (; d < ndepth; d++) {

			/* transform to wavenumber domain */

			cc1fft(locdat, nkx, 1);
			cc1fft(locsrc, nkx, 1);

			/* Extrapolation of data and source */

			c = 0.0;
			for (ix = ixmin; ix <= ixmax; ix++) c += velmod[d*nxm+ix];
			k = nx*om/c;
			k2 = k*k;

			/* kx = 0 */
			ez.r = cos(k*dz);
			ez.i = -sin(k*dz);

			tmp.r  = ez.r*locdat[0].r;
			tmp.r += ez.i*locdat[0].i;
			tmp.i  = ez.r*locdat[0].i;
			tmp.i -= ez.i*locdat[0].r;
			locdat[0] = tmp;
	
			tmp.r  = ez.r*locsrc[0].r;
			tmp.r -= ez.i*locsrc[0].i;
			tmp.i  = ez.r*locsrc[0].i;
			tmp.i += ez.i*locsrc[0].r;
			locsrc[0] = tmp;
	
			/* kx != 0 */
			endkx = MIN((int)(k/dkx),nkx/2);
        	for (ikx = 1; ikx <= endkx; ikx++) {
                kx  = ikx*dkx;
               	kx2 = kx*kx;
               	kz2 = k2 - kx2;

				ez.r = cos(sqrt(kz2)*dz);
				ez.i = -sin(sqrt(kz2)*dz);

				tmp.r  = ez.r*locdat[ikx].r;
				tmp.r += ez.i*locdat[ikx].i;
				tmp.i  = ez.r*locdat[ikx].i;
				tmp.i -= ez.i*locdat[ikx].r;
				locdat[ikx] = tmp;

				tmp.r  = ez.r*locdat[nkx-ikx].r;
				tmp.r += ez.i*locdat[nkx-ikx].i;
				tmp.i  = ez.r*locdat[nkx-ikx].i;
				tmp.i -= ez.i*locdat[nkx-ikx].r;
				locdat[nkx-ikx] = tmp;

				tmp.r  = ez.r*locsrc[ikx].r;
				tmp.r -= ez.i*locsrc[ikx].i;
				tmp.i  = ez.r*locsrc[ikx].i;
				tmp.i += ez.i*locsrc[ikx].r;
				locsrc[ikx] = tmp;

				tmp.r  = ez.r*locsrc[nkx-ikx].r;
				tmp.r -= ez.i*locsrc[nkx-ikx].i;
				tmp.i  = ez.r*locsrc[nkx-ikx].i;
				tmp.i += ez.i*locsrc[nkx-ikx].r;
				locsrc[nkx-ikx] = tmp;
       		}

			/* transform data back to space domain */

			cc1fft(locdat, nkx, -1);
			cc1fft(locsrc, nkx, -1);

			for (j = 0; j < nkx; j++) {
				locdat[j].r *= scl2;
				locdat[j].i *= scl2;
				locsrc[j].r *= scl2;
				locsrc[j].i *= scl2;
			}

			/* imaging condition */

			for (ix = ixmin; ix <= ixmax; ix++) {
					wa = locdat[ntap+ix];
					da = locsrc[ntap+ix];
		
					sr = (da.r*wa.r+da.i*wa.i);
					if (imc == 0) {
						locima[d*nx+ix] += sr;
					}
					else if (imc == 1) {
						locima[d*nx+ix] += sr/(da.r*da.r+da.i*da.i+eps);
					}
					else if (imc == 2) {
						locima[d*nx+ix] += sr;
						locima2[d*nx+ix] += da.r*da.r+da.i*da.i;
					}
			}

			/* save extrapolated field at requested depth */

			if (d==ndepthex-1) {
				if ( cexsrc ) {
					for (ix = 0; ix < nxm; ix++) {
						cexsrc[iom*nxm+ix] = locsrc[ntap+ix];
					}
				}
				if ( cexrcv ) {
					for (ix = 0; ix < nxm; ix++) {
						cexrcv[iom*nxm+ix] = locdat[ntap+ix];
					}
				}
			}

			/* taper extrapolated data at edges */

			for (j = 0; j < ntap; j++) {
				locdat[j].r *= taper[j];
				locdat[j].i *= taper[j];
				locsrc[j].r *= taper[j];
				locsrc[j].i *= taper[j];
				locdat[nkx-j-1].r *= taper[j];
				locdat[nkx-j-1].i *= taper[j];
				locsrc[nkx-j-1].r *= taper[j];
				locsrc[nkx-j-1].i *= taper[j];
			}

		} /* end of depth loop */

	} /* end of iom loop */

	if (imc < 2) {
		for (d = 0; d < ndepth; d++) {
			for (ix = ixmin; ix <= ixmax; ix++) {
				image[ix*nzm+d+1] += scl*locima[d*nx+ix];
			}
		}
	}
	else if (imc == 2) {
		for (d = 0; d < ndepth; d++) {
			for (ix = ixmin; ix <= ixmax; ix++) {
				tmpim[d*nx+ix] += locima[d*nx+ix];
				tmpim2[d*nx+ix] += locima2[d*nx+ix];
			}
		}
	}

	free(locdat);
	free(locsrc);
	free(locima);

	if (imc == 2) {
		for (d = 0; d < ndepth; d++) {
			for (ix = ixmin; ix <= ixmax; ix++) {
				image[ix*nzm+d+1] += scl*tmpim[d*nx+ix]/(tmpim2[d*nx+ix]+eps);
			}
		}
		free(tmpim);
		free(tmpim2);
		free(locima2);
	}
	t1 = wallclock_time();
	vmess("kwMigr took: %f seconds", t1-t0);

	if(exsrc) wx2xt(cexsrc, exsrc, optn, nxm, nxm, optn);
	if(exrcv) wx2xt(cexrcv, exrcv, optn, nxm, nxm, optn);

	free(cdata);
	free(csrc);
	free(taper);

	return;
}
