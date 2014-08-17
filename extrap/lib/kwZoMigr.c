#include "optim.h"
#include "genfft.h"
#include "par.h"

void kwZoMigr(float *data, int nx, int nt, float dt, float *velmod, int nxm, int nzm, int ixa, int ixb, float fmin, float fmax, float *xrcv, int izrcv, float ox, float dxm, float dz, int ntap, int conjg, int ndepth, float *image, int verbose, float *exrcv, int ndepthex)
{
	int     iomin, iomax, iom, ix, d, i, j;
	int     nfreq, optn, nkx, sign, endkx;
	int     ixrcv, ixmin, ixmax, ixo, ixn, ikx;
	float   k, k2, kz2, kx, kx2;
	float   dom, om, c, dkx, df, sr;
	float   *taper, scl, scl2, *pdata;
	float 	*trace;
	complex *ctrace;
	float   t0, t1;
	complex *cdata, tmp, ez;
	complex wa, *ctmp, da, *locdat;
	complex *cexrcv=(complex *) exrcv;

/* define some constants */

	optn  = optncr(nt);
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

	if (conjg) scl = -1.0;
	else scl = 1.0;

    sign  = -1;
	for (ix = 0; ix < nx; ix++) {

		memcpy(trace,&data[ix*nt],nt*sizeof(float));
        if (optn > nt) 
			memset( &trace[nt], 0, sizeof(float)*(optn-nt) );
            
        rc1fft(trace,ctrace,optn,sign);

   		ixrcv = NINT((xrcv[ix]-ox)/dxm);
		if (ixrcv < 0 || ixrcv > nxm-1) {
			fprintf(stderr,"kwZoMigr: ixrcv %f (%d) outside model\n", xrcv[ix], ixrcv);
			continue;
		}
        for (iom=0; iom<nfreq; iom++) {
			/* positioning of shot record into velocity model */
			cdata[iom*nxm+ixrcv].r = ctrace[iom].r;
			cdata[iom*nxm+ixrcv].i = ctrace[iom].i*scl;
        }
	}

/* determine aperture to be calculated */

    ixo = nxm;
    ixn = 0;
	for (ix = 0; ix < nx; ix++) {
		ixrcv = NINT((xrcv[ix]-ox)/dxm);
		if (ixrcv < ixo) ixo = ixrcv;
		if (ixrcv > ixn) ixn = ixrcv;
	}
	ixmin = MAX(0, ixo-ixb-1);
	ixmax = MIN(ixn+ixa+1, nxm-1);
	nx    = (ixmax-ixmin)+1;

	if (verbose>=2) {
		vmess("kwZoMigr: calculation aperture: %.2f (%d) <--> %.2f (%d) (%d positions)", ixmin*dxm+ox, ixmin, ixmax*dxm+ox, ixmax, nx);
	}

/* define some constants */

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

	/* calculate image at depth = 0 */

	if(izrcv==0 ) {
		for (ix = ixmin; ix <= ixmax; ix++) {
			for (iom = iomin; iom <= iomax; iom++) {
				image[ix*nzm+0] += scl*cdata[iom*nxm+ix].r;
			}
		}
	}	

	t0 = wallclock_time();

	locdat  = (complex *)malloc(nkx*sizeof(complex));

	/* start extrapolation for all frequencies, depths and x-positions */

	for (iom = iomin; iom <= iomax; iom++) {
		memset(locdat,0,nkx*sizeof(complex));

		for (ix = ixmin; ix <= ixmax; ix++) {
			locdat[ntap+ix] = cdata[iom*nxm+ix];
		}
		om = iom*dom;
		d  = izrcv;

		/* start extrapolation of receiver arrays */

		for (; d < ndepth; d++) {

			/* transform to wavenumber domain */

			cc1fft(locdat, nkx, 1);

			/* Extrapolation of data */

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

			/* transform data back to space domain */

			cc1fft(locdat, nkx, -1);

			for (j = 0; j < nkx; j++) {
				locdat[j].r *= scl2;
				locdat[j].i *= scl2;
			}

			/* imaging condition */

			for (ix = ixmin; ix <= ixmax; ix++) {
				image[ix*nzm+d+1]+= scl*locdat[ntap+ix].r;
			}

			/* save extrapolated field at requested depth */

			if	(d==ndepthex-1) {
				if ( cexrcv != NULL ) {
					for (ix = 0; ix < nxm; ix++) {
						cexrcv[iom*nxm+ix] = locdat[ntap+ix];
					}
				}
			}

			/* taper extrapolated data at edges */

			for (j = 0; j < ntap; j++) {
				locdat[j].r *= taper[j];
				locdat[j].i *= taper[j];
				locdat[nkx-j-1].r *= taper[j];
				locdat[nkx-j-1].i *= taper[j];
			}

		} /* end of depth loop */

	} /* end of iom loop */


	if(exrcv) wx2xt(cexrcv, exrcv, optn, nxm, nxm, optn);

	free(locdat);
	free(cdata);
	free(taper);

	t1 = wallclock_time();
	vmess("kwZoMigr took: %f seconds", t1-t0);

	return;
}
