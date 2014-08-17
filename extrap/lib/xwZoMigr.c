#include "optim.h"
#include "genfft.h"
#include "par.h"

void readtable_opt(complex *oper, float k, int *hopl);

void xwZoMigr(float *data, int nx, int nt, float dt, float *velmod, int nxm, int nzm, int ixa, int ixb, float fmin, float fmax, float *xrcv, int izrcv, float ox, float dxm, int opl_max, int ntap, int conjg, int ndepth, float *image, int verbose, float *exrcv, int ndepthex)
{
	int     iomin, iomax, iom, niom, ix, jx, d, hopl, hopl2, i, j;
	int     index1, i1, i2, nfreq, optn, lenx, hoplen, sign;
	int     ixrcv, ixsrc, ixmin, ixmax, ixo, ixn;
	float   dom, om, c, cprev, df, sr, max_e, eps;
	float   *taper, scl, *pdata;
	float 	*trace;
	complex *ctrace;
	float   t0, t1;
	complex *opx, *cdata, *tmp1;
	complex wa, *ctmp, *locdat;
	complex *cexrcv=(complex *) exrcv;

/* define some constants */

	optn  = optncr(nt);
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
			if (verbose>1) fprintf(stderr,"xwZoMigr: ixrcv %f (%d) outside model\n", xrcv[ix], ixrcv);
			continue;
		}
        for (iom=0; iom<nfreq; iom++) {
			/* positioning of shot record into velocity model */
			cdata[iom*nxm+ixrcv].r = ctrace[iom].r;
			cdata[iom*nxm+ixrcv].i = ctrace[iom].i*scl;
        }
	}
    free(ctrace);
    free(trace);

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
	lenx  = 2*hopl2+nx;
	scl   = 2.0/nfreq;
	scl   = 1.0/(dt*dxm*dxm);

	if (verbose) {
		vmess("xwZoMigr: calculation aperture: %.2f (%d) <--> %.2f (%d) (%d positions)\n", ixmin*dxm+ox, ixmin, ixmax*dxm+ox, ixmax, nx);
	}

	/* allocate taper and local image arrays */

	if (ntap) {
		taper = (float *)malloc(ntap*sizeof(float));
		for (ix = 0; ix < ntap; ix++) {
			taper[ix] = exp(-1.0*(pow((2.0*(ntap-ix)/ntap), 2)));
		}
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

	tmp1  = (complex *)calloc(lenx, sizeof(complex));
	opx   = (complex *)calloc(hopl, sizeof(complex));
	locdat  = (complex *)calloc(lenx, sizeof(complex));

/* start extrapolation for all frequencies, depths and x-positions */

	for (iom = iomin; iom <= iomax; iom++) {
		for (j = 0, ix = ixmin; j < nx; j++, ix++) {
			locdat[hopl2+j] = cdata[iom*nxm+ix];
		}
		om    = iom*dom;
		d  = izrcv;
		cprev = 0;

		/* start extrapolation of receiver arrays */

		for (; d < ndepth; d++) {

			for (ix = ixmin, jx = 0; ix <= ixmax; ix++, jx++) {

				/* read operator from operator table (static array)*/

				c = velmod[d*nxm+ix];
				if (c != cprev) {
					readtable_opt(opx, om/c, &hoplen);
					cprev = c;
				}

				/* Extrapolation of data and source */

				wa.r = wa.i = 0.0;
				index1 = jx + hopl2;

    			for (j = 0; j < hoplen; j++) {
        			i1 = index1+j;
        			wa.r += locdat[i1].r*opx[j].r;
        			wa.r += locdat[i1].i*opx[j].i;
        			wa.i += locdat[i1].i*opx[j].r;
        			wa.i -= locdat[i1].r*opx[j].i;
			
        			i2 = index1-j;
        			wa.r += locdat[i2].r*opx[j].r;
        			wa.r += locdat[i2].i*opx[j].i;
        			wa.i += locdat[i2].i*opx[j].r;
        			wa.i -= locdat[i2].r*opx[j].i;
    			}

				tmp1[index1] = wa;
		
				/* imaging condition */
				image[ix*nzm+d+1] += wa.r;

			}

			for (j = 0; j < lenx; j++) {
				locdat[j] = tmp1[j];
			}

			/* taper edges of extrapolated data */

			for (j = 0; j < ntap; j++) {
				locdat[j].r *= taper[j];
				locdat[j].i *= taper[j];
				locdat[lenx-j-1].r *= taper[j];
				locdat[lenx-j-1].i *= taper[j];
			}

			if(d==ndepthex-1) {
				if ( cexrcv != NULL) {
					for (ix = 0; ix < nx; ix++) 
						cexrcv[iom*nxm+ix+ixmin] = locdat[ix];
				}
			}

		}  /* end of depth loop */
	} /* end of iom loop */


	free(opx);
	free(tmp1);
	free(locdat);

	if(exrcv) wx2xt(cexrcv, exrcv, optn, nxm, nxm, optn);

	free(cdata);
	if (ntap) free(taper);

	t1 = wallclock_time();
	vmess("xwZoMigr took: %f seconds", t1-t0);

	return;
}


