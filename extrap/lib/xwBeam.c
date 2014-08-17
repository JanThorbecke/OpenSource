#include "optim.h"
#include "genfft.h"
#include "par.h"

void readtable_opt(complex *oper, float k, int *hopl);

void xwBeam(float *data, int nx, int nt, float dt, float *velmod, int ndepth, float fmin, float fmax, int opl, int ntap, int conjg, float *beams, int nxm, float ox, float dxm, float *xrcv, int verbose)
{
	int    	iomin, iomax, iom, ix, d, hopl, hopl2, i, j, ixrcv;
	int	index1, nfreq, optn, lenx, i1, i2, hoplen;
	float  	dom, om, c, cprev, df, scale, scl;
	float	*taper, *locbeam, *pdata;
	complex	*opx, *cdata, *tmp1, *tmp2, wa;

	optn  = optncr(nt);
	nfreq = optn/2 + 1;
	tmp1 = (complex *)calloc(nx*nfreq,sizeof(complex)); 

	/* pad zero's to get fourier length */
	if( nt != optn ) {
		if (verbose >1) vmess("xwBeam: padding zeros to data from %d to %d",nt,optn);
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

	if (conjg) scl = -1.0;
	else scl = 1.0;
	df    = 1.0/(optn*dt);
	dom   = 2.*PI*df;
	iomin = (int)MIN((fmin*dt*optn), (nfreq-1));
	iomin = MAX(iomin, 1);
	iomax = MIN((int)(fmax*dt*optn), (nfreq-1));
	hopl  = (opl+1)/2;
	hopl2 = hopl-1;
	lenx  = 2*hopl2+nxm;
	scale = 1.0/(float)(iomax-iomin+1);

/* plane wave correction (rotation) at the surface 

    for (iom = 0; iom < nfreq; iom++) {
        for (ix = 0; ix < nx; ix++) {
            xoff = -1500.0+ix*15.0;
            shift = sin(15.0*M_PI/180)*xoff/2000.;
            if (iom ==0 ) fprintf(stderr,"xoff = %f shift = %f \n",xoff, shift);
            om = iom*dom*shift;
            wa = cdata[iom*nx+ix];
            cdata[iom*nx+ix].r = wa.r*cos(-om) - wa.i*sin(-om);
            cdata[iom*nx+ix].i = wa.i*cos(-om) + wa.r*sin(-om);
        }
    }
*/

	taper = (float *)malloc(nxm*sizeof(float));
	for (ix = 0; ix < nxm; ix++) taper[ix] = 1.0;
	for (ix = 0; ix < ntap; ix++) {
		taper[ix] = exp(-1.0*(pow((0.4*(ntap-ix)/ntap), 2)));
		taper[nxm-1-ix] = taper[ix];
	}

	cdata = (complex *)calloc(nxm*nfreq,sizeof(complex));
	for (iom = iomin; iom <= iomax; iom++) {
		for (ix = 0; ix < nx; ix++) {
			ixrcv = NINT((xrcv[ix]-ox)/dxm);
			cdata[iom*nxm+ixrcv].r += tmp1[iom*nx+ix].r;
			cdata[iom*nxm+ixrcv].i += tmp1[iom*nx+ix].i*scl;
		}
	}
	free(tmp1);

	if(verbose) {
		fprintf(stderr,"    xwBeam: number of depth steps = %d\n", ndepth);
		fprintf(stderr,"    xwBeam: number of frequencies = %d\n", iomax-iomin+1);
	}


	tmp1  = (complex *)calloc(lenx, sizeof(complex));
	tmp2  = (complex *)calloc(lenx, sizeof(complex));
	opx   = (complex *)calloc(hopl, sizeof(complex));
	locbeam = (float *)calloc(nxm*ndepth, sizeof(float));

	for (iom = iomin; iom <= iomax; iom++) {
		om = iom*dom;
		cprev = 0;
		for (j = 0; j < nxm; j++) {
			tmp1[hopl2+j].r = cdata[iom*nxm+j].r*taper[j];
			tmp1[hopl2+j].i = cdata[iom*nxm+j].i*taper[j];
		}

		for (d = 0; d < ndepth; d++) {

			for (ix = 0; ix < nxm; ix++) {
				locbeam[d*nxm+ix] += sqrt(tmp1[hopl2+ix].r*tmp1[hopl2+ix].r+
						tmp1[hopl2+ix].i*tmp1[hopl2+ix].i)*scale;
			}

			for (ix = 0; ix < nxm; ix++) {
				c = velmod[d*nxm+ix];
				if (c != cprev && c!=0.0) {
					readtable_opt(opx, om/c, &hoplen);
					cprev = c;
				}
				if (c==0.0) hoplen=0;

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
				if (hoplen != 0) tmp2[index1] = wa;
				else tmp2[index1] = tmp1[index1];
			}
			for (j = 0; j < lenx; j++) tmp1[j] = tmp2[j];
		}
	}

	for (d = 0; d < ndepth; d++) {
		for (ix = 0; ix < nxm; ix++) {
			beams[ix*ndepth+d] += locbeam[d*nxm+ix];
		}
	}

	free(opx);
	free(tmp1);
	free(tmp2);
	free(locbeam);
	free(cdata);
	free(taper);

	return;
}
