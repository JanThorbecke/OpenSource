#include "optim.h"
#include "genfft.h"
#include "par.h"

void readtable_opt(complex *oper, float k, int *hopl);

int extrapEdge(complex *data, complex *opx, complex *tmp, int ntap, int hoplen, int mode, int i0);

void xwExtr(float *data, int nx, int nt, float dt, float *velmod, 
	int id0, int id1, int ds, float fmin, float fmax, int *xi, int *zi, int nrec, 
	int opl_max, int ntap, int conjg, float *extr, int nxm, 
	float ox, float dxm, float *xrcv, int verbose)
{
	int    	iomin, iomax, iom, ix, d, hopl, hopl2, i, j;
	int	index1, i1, i2, nfreq, optn, lenx, ixrcv, hoplen;
	float  	dom, om, c, cprev, df, scl;
	float	*rdata, *taper, *pdata;
	complex	*opx, *cdata, *cextr, *tmp1, *tmp2, wa;

	optn  = optncr(nt);
	nfreq = optn/2 + 1;
	tmp1 = (complex *)malloc(nx*nfreq*sizeof(complex));

	/* pad zero's to get fourier length */
	if( nt != optn ) {
		if (verbose >1) vmess("xwExtr: padding zeros to data from %d to %d",nt,optn);
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

/* plane wave correction (rotation) at the surface 

      for (iom = 0; iom < nfreq; iom++) {
              for (ix = 0; ix < nx; ix++) {
                      xoff = NINT((xrcv[ix]-xrcv[nx/2])/dxm);
            shift = sin(45.0*M_PI/180)*xoff/2000.;
                      if (iom ==0 ) fprintf(stderr,"xoff = %f shift = %f \n",xoff, shift);
            om = dom*iom;
            tom = om*shift;
                      tmpt = cdata[iom*nxm+ixrcv];
            cdata[iom*nxm+ixrcv].r = tmpt.r*cos(-tom) - tmpt.i*sin(-tom);
            cdata[iom*nxm+ixrcv].i = tmpt.i*cos(-tom) + tmpt.r*sin(-tom);
              }
      }
*/

/* define some constants */

	nx    = nxm;
	df    = 1.0/(optn*dt);
	dom   = 2.*PI*df;
	iomin = (int)MIN((fmin*dt*optn), (nfreq-1));
	iomin = MAX(iomin, 1);
	iomax = MIN((int)(fmax*dt*optn), (nfreq-1));
	hopl  = (opl_max+1)/2;
	hopl2 = hopl-1;
    lenx  = 2*hopl2+nx+2*ntap;

	if(verbose) {
		fprintf(stderr,"    xwExtr: number of depth steps = %d\n", abs(id1-id0));
		fprintf(stderr,"    xwExtr: number of frequencies = %d\n", iomax-iomin+1);
	}

	if (ntap) {
		taper = (float *)malloc(ntap*sizeof(float));
		for (ix = 0; ix < ntap; ix++) {
			taper[ix] = exp(-1.0*(pow((2.0*(ntap-ix)/ntap), 2)));
		}
	}

	cextr = (complex *)calloc(nrec*nfreq, sizeof(complex));
	tmp1  = (complex *)calloc(lenx, sizeof(complex));
	tmp2  = (complex *)calloc(lenx, sizeof(complex));
	opx   = (complex *)calloc(hopl, sizeof(complex));

	for (iom = iomin; iom <= iomax; iom++) {
		om = iom*dom;
		cprev = 0;
		for (j = 0; j < nx; j++) {
			tmp1[hopl2+ntap+j] = cdata[iom*nx+j];
		}

/* experimental: taper edges of data */

		for (j = 0; j < ntap; j++) {
			tmp1[hopl2+ntap+j].r *= taper[j];
			tmp1[hopl2+ntap+j].i *= taper[j];
	
			tmp1[hopl2+ntap+ixrcv-j].r *= taper[j];
			tmp1[hopl2+ntap+ixrcv-j].i *= taper[j];
		}

		for (j=0; j<nrec; j++) 
			if (zi[j] == id0) cextr[iom*nrec+j] = tmp1[hopl2+ntap+xi[j]];

        for (d = id0; d != id1; d += ds) {

/* Extrapolation of data at the ntap edges of the model */

			c = velmod[d*nxm+0];
			readtable_opt(opx, om/c, &hoplen);
			if (ntap && c!=0.0) {
				extrapEdge(tmp1, opx, tmp2, ntap, hoplen, 1, hopl2);
			}

			for (ix = 0; ix < nx; ix++) {
				c = velmod[d*nxm+ix];
				if (c != cprev && c!=0.0) {
					readtable_opt(opx, om/c, &hoplen);
					cprev = c;
				}
				if (c==0.0) hoplen=0;

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
				if (hoplen != 0) tmp2[index1] = wa;
				else tmp2[index1] = tmp1[index1];
			}

/* Extrapolation of data at the ntap edges of the model */

			if (ntap) {
				extrapEdge(tmp1, opx, tmp2, ntap, hoplen, 1, hopl2+ntap+nx);
			}
			for (ix = 0; ix < lenx; ix++)  tmp1[ix] = tmp2[ix];

			for (j=0; j<nrec; j++) 
				if (zi[j] == (d+ds)) cextr[iom*nrec+j] = tmp2[hopl2+ntap+xi[j]];

			for (j = 0; j < ntap; j++) {
				tmp1[hopl2+j].r *= taper[j];
				tmp1[hopl2+j].i *= taper[j];

				tmp1[lenx-j-hopl2-1].r *= taper[j];
				tmp1[lenx-j-hopl2-1].i *= taper[j];
			}

		}
	}

	free(opx);
	free(tmp1);
	free(tmp2);

	rdata = (float *)malloc(nrec*optn*sizeof(float));
	wx2xt(cextr, rdata, optn, nrec, nrec, optn);

	for (ix = 0; ix < nrec; ix++) {
		for (j = 0; j < nt; j++) extr[ix*nt+j] = rdata[ix*optn+j];
	}

	free(cdata);
	free(cextr);
	if (ntap) free(taper);
	free(rdata);

	return;
}
