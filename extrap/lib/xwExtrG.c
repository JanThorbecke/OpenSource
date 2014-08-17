#include "optim.h"
#include "genfft.h"
#include "par.h"

void xwExtrG(float *data, int nx, int nt, float dt, float *velmod, 
	int ndepth, float fmin, float fmax, int *xi, int *zi, int nrec, 
	int opl_max, int ntap, int conjg, float *extr, int nxm, 
	float ox, float dxm, float dzm, float *xrcv, int verbose)
{
	int    	iomin, iomax, iom, ix, d, hopl, hopl2, i, j;
	int	index1, i1, i2, nfreq, optn, lenx, ixrcv,hoplen;
	float  	dom, om, c, cprev, df, scl;
	float	*rdata, *taper, *pdata;
	float   zr, xr, r, dz, pi2, pi4, cosphi, k1r, k, vel, *sqa, phi;
	complex	*opx, *cdata, *cextr, *tmp1, *tmp2, wa, tmp;

	pi4  = PI*0.25;
	pi2  = 1.0/(2.0*PI);
	dz = dzm;

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

/* define some constants */

	nx    = nxm;
	df    = 1.0/(optn*dt);
	dom   = 2.*PI*df;
	iomin = (int)MIN((fmin*dt*optn), (nfreq-1));
	iomin = MAX(iomin, 1);
	iomax = MIN((int)(fmax*dt*optn), (nfreq-1));
	hopl  = (opl_max+1)/2;
	hopl2 = hopl-1;
	lenx  = 2*hopl2+nx;

	if(verbose) {
		fprintf(stderr,"    xwExtrG: number of depth steps = %d\n", ndepth);
		fprintf(stderr,"    xwExtrG: number of frequencies = %d\n", iomax-iomin+1);
	}

	cextr = (complex *)calloc(nrec*nfreq, sizeof(complex));
	taper = (float *)malloc(lenx*sizeof(float));
	for (ix = 0; ix < lenx; ix++) taper[ix] = 1.0;
	for (ix = 0; ix < ntap; ix++) {
		taper[ix] = exp(-1.0*(pow((0.4*(ntap-ix)/ntap), 2)));
		taper[lenx-1-ix] = taper[ix];
	}

	tmp1  = (complex *)calloc(lenx, sizeof(complex));
	tmp2  = (complex *)calloc(lenx, sizeof(complex));
	opx   = (complex *)calloc(hopl, sizeof(complex));

	for (ix = 0; ix < nx; ix++) {
		xr = ix*dxm-ox;
		zr = dz;
		r  = sqrt(xr*xr + zr*zr);
	}
	sqa   = (float *)malloc(nfreq*sizeof(float));
    for (iom = 0; iom < nfreq; iom++) {
        sqa[iom] = sqrt(iom*dom*pi2/2000.0);
    }

	for (iom = iomin; iom <= iomax; iom++) {
		om = iom*dom;
		cprev = 0;
		for (j = 0; j < nx; j++) {
			tmp1[hopl2+j] = cdata[iom*nx+j];
		}
		for (j=0; j<nrec; j++) 
			if (zi[j] == 0) cextr[iom*nrec+j] = tmp1[hopl2+xi[j]];

		for (d = 0; d < ndepth; d++) {

			for (ix = 0; ix < nx; ix++) {
				c = velmod[d*nxm+ix];

				wa.r = wa.i = 0.0;
				index1 = ix + hopl2;
				for (j = 0; j < nx; j++) {
/*
					xr = (ix-j)*dxm;
					r  = sqrt(xr*xr + dz*dz);
					cosphi = dz/r*sqrt(r);
					vel = 0.5*(c+velmod[(d+1)*nxm+j]);
					k = om/vel;
					k1r = k*r;
					tmp.r  = sqrt(om*pi2/vel)*cosphi*cos(k1r-pi4);
					tmp.i  = -sqrt(om*pi2/vel)*cosphi*sin(k1r-pi4);
*/

					vel = 0.5*(c+velmod[(d+1)*nxm+j]);
					k = om/vel;

					xr = (ix-j)*dxm;
				        r   = sqrt(xr*xr + dz*dz);
			            phi = acos(dz/r);
                        k1r = k*r;
                        
                        cosphi = cos(phi)/sqrt(r);
				        tmp.r  = sqa[iom]*cosphi*cos(k1r-pi4);
				        tmp.i  = -sqa[iom]*cosphi*sin(k1r-pi4);

					wa.r += tmp.r*tmp1[hopl2+j].r;
					wa.r -= tmp.i*tmp1[hopl2+j].i;
					wa.i += tmp.r*tmp1[hopl2+j].i;
					wa.i += tmp.i*tmp1[hopl2+j].r;
				}
				tmp2[index1].r = wa.r;
				tmp2[index1].i = wa.i;

//				fprintf(stderr,"wa = %e %e\n", wa.r, wa.i);
			}
			for (j=0; j<nrec; j++) 
				if (zi[j] == (d+1)) cextr[iom*nrec+j] = tmp2[hopl2+xi[j]];

			for (ix = 0; ix < lenx; ix++) {
				tmp1[ix].r = tmp2[ix].r*taper[ix];	
				tmp1[ix].i = tmp2[ix].i*taper[ix];	
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
	free(taper);
	free(rdata);

	return;
}
