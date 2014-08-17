#include "optim.h"
#include "genfft.h"
#include "par.h"

void readtable_opt(complex *oper, float k, int *hopl);

int extrapEdge(complex *data, complex *opx, complex *tmp, int ntap, int hoplen, int mode, int i0);

void xwSnap(float *data, int nx, int nt, float dt, float *velmod, 
	int ndepth, float fmin, float fmax, int opl, int ntap, int conjg, 
	float *snaps, int nxm, float ox, float dxm, float *xrcv, int verbose)
{
	int    	iomin, iomax, iom, ix, d, hopl, hopl2, i, j, ixrcv;
	int	index1, i1, i2, nfreq, optn, lenx, *it, nsnap, hoplen;
	int     reverse, id, id0, id1, dstep, rnd;
	float  	dom, om, c, cprev, df, tom, shift;
	float	*taper, *rdata, *pdata;
	float	tsnap1, tsnap2, dtsnap, t, scl;
	complex	*opx, *cdata, *tmp1, *tmp2, wa;

	if(!getparint("reverse", &reverse)) reverse = 0;
    if(!getparint("rnd", &rnd)) rnd = 0;

	optn  = optncr(nt);
	nfreq = (optn+2)/2;
	tmp1 = (complex *)calloc(nx*nfreq,sizeof(complex));

	/* pad zero's to get fourier length */
	if( nt != optn ) {
		if (verbose >1) vmess("xwSnap: padding zeros to data from %d to %d",nt,optn);
		pdata = (float *) malloc(optn*nx*sizeof(float));
		for (i=0; i<nx; i++) {
			for (j=0; j<nt; j++) pdata[i*optn+j] = data[i*nt+j];
			for (; j<optn; j++) pdata[i*optn+j] = 0.0;
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
	df    = 1.0/(optn*dt);
	dom   = 2.*PI*df;
	iomin = (int)MIN((fmin*dt*optn), (nfreq-1));
	iomin = MAX(iomin, 1);
	iomax = MIN((int)(fmax*dt*optn), (nfreq-1));
	hopl  = (opl+1)/2;
	hopl2 = hopl-1;
	lenx  = 2*hopl2+nxm+2*ntap;

	cdata = (complex *)calloc(nxm*nfreq, sizeof(complex));
    if (rnd) {
        srand(10);
        for (ix = 0; ix < nx; ix++) {
            ixrcv = NINT((xrcv[ix]-ox)/dxm);
            shift = ( (float)(1.0*rand()/((float)RAND_MAX)));
            fprintf(stderr,"random %d = %f %d\n",ix, shift, ixrcv);
            for (iom = 0; iom < nfreq; iom++) {
                om = dom*iom;
                tom = om*shift;
                cdata[iom*nxm+ixrcv].r += tmp1[iom*nx+ix].r*cos(-tom) - tmp1[iom*nx+ix].i*sin(-tom);
                cdata[iom*nxm+ixrcv].i += (tmp1[iom*nx+ix].i*cos(-tom) + tmp1[iom*nx+ix].r*sin(-tom))*scl;
            }
        }
    }
    else {
        for (iom = 0; iom < nfreq; iom++) {
            for (ix = 0; ix < nx; ix++) {
                ixrcv = NINT((xrcv[ix]-ox)/dxm);
                cdata[iom*nxm+ixrcv].r += tmp1[iom*nx+ix].r;
                cdata[iom*nxm+ixrcv].i += tmp1[iom*nx+ix].i*scl;
            }
        }
    }
	free(tmp1);

	if(verbose) {
		fprintf(stderr,"    xwSnap: number of depth steps = %d\n", ndepth);
		fprintf(stderr,"    xwSnap: number of frequencies = %d\n", iomax-iomin+1);
	}

	tmp1  = (complex *)calloc(lenx, sizeof(complex));
	tmp2  = (complex *)calloc(lenx, sizeof(complex));
	opx   = (complex *)calloc(hopl, sizeof(complex));
	rdata = (float *)malloc(nxm*optn*sizeof(float));

	if (ntap) {
		taper = (float *)malloc(ntap*sizeof(float));
		for (ix = 0; ix < ntap; ix++) {
			taper[ix] = exp(-1.0*(pow((2.0*(ntap-ix)/ntap), 2)));
		}
	}

	if(!getparfloat("dtsnap", &dtsnap)) dtsnap = 25.0*dt;
	if(!getparfloat("tsnap1", &tsnap1)) 
		tsnap1 = -floor(0.5*nt*dt/dtsnap)*dtsnap;
	if(!getparfloat("tsnap2", &tsnap2)) 
		tsnap2 = floor(0.5*nt*dt/dtsnap)*dtsnap;
	nsnap = 1+NINT((tsnap2-tsnap1)/dtsnap);
	if (verbose) 
		fprintf(stderr,"    xwSnap: number of snapshots = %d\n", nsnap);

	it = (int *)malloc(nsnap*sizeof(int));
	for (j = 0 ; j < nsnap ; j++) {
		t = (tsnap1+j*dtsnap);
		if (t >= 0) it[j] = NINT(t/dt);
		else it[j] = optn+NINT(t/dt);
		if (verbose >= 2) 
			fprintf(stderr,"    xwSnap: snapshot[%d] at time = %f[%d]\n",j,t,it[j]);
	}

	for (iom = 0; iom < iomin; iom++) {
		for (j = 0; j < nxm; j++) {
			cdata[iom*nxm+j].r = 0.0;
			cdata[iom*nxm+j].i = 0.0;
		}
	}
	for (iom = iomax+1; iom < nfreq; iom++) {
		for (j = 0; j < nxm; j++) {
			cdata[iom*nxm+j].r = 0.0;
			cdata[iom*nxm+j].i = 0.0;
		}
	}

	if (reverse) {
		id0=-1*(ndepth-1);
		id1=0;
		dstep=-1;
	}
	else {
		id0=0;
		id1=ndepth-1;
		dstep=1;
	}

	for (id = id0; id <= id1; id++) {
		d = id*dstep;

		for (j = 0; j < nxm*optn; j++) rdata[j] = 0.0;
		wx2xt(cdata, rdata, optn, nxm, nxm, optn);

		for (ix = 0; ix < nxm; ix++) {
			for (j = 0 ; j < nsnap ; j++) {
				snaps[ix*ndepth+d] += rdata[ix*optn+it[j]];
			}
		}

		for (iom = iomin; iom <= iomax; iom++) {
			om = iom*dom;
			cprev = 0;
			for (j = 0; j < nxm; j++) {
				tmp1[hopl2+ntap+j] = cdata[iom*nxm+j];
			}

			c = velmod[d*nxm+0];
			readtable_opt(opx, om/c, &hoplen);
			if (ntap && c!=0.0) {
				extrapEdge(tmp1, opx, tmp2, ntap, hoplen, 1, hopl2);
			} 

			for (ix = 0; ix < nxm; ix++) {
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

			if (ntap) {
				extrapEdge(tmp1, opx, tmp2, ntap, hoplen, 1, hopl2+ntap+nxm);
			}

            for (j = 0; j < ntap; j++) {
                tmp1[hopl2+j].r *= taper[j];
                tmp1[hopl2+j].i *= taper[j];

                tmp1[lenx-j-hopl2-1].r *= taper[j];
                tmp1[lenx-j-hopl2-1].i *= taper[j];
            }

			for (ix = 0; ix < nxm; ix++) {
				cdata[iom*nxm+ix] = tmp2[hopl2+ntap+ix];
			}

		}
	}

	free(opx);
	free(tmp1);
	free(tmp2);
	free(cdata);
	free(rdata);
	if (ntap) free(taper);
	free(it);

	return;
}
