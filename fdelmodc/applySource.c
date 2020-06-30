#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc.h"

void vmess(char *fmt, ...);

#define c1 (9.0/8.0)
#define c2 (-1.0/24.0)

/*********************************************************************
 * 
 * Add's the source amplitude(s) to the grid.
 * 
 * For the acoustic schemes, the source-type must not be txx tzz or txz.
 *
 *   AUTHOR:
 *           Jan Thorbecke (janth@xs4all.nl)
 *           The Netherlands 
 *
 **********************************************************************/

int applySource(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float **src_nwav, int verbose)
{
	int is0, ibndz, ibndx;
	int isrc, ix, iz, n1;
	int id1, id2;
	float src_ampl, time, scl, dt, sdx;
	float Mxx, Mzz, Mxz;
	static int first=1;

	if (src.type==6) {
    	ibndz = mod.ioXz;
    	ibndx = mod.ioXx;
	}
	else if (src.type==7) {
    	ibndz = mod.ioZz;
    	ibndx = mod.ioZx;
	}
	else if (src.type==2) {
    	ibndz = mod.ioTz;
    	ibndx = mod.ioTx;
    	if (bnd.lef==4 || bnd.lef==2) ibndx += bnd.ntap;
    	if (bnd.top==4 || bnd.top==2) ibndz += bnd.ntap;
	}
	else {	
    	ibndz = mod.ioPz;
    	ibndx = mod.ioPx;
    	if (bnd.lef==4 || bnd.lef==2) ibndx += bnd.ntap;
    	if (bnd.top==4 || bnd.top==2) ibndz += bnd.ntap;
	}

	n1   = mod.naz;
	dt   = mod.dt;
	sdx  = 1.0/mod.dx;

	/* special txz source activated? */

	if ((bnd.top==1) && (src.type==2)) {
		iz = izsrc + ibndz;
		if (iz==ibndz) {
            if (src.orient != 1) {
				if (first) {
					vmess("Only monopole Txz source allowed at surface. Reset to monopole");
					first = 0;
				}
				src.orient=1;
			}
		}
	}
             
/*
* for plane wave sources the sources are placed 
* around the central shot position 
* the first source position has an offset in x of is0
*
* itime = 0 corresponds with time=0
* itime = 1 corresponds with time=dt
* src[0] (the first sample) corresponds with time = 0
*/

	is0 = -1*floor((src.n-1)/2);
#pragma omp	for private (isrc, src_ampl, ix, iz, time, id1, id2, scl) 
	for (isrc=0; isrc<src.n; isrc++) {
		src_ampl=0.0;
		/* calculate the source position */
		if (src.random || src.multiwav) {
			ix = src.x[isrc] + ibndx;
			iz = src.z[isrc] + ibndz;
		}
        else if (src.plane) {/* plane wave sources */
            ix = ixsrc + ibndx + src.x[isrc];
            iz = izsrc + ibndz + src.z[isrc];
		}
		else { /* point sources */
            ix = ixsrc + ibndx + is0 + isrc;
            iz = izsrc + ibndz;
		}
		time = itime*dt - src.tbeg[isrc];
		id1 = floor(time/dt);
		id2 = id1+1;
        
		/* delay not reached or no samples left in source wavelet? */
		if ( (time < 0.0) || ( (itime*dt) >= src.tend[isrc]) ) continue;

//		fprintf(stderr,"isrc=%d ix=%d iz=%d src.x=%d src.z=%d\n", isrc, ix, iz, src.x[isrc], src.z[isrc]);

		if (!src.multiwav) { /* only one wavelet for all sources */
			src_ampl = src_nwav[0][id1]*(id2-time/dt) + src_nwav[0][id2]*(time/dt-id1);
		}
		else { /* multi-wavelet sources */
			src_ampl = src_nwav[isrc][id1]*(id2-time/dt) + src_nwav[isrc][id2]*(time/dt-id1);
		}

		if (src_ampl==0.0) continue;
		if ( ((ix-ibndx)<0) || ((ix-ibndx)>mod.nx) ) continue; /* source outside grid */

		if (verbose>=4 && itime==0) {
			vmess("Source %d positioned at grid ix=%d iz=%d",isrc, ix, iz);
		}

		/* cosine squared windowing to reduce edge effects on shot arrays */
		if ( (src.n>1) && src.window) {
            scl = 1.0;
			if (isrc < src.window) {
				scl = cos(0.5*M_PI*(src.window - isrc)/src.window);
			}
			else if (isrc > src.n-src.window+1) {
				scl = cos(0.5*M_PI*(src.window - (src.n-isrc+1))/src.window);
			}
			src_ampl *= scl*scl;
		}

		/* source scaling factor to compensate for discretisation */

		/* old amplitude setting does not obey reciprocity */
		// src_ampl *= rox[ix*n1+iz]*l2m[ix*n1+iz]/(dt);

/* in older version added factor 2.0 to be compliant with defined Green's functions in Marchenko algorithm */
/* this is now set to 1.0 */
		src_ampl *= (1.0/mod.dx)*l2m[ix*n1+iz];

		if (verbose>5) {
			vmess("Source %d at grid [ix=%d,iz=%d] at itime %d has value %e",isrc, ix,iz, itime, src_ampl);
		}

		/* Force source */

#pragma omp critical 
{
		if (src.type == 6) {
			vx[ix*n1+iz] += src_ampl*rox[ix*n1+iz]/(l2m[ix*n1+iz]);
			/* stable implementation from "Numerical Techniques for Conservation Laws with Source Terms" by Justin Hudson */
			//vx[ix*n1+iz] = 0.5*(vx[(ix+1)*n1+iz]+vx[(ix-1)*n1+iz])+src_ampl*rox[ix*n1+iz]/(l2m[ix*n1+iz]);
		}
		else if (src.type == 7) {
			vz[ix*n1+iz] += src_ampl*roz[ix*n1+iz]/(l2m[ix*n1+iz]);
			/* stable implementation from "Numerical Techniques for Conservation Laws with Source Terms" by Justin Hudson */
			/* stable implementation changes amplitude and more work is needed */
			//vz[ix*n1+iz] = 0.5*(vz[ix*n1+iz-1]+vz[ix*n1+iz+1])+src_ampl*roz[ix*n1+iz]/(l2m[ix*n1+iz]);
			//vz[ix*n1+iz] = 0.25*(vz[ix*n1+iz-2]+vz[ix*n1+iz-1]+vz[ix*n1+iz]+vz[ix*n1+iz+1])+src_ampl*roz[ix*n1+iz]/(l2m[ix*n1+iz]);
        } 
		else if (src.type == 10) { /* scale with 1/(ro*2dx) note that roz=dt/(ro*dx) */
		    tzz[ix*n1+iz-1] -= src_ampl*roz[ix*n1+iz]/(2.0*mod.dt);
		    tzz[ix*n1+iz+1] += src_ampl*roz[ix*n1+iz]/(2.0*mod.dt);
        } /* src.type */

        
		/* Stress source */

		if (mod.ischeme <= 2) { /* Acoustic scheme */
			/* Compressional source */
			if (src.type == 1) {
				if (src.orient != 1) src_ampl=src_ampl/mod.dx;

				if (src.orient==1) { /* monopole */
					tzz[ix*n1+iz] += src_ampl;
				}
				else if (src.orient==2) { /* dipole +/- */
					tzz[ix*n1+iz] += src_ampl;
					tzz[ix*n1+iz+1] -= src_ampl;
				}
				else if (src.orient==3) { /* dipole - + */
					tzz[ix*n1+iz] += src_ampl;
					tzz[(ix-1)*n1+iz] -= src_ampl;
				}
				else if (src.orient==4) { /* dipole +/0/- */
					if (iz > ibndz) 
						tzz[ix*n1+iz-1]+= 0.5*src_ampl;
					if (iz < mod.nz+ibndz-1) 
						tzz[ix*n1+iz+1] -= 0.5*src_ampl;
				}
				else if (src.orient==5) { /* dipole + - */
					tzz[ix*n1+iz] += src_ampl;
					tzz[(ix+1)*n1+iz] -= src_ampl;
				}
			}
		}
		else { /* Elastic scheme */
			/* Compressional source */
			if (src.type == 1) {
				if (src.orient != 1) src_ampl=src_ampl/mod.dx;
				if (src.orient==1) { /* monopole */
					txx[ix*n1+iz] += src_ampl;
					tzz[ix*n1+iz] += src_ampl;
				}
				else if (src.orient==2) { /* dipole +/- */
					txx[ix*n1+iz] += src_ampl;
					tzz[ix*n1+iz] += src_ampl;
					txx[ix*n1+iz+1] -= src_ampl;
					tzz[ix*n1+iz+1] -= src_ampl;
				}
				else if (src.orient==3) { /* dipole - + */
					txx[ix*n1+iz] += src_ampl;
					tzz[ix*n1+iz] += src_ampl;
					txx[(ix-1)*n1+iz] -= src_ampl;
					tzz[(ix-1)*n1+iz] -= src_ampl;
				}
				else if (src.orient==4) { /* dipole +/0/- */
					if (iz > ibndz) {
						txx[ix*n1+iz-1]+= 0.5*src_ampl;
						tzz[ix*n1+iz-1]+= 0.5*src_ampl;
					}
					if (iz < mod.nz+ibndz-1) {
						txx[ix*n1+iz+1] -= 0.5*src_ampl;
						tzz[ix*n1+iz+1] -= 0.5*src_ampl;
					}
				}
				else if (src.orient==5) { /* dipole + - */
					txx[ix*n1+iz] += src_ampl;
					tzz[ix*n1+iz] += src_ampl;
					txx[(ix+1)*n1+iz] -= src_ampl;
					tzz[(ix+1)*n1+iz] -= src_ampl;
				}
			}
			else if (src.type == 2) {
				/* Txz source */
				if ((iz == ibndz) && bnd.top==1) {
					txz[(ix-1)*n1+iz-1] += src_ampl;
					txz[ix*n1+iz-1] += src_ampl;
				}
				else {
					txz[ix*n1+iz] += src_ampl;
				}
				/* possible dipole orientations for a txz source */
				if (src.orient == 2) { /* dipole +/- */
					txz[ix*n1+iz+1] -= src_ampl;
				}
				else if (src.orient == 3) { /* dipole - + */
					txz[(ix-1)*n1+iz] -= src_ampl;
				}
				else if (src.orient == 4) { /*  dipole +/O/- */
					/* correction: subtrace previous value to prevent z-1 values. */
					txz[ix*n1+iz] -= 2.0*src_ampl;
					txz[ix*n1+iz+1] += src_ampl;
				}
				else if (src.orient == 5) { /* dipole + - */
					txz[(ix+1)*n1+iz] -= src_ampl;
				}
			}
			/* Tzz source */
			else if(src.type == 3) {
				tzz[ix*n1+iz] += src_ampl;
			} 
			/* Txx source */
			else if(src.type == 4) {
				txx[ix*n1+iz] += src_ampl;
			} 

/***********************************************************************
* pure potential shear S source (experimental)
* Curl S-pot = CURL(F) = dF_x/dz - dF_z/dx
***********************************************************************/
			else if(src.type == 5) {
				src_ampl = src_ampl*rox[ix*n1+iz]/(l2m[ix*n1+iz]);
				if (src.orient == 3) src_ampl = -src_ampl;
                /* first order derivatives */
				vx[ix*n1+iz]     += src_ampl*sdx;
				vx[ix*n1+iz-1]   -= src_ampl*sdx;
				vz[ix*n1+iz]     -= src_ampl*sdx;
				vz[(ix-1)*n1+iz] += src_ampl*sdx;
                
                /* second order derivatives */
                /*
				vx[ix*n1+iz]     += c1*src_ampl*sdx;
                vx[ix*n1+iz-1]   -= c1*src_ampl*sdx;
				vx[ix*n1+iz+1]   += c2*src_ampl*sdx;
                vx[ix*n1+iz-2]   -= c2*src_ampl*sdx;

                vz[ix*n1+iz]     -= c1*src_ampl*sdx;
				vz[(ix-1)*n1+iz] += c1*src_ampl*sdx;
				vz[(ix+1)*n1+iz] -= c2*src_ampl*sdx;
				vz[(ix-2)*n1+iz] += c2*src_ampl*sdx;
                 */

				/* determine second position of dipole */
				if (src.orient == 2) { /* dipole +/- vertical */
					iz += 1;
                    vx[ix*n1+iz]     -= src_ampl*sdx;
                    vx[ix*n1+iz-1]   += src_ampl*sdx;
                    vz[ix*n1+iz]     += src_ampl*sdx;
                    vz[(ix-1)*n1+iz] -= src_ampl*sdx;
				}
				else if (src.orient == 3) { /* dipole - + horizontal */
					ix += 1;
                    vx[ix*n1+iz]     -= src_ampl*sdx;
                    vx[ix*n1+iz-1]   += src_ampl*sdx;
                    vz[ix*n1+iz]     += src_ampl*sdx;
                    vz[(ix-1)*n1+iz] -= src_ampl*sdx;
				}
            }
/***********************************************************************
* pure potential pressure P source (experimental)
* Divergence P-pot = DIV(F) = dF_x/dx + dF_z/dz
***********************************************************************/
            else if(src.type == 8) {
			    src_ampl = src_ampl*rox[ix*n1+iz]/(l2m[ix*n1+iz]);
                if (src.orient == 3) src_ampl = -src_ampl;
                vx[(ix+1)*n1+iz] += src_ampl*sdx;
                vx[ix*n1+iz]     -= src_ampl*sdx;
                vz[ix*n1+iz+1]   += src_ampl*sdx;
                vz[ix*n1+iz]     -= src_ampl*sdx;
                /* determine second position of dipole */
                if (src.orient == 2) { /* dipole +/- */
                    iz += 1;
                    vx[(ix+1)*n1+iz] -= src_ampl*sdx;
                    vx[ix*n1+iz]     += src_ampl*sdx;
                    vz[ix*n1+iz+1]   -= src_ampl*sdx;
                    vz[ix*n1+iz]     += src_ampl*sdx;
                }
                else if (src.orient == 3) { /* dipole - + */
                    ix += 1;
                    vx[(ix+1)*n1+iz] -= src_ampl*sdx;
                    vx[ix*n1+iz]     += src_ampl*sdx;
                    vz[ix*n1+iz+1]   -= src_ampl*sdx;
                    vz[ix*n1+iz]     += src_ampl*sdx;
                }
			}
            else if(src.type == 9 || src.type == 11) {
				txx[ix*n1+iz] -= src.Mxx*src_ampl;
				tzz[ix*n1+iz] -= src.Mzz*src_ampl;
				txz[ix*n1+iz] -= src.Mxz*src_ampl;
			} /* src.type */
		} /* ischeme */
}
	} /* loop over isrc */

	return 0;
}
