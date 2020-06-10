#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc3D.h"

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
long applySource3D(modPar mod, srcPar src, wavPar wav, bndPar bnd, long itime, long ixsrc, long iysrc, long izsrc,
	float *vx, float *vy, float *vz, float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz,
	float ***rox, float ***roy, float ***roz, float ***l2m, float **src_nwav, long verbose)
{
	long is0, ibndz, ibndy, ibndx;
	long isrc, ix, iy, iz, n1, n2, ix0, iy0, ixe, iye;
	long id1, id2, id3;
	float src_ampl, time, scl, dt, sdx;
	float Mxx, Myy, Mzz, Mxz, Myz, Mxy;
	static long first=1;

	if (src.type==6) {
    	ibndz = mod.ioXz;
    	ibndy = mod.ioXy;
    	ibndx = mod.ioXx;
	}
	else if (src.type==7) {
    	ibndz = mod.ioZz;
    	ibndy = mod.ioZy;
    	ibndx = mod.ioZx;
	}
	else if (src.type==2) {
    	ibndz = mod.ioTz;
    	ibndy = mod.ioTy;
    	ibndx = mod.ioTx;
    	if (bnd.lef==4 || bnd.lef==2) ibndx += bnd.ntap;
    	if (bnd.bac==4 || bnd.fro==2) ibndy += bnd.ntap;
    	if (bnd.top==4 || bnd.top==2) ibndz += bnd.ntap;
	}
	else {	
    	ibndz = mod.ioPz;
    	ibndy = mod.ioPy;
    	ibndx = mod.ioPx;
    	if (bnd.lef==4 || bnd.lef==2) ibndx += bnd.ntap;
    	if (bnd.bac==4 || bnd.fro==2) ibndy += bnd.ntap;
    	if (bnd.top==4 || bnd.top==2) ibndz += bnd.ntap;
	}

	n1   = mod.naz;
    n2   = mod.nax;
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
#pragma omp	for private (isrc, src_ampl, ix, iy, iz, time, id1, id2, id3, scl) 
	for (isrc=0; isrc<src.n; isrc++) {
		src_ampl=0.0;
		/* calculate the source position */
		if (src.random || src.multiwav) {
			ix = src.x[isrc] + ibndx;
			iy = src.y[isrc] + ibndy;
			iz = src.z[isrc] + ibndz;
		}
		else if (src.plane) {/* plane wave sources */
            ix = ixsrc + ibndx + src.x[isrc];
            iy = iysrc + ibndy + src.y[isrc];
            iz = izsrc + ibndz + src.z[isrc];
			ix0 = ixsrc + ibndx + src.x[0];
			ixe = ixsrc + ibndx + src.x[src.n-1];
			iy0 = iysrc + ibndy + src.y[0];
			iye = iysrc + ibndy + src.y[src.n-1];
			// vmess("ix0 %li ixe %li iy0 %li iye %li",ix0,ixe,iy0,iye);
		}
		else { /* point sources */
            ix = ixsrc + ibndx + isrc;
            iy = iysrc + ibndy + isrc;
            iz = izsrc + ibndz;
		}
		time = itime*dt - src.tbeg[isrc];
		id1 = floor(time/dt);
		id2 = id1+1;


		/* delay not reached or no samples left in source wavelet? */
		if ( (time < 0.0) || ( (itime*dt) >= src.tend[isrc]) ) continue;

		// fprintf(stderr,"isrc=%li ix=%li iy=%li iz=%li src.x=%li src.y=%li src.z=%li\n", isrc, ix, iy, iz, src.x[isrc], src.y[isrc], src.z[isrc]);

		if (!src.multiwav) { /* only one wavelet for all sources */
			src_ampl = src_nwav[0][id1]*(id2-time/dt) + src_nwav[0][id2]*(time/dt-id1);
		}
		else { /* multi-wavelet sources */
			src_ampl = src_nwav[isrc][id1]*(id2-time/dt) + src_nwav[isrc][id2]*(time/dt-id1);
		}
		if (src_ampl==0.0) continue;
		if ( ((ix-ibndx)<0) || ((ix-ibndx)>mod.nx) ) continue; /* source outside grid */
        if ( ((iy-ibndy)<0) || ((iy-ibndy)>mod.ny) ) continue; /* source outside grid */

		if (verbose>=4 && itime==0) {
			vmess("Source %li positioned at grid ix=%li iy=%li iz=%li",isrc, ix, iy, iz);
		}

		/* cosine squared windowing to reduce edge effects on shot arrays */
		if (src.plane) {
			if (src.nxwindow > 0){
				scl = 1.0;
				if (ix-ix0 < src.nxwindow) scl = cos(0.5*M_PI*(src.nxwindow - (ix-ix0))/src.nxwindow);
				else if (ixe-ix < src.nxwindow) scl = cos(0.5*M_PI*(src.nxwindow - (ixe-ix))/src.nxwindow);
				src_ampl *= scl*scl;
			}
			if (src.nywindow > 0){
				scl = 1.0;
				if (iy-iy0 < src.nywindow) scl = cos(0.5*M_PI*(src.nywindow - (iy-iy0))/src.nywindow);
				else if (iye-iy < src.nywindow) scl = cos(0.5*M_PI*(src.nywindow - (iye-iy))/src.nywindow);
				src_ampl *= scl*scl;
			}
		}

		/* source scaling factor to compensate for discretisation */

		/* old amplitude setting does not obey reciprocity */
		// src_ampl *= rox[ix*n1+iz]*l2m[ix*n1+iz]/(dt);

/* in older version added factor 2.0 to be compliant with defined Green's functions in Marchenko algorithm */
/* this is now set to 1.0 */
		src_ampl *= (1.0/(mod.dx*mod.dx))*l2m[iy][ix][iz];

		if (verbose>5) {
			vmess("Source %li at grid [ix=%li,iy=%li,iz=%li] at itime %li has value %e",isrc, ix, iy, iz, itime, src_ampl);
		}

		/* Force source */

		if (src.type == 6) {
			vx[iy*n1*n2+ix*n1+iz] += src_ampl*rox[iy][ix][iz]/(l2m[iy][ix][iz]);
			/* stable implementation from "Numerical Techniques for Conservation Laws with Source Terms" by Justin Hudson */
			//vx[ix*n1+iz] = 0.5*(vx[(ix+1)*n1+iz]+vx[(ix-1)*n1+iz])+src_ampl*rox[ix*n1+iz]/(l2m[ix*n1+iz]);
		}
		else if (src.type == 7) {
			vz[iy*n1*n2+ix*n1+iz] += src_ampl*roz[iy][ix][iz]/(l2m[iy][ix][iz]);
			/* stable implementation from "Numerical Techniques for Conservation Laws with Source Terms" by Justin Hudson */
			/* stable implementation changes amplitude and more work is needed */
			//vz[ix*n1+iz] = 0.5*(vz[ix*n1+iz-1]+vz[ix*n1+iz+1])+src_ampl*roz[ix*n1+iz]/(l2m[ix*n1+iz]);
			//vz[ix*n1+iz] = 0.25*(vz[ix*n1+iz-2]+vz[ix*n1+iz-1]+vz[ix*n1+iz]+vz[ix*n1+iz+1])+src_ampl*roz[ix*n1+iz]/(l2m[ix*n1+iz]);
        } /* src.type */

        
		/* Stress source */

		if (mod.ischeme <= 2) { /* Acoustic scheme */
			/* Compressional source */
			if (src.type == 1) {
				if (src.orient != 1) src_ampl=src_ampl/mod.dx;

				if (src.orient==1) { /* monopole */
					tzz[iy*n1*n2+ix*n1+iz] += src_ampl;
				}
				else if (src.orient==2) { /* dipole +/- */
					tzz[iy*n1*n2+ix*n1+iz] += src_ampl;
					tzz[iy*n1*n2+ix*n1+iz+1] -= src_ampl;
				}
				else if (src.orient==3) { /* dipole - + */
					tzz[iy*n1*n2+ix*n1+iz] += src_ampl;
					tzz[iy*n1*n2+(ix-1)*n1+iz] -= src_ampl;
				}
				else if (src.orient==4) { /* dipole +/0/- */
					if (iz > ibndz) 
						tzz[iy*n1*n2+ix*n1+iz-1]+= 0.5*src_ampl;
					if (iz < mod.nz+ibndz-1) 
						tzz[iy*n1*n2+ix*n1+iz+1] -= 0.5*src_ampl;
				}
				else if (src.orient==5) { /* dipole + - */
					tzz[iy*n1*n2+ix*n1+iz] += src_ampl;
					tzz[iy*n1*n2+(ix+1)*n1+iz] -= src_ampl;
				}
			}
		}
		else { /* Elastic scheme */
			/* Compressional source */
			if (src.type == 1) {
				if (src.orient==1) { /* monopole */
					txx[iy*n1*n2+ix*n1+iz] += src_ampl;
					tyy[iy*n1*n2+ix*n1+iz] += src_ampl;
					tzz[iy*n1*n2+ix*n1+iz] += src_ampl;
				}
				else if (src.orient==2) { /* dipole +/- */
					txx[iy*n1*n2+ix*n1+iz] += src_ampl;
					tyy[iy*n1*n2+ix*n1+iz] += src_ampl;
					tzz[iy*n1*n2+ix*n1+iz] += src_ampl;
					txx[iy*n1*n2+ix*n1+iz+1] -= src_ampl;
					tyy[iy*n1*n2+ix*n1+iz+1] -= src_ampl;
					tzz[iy*n1*n2+ix*n1+iz+1] -= src_ampl;
				}
				else if (src.orient==3) { /* dipole - + */
					txx[iy*n1*n2+ix*n1+iz] += src_ampl;
					tyy[iy*n1*n2+ix*n1+iz] += src_ampl;
					tzz[iy*n1*n2+ix*n1+iz] += src_ampl;
					txx[iy*n1*n2+(ix-1)*n1+iz] -= src_ampl;
					tyy[iy*n1*n2+(ix-1)*n1+iz] -= src_ampl;
					tzz[iy*n1*n2+(ix-1)*n1+iz] -= src_ampl;
				}
				else if (src.orient==4) { /* dipole +/0/- */
					if (iz > ibndz) {
						txx[iy*n1*n2+ix*n1+iz-1]+= 0.5*src_ampl;
						tyy[iy*n1*n2+ix*n1+iz-1]+= 0.5*src_ampl;
						tzz[iy*n1*n2+ix*n1+iz-1]+= 0.5*src_ampl;
					}
					if (iz < mod.nz+ibndz-1) {
						txx[iy*n1*n2+ix*n1+iz+1] -= 0.5*src_ampl;
						tyy[iy*n1*n2+ix*n1+iz+1] -= 0.5*src_ampl;
						tzz[iy*n1*n2+ix*n1+iz+1] -= 0.5*src_ampl;
					}
				}
				else if (src.orient==5) { /* dipole + - */
					txx[iy*n1*n2+ix*n1+iz] += src_ampl;
					tyy[iy*n1*n2+ix*n1+iz] += src_ampl;
					tzz[iy*n1*n2+ix*n1+iz] += src_ampl;
					txx[iy*n1*n2+(ix+1)*n1+iz] -= src_ampl;
					tyy[iy*n1*n2+(ix+1)*n1+iz] -= src_ampl;
					tzz[iy*n1*n2+(ix+1)*n1+iz] -= src_ampl;
				}
			}
			else if (src.type == 2) {
				/* Txz source */
				if ((iz == ibndz) && bnd.top==1) {
					txz[iy*n1*n2+(ix-1)*n1+iz-1] += src_ampl;
					txz[iy*n1*n2+ix*n1+iz-1] += src_ampl;
				}
				else {
					txz[iy*n1*n2+ix*n1+iz] += src_ampl;
				}
				/* possible dipole orientations for a txz source */
				if (src.orient == 2) { /* dipole +/- */
					txz[iy*n1*n2+ix*n1+iz+1] -= src_ampl;
				}
				else if (src.orient == 3) { /* dipole - + */
					txz[iy*n1*n2+(ix-1)*n1+iz] -= src_ampl;
				}
				else if (src.orient == 4) { /*  dipole +/O/- */
					/* correction: subtrace previous value to prevent z-1 values. */
					txz[iy*n1*n2+ix*n1+iz] -= 2.0*src_ampl;
					txz[iy*n1*n2+ix*n1+iz+1] += src_ampl;
				}
				else if (src.orient == 5) { /* dipole + - */
					txz[iy*n1*n2+(ix+1)*n1+iz] -= src_ampl;
				}
			}
			/* Tzz source */
			else if(src.type == 3) {
				tzz[iy*n1*n2+ix*n1+iz] += src_ampl;
			} 
			/* Txx source */
			else if(src.type == 4) {
				txx[iy*n1*n2+ix*n1+iz] += src_ampl;
			} 

/***********************************************************************
* pure potential shear S source (experimental)
* Curl S-pot = CURL(F) = dF_x/dz - dF_z/dx
***********************************************************************/
			else if(src.type == 5) {
				src_ampl = src_ampl*rox[iy][ix][iz]/(l2m[iy][ix][iz]);
				if (src.orient == 3) src_ampl = -src_ampl;
                /* first order derivatives */
				vx[iy*n1*n2+ix*n1+iz]         += src_ampl*sdx;
				vx[(iy-1)*n1*n2+ix*n1+iz-1]   -= src_ampl*sdx;
				vy[iy*n1*n2+ix*n1+iz]         += src_ampl*sdx;
				vy[iy*n1*n2+(ix-1)*n1+iz-1]   -= src_ampl*sdx;
				vz[iy*n1*n2+ix*n1+iz]         -= src_ampl*sdx;
				vz[(iy-1)*n1*n2+(ix-1)*n1+iz] += src_ampl*sdx;
                
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
                    vx[iy*n1*n2+ix*n1+iz]         -= src_ampl*sdx;
                    vx[(iy-1)*n1*n2+ix*n1+iz-1]   += src_ampl*sdx;
                    vy[iy*n1*n2+ix*n1+iz]         -= src_ampl*sdx;
                    vy[iy*n1*n2+(ix-1)*n1+iz-1]   += src_ampl*sdx;
                    vz[iy*n1*n2+ix*n1+iz]         += src_ampl*sdx;
                    vz[(iy-1)*n1*n2+(ix-1)*n1+iz] -= src_ampl*sdx;
				}
				else if (src.orient == 3) { /* dipole - + horizontal */
					ix += 1;
                    vx[iy*n1*n2+ix*n1+iz]         -= src_ampl*sdx;
                    vx[(iy-1)*n1*n2+ix*n1+iz-1]   += src_ampl*sdx;
                    vy[iy*n1*n2+ix*n1+iz]         -= src_ampl*sdx;
                    vy[iy*n1*n2+(ix-1)*n1+iz-1]   += src_ampl*sdx;
                    vz[iy*n1*n2+ix*n1+iz]         += src_ampl*sdx;
                    vz[(iy-1)*n1*n2+(ix-1)*n1+iz] -= src_ampl*sdx;
				}
            }
/***********************************************************************
* pure potential pressure P source (experimental)
* Divergence P-pot = DIV(F) = dF_x/dx + dF_z/dz
***********************************************************************/
            else if(src.type == 8) {
			    src_ampl = src_ampl*rox[iy][ix][iz]/(l2m[iy][ix][iz]);
                if (src.orient == 3) src_ampl = -src_ampl;
                vx[iy*n1*n2+(ix+1)*n1+iz] += src_ampl*sdx;
                vx[iy*n1*n2+ix*n1+iz]     -= src_ampl*sdx;
                vy[(iy+1)*n1*n2+ix*n1+iz] += src_ampl*sdx;
                vy[iy*n1*n2+ix*n1+iz]     -= src_ampl*sdx;
                vz[iy*n1*n2+ix*n1+iz+1]   += src_ampl*sdx;
                vz[iy*n1*n2+ix*n1+iz]     -= src_ampl*sdx;
                /* determine second position of dipole */
                if (src.orient == 2) { /* dipole +/- */
                    iz += 1;
                    vx[iy*n1*n2+(ix+1)*n1+iz] -= src_ampl*sdx;
                    vx[iy*n1*n2+ix*n1+iz]     += src_ampl*sdx;
                    vy[(iy+1)*n1*n2+ix*n1+iz] -= src_ampl*sdx;
                    vy[iy*n1*n2+ix*n1+iz]     += src_ampl*sdx;
                    vz[iy*n1*n2+ix*n1+iz+1]   -= src_ampl*sdx;
                    vz[iy*n1*n2+ix*n1+iz]     += src_ampl*sdx;
                }
                else if (src.orient == 3) { /* dipole - + */
                    ix += 1;
                    vx[iy*n1*n2+(ix+1)*n1+iz] -= src_ampl*sdx;
                    vx[iy*n1*n2+ix*n1+iz]     += src_ampl*sdx;
                    vy[(iy+1)*n1*n2+ix*n1+iz] -= src_ampl*sdx;
                    vy[iy*n1*n2+ix*n1+iz]     += src_ampl*sdx;
                    vz[iy*n1*n2+ix*n1+iz+1]   -= src_ampl*sdx;
                    vz[iy*n1*n2+ix*n1+iz]     += src_ampl*sdx;
                }
			}
            else if(src.type == 9) {
				Mxx = -1.0*(sin(src.dip)*cos(src.rake)*sin(2.0*src.strike)+sin(src.dip*2.0)*sin(src.rake)*sin(src.strike)*sin(src.strike));
				Myy = sin(src.dip)*cos(src.rake)*sin(2.0*src.strike)-sin(src.dip*2.0)*sin(src.rake)*cos(src.strike)*cos(src.strike);
				Mzz = sin(src.dip*2.0)*sin(src.rake);
				Mxz = -1.0*(cos(src.dip)*cos(src.rake)*cos(src.strike)+cos(src.dip*2.0)*sin(src.rake)*sin(src.strike));
				Mxy = sin(src.dip)*cos(src.rake)*cos(src.strike*2.0)+0.5*(sin(src.dip*2.0)*sin(src.rake)*sin(src.strike*2.0));
				Myz = -1.0*(cos(src.dip)*cos(src.rake)*sin(src.strike)-cos(src.dip*2.0)*sin(src.rake)*cos(src.strike));

				txx[iy*n1*n2+ix*n1+iz] -= Mxx*src_ampl;
				tyy[iy*n1*n2+ix*n1+iz] -= Myy*src_ampl;
				tzz[iy*n1*n2+ix*n1+iz] -= Mzz*src_ampl;
				txz[iy*n1*n2+ix*n1+iz] -= Mxz*src_ampl;
				txy[iy*n1*n2+ix*n1+iz] -= Mxy*src_ampl;
				tyz[iy*n1*n2+ix*n1+iz] -= Myz*src_ampl;
			} /* src.type */
		} /* ischeme */
	} /* loop over isrc */

	return 0;
}
