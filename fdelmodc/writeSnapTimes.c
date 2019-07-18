#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#define ISODD(n) ((n) & 01)

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include "par.h"
#include "segy.h"
#include "fdelmodc.h"

/**
*  Writes gridded wavefield(s) at a desired time to output file(s) 
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/


FILE *fileOpen(char *file, char *ext, int append);
int traceWrite(segy *hdr, float *data, int n, FILE *fp);

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int writeSnapTimes(modPar mod, snaPar sna, bndPar bnd, wavPar wav, int ixsrc, int izsrc, int itime, float *vx, float *vz, float *tzz, float *txx, float *txz, int verbose)
{
	FILE    *fpvx, *fpvz, *fptxx, *fptzz, *fptxz, *fpp, *fppp, *fpss;
	int append, isnap;
	static int first=1;
	int n1, ibndx, ibndz, ixs, izs, ize, i, j;
	int ix, iz, ix2;
	float *snap, sdx, stime;
	segy hdr;

	if (sna.nsnap==0) return 0;

    ibndx = mod.ioXx;
    ibndz = mod.ioXz;
	n1    = mod.naz;
	sdx   = 1.0/mod.dx;

	if (sna.withbnd) {
		sna.nz=mod.naz;
		sna.z1=0;
		sna.z2=mod.naz-1;
		sna.skipdz=1;

		sna.nx=mod.nax;
		sna.x1=0;
		sna.x2=mod.nax-1;
		sna.skipdx=1;
	}

	/* check if this itime is a desired snapshot time */
	if ( (((itime-sna.delay) % sna.skipdt)==0) && 
		  (itime >= sna.delay) &&
		  (itime <= sna.delay+(sna.nsnap-1)*sna.skipdt) ) {

		isnap = NINT((itime-sna.delay)/sna.skipdt);

        if (mod.grid_dir) stime = (-wav.nt+1+itime+1)*mod.dt;  /* reverse time modeling */
        else  stime = itime*mod.dt;
		if (verbose>1) vmess("Writing snapshot(%d) at time=%.4f", isnap+1, stime);
	
		if (first) {
			append=0;
			first=0;
		}
		else {
			append=1;
		}

		if (sna.type.vx)  fpvx  = fileOpen(sna.file_snap, "_svx", append);
		if (sna.type.vz)  fpvz  = fileOpen(sna.file_snap, "_svz", append);
		if (sna.type.p)   fpp   = fileOpen(sna.file_snap, "_sp", append);
		if (sna.type.txx) fptxx = fileOpen(sna.file_snap, "_stxx", append);
		if (sna.type.tzz) fptzz = fileOpen(sna.file_snap, "_stzz", append);
		if (sna.type.txz) fptxz = fileOpen(sna.file_snap, "_stxz", append);
		if (sna.type.pp)  fppp  = fileOpen(sna.file_snap, "_spp", append);
		if (sna.type.ss)  fpss  = fileOpen(sna.file_snap, "_sss", append);
	
		memset(&hdr,0,TRCBYTES);
		hdr.dt     = 1000000*(sna.skipdt*mod.dt);
		hdr.ungpow  = (sna.delay*mod.dt);
		hdr.scalco = -1000;
		hdr.scalel = -1000;
		hdr.sx     = 1000*(mod.x0+ixsrc*mod.dx);
		hdr.sdepth = 1000*(mod.z0+izsrc*mod.dz);
		hdr.fldr   = isnap+1;
		hdr.trid   = 1;
		hdr.ns     = sna.nz;
		hdr.trwf   = sna.nx;
		hdr.ntr    = (isnap+1)*sna.nx;
		hdr.f1     = sna.z1*mod.dz+mod.z0;
		hdr.f2     = sna.x1*mod.dx+mod.x0;
		hdr.d1     = mod.dz*sna.skipdz;
		hdr.d2     = mod.dx*sna.skipdx;
		if (sna.withbnd) {
        	if ( !ISODD(bnd.top)) hdr.f1 = mod.z0 - bnd.ntap*mod.dz;
        	if ( !ISODD(bnd.lef)) hdr.f2 = mod.x0 - bnd.ntap*mod.dx;
        	//if ( !ISODD(bnd.rig)) ;
        	//if ( !ISODD(bnd.bot)) store=1;
		}

/***********************************************************************
* vx velocities have one sample less in x-direction
* vz velocities have one sample less in z-direction
* txz stresses have one sample less in z-direction and x-direction
***********************************************************************/

		snap = (float *)malloc(sna.nz*sizeof(float));

		/* Decimate, with skipdx and skipdz, the number of gridpoints written to file 
		   and write to file. */
		for (ixs=sna.x1, i=0; ixs<=sna.x2; ixs+=sna.skipdx, i++) {
			hdr.tracf  = i+1;
			hdr.tracl  = isnap*sna.nx+i+1;
			hdr.gx     = 1000*(mod.x0+ixs*mod.dx);
			ix = ixs+ibndx;
			ix2 = ix+1;

			izs = sna.z1+ibndz;
			ize = sna.z2+ibndz;

			if (sna.withbnd) {
				izs = 0;
				ize = sna.z2;
				ix = ixs;
				ix2 = ix;
				if (sna.type.vz || sna.type.txz) izs = -1;
        		if ( !ISODD(bnd.lef)) hdr.gx = 1000*(mod.x0 - bnd.ntap*mod.dx);
			}

			if (sna.type.vx) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					snap[j] = vx[ix2*n1+iz];
				}
				traceWrite(&hdr, snap, sna.nz, fpvx);
			}
			if (sna.type.vz) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					snap[j] = vz[ix*n1+iz+1];
				}
				traceWrite(&hdr, snap, sna.nz, fpvz);
			}
			if (sna.type.p) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					snap[j] = tzz[ix*n1+iz];
				}
				traceWrite(&hdr, snap, sna.nz, fpp);
			}
			if (sna.type.tzz) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					snap[j] = tzz[ix*n1+iz];
				}
				traceWrite(&hdr, snap, sna.nz, fptzz);
			}
			if (sna.type.txx) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					snap[j] = txx[ix*n1+iz];
				}
				traceWrite(&hdr, snap, sna.nz, fptxx);
			}
			if (sna.type.txz) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					snap[j] = txz[ix2*n1+iz+1];
				}
				traceWrite(&hdr, snap, sna.nz, fptxz);
			}
			/* calculate divergence of velocity field */
			if (sna.type.pp) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					snap[j] = sdx*((vx[(ix+1)*n1+iz]-vx[ix*n1+iz])+
									(vz[ix*n1+iz+1]-vz[ix*n1+iz]));
				}
				traceWrite(&hdr, snap, sna.nz, fppp);
			}
			/* calculate rotation of velocity field */
			if (sna.type.ss) {
				for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
					snap[j] = sdx*((vx[ix*n1+iz]-vx[ix*n1+iz-1])-
									(vz[ix*n1+iz]-vz[(ix-1)*n1+iz]));
				}
				traceWrite(&hdr, snap, sna.nz, fpss);
			}

		}

		if (sna.type.vx) fclose(fpvx);
		if (sna.type.vz) fclose(fpvz);
		if (sna.type.p) fclose(fpp);
		if (sna.type.txx) fclose(fptxx);
		if (sna.type.tzz) fclose(fptzz);
		if (sna.type.txz) fclose(fptxz);
		if (sna.type.pp) fclose(fppp);
		if (sna.type.ss) fclose(fpss);

		free(snap);
	}

	return 0;
}

