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
#include "fdelmodc3D.h"

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
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

long writeSnapTimes3D(modPar mod, snaPar sna, bndPar bnd, wavPar wav, long ixsrc, long iysrc, long izsrc, long itime, float *vx, float *vy, float *vz, float *tzz, float *tyy, float *txx, float *txz, float *tyz, float *txy, long verbose)
{
	FILE    *fpvx, *fpvy, *fpvz, *fptxx, *fptyy, *fptzz, *fptxz, *fptyz, *fptxy, *fpp, *fppp, *fpss;
	long append, isnap;
	static long first=1;
	long n1, n2, ibndx, ibndy, ibndz, ixs, iys, izs, ize, i, j, l;
	long ix, iy, iz, ix2, iy2;
	float *snap, sdx, stime;
	segy hdr;

	if (sna.nsnap==0) return 0;

    ibndx = mod.ioXx;
    ibndy = mod.ioXy;
    ibndz = mod.ioXz;
	n1    = mod.naz;
	n2    = mod.nax;
	sdx   = 1.0/mod.dx;

	if (sna.withbnd) {
		sna.nz=mod.naz;
		sna.z1=0;
		sna.z2=mod.naz-1;
		sna.skipdz=1;

		sna.ny=mod.nax;
		sna.y1=0;
		sna.y2=mod.nay-1;
		sna.skipdy=1;

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
		if (verbose) vmess("Writing snapshot(%li) at time=%.4f", isnap+1, stime);
	
		if (first) {
			append=0;
			first=0;
		}
		else {
			append=1;
		}

		if (sna.type.vx)  fpvx  = fileOpen(sna.file_snap, "_svx", (int)append);
		if (sna.type.vy)  fpvy  = fileOpen(sna.file_snap, "_svy", (int)append);
		if (sna.type.vz)  fpvz  = fileOpen(sna.file_snap, "_svz", (int)append);
		if (sna.type.p)   fpp   = fileOpen(sna.file_snap, "_sp", (int)append);
		if (sna.type.txx) fptxx = fileOpen(sna.file_snap, "_stxx", (int)append);
		if (sna.type.tyy) fptyy = fileOpen(sna.file_snap, "_styy", (int)append);
		if (sna.type.tzz) fptzz = fileOpen(sna.file_snap, "_stzz", (int)append);
		if (sna.type.txz) fptxz = fileOpen(sna.file_snap, "_stxz", (int)append);
		if (sna.type.tyz) fptyz = fileOpen(sna.file_snap, "_styz", (int)append);
		if (sna.type.txy) fptxy = fileOpen(sna.file_snap, "_stxy", (int)append);
		if (sna.type.pp)  fppp  = fileOpen(sna.file_snap, "_spp", (int)append);
		if (sna.type.ss)  fpss  = fileOpen(sna.file_snap, "_sss", (int)append);
	
		memset(&hdr,0,TRCBYTES);
		hdr.dt     = 1000000*(sna.skipdt*mod.dt);
		hdr.ungpow  = (sna.delay*mod.dt);
		hdr.scalco = -1000;
		hdr.scalel = -1000;
		hdr.sx     = 1000*(mod.x0+ixsrc*mod.dx);
		hdr.sy     = 1000*(mod.y0+iysrc*mod.dy);
		hdr.sdepth = 1000*(mod.z0+izsrc*mod.dz);
		hdr.fldr   = isnap+1;
		hdr.trid   = 1;
		hdr.ns     = sna.nz;
		hdr.trwf   = sna.nx*sna.nx;
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
		for (iys=sna.y1, l=0; iys<=sna.x2; iys+=sna.skipdx, l++) {
			for (ixs=sna.x1, i=0; ixs<=sna.x2; ixs+=sna.skipdx, i++) {
				hdr.tracf  = l*sna.nx+i+1;
				hdr.tracl  = isnap*sna.nx*sna.ny+l*sna.nx+i+1;
				hdr.gx     = 1000*(mod.x0+ixs*mod.dx);
				hdr.gy     = 1000*(mod.y0+ixs*mod.dy);
				ix  = ixs+ibndx;
				ix2 = ix+1;
				iy  = iys+ibndy;
				iy2 = iy+1;

				izs = sna.z1+ibndz;
				ize = sna.z2+ibndz;

				if (sna.withbnd) {
					izs = 0;
					ize = sna.z2;
					ix  = ixs;
					ix2 = ix;
					iy  = iys;
					iy2 = iy;
					if (sna.type.vz || sna.type.txz || sna.type.tyz) izs = -1;
					if ( !ISODD(bnd.lef)) hdr.gx = 1000*(mod.x0 - bnd.ntap*mod.dx);
					if ( !ISODD(bnd.fro)) hdr.gy = 1000*(mod.y0 - bnd.ntap*mod.dy);
				}

				if (sna.type.vx) {
					for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
						snap[j] = vx[iy*n1*n2+ix2*n1+iz];
					}
					traceWrite(&hdr, snap, (int)sna.nz, fpvx);
				}
				if (sna.type.vy) {
					for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
						snap[j] = vy[iy2*n1*n2+ix*n1+iz];
					}
					traceWrite(&hdr, snap, (int)sna.nz, fpvy);
				}
				if (sna.type.vz) {
					for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
						snap[j] = vz[iy*n1*n2+ix*n1+iz+1];
					}
					traceWrite(&hdr, snap, (int)sna.nz, fpvz);
				}
				if (sna.type.p) {
					for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
						snap[j] = tzz[iy*n1*n2+ix*n1+iz];
					}
					traceWrite(&hdr, snap, (int)sna.nz, fpp);
				}
				if (sna.type.tzz) {
					for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
						snap[j] = tzz[iy*n1*n2+ix*n1+iz];
					}
					traceWrite(&hdr, snap, (int)sna.nz, fptzz);
				}
				if (sna.type.tyy) {
					for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
						snap[j] = tyy[iy*n1*n2+ix*n1+iz];
					}
					traceWrite(&hdr, snap, (int)sna.nz, fptyy);
				}
				if (sna.type.txx) {
					for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
						snap[j] = txx[iy*n1*n2+ix*n1+iz];
					}
					traceWrite(&hdr, snap, (int)sna.nz, fptxx);
				}
				if (sna.type.txz) {
					for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
						snap[j] = txz[iy*n1*n2+ix2*n1+iz+1];
					}
					traceWrite(&hdr, snap, (int)sna.nz, fptxz);
				}
				if (sna.type.txy) {
					for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
						snap[j] = txy[iy2*n1*n2+ix2*n1+iz];
					}
					traceWrite(&hdr, snap, (int)sna.nz, fptxy);
				}
				if (sna.type.tyz) {
					for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
						snap[j] = tyz[iy2*n1*n2+ix*n1+iz+1];
					}
					traceWrite(&hdr, snap, (int)sna.nz, fptyz);
				}
				/* calculate divergence of velocity field */
				if (sna.type.pp) {
					for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
						snap[j] =  sdx*((vx[iy*n1*n2+(ix+1)*n1+iz]-vx[iy*n1*n2+ix*n1+iz])+
										(vy[(iy+1)*n1*n2+ix*n1+iz]-vy[iy*n1*n2+ix*n1+iz])+
										(vz[iy*n1*n2+ix*n1+iz+1]-vz[iy*n1*n2+ix*n1+iz]));
					}
					traceWrite(&hdr, snap, (int)sna.nz, fppp);
				}
				/* calculate rotation of velocity field */
				if (sna.type.ss) {
					for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
						snap[j] =  sdx*((vx[iy*n1*n2+ix*n1+iz]-vx[(iy-1)*n1*n2+ix*n1+iz-1])-
										(vy[iy*n1*n2+ix*n1+iz]-vy[iy*n1*n2+(ix-1)*n1+iz-1])-
										(vz[iy*n1*n2+ix*n1+iz]-vz[(iy-1)*n1*n2+(ix-1)*n1+iz]));
					}
					traceWrite(&hdr, snap, (int)sna.nz, fpss);
				}

			}
		}

		if (sna.type.vx) fclose(fpvx);
		if (sna.type.vy) fclose(fpvy);
		if (sna.type.vz) fclose(fpvz);
		if (sna.type.p) fclose(fpp);
		if (sna.type.txx) fclose(fptxx);
		if (sna.type.tyy) fclose(fptyy);
		if (sna.type.tzz) fclose(fptzz);
		if (sna.type.txz) fclose(fptxz);
		if (sna.type.tyz) fclose(fptyz);
		if (sna.type.txy) fclose(fptxy);
		if (sna.type.pp) fclose(fppp);
		if (sna.type.ss) fclose(fpss);

		free(snap);
	}

	return 0;
}

