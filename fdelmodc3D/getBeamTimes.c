#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include "segy.h"
#include "fdelmodc.h"

/**
*  getBeamTimes: stores energy fields (beams) in arrays at certain time steps 
*  writeBeams: writes the stored fields to output file(s) 
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/


FILE *fileOpen(char *file, char *ext, int append);
int traceWrite(segy *hdr, float *data, int n, FILE *fp);
void name_ext(char *filename, char *extension);
void vmess(char *fmt, ...);

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int getBeamTimes(modPar mod, snaPar sna, float *vx, float *vz, float *tzz, float *txx, float *txz, 
				 float *beam_vx, float *beam_vz, float *beam_txx, float *beam_tzz, float *beam_txz, 
				 float *beam_p, float *beam_pp, float *beam_ss, int verbose)
{
	int n1, ibndx, ibndz, ixs, izs, ize, i, j;
	int ix, iz, ix2, iz2;
	float sdx, s, p;

    ibndx = mod.ioPx;
    ibndz = mod.ioPz;
	n1   = mod.naz;
	sdx  = 1.0/mod.dx;
	izs = sna.z1+ibndx;
	ize = sna.z2+ibndz;

	for (ixs=sna.x1, i=0; ixs<=sna.x2; ixs+=sna.skipdx, i++) {
		ix  = ixs+ibndx;
		ix2 = ix+1;

		if (sna.type.vx) {
			for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
				beam_vx[i*sna.nz+j] += sqrt(vx[ix2*n1+iz]*vx[ix2*n1+iz]);
			}
		}
		if (sna.type.vz) {
			for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
				beam_vz[i*sna.nz+j] += sqrt(vz[ix*n1+iz+1]*vz[ix*n1+iz+1]);
			}
		}
		if (sna.type.p) {
			for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
				beam_p[i*sna.nz+j] += sqrt(tzz[ix*n1+iz]*tzz[ix*n1+iz]);
			}
		}
		if (sna.type.tzz) {
			for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
				beam_tzz[i*sna.nz+j] += sqrt(tzz[ix*n1+iz]*tzz[ix*n1+iz]);
			}
		}
		if (sna.type.txx) {
			for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
				beam_txx[i*sna.nz+j] += sqrt(txx[ix*n1+iz]*txx[ix*n1+iz]);
			}
		}
		if (sna.type.txz) {
			for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
				beam_txz[i*sna.nz+j] += sqrt(txz[ix2*n1+iz+1]*txz[ix2*n1+iz+1]);
			}
		}
		/* calculate divergence of velocity field */
		if (sna.type.pp) {
			for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
				iz2 = iz+1;
				p = sdx*((vx[ix2*n1+iz]-vx[ix*n1+iz])+
							   (vz[ix*n1+iz2]-vz[ix*n1+iz]));
				beam_pp[i*sna.nz+j] += sqrt(p*p);
			}
		}
		/* calculate rotation of velocity field */
		if (sna.type.ss) {
			for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
				iz2 = iz+1;
				s = sdx*((vx[ix2*n1+iz2]-vx[ix2*n1+iz])-
							   (vz[ix2*n1+iz2]-vz[ix*n1+iz2]));
				beam_ss[i*sna.nz+j] += sqrt(s*s);
			}
		}
	}
	return 0;
}


int writeBeams(modPar mod, snaPar sna, int ixsrc, int izsrc, int ishot, int fileno, 
			   float *beam_vx, float *beam_vz, float *beam_txx, float *beam_tzz, float *beam_txz, 
			   float *beam_p, float *beam_pp, float *beam_ss, int verbose)
{
	FILE    *fpvx, *fpvz, *fptxx, *fptzz, *fptxz, *fpp, *fppp, *fpss;
	int append;
	int ix;
	char number[16], filename[1024];
	segy hdr;

	if (sna.beam==0) return 0;
	/* all beam snapshots are written to the same output file(s) */
	if (ishot) append=1;
	else append=0;
	
	strcpy(filename, sna.file_beam);
	if (fileno) {
		sprintf(number,"_%03d",fileno);
		name_ext(filename, number);
	}
	if (verbose>2) vmess("Writing beam data to file %s", filename);


	if (sna.type.vx)  fpvx  = fileOpen(filename, "_bvx", append);
	if (sna.type.vz)  fpvz  = fileOpen(filename, "_bvz", append);
	if (sna.type.p)   fpp   = fileOpen(filename, "_bp", append);
	if (sna.type.txx) fptxx = fileOpen(filename, "_btxx", append);
	if (sna.type.tzz) fptzz = fileOpen(filename, "_btzz", append);
	if (sna.type.txz) fptxz = fileOpen(filename, "_btxz", append);
	if (sna.type.pp)  fppp  = fileOpen(filename, "_bpp", append);
	if (sna.type.ss)  fpss  = fileOpen(filename, "_bss", append);
	
	memset(&hdr,0,TRCBYTES);
	hdr.dt     = 1000000*(mod.dt);
	hdr.scalco = -1000;
	hdr.scalel = -1000;
	hdr.sx     = 1000*(mod.x0+ixsrc*mod.dx);
	hdr.sdepth = 1000*(mod.z0+izsrc*mod.dz);
	hdr.fldr   = ishot+1;
	hdr.trid   = 1;
	hdr.ns     = sna.nz;
	hdr.trwf   = sna.nx;
	hdr.ntr    = sna.nx;
	hdr.f1     = sna.z1*mod.dz+mod.z0;
	hdr.f2     = sna.x1*mod.dx+mod.x0;
	hdr.d1     = mod.dz*sna.skipdz;
	hdr.d2     = mod.dx*sna.skipdx;

	for (ix=0; ix<sna.nx; ix++) {
		hdr.tracf  = ix+1;
		hdr.tracl  = ix+1;
		hdr.gx     = 1000*(mod.x0+(sna.x1+ix)*mod.dx);

		if (sna.type.vx) {
			traceWrite( &hdr, &beam_vx[ix*sna.nz], sna.nz, fpvx) ;
		}
		if (sna.type.vz) {
			traceWrite( &hdr, &beam_vz[ix*sna.nz], sna.nz, fpvz) ;
		}
		if (sna.type.p) {
			traceWrite( &hdr, &beam_p[ix*sna.nz], sna.nz, fpp) ;
		}
		if (sna.type.tzz) {
			traceWrite( &hdr, &beam_tzz[ix*sna.nz], sna.nz, fptzz) ;
		}
		if (sna.type.txx) {
			traceWrite( &hdr, &beam_txx[ix*sna.nz], sna.nz, fptxx) ;
		}
		if (sna.type.txz) {
			traceWrite( &hdr, &beam_txz[ix*sna.nz], sna.nz, fptxz) ;
		}
		if (sna.type.pp) {
			traceWrite( &hdr, &beam_pp[ix*sna.nz], sna.nz, fppp) ;
		}
		if (sna.type.ss) {
			traceWrite( &hdr, &beam_ss[ix*sna.nz], sna.nz, fpss) ;
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

	return 0;
}

