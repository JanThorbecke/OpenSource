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
#include "fdelmodc3D.h"

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
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

long getBeamTimes3D(modPar mod, snaPar sna, float *vx, float *vy, float *vz,
    float *tzz, float *tyy, float *txx, float *txz, float *tyz, float *txy,
	float *beam_vx, float *beam_vy, float *beam_vz,
    float *beam_txx, float *beam_tyy, float *beam_tzz,
    float *beam_txz, float *beam_tyz, float *beam_txy,
	float *beam_p, float *beam_pp, float *beam_ss, long verbose)
{
	long n1, n2, ibndx, ibndy, ibndz, ixs, iys, izs, ize, i, j, l;
	long ix, iy, iz, ix2, iy2, iz2;
	float sdx, s, p;

    ibndx = mod.ioPx;
    ibndy = mod.ioPy;
    ibndz = mod.ioPz;
	n1   = mod.naz;
	n2   = mod.nax;
	sdx  = 1.0/mod.dx;
	izs = sna.z1+ibndx;
	ize = sna.z2+ibndz;

    for (iys=sna.y1, l=0; iys<=sna.y2; iys+=sna.skipdy, l++) {
        iy  = iys+ibndy;
        iy2 = iy+1;
        for (ixs=sna.x1, i=0; ixs<=sna.x2; ixs+=sna.skipdx, i++) {
            ix  = ixs+ibndx;
            ix2 = ix+1;

            if (sna.type.vx) {
                for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
                    beam_vx[l*sna.nz*sna.nz+i*sna.nz+j] += sqrt(vx[iy*n1*n2+ix2*n1+iz]*vx[iy*n1*n2+ix2*n1+iz]);
                }
            }
            if (sna.type.vy) {
                for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
                    beam_vy[l*sna.nz*sna.nz+i*sna.nz+j] += sqrt(vy[iy2*n1*n2+ix*n1+iz]*vx[iy2*n1*n2+ix*n1+iz]);
                }
            }
            if (sna.type.vz) {
                for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
                    beam_vz[l*sna.nz*sna.nz+i*sna.nz+j] += sqrt(vz[iy*n1*n2+ix*n1+iz+1]*vz[iy*n1*n2+ix*n1+iz+1]);
                }
            }
            if (sna.type.p) {
                for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
                    beam_p[l*sna.nz*sna.nz+i*sna.nz+j] += sqrt(tzz[iy*n1*n2+ix*n1+iz]*tzz[iy*n1*n2+ix*n1+iz]);
                }
            }
            if (sna.type.tzz) {
                for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
                    beam_tzz[l*sna.nz*sna.nz+i*sna.nz+j] += sqrt(tzz[iy*n1*n2+ix*n1+iz]*tzz[iy*n1*n2+ix*n1+iz]);
                }
            }
            if (sna.type.tyy) {
                for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
                    beam_tyy[l*sna.nz*sna.nz+i*sna.nz+j] += sqrt(tyy[iy*n1*n2+ix*n1+iz]*tyy[iy*n1*n2+ix*n1+iz]);
                }
            }
            if (sna.type.txx) {
                for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
                    beam_txx[l*sna.nz*sna.nz+i*sna.nz+j] += sqrt(txx[iy*n1*n2+ix*n1+iz]*txx[iy*n1*n2+ix*n1+iz]);
                }
            }
            if (sna.type.txz) {
                for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
                    beam_txz[l*sna.nz*sna.nz+i*sna.nz+j] += sqrt(txz[iy*n1*n2+ix2*n1+iz+1]*txz[iy*n1*n2+ix2*n1+iz+1]);
                }
            }
            if (sna.type.tyz) {
                for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
                    beam_tyz[l*sna.nz*sna.nz+i*sna.nz+j] += sqrt(tyz[iy2*n1*n2+ix*n1+iz+1]*tyz[iy2*n1*n2+ix*n1+iz+1]);
                }
            }
            if (sna.type.txz) {
                for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
                    beam_txy[l*sna.nz*sna.nz+i*sna.nz+j] += sqrt(txy[iy2*n1*n2+ix2*n1+iz]*txy[iy2*n1*n2+ix2*n1+iz]);
                }
            }
            /* calculate divergence of velocity field */
            if (sna.type.pp) {
                for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
                    iz2 = iz+1;
                    p =    sdx*((vx[iy*n1*n2+ix2*n1+iz]-vx[iy*n1*n2+ix*n1+iz])+
                                (vy[iy2*n1*n2+ix*n1+iz]-vy[iy*n1*n2+ix*n1+iz])+
                                (vz[iy*n1*n2+ix*n1+iz2]-vz[iy*n1*n2+ix*n1+iz]));
                    beam_pp[l*sna.nz*sna.nz+i*sna.nz+j] += sqrt(p*p);
                }
            }
            /* calculate rotation of velocity field */
            if (sna.type.ss) {
                for (iz=izs, j=0; iz<=ize; iz+=sna.skipdz, j++) {
                    iz2 = iz+1;
                    s =    sdx*((vx[iy2*n1*n2+ix2*n1+iz2]-vx[iy*n1*n2+ix2*n1+iz])-
                                (vy[iy2*n1*n2+ix2*n1+iz2]-vy[iy2*n1*n2+ix*n1+iz])-
                                (vz[iy2*n1*n2+ix2*n1+iz2]-vz[iy*n1*n2+ix*n1+iz2]));
                    beam_ss[l*sna.nz*sna.nz+i*sna.nz+j] += sqrt(s*s);
                }
            }
        }
    }
	return 0;
}


long writeBeams3D(modPar mod, snaPar sna, long ixsrc, long iysrc, long izsrc,
    long ishot, long fileno, float *beam_vx, float *beam_vy, float *beam_vz,
    float *beam_txx, float *beam_tyy, float *beam_tzz,
    float *beam_txz, float *beam_tyz, float *beam_txy, 
	float *beam_p, float *beam_pp, float *beam_ss, long verbose)
{
	FILE    *fpvx, *fpvy, *fpvz, *fptxx, *fptyy, *fptzz, *fptxz, *fptyz, *fptxy, *fpp, *fppp, *fpss;
	long append;
	long ix, iy;
	char number[16], filename[1024];
	segy hdr;

	if (sna.beam==0) return 0;
	/* all beam snapshots are written to the same output file(s) */
	if (ishot) append=1;
	else append=0;
	
	strcpy(filename, sna.file_beam);
	if (fileno) {
		sprintf(number,"_%03ld",fileno);
		name_ext(filename, number);
	}
	if (verbose>2) vmess("Writing beam data to file %s", filename);


	if (sna.type.vx)  fpvx  = fileOpen(filename, "_bvx", (int)append);
	if (sna.type.vy)  fpvy  = fileOpen(filename, "_bvy", (int)append);
	if (sna.type.vz)  fpvz  = fileOpen(filename, "_bvz", (int)append);
	if (sna.type.p)   fpp   = fileOpen(filename, "_bp", (int)append);
	if (sna.type.txx) fptxx = fileOpen(filename, "_btxx", (int)append);
	if (sna.type.tyy) fptyy = fileOpen(filename, "_btyy", (int)append);
	if (sna.type.tzz) fptzz = fileOpen(filename, "_btzz", (int)append);
	if (sna.type.txz) fptxz = fileOpen(filename, "_btxz", (int)append);
	if (sna.type.tyz) fptyz = fileOpen(filename, "_btyz", (int)append);
	if (sna.type.txy) fptxy = fileOpen(filename, "_btxy", (int)append);
	if (sna.type.pp)  fppp  = fileOpen(filename, "_bpp", (int)append);
	if (sna.type.ss)  fpss  = fileOpen(filename, "_bss", (int)append);
	
	memset(&hdr,0,TRCBYTES);
	hdr.dt     = 1000000*(mod.dt);
	hdr.scalco = -1000;
	hdr.scalel = -1000;
	hdr.sx     = 1000*(mod.x0+ixsrc*mod.dx);
	hdr.sy     = 1000*(mod.y0+iysrc*mod.dy);
	hdr.sdepth = 1000*(mod.z0+izsrc*mod.dz);
	hdr.fldr   = ishot+1;
	hdr.trid   = 1;
	hdr.ns     = sna.nz;
	hdr.trwf   = sna.nx*sna.ny;
	hdr.ntr    = sna.nx*sna.ny;
	hdr.f1     = sna.z1*mod.dz+mod.z0;
	hdr.f2     = sna.x1*mod.dx+mod.x0;
	hdr.d1     = mod.dz*sna.skipdz;
	hdr.d2     = mod.dx*sna.skipdx;

    for (iy=0; iy<sna.ny; iy++) {
        for (ix=0; ix<sna.nx; ix++) {
            hdr.tracf  = iy*sna.nx+ix+1;
            hdr.tracl  = iy*sna.nx+ix+1;
            hdr.gx     = 1000*(mod.x0+(sna.x1+ix)*mod.dx);
            hdr.gy     = 1000*(mod.y0+(sna.y1+ix)*mod.dy);

            if (sna.type.vx) {
                traceWrite( &hdr, &beam_vx[iy*sna.nx*sna.nz+ix*sna.nz], (int)sna.nz, fpvx) ;
            }
            if (sna.type.vy) {
                traceWrite( &hdr, &beam_vy[iy*sna.nx*sna.nz+ix*sna.nz], (int)sna.nz, fpvy) ;
            }
            if (sna.type.vz) {
                traceWrite( &hdr, &beam_vz[iy*sna.nx*sna.nz+ix*sna.nz], (int)sna.nz, fpvz) ;
            }
            if (sna.type.p) {
                traceWrite( &hdr, &beam_p[iy*sna.nx*sna.nz+ix*sna.nz], (int)sna.nz, fpp) ;
            }
            if (sna.type.tzz) {
                traceWrite( &hdr, &beam_tzz[iy*sna.nx*sna.nz+ix*sna.nz], (int)sna.nz, fptzz) ;
            }
            if (sna.type.tyy) {
                traceWrite( &hdr, &beam_tyy[iy*sna.nx*sna.nz+ix*sna.nz], (int)sna.nz, fptyy) ;
            }
            if (sna.type.txx) {
                traceWrite( &hdr, &beam_txx[iy*sna.nx*sna.nz+ix*sna.nz], (int)sna.nz, fptxx) ;
            }
            if (sna.type.txz) {
                traceWrite( &hdr, &beam_txz[iy*sna.nx*sna.nz+ix*sna.nz], (int)sna.nz, fptxz) ;
            }
            if (sna.type.txy) {
                traceWrite( &hdr, &beam_txy[iy*sna.nx*sna.nz+ix*sna.nz], (int)sna.nz, fptxy) ;
            }
            if (sna.type.tyz) {
                traceWrite( &hdr, &beam_tyz[iy*sna.nx*sna.nz+ix*sna.nz], (int)sna.nz, fptyz) ;
            }
            if (sna.type.pp) {
                traceWrite( &hdr, &beam_pp[iy*sna.nx*sna.nz+ix*sna.nz], (int)sna.nz, fppp) ;
            }
            if (sna.type.ss) {
                traceWrite( &hdr, &beam_ss[iy*sna.nx*sna.nz+ix*sna.nz], (int)sna.nz, fpss) ;
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

	return 0;
}
