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
#include <genfft.h>

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
typedef struct _dcomplexStruct { /* complex number */
    double r,i;
} dcomplex;
#endif/* complex */

/**
*  Writes the receiver array(s) to output file(s)
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

FILE *fileOpen(char *file, char *ext, int append);
int traceWrite(segy *hdr, float *data, int n, FILE *fp) ;
void name_ext(char *filename, char *extension);
void kxwdecomp(complex *rp, complex *rvz, complex *up, complex *down,
               int nkx, float dx, int nt, float dt, float fmin, float fmax,
               float cp, float rho, int vznorm, int verbose);


#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

long writeRec3D(recPar rec, modPar mod, bndPar bnd, wavPar wav, long ixsrc, long iysrc, long izsrc,
    long nsam, long ishot, long fileno, float *rec_vx, float *rec_vy, float *rec_vz,
    float *rec_txx, float *rec_tyy, float *rec_tzz, float *rec_txz,  float *rec_tyz, 
    float *rec_txy, float *rec_p, float *rec_pp, float *rec_ss,
    float *rec_udp, float *rec_udvz, long verbose)
{
    FILE    *fpvx, *fpvy, *fpvz, *fptxx, *fptyy, *fptzz, *fptxz, *fptyz, *fptxy, *fpp, *fppp, *fpss, *fpup, *fpdown;
    float *rec_up, *rec_down, *trace, *rec_vze, *rec_pe;
    float dx, dy, dt, cp, rho, fmin, fmax;
    complex *crec_vz, *crec_p, *crec_up, *crec_dw;
    long irec, ntfft, nfreq, nkx, nky, xorig, yorig, ix, iy, iz, it, ibndx, ibndy;
    long append, vznorm, sx, sy;
    double ddt;
    char number[16], filename[1024];
    segy hdr;
    
    if (!rec.n) return 0;
    if (ishot) append=1;
    else append=0;

    /* if the total number of samples exceeds rec_ntsam then a new (numbered) file is opened */
    /* fileno has a non-zero value (from fdelmodc.c) if the number of samples exceeds rec_ntsam. */
    strcpy(filename, rec.file_rcv);
    if (fileno) {
        sprintf(number,"_%03d",fileno);
        name_ext(filename, number);
    }
#ifdef MPI
    sx = (long)mod.x0+ixsrc*mod.dx;
    sprintf(number,"_%06d",sx);
    name_ext(filename, number);
#endif

    if (verbose>2) vmess("Writing receiver data to file %s", filename);
    if (nsam != rec.nt && verbose) vmess("Number of samples written to last file = %li",nsam);

    memset(&hdr,0,TRCBYTES);
    ddt = (double)mod.dt;/* to avoid rounding in 32 bit precision */
    dt  = (float)ddt*rec.skipdt;
    dx  = (rec.x[1]-rec.x[0])*mod.dx;
    dy  = (rec.y[1]-rec.y[0])*mod.dy;
    hdr.dt     = (unsigned short)lround((((double)1.0e6*ddt*rec.skipdt)));
    hdr.scalco = -1000;
    hdr.scalel = -1000;
    hdr.sx     = 1000*(mod.x0+ixsrc*mod.dx);
    hdr.sy     = 1000*(mod.y0+ixsrc*mod.dy);
    hdr.sdepth = 1000*(mod.z0+izsrc*mod.dz);
    hdr.selev  = (int)(-1000.0*(mod.z0+izsrc*mod.dz));
    hdr.fldr   = ishot+1;
    hdr.trid   = 1;
    hdr.ns     = nsam;
    hdr.trwf   = rec.n;
    hdr.ntr    = rec.n;
    if (mod.grid_dir) { /* reverse time modeling */
        hdr.f1 = (-mod.nt+1)*mod.dt;
    }
    else {
        hdr.f1 = 0.0;
    }
    hdr.d1     = mod.dt*rec.skipdt;
    hdr.d2     = (rec.x[1]-rec.x[0])*mod.dx;
    hdr.f2     = mod.x0+rec.x[0]*mod.dx;

    if (rec.type.vx)  fpvx  = fileOpen(filename, "_rvx", (int)append);
    if (rec.type.vy)  fpvy  = fileOpen(filename, "_rvy", (int)append);
    if (rec.type.vz)  fpvz  = fileOpen(filename, "_rvz", (int)append);
    if (rec.type.p)   fpp   = fileOpen(filename, "_rp", (int)append);
    if (rec.type.txx) fptxx = fileOpen(filename, "_rtxx", (int)append);
    if (rec.type.tyy) fptyy = fileOpen(filename, "_rtyy", (int)append);
    if (rec.type.tzz) fptzz = fileOpen(filename, "_rtzz", (int)append);
    if (rec.type.txz) fptxz = fileOpen(filename, "_rtxz", (int)append);
    if (rec.type.tyz) fptyz = fileOpen(filename, "_rtyz", (int)append);
    if (rec.type.txy) fptxy = fileOpen(filename, "_rtxy", (int)append);
    if (rec.type.pp)  fppp  = fileOpen(filename, "_rpp", (int)append);
    if (rec.type.ss)  fpss  = fileOpen(filename, "_rss", (int)append);

    /* decomposed wavefield */
    if (rec.type.ud && (mod.ischeme==1 || mod.ischeme==2) )  {
        fpup   = fileOpen(filename, "_ru", (int)append);
        fpdown = fileOpen(filename, "_rd", (int)append);
        ntfft = optncr(nsam);
        nfreq = ntfft/2+1;
        fmin = 0.0;
        fmax = wav.fmax;
        nkx = optncc(2*mod.nax);
        nky = optncc(2*mod.nay);
        ibndx = mod.ioPx;
        ibndy = mod.ioPy;
        if (bnd.lef==4 || bnd.lef==2) ibndx += bnd.ntap;
        if (bnd.fro==4 || bnd.fro==2) ibndy += bnd.ntap;
        cp  = rec.cp;
        rho = rec.rho;
		if (rec.type.ud==2) vznorm=1;
		else vznorm=0;
        if (verbose) vmess("Decomposition array at z=%.2f with cp=%.2f rho=%.2f", rec.zr[0]+mod.z0, cp, rho);
        rec_up  = (float *)calloc(ntfft*nkx*nky,sizeof(float));
        rec_down= (float *)calloc(ntfft*nkx*nky,sizeof(float));
        crec_vz = (complex *)malloc(nfreq*nkx*nky*sizeof(complex));
        crec_p  = (complex *)malloc(nfreq*nkx*nky*sizeof(complex));
        crec_up = (complex *)malloc(nfreq*nkx*nky*sizeof(complex));
        crec_dw = (complex *)malloc(nfreq*nkx*nky*sizeof(complex));

        rec_vze = rec_up;
        rec_pe  = rec_down;
        /* copy input data into extended arrays with padded zeroes */
        for (iy=0; iy<mod.nay; iy++) {
            for (ix=0; ix<mod.nax; ix++) {
                memcpy(&rec_vze[iy*mod.nax*ntfft+ix*ntfft],&rec_udvz[iy*mod.nax*rec.nt+ix*rec.nt],nsam*sizeof(float));
                memcpy(&rec_pe[iy*mod.nax*ntfft+ix*ntfft], &rec_udp[iy*mod.nax*rec.nt+ix*rec.nt], nsam*sizeof(float));
            }
        }

        /* transform from t-x to kx-w */
        xorig = ixsrc+ibndx;
        xt2wkx(rec_vze, crec_vz, ntfft, nkx, ntfft, nkx, xorig);
        xt2wkx(rec_pe, crec_p, ntfft, nkx, ntfft, nkx, xorig);

        /* apply decomposition operators */
        kxwdecomp(crec_p, crec_vz, crec_up, crec_dw,
            nkx, mod.dx, nsam, dt, fmin, fmax, cp, rho, vznorm, verbose);

        /* transform back to t-x */
        wkx2xt(crec_up, rec_up, ntfft, nkx, nkx, ntfft, xorig);
        wkx2xt(crec_dw, rec_down, ntfft, nkx, nkx, ntfft, xorig);

        /* reduce array to rec.nt samples rec.n traces */
        for (irec=0; irec<rec.n; irec++) {
            ix = rec.x[irec]+ibndx;
            iy = rec.y[irec]+ibndy;
            for (it=0; it<rec.nt; it++) {
                rec_up[irec*rec.nt+it]   = rec_up[iy*nkx*ntfft+ix*ntfft+it];
                rec_down[irec*rec.nt+it] = rec_down[iy*nkx*ntfft+ix*ntfft+it];
            }
        }
        free(crec_vz);
        free(crec_p);
        free(crec_up);
        free(crec_dw);
    }
    if (rec.type.ud && (mod.ischeme==3 || mod.ischeme==4) )  {
    }

    for (irec=0; irec<rec.n; irec++) {
        hdr.tracf  = irec+1;
        hdr.tracl  = ishot*rec.n+irec+1;
        hdr.gx     = 1000*(mod.x0+rec.x[irec]*mod.dx);
        hdr.gy     = 1000*(mod.y0+rec.y[irec]*mod.dy);
        hdr.offset = (rec.x[irec]-ixsrc)*mod.dx;
        hdr.gelev  = (int)(-1000*(mod.z0+rec.z[irec]*mod.dz));

        if (rec.type.vx) {
            traceWrite( &hdr, &rec_vx[irec*rec.nt], (int)nsam, fpvx) ;
        }
        if (rec.type.vy) {
            traceWrite( &hdr, &rec_vy[irec*rec.nt], (int)nsam, fpvy) ;
        }
        if (rec.type.vz) {
            traceWrite( &hdr, &rec_vz[irec*rec.nt], (int)nsam, fpvz) ;
        }
        if (rec.type.p) {
            traceWrite( &hdr, &rec_p[irec*rec.nt], (int)nsam, fpp) ;
        }
        if (rec.type.txx) {
            traceWrite( &hdr, &rec_txx[irec*rec.nt], (int)nsam, fptxx) ;
        }
        if (rec.type.tyy) {
            traceWrite( &hdr, &rec_tyy[irec*rec.nt], (int)nsam, fptyy) ;
        }
        if (rec.type.tzz) {
            traceWrite( &hdr, &rec_tzz[irec*rec.nt], (int)nsam, fptzz) ;
        }
        if (rec.type.txz) {
            traceWrite( &hdr, &rec_txz[irec*rec.nt], (int)nsam, fptxz) ;
        }
        if (rec.type.txy) {
            traceWrite( &hdr, &rec_txy[irec*rec.nt], (int)nsam, fptxy) ;
        }
        if (rec.type.tyz) {
            traceWrite( &hdr, &rec_tyz[irec*rec.nt], (int)nsam, fptyz) ;
        }
        if (rec.type.pp) {
            traceWrite( &hdr, &rec_pp[irec*rec.nt], (int)nsam, fppp) ;
        }
        if (rec.type.ss) {
            traceWrite( &hdr, &rec_ss[irec*rec.nt], (int)nsam, fpss) ;
        }
        if (rec.type.ud && mod.ischeme==1)  {
            traceWrite( &hdr, &rec_up[irec*rec.nt], (int)nsam, fpup) ;
            traceWrite( &hdr, &rec_down[irec*rec.nt], (int)nsam, fpdown) ;
        }
    }

    if (rec.type.vx) fclose(fpvx);
    if (rec.type.vy) fclose(fpvy);
    if (rec.type.vz) fclose(fpvz);
    if (rec.type.p) fclose(fpp);
    if (rec.type.txx) fclose(fptxx);
    if (rec.type.tyy) fclose(fptyy);
    if (rec.type.tzz) fclose(fptzz);
    if (rec.type.txz) fclose(fptxz);
    if (rec.type.txy) fclose(fptxy);
    if (rec.type.tyz) fclose(fptyz);
    if (rec.type.pp) fclose(fppp);
    if (rec.type.ss) fclose(fpss);
    if (rec.type.ud) {
        fclose(fpup);
        fclose(fpdown);
        free(rec_up);
        free(rec_down);
    }

    return 0;
}

