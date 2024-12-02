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
#include <genfft.h>
#ifdef MPI
#include<mpi.h>
#endif

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

#ifdef MPI
static MPI_File fpvx, fpvz, fptxx, fptzz, fptxz, fpp, fppp, fpss, fpup, fpdown, fpdxvx, fpdzvz;
MPI_File fileOpen(char *file, char *ext, int append);
int traceWrite(segy *hdr, float *data, int n, long long offset, MPI_File fh);
void fileClose(MPI_File fh);
static int opened;
#else
static FILE *fpvx, *fpvz, *fptxx, *fptzz, *fptxz, *fpp, *fppp, *fpss, *fpup, *fpdown, *fpdxvx, *fpdzvz, *fpq;
FILE *fileOpen(char *file, char *ext, int append);
int traceWrite(segy *hdr, float *data, int n, long long offset, FILE *fp);
void fileClose(FILE *fp);
#endif
void name_ext(char *filename, char *extension);
void kxwdecomp(complex *rp, complex *rvz, complex *up, complex *down,
               int nkx, float dx, int nt, float dt, float fmin, float fmax,
               float cp, float rho, int vznorm, int verbose);

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int writeRec(recPar rec, modPar mod, bndPar bnd, wavPar wav, int ixsrc, int izsrc, int nsam, int ishot, int nshots, int fileno, 
             float *rec_vx, float *rec_vz, float *rec_txx, float *rec_tzz, float *rec_txz, 
             float *rec_p, float *rec_pp, float *rec_ss, float *rec_q, float *rec_udp, float *rec_udvz, float *rec_dxvx, float *rec_dzvz, int verbose)
{
    float *rec_up, *rec_down, *trace, *rec_vze, *rec_pe;
    float dx, dt, cp, rho, fmin, fmax;
    complex *crec_vz, *crec_p, *crec_up, *crec_dw;
    int irec, ntfft, nfreq, nkx, xorig, ix, iz, it, ibndx;
    int append, vznorm, sx;
    long long offset;
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

    if (verbose>2) vmess("Writing receiver data to file %s", filename);
    if (nsam != rec.nt && verbose) vmess("Number of samples written to last file = %d",nsam);

#ifdef MPI
    int pe;
    MPI_Comm_rank( MPI_COMM_WORLD, &pe );
    if (verbose>2) vmess("PE %d writes to file %s", pe, filename);
#endif
    memset(&hdr,0,TRCBYTES);
    ddt = (double)mod.dt;/* to avoid rounding in 32 bit precision */
    dt  = (float)ddt*rec.skipdt;
    dx  = (rec.x[1]-rec.x[0])*mod.dx;
    hdr.dt     = (unsigned short)lround((((double)1.0e6*ddt*rec.skipdt)));
    hdr.scalco = -1000;
    hdr.scalel = -1000;
    hdr.sx     = 1000*(mod.x0+ixsrc*mod.dx);
    hdr.sdepth = 1000*(mod.z0+izsrc*mod.dz);
    hdr.selev  = (int)(-1000.0*(mod.z0+izsrc*mod.dz));
    hdr.fldr   = ishot+1;
    hdr.trid   = 1;
    hdr.ns     = nsam;
    hdr.trwf   = 0;
    hdr.ntr    = nshots*rec.n;
    if (mod.grid_dir) { /* reverse time modeling */
        hdr.f1 = (-mod.nt+1)*mod.dt;
    }
    else {
        hdr.f1 = 0.0;
    }
    hdr.d1     = mod.dt*rec.skipdt;
    hdr.d2     = (rec.x[1]-rec.x[0])*mod.dx;
    hdr.f2     = mod.x0+rec.x[0]*mod.dx;

#ifdef MPI
    if (!opened) {
#endif 
    if (rec.type.vx)  fpvx  = fileOpen(filename, "_rvx", append);
    if (rec.type.vz)  fpvz  = fileOpen(filename, "_rvz", append);
    if (rec.type.p)   fpp   = fileOpen(filename, "_rp", append);
    if (rec.type.txx) fptxx = fileOpen(filename, "_rtxx", append);
    if (rec.type.tzz) fptzz = fileOpen(filename, "_rtzz", append);
    if (rec.type.txz) fptxz = fileOpen(filename, "_rtxz", append);
    if (rec.type.dxvx) fpdxvx = fileOpen(filename, "_rdxvx", append);
    if (rec.type.dzvz) fpdzvz = fileOpen(filename, "_rdzvz", append);
    if (rec.type.pp)  fppp  = fileOpen(filename, "_rpp", append);
    if (rec.type.ss)  fpss  = fileOpen(filename, "_rss", append);
    if (rec.type.q)  fpq  = fileOpen(filename, "_rq", append);
    if (rec.type.ud && (mod.ischeme==1 || mod.ischeme==2) )  {
        fpup   = fileOpen(filename, "_ru", append);
        fpdown = fileOpen(filename, "_rd", append);
    }
#ifdef MPI
    opened=1;
    }
#endif 

    /* decomposed wavefield */
    if (rec.type.ud && (mod.ischeme==1 || mod.ischeme==2) )  {
        ntfft = optncr(nsam);
        nfreq = ntfft/2+1;
        fmin = 0.0;
        fmax = wav.fmax;
        nkx = optncc(2*mod.nax);
        ibndx = mod.ioPx;
        if (bnd.lef==4 || bnd.lef==2) ibndx += bnd.ntap;
        cp  = rec.cp;
        rho = rec.rho;
		if (rec.type.ud==1) vznorm=0;
		else if (rec.type.ud==2) vznorm=1;
		else if (rec.type.ud==3) vznorm=-1;
        if (verbose) vmess("Decomposition array at z=%.2f with cp=%.2f rho=%.2f", rec.zr[0]+mod.z0, cp, rho);
        rec_up  = (float *)calloc(ntfft*nkx,sizeof(float));
        rec_down= (float *)calloc(ntfft*nkx,sizeof(float));
        crec_vz = (complex *)malloc(nfreq*nkx*sizeof(complex));
        crec_p  = (complex *)malloc(nfreq*nkx*sizeof(complex));
        crec_up = (complex *)malloc(nfreq*nkx*sizeof(complex));
        crec_dw = (complex *)malloc(nfreq*nkx*sizeof(complex));

        rec_vze = rec_up;
        rec_pe  = rec_down;
        /* copy input data into extended arrays with padded zeroes */
        for (ix=0; ix<mod.nax; ix++) {
            memcpy(&rec_vze[ix*ntfft],&rec_udvz[ix*rec.nt],nsam*sizeof(float));
            memcpy(&rec_pe[ix*ntfft], &rec_udp[ix*rec.nt], nsam*sizeof(float));
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
            for (it=0; it<rec.nt; it++) {
                rec_up[irec*rec.nt+it]   = rec_up[ix*ntfft+it];
                rec_down[irec*rec.nt+it] = rec_down[ix*ntfft+it];
            }
        }
        free(crec_vz);
        free(crec_p);
        free(crec_up);
        free(crec_dw);
    }
    if (rec.type.ud && (mod.ischeme==3 || mod.ischeme==4) )  {
// Not yet implemented
    }

    for (irec=0; irec<rec.n; irec++) {
        hdr.tracf  = irec+1;
        hdr.tracl  = ishot*rec.n+irec+1;
        hdr.gx     = 1000*(mod.x0+rec.x[irec]*mod.dx);
        hdr.offset = (rec.x[irec]-ixsrc)*mod.dx;
        hdr.cdp    = (hdr.gx + hdr.sx)/2;
        hdr.gelev  = (int)(-1000*(mod.z0+rec.z[irec]*mod.dz));
        offset     = (ishot*rec.n+irec)*(TRCBYTES+rec.nt*sizeof(float));

        if (rec.type.vx) {
            traceWrite( &hdr, &rec_vx[irec*rec.nt], nsam, offset, fpvx) ;
        }
        if (rec.type.vz) {
            traceWrite( &hdr, &rec_vz[irec*rec.nt], nsam, offset, fpvz) ;
        }
        if (rec.type.p) {
            traceWrite( &hdr, &rec_p[irec*rec.nt], nsam, offset, fpp) ;
        }
        if (rec.type.txx) {
            traceWrite( &hdr, &rec_txx[irec*rec.nt], nsam, offset, fptxx) ;
        }
        if (rec.type.tzz) {
            traceWrite( &hdr, &rec_tzz[irec*rec.nt], nsam, offset, fptzz) ;
        }
        if (rec.type.txz) {
            traceWrite( &hdr, &rec_txz[irec*rec.nt], nsam, offset, fptxz) ;
        }
        if (rec.type.dxvx) {
            traceWrite( &hdr, &rec_dxvx[irec*rec.nt], nsam, offset, fpdxvx) ;
        }
        if (rec.type.dzvz) {
            traceWrite( &hdr, &rec_dzvz[irec*rec.nt], nsam, offset, fpdzvz) ;
        }
        if (rec.type.pp) {
            traceWrite( &hdr, &rec_pp[irec*rec.nt], nsam, offset, fppp) ;
        }
        if (rec.type.ss) {
            traceWrite( &hdr, &rec_ss[irec*rec.nt], nsam, offset, fpss) ;
        }
        if (rec.type.q) {
            traceWrite( &hdr, &rec_q[irec*rec.nt], nsam, offset, fpq) ;
        }
        if (rec.type.ud && mod.ischeme==1)  {
            traceWrite( &hdr, &rec_up[irec*rec.nt], nsam, offset, fpup) ;
            traceWrite( &hdr, &rec_down[irec*rec.nt], nsam, offset, fpdown) ;
        }
    }

//#ifdef MPI
//	    fprintf(stderr,"PE %d has writen to file %s\n", pe, filename);
//	    fflush(stderr);
//#endif

#ifndef MPI
    if (rec.type.vx) fileClose(fpvx);
    if (rec.type.vz) fileClose(fpvz);
    if (rec.type.p) fileClose(fpp);
    if (rec.type.txx) fileClose(fptxx);
    if (rec.type.tzz) fileClose(fptzz);
    if (rec.type.txz) fileClose(fptxz);
    if (rec.type.dxvx) fileClose(fpdxvx);
    if (rec.type.dzvz) fileClose(fpdzvz);
    if (rec.type.pp) fileClose(fppp);
    if (rec.type.ss) fileClose(fpss);
    if (rec.type.ud) {
        fileClose(fpup);
        fileClose(fpdown);
        free(rec_up);
        free(rec_down);
    }
#else
    if (rec.type.ud) {
        free(rec_up);
        free(rec_down);
    }
#endif

    return 0;
}

#ifdef MPI
int closeRec(recPar rec) 
{
    if (rec.type.vx) fileClose(fpvx);
    if (rec.type.vz) fileClose(fpvz);
    if (rec.type.p) fileClose(fpp);
    if (rec.type.txx) fileClose(fptxx);
    if (rec.type.tzz) fileClose(fptzz);
    if (rec.type.txz) fileClose(fptxz);
    if (rec.type.dxvx) fileClose(fpdxvx);
    if (rec.type.dzvz) fileClose(fpdzvz);
    if (rec.type.pp) fileClose(fppp);
    if (rec.type.ss) fileClose(fpss);
    if (rec.type.ud) {
        fileClose(fpup);
        fileClose(fpdown);
    }
    return 0;
}
#endif
