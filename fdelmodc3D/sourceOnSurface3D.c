#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc3D.h"
#define ISODD(n) ((n) & 01)
static float *src1x, *src1y, *src1z, *src2x, *src2y, *src2z;
static long first=1;
void vmess(char *fmt, ...);

long allocStoreSourceOnSurface3D(srcPar src)
{
    /* allocated 2x size for dipole Potential sources */
    src1x = (float *)calloc(2*src.n, sizeof(float));
    src1y = (float *)calloc(2*src.n, sizeof(float));
    src1z = (float *)calloc(2*src.n, sizeof(float));
    src2x = (float *)calloc(2*src.n, sizeof(float));
    src2y = (float *)calloc(2*src.n, sizeof(float));
    src2z = (float *)calloc(2*src.n, sizeof(float));
    first = 0;
    return 0;
}

long freeStoreSourceOnSurface3D(void)
{
    free(src1x);
    free(src1y);
    free(src1z);
    free(src2x);
    free(src2y);
    free(src2z);
    first = 1;
    return 0;
}

long storeSourceOnSurface3D(modPar mod, srcPar src, bndPar bnd,
    long ixsrc, long iysrc, long izsrc, float *vx, float *vy, float *vz,
    float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz, long verbose)
{
/**********************************************************************

   AUTHOR:
           Jan Thorbecke (janth@xs4all.nl)
           The Netherlands 

***********************************************************************/

	long   ixs, iys, izs, isrc, is0;
    long   ibndz, ibndy, ibndx, store;
	long   nx, ny, nz, n1, n2;

	nx  = mod.nx;
	ny  = mod.ny;
	nz  = mod.nz;
	n1  = mod.naz;
    n2  = mod.nax;

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
        if (bnd.fro==4 || bnd.fro==2) ibndy += bnd.ntap;
	}
	else {	
    	ibndz = mod.ioPz;
    	ibndy = mod.ioPy;
    	ibndx = mod.ioPx;
    	if (bnd.lef==4 || bnd.lef==2) ibndx += bnd.ntap;
    	if (bnd.fro==4 || bnd.fro==2) ibndy += bnd.ntap;
	}

/* check if there are sources placed on the boundaries */
	is0 = -1*floor((src.n-1)/2);
#pragma omp	for private (isrc, ixs, iys, izs, store) 
	for (isrc=0; isrc<src.n; isrc++) {
		/* calculate the source position */
		if (src.random || src.multiwav) {
			ixs = src.x[isrc] + ibndx;
			iys = src.y[isrc] + ibndy;
			izs = src.z[isrc] + ibndz;
		}
		else if (src.plane) {/* plane wave sources */
            ixs = ixsrc + ibndx + src.x[isrc];
            iys = iysrc + ibndy + src.y[isrc];
            izs = izsrc + ibndz + src.z[isrc];
		}
		else { /* point sources */
            ixs = ixsrc + ibndx + isrc;
            iys = iysrc + ibndy + isrc;
            izs = izsrc + ibndz;
		}

/* check if there are source positions placed on the boundaries. 
 * In that case store them and reapply them after the boundary
 * conditions have been set */

        store=0;
		if ( (ixs <= ibndx+1)  && ISODD(bnd.lef)) store=1;
		if ( (ixs >= nx+ibndx) && ISODD(bnd.rig)) store=1;
		if ( (iys <= ibndy+1)  && ISODD(bnd.fro)) store=1;
		if ( (iys >= ny+ibndy) && ISODD(bnd.bac)) store=1;
		if ( (izs <= ibndz+1)  && ISODD(bnd.top)) store=1;
		if ( (izs >= nz+ibndz) && ISODD(bnd.bot)) store=1;

		if (mod.ischeme <= 2) { /* Acoustic scheme */
            
            if (store) {
                if (verbose>=5) vmess("source at x=%li y=%li z=%li stored before free surface", ixs, iys, izs);

                /* Compressional source */
                if (src.type == 1) {
                
                    if (src.orient==1) { /* monopole */
                        src1z[isrc] = tzz[iys*n1*n2+ixs*n1+izs];
                    }
                    else if (src.orient==2) { /* dipole +/- */
                        src1z[isrc] = tzz[iys*n1*n2+ixs*n1+izs];
                        src2z[isrc] = tzz[iys*n1*n2+ixs*n1+izs+1];
                    }
                    else if (src.orient==3) { /* dipole - + */
                        src1z[isrc] = tzz[iys*n1*n2+ixs*n1+izs];
                        src2z[isrc] = tzz[iys*n1*n2+(ixs-1)*n1+izs];
                    }
                    else if (src.orient==4) { /* dipole +/0/- */
                        if (izs > ibndz) 
                            src1z[isrc] = tzz[iys*n1*n2+ixs*n1+izs-1];
                        if (izs < mod.nz+ibndz-1) 
                            src2z[isrc] = tzz[iys*n1*n2+ixs*n1+izs+1];
                    }
                    else if (src.orient==5) { /* dipole + - */
                        src1z[isrc] = tzz[iys*n1*n2+ixs*n1+izs];
                        src2z[isrc] = tzz[iys*n1*n2+(ixs+1)*n1+izs];
                    }
                }
                else if (src.type==6) {
                    src1x[isrc] = vx[iys*n1*n2+ixs*n1+izs];
                }
                else if (src.type==7) {
                    src1z[isrc] = vz[iys*n1*n2+ixs*n1+izs];
                }
                
            }
        }
        else { /* Elastic scheme */

          	if (store) {
                if (verbose>=5) vmess("source at x=%li y=%li z=%li stored before free surface", ixs, iys, izs);

              	if (src.type==1) {
                    if (src.orient==1) { /* monopole */
                        src1x[isrc] = txx[iys*n1*n2+ixs*n1+izs];
                        src1y[isrc] = tyy[iys*n1*n2+ixs*n1+izs];
                        src1z[isrc] = tzz[iys*n1*n2+ixs*n1+izs];
                    }
                    else if (src.orient==2) { /* dipole +/- */
                        src1x[isrc] = txx[iys*n1*n2+ixs*n1+izs];
                        src1y[isrc] = tyy[iys*n1*n2+ixs*n1+izs];
                        src1z[isrc] = tzz[iys*n1*n2+ixs*n1+izs];
                        src2x[isrc] = txx[iys*n1*n2+ixs*n1+izs+1];
                        src2y[isrc] = tyy[iys*n1*n2+ixs*n1+izs+1];
                        src2z[isrc] = tzz[iys*n1*n2+ixs*n1+izs+1];
                    }
                    else if (src.orient==3) { /* dipole - + */
                        src1x[isrc] = txx[iys*n1*n2+ixs*n1+izs];
                        src1y[isrc] = tyy[iys*n1*n2+ixs*n1+izs];
                        src1z[isrc] = tzz[iys*n1*n2+ixs*n1+izs];
                        src2x[isrc] = txx[iys*n1*n2+(ixs-1)*n1+izs];
                        src2y[isrc] = tyy[iys*n1*n2+(ixs-1)*n1+izs];
                        src2z[isrc] = tzz[iys*n1*n2+(ixs-1)*n1+izs];
                    }
                    else if (src.orient==4) { /* dipole +/0/- */
                        if (izs > ibndz) {
                            src1x[isrc] = txx[iys*n1*n2+ixs*n1+izs-1];
                            src1y[isrc] = tyy[iys*n1*n2+ixs*n1+izs-1];
                            src1z[isrc] = tzz[iys*n1*n2+ixs*n1+izs-1];
                        }
                        if (izs < mod.nz+ibndz-1) {
                            src1x[isrc] = txx[iys*n1*n2+ixs*n1+izs+1];
                            src1y[isrc] = tyy[iys*n1*n2+ixs*n1+izs+1];
                            src1z[isrc] = tzz[iys*n1*n2+ixs*n1+izs+1];
                        }
                    }
                    else if (src.orient==5) { /* dipole + - */
                        src1x[isrc] = txx[iys*n1*n2+ixs*n1+izs];
                        src1y[isrc] = tyy[iys*n1*n2+ixs*n1+izs];
                        src1z[isrc] = tzz[iys*n1*n2+ixs*n1+izs];
                        src2x[isrc] = txx[iys*n1*n2+(ixs+1)*n1+izs];
                        src2y[isrc] = tyy[iys*n1*n2+(ixs+1)*n1+izs];
                        src2z[isrc] = tzz[iys*n1*n2+(ixs+1)*n1+izs];
                    }
              	}
              	else if (src.type==2) {
                    
                    /* Txz source */
                    if ((izs == ibndz) && bnd.top==1) {
                        src1x[isrc] = txz[iys*n1*n2+(ixs-1)*n1+izs-1];
                        src2x[isrc] = txz[iys*n1*n2+ixs*n1+izs-1];
                    }
                    else {
                        src1x[isrc] = txz[iys*n1*n2+ixs*n1+izs];
                    }
                    /* possible dipole orientations for a txz source */
                    if (src.orient == 2) { /* dipole +/- */
                        src2x[isrc] = txz[iys*n1*n2+ixs*n1+izs+1];
                    }
                    else if (src.orient == 3) { /* dipole - + */
                        src2x[isrc] = txz[iys*n1*n2+(ixs-1)*n1+izs];
                    }
                    else if (src.orient == 4) { /*  dipole +/O/- */
                        /* correction: subtrace previous value to prevent z-1 values. */
                        src1x[isrc] = txz[iys*n1*n2+ixs*n1+izs];
                        src2x[isrc] = txz[iys*n1*n2+ixs*n1+izs+1];
                    }
                    else if (src.orient == 5) { /* dipole + - */
                        src2x[isrc] = txz[iys*n1*n2+(ixs+1)*n1+izs];
                    }

              	}
               	else if (src.type==3) {
                   	src1z[isrc] = tzz[iys*n1*n2+ixs*n1+izs];
               	}
               	else if (src.type==4) {
                   	src1x[isrc] = txx[iys*n1*n2+ixs*n1+izs];
               	}
               	else if (src.type==5) {
                                        
                    src1x[isrc] = vx[iys*n1*n2+ixs*n1+izs];
                    src1y[isrc] = vy[iys*n1*n2+ixs*n1+izs];
                    src1z[isrc] = vz[iys*n1*n2+ixs*n1+izs];
                    src2x[isrc] = vx[iys*n1*n2+ixs*n1+izs-1];
                    src2y[isrc] = vy[(iys-1)*n1*n2+ixs*n1+izs];
                    src2z[isrc] = vz[iys*n1*n2+(ixs-1)*n1+izs];

                    /* determine second position of dipole */
                    if (src.orient == 2) { /* dipole +/- vertical */
                        izs += 1;
                        src1x[isrc+src.n] = vx[iys*n1*n2+ixs*n1+izs];
                        src1y[isrc+src.n] = vy[iys*n1*n2+ixs*n1+izs];
                        src1z[isrc+src.n] = vz[iys*n1*n2+ixs*n1+izs];
                        src2x[isrc+src.n] = vx[iys*n1*n2+ixs*n1+izs-1];
                        src2y[isrc+src.n] = vz[(iys-1)*n1*n2+ixs*n1+izs];
                        src2z[isrc+src.n] = vz[iys*n1*n2+(ixs-1)*n1+izs];
                    }
                    else if (src.orient == 3) { /* dipole - + horizontal */
                        ixs += 1;
                        src1x[isrc+src.n] = vx[iys*n1*n2+ixs*n1+izs];
                        src1y[isrc+src.n] = vy[iys*n1*n2+ixs*n1+izs];
                        src1z[isrc+src.n] = vz[iys*n1*n2+ixs*n1+izs];
                        src2x[isrc+src.n] = vx[iys*n1*n2+ixs*n1+izs-1];
                        src2y[isrc+src.n] = vz[(iys-1)*n1*n2+ixs*n1+izs];
                        src2z[isrc+src.n] = vz[iys*n1*n2+(ixs-1)*n1+izs];
                    }

				}
               	else if (src.type==6) {
                   	src1x[isrc] = vx[iys*n1*n2+ixs*n1+izs];
               	}
               	else if (src.type==7) {
                   	src1z[isrc] = vz[iys*n1*n2+ixs*n1+izs];
               	}
               	else if (src.type==8) {
                    
                    src1x[isrc] = vx[iys*n1*n2+(ixs+1)*n1+izs];
                    src1y[isrc] = vy[(iys+1)*n1*n2+ixs*n1+izs];
                    src1z[isrc] = vz[iys*n1*n2+ixs*n1+izs+1];
                    src2x[isrc] = vx[iys*n1*n2+ixs*n1+izs];
                    src2y[isrc] = vy[iys*n1*n2+ixs*n1+izs];
                    src2z[isrc] = vz[iys*n1*n2+ixs*n1+izs];
                    
                    /* determine second position of dipole */
                    if (src.orient == 2) { /* dipole +/- vertical */
                        izs += 1;
                        src1x[isrc+src.n] = vx[iys*n1*n2+(ixs+1)*n1+izs];
                        src1y[isrc]       = vy[(iys+1)*n1*n2+ixs*n1+izs];
                        src1z[isrc+src.n] = vz[iys*n1*n2+ixs*n1+izs+1];
                        src2x[isrc+src.n] = vx[iys*n1*n2+ixs*n1+izs];
                        src2y[isrc]       = vy[iys*n1*n2+ixs*n1+izs];
                        src2z[isrc+src.n] = vz[iys*n1*n2+ixs*n1+izs];
                    }
                    else if (src.orient == 3) { /* dipole - + horizontal */
                        ixs += 1;
                        src1x[isrc+src.n] = vx[iys*n1*n2+(ixs+1)*n1+izs];
                        src1y[isrc]       = vy[(iys+1)*n1*n2+ixs*n1+izs];
                        src1z[isrc+src.n] = vz[iys*n1*n2+ixs*n1+izs+1];
                        src2x[isrc+src.n] = vx[iys*n1*n2+ixs*n1+izs];
                        src2y[isrc]       = vy[iys*n1*n2+ixs*n1+izs];
                        src2z[isrc+src.n] = vz[iys*n1*n2+ixs*n1+izs];
                    }

               	} /* end of source.type */
           	}
		}
    }
    
    return 0;
}

    
    
long reStoreSourceOnSurface3D(modPar mod, srcPar src, bndPar bnd,
    long ixsrc, long iysrc, long izsrc, float *vx, float *vy, float *vz,
    float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz, long verbose)
{
    /**********************************************************************
     
     AUTHOR:
     Jan Thorbecke (janth@xs4all.nl)
     The Netherlands 
     
     ***********************************************************************/
    
	long   ixs, iys, izs, isrc, is0;
    long   ibndz, ibndy, ibndx, store;
	long   nx, ny, nz, n1, n2;
    
	nx  = mod.nx;
    ny  = mod.ny;
	nz  = mod.nz;
	n1  = mod.naz;
    n2  = mod.nax;
    
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
        if (bnd.fro==4 || bnd.fro==2) ibndy += bnd.ntap;
	}
	else {	
    	ibndz = mod.ioPz;
    	ibndy = mod.ioPy;
    	ibndx = mod.ioPx;
    	if (bnd.lef==4 || bnd.lef==2) ibndx += bnd.ntap;
        if (bnd.fro==4 || bnd.fro==2) ibndy += bnd.ntap;
	}

	/* restore source positions on the edge */
	is0 = -1*floor((src.n-1)/2);
#pragma omp	for private (isrc, ixs, iys, izs, store) 
	for (isrc=0; isrc<src.n; isrc++) {
		/* calculate the source position */
		if (src.random || src.multiwav) {
			ixs = src.x[isrc] + ibndx;
			iys = src.y[isrc] + ibndy;
			izs = src.z[isrc] + ibndz;
		}
		else { /* plane wave and point sources */
			ixs = ixsrc + ibndx + is0 + isrc;
			iys = iysrc + ibndy + is0 + isrc;
			izs = izsrc + ibndz;
		}
        
        store=0;
		if ( (ixs <= ibndx+1)  && ISODD(bnd.lef)) store=1;
		if ( (ixs >= nx+ibndx) && ISODD(bnd.rig)) store=1;
		if ( (iys <= ibndy+1)  && ISODD(bnd.fro)) store=1;
		if ( (iys >= ny+ibndy) && ISODD(bnd.bac)) store=1;
		if ( (izs <= ibndz+1)  && ISODD(bnd.top)) store=1;
		if ( (izs >= nz+ibndz) && ISODD(bnd.bot)) store=1;
        
		if (mod.ischeme <= 2) { /* Acoustic scheme */
            
            if (store) {
                if (verbose>=5) vmess("source at x=%li y=%li z=%li restored at free surface", ixs, iys, izs);

                /* Compressional source */
                if (src.type == 1) {
                    
                    if (src.orient==1) { /* monopole */
                        tzz[iys*n1*n2+ixs*n1+izs]= src1z[isrc];
                    }
                    else if (src.orient==2) { /* dipole +/- */
                        tzz[iys*n1*n2+ixs*n1+izs] = src1z[isrc];
                        tzz[iys*n1*n2+ixs*n1+izs+1] = src2z[isrc];
                    }
                    else if (src.orient==3) { /* dipole - + */
                        tzz[iys*n1*n2+ixs*n1+izs] = src1z[isrc];
                        tzz[iys*n1*n2+(ixs-1)*n1+izs] = src2z[isrc];
                    }
                    else if (src.orient==4) { /* dipole +/0/- */
                        if (izs > ibndz) 
                            tzz[iys*n1*n2+ixs*n1+izs-1] = src1z[isrc];
                        if (izs < mod.nz+ibndz-1) 
                            tzz[iys*n1*n2+ixs*n1+izs+1] = src2z[isrc];
                    }
                    else if (src.orient==5) { /* dipole + - */
                        tzz[iys*n1*n2+ixs*n1+izs] = src1z[isrc];
                        tzz[iys*n1*n2+(ixs+1)*n1+izs] = src2z[isrc];
                    }
                }
                else if (src.type==6) {
                    vx[iys*n1*n2+ixs*n1+izs] = src1x[isrc];
                }
                else if (src.type==7) {
                    vz[iys*n1*n2+ixs*n1+izs] = src1z[isrc];
                }
                
            }
            
        }
        else { /* Elastic scheme */
            
          	if (store) {
                if (verbose>=5) vmess("source at x=%li y=%li z=%li restored at free surface", ixs, iys, izs);

              	if (src.type==1) {
                    if (src.orient==1) { /* monopole */
                        txx[iys*n1*n2+ixs*n1+izs] = src1x[isrc];
                        tyy[iys*n1*n2+ixs*n1+izs] = src1y[isrc];
                        tzz[iys*n1*n2+ixs*n1+izs] = src1z[isrc];
                    }
                    else if (src.orient==2) { /* dipole +/- */
                        txx[iys*n1*n2+ixs*n1+izs] = src1x[isrc];
                        tyy[iys*n1*n2+ixs*n1+izs] = src1y[isrc];
                        tzz[iys*n1*n2+ixs*n1+izs] = src1z[isrc];
                        txx[iys*n1*n2+ixs*n1+izs+1] = src2x[isrc];
                        tyy[iys*n1*n2+ixs*n1+izs+1] = src2y[isrc];
                        tzz[iys*n1*n2+ixs*n1+izs+1] = src2z[isrc];
                    }
                    else if (src.orient==3) { /* dipole - + */
                        txx[iys*n1*n2+ixs*n1+izs] = src1x[isrc];
                        tyy[iys*n1*n2+ixs*n1+izs] = src1y[isrc];
                        tzz[iys*n1*n2+ixs*n1+izs] = src1z[isrc];
                        txx[iys*n1*n2+(ixs-1)*n1+izs] = src2x[isrc];
                        tyy[iys*n1*n2+(ixs-1)*n1+izs] = src2y[isrc];
                        tzz[iys*n1*n2+(ixs-1)*n1+izs] = src2z[isrc];
                    }
                    else if (src.orient==4) { /* dipole +/0/- */
                        if (izs > ibndz) {
                            txx[iys*n1*n2+ixs*n1+izs-1] = src1x[isrc];
                            tyy[iys*n1*n2+ixs*n1+izs-1] = src1y[isrc];
                            tzz[iys*n1*n2+ixs*n1+izs-1] = src1z[isrc];
                        }
                        if (izs < mod.nz+ibndz-1) {
                            txx[iys*n1*n2+ixs*n1+izs+1] = src1x[isrc];
                            tyy[iys*n1*n2+ixs*n1+izs+1] = src1y[isrc];
                            tzz[iys*n1*n2+ixs*n1+izs+1] = src1z[isrc];
                        }
                    }
                    else if (src.orient==5) { /* dipole + - */
                        txx[iys*n1*n2+ixs*n1+izs] = src1x[isrc];
                        tyy[iys*n1*n2+ixs*n1+izs] = src1y[isrc];
                        tzz[iys*n1*n2+ixs*n1+izs] = src1z[isrc];
                        txx[iys*n1*n2+(ixs+1)*n1+izs] = src2x[isrc];
                        tyy[iys*n1*n2+(ixs+1)*n1+izs] = src2y[isrc];
                        tzz[iys*n1*n2+(ixs+1)*n1+izs] = src2z[isrc];
                    }
              	}
              	else if (src.type==2) {
                    
                    /* Txz source */
                    if ((izs == ibndz) && bnd.top==1) {
                        txz[iys*n1*n2+(ixs-1)*n1+izs-1] = src1x[isrc];
                        txz[iys*n1*n2+ixs*n1+izs-1] = src2x[isrc];
                    }
                    else {
                        txz[iys*n1*n2+ixs*n1+izs] = src1x[isrc];
                    }
                    /* possible dipole orientations for a txz source */
                    if (src.orient == 2) { /* dipole +/- */
                        txz[iys*n1*n2+ixs*n1+izs+1] = src2x[isrc];
                    }
                    else if (src.orient == 3) { /* dipole - + */
                        txz[iys*n1*n2+(ixs-1)*n1+izs] = src2x[isrc];
                    }
                    else if (src.orient == 4) { /*  dipole +/O/- */
                        /* correction: subtrace previous value to prevent z-1 values. */
                        txz[iys*n1*n2+ixs*n1+izs] = src1x[isrc];
                        txz[iys*n1*n2+ixs*n1+izs+1] = src2x[isrc];
                    }
                    else if (src.orient == 5) { /* dipole + - */
                        txz[iys*n1*n2+(ixs+1)*n1+izs] = src2x[isrc];
                    }
                    
              	}
               	else if (src.type==3) {
                   	tzz[iys*n1*n2+ixs*n1+izs] = src1z[isrc];
               	}
               	else if (src.type==4) {
                   	txx[iys*n1*n2+ixs*n1+izs] = src1x[isrc];
               	}
               	else if (src.type==5) {
                    
                    vx[iys*n1*n2+ixs*n1+izs]= src1x[isrc];
                    vy[iys*n1*n2+ixs*n1+izs] = src1y[isrc];
                    vz[iys*n1*n2+ixs*n1+izs] = src1z[isrc];
                    vx[iys*n1*n2+ixs*n1+izs-1] = src2x[isrc];
                    vy[(iys-1)*n1*n2+ixs*n1+izs] = src2y[isrc];
                    vz[iys*n1*n2+(ixs-1)*n1+izs] = src2z[isrc];
                    
                    /* determine second position of dipole */
                    if (src.orient == 2) { /* dipole +/- vertical */
                        izs += 1;
                        vx[iys*n1*n2+ixs*n1+izs] = src1x[isrc+src.n];
                        vy[iys*n1*n2+ixs*n1+izs] = src1y[isrc+src.n];
                        vz[iys*n1*n2+ixs*n1+izs] = src1z[isrc+src.n];
                        vx[iys*n1*n2+ixs*n1+izs-1] = src2x[isrc+src.n];
                        vy[(iys-1)*n1*n2+ixs*n1+izs] = src2y[isrc+src.n];
                        vz[iys*n1*n2+(ixs-1)*n1+izs] = src2z[isrc+src.n];
                    }
                    else if (src.orient == 3) { /* dipole - + horizontal */
                        ixs += 1;
                        vx[iys*n1*n2+ixs*n1+izs] = src1x[isrc+src.n];
                        vy[iys*n1*n2+ixs*n1+izs] = src1y[isrc+src.n];
                        vz[iys*n1*n2+ixs*n1+izs] = src1z[isrc+src.n];
                        vx[iys*n1*n2+ixs*n1+izs-1] = src2x[isrc+src.n];
                        vy[(iys-1)*n1*n2+ixs*n1+izs] = src2y[isrc+src.n];
                        vz[iys*n1*n2+(ixs-1)*n1+izs] = src2z[isrc+src.n];
                    }
                    
				}
               	else if (src.type==6) {
                   	vx[iys*n1*n2+ixs*n1+izs] = src1x[isrc];
               	}
               	else if (src.type==7) {
                   	vz[iys*n1*n2+ixs*n1+izs] = src1z[isrc];
               	}
               	else if (src.type==8) {
                    
                    vx[iys*n1*n2+(ixs+1)*n1+izs] = src1x[isrc];
                    vy[(iys+1)*n1*n2+ixs*n1+izs] = src1z[isrc];
                    vz[iys*n1*n2+ixs*n1+izs+1] = src1z[isrc];
                    vx[iys*n1*n2+ixs*n1+izs] = src2x[isrc];
                    vy[iys*n1*n2+ixs*n1+izs] = src2y[isrc];
                    vz[iys*n1*n2+ixs*n1+izs] = src2z[isrc];
                    
                    /* determine second position of dipole */
                    if (src.orient == 2) { /* dipole +/- vertical */
                        izs += 1;
                        vx[iys*n1*n2+(ixs+1)*n1+izs] = src1x[isrc+src.n];
                        vy[(iys+1)*n1*n2+ixs*n1+izs] = src1y[isrc+src.n];
                        vz[iys*n1*n2+ixs*n1+izs+1] = src1z[isrc+src.n];
                        vx[iys*n1*n2+ixs*n1+izs] = src2x[isrc+src.n];
                        vy[iys*n1*n2+ixs*n1+izs] = src2y[isrc+src.n];
                        vz[iys*n1*n2+ixs*n1+izs] = src2z[isrc+src.n];
                    }
                    else if (src.orient == 3) { /* dipole - + horizontal */
                        ixs += 1;
                        vx[iys*n1*n2+(ixs+1)*n1+izs] = src1x[isrc+src.n];
                        vy[(iys+1)*n1*n2+ixs*n1+izs] = src1y[isrc+src.n];
                        vz[iys*n1*n2+ixs*n1+izs+1] = src1z[isrc+src.n];
                        vx[iys*n1*n2+ixs*n1+izs] = src2x[isrc+src.n];
                        vy[iys*n1*n2+ixs*n1+izs] = src2y[isrc+src.n];
                        vz[iys*n1*n2+ixs*n1+izs] = src2z[isrc+src.n];
                    }
                    
               	}
           	}
		}
    }

    return 0;
}
