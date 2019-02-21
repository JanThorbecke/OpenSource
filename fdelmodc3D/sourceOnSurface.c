#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc.h"
#define ISODD(n) ((n) & 01)
static float *src1x, *src1z, *src2x, *src2z;
static int first=1;
void vmess(char *fmt, ...);

int allocStoreSourceOnSurface(srcPar src)
{
    /* allocated 2x size for dipole Potential sources */
    src1x = (float *)calloc(2*src.n, sizeof(float));
    src1z = (float *)calloc(2*src.n, sizeof(float));
    src2x = (float *)calloc(2*src.n, sizeof(float));
    src2z = (float *)calloc(2*src.n, sizeof(float));
    first = 0;
    return 0;
}

int freeStoreSourceOnSurface(void)
{
    free(src1x);
    free(src1z);
    free(src2x);
    free(src2z);
    first = 1;
    return 0;
}

int storeSourceOnSurface(modPar mod, srcPar src, bndPar bnd, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, int verbose)
{
/**********************************************************************

   AUTHOR:
           Jan Thorbecke (janth@xs4all.nl)
           The Netherlands 

***********************************************************************/

	int   ixs, izs, isrc, is0;
    int   ibndz, ibndx, store;
	int   nx, nz, n1;

	nx  = mod.nx;
	nz  = mod.nz;
	n1  = mod.naz;

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
	}
	else {	
    	ibndz = mod.ioPz;
    	ibndx = mod.ioPx;
    	if (bnd.lef==4 || bnd.lef==2) ibndx += bnd.ntap;
	}

/* check if there are sources placed on the boundaries */
	is0 = -1*floor((src.n-1)/2);
#pragma omp	for private (isrc, ixs, izs, store) 
	for (isrc=0; isrc<src.n; isrc++) {
		/* calculate the source position */
		if (src.random || src.multiwav) {
			ixs = src.x[isrc] + ibndx;
			izs = src.z[isrc] + ibndz;
		}
		else { /* plane wave and point sources */
			ixs = ixsrc + ibndx + is0 + isrc;
			izs = izsrc + ibndz;
		}

//        vmess("source at x=%d bounds at %d %d : %d %d ", ixs, ibndx+1, nx+ibndx, mod.ioXz, mod.ieXz);
//        vmess("source at z=%d bounds at %d %d : %d %d ", izs, ibndz+1, nz+ibndz, mod.ioXx, mod.ieXx);

/* check if there are source positions placed on the boundaries. 
 * In that case store them and reapply them after the boundary
 * conditions have been set */

        store=0;
		if ( (ixs <= ibndx+1)  && ISODD(bnd.lef)) store=1;
		if ( (ixs >= nx+ibndx) && ISODD(bnd.rig)) store=1;
		if ( (izs <= ibndz+1)  && ISODD(bnd.top)) store=1;
		if ( (izs >= nz+ibndz) && ISODD(bnd.bot)) store=1;

		if (mod.ischeme <= 2) { /* Acoustic scheme */
            
            if (store) {
                if (verbose>=5) vmess("source at x=%d z=%d stored before free surface", ixs, izs);

                /* Compressional source */
                if (src.type == 1) {
                
                    if (src.orient==1) { /* monopole */
                        src1z[isrc] = tzz[ixs*n1+izs];
                    }
                    else if (src.orient==2) { /* dipole +/- */
                        src1z[isrc] = tzz[ixs*n1+izs];
                        src2z[isrc] = tzz[ixs*n1+izs+1];
                    }
                    else if (src.orient==3) { /* dipole - + */
                        src1z[isrc] = tzz[ixs*n1+izs];
                        src2z[isrc] = tzz[(ixs-1)*n1+izs];
                    }
                    else if (src.orient==4) { /* dipole +/0/- */
                        if (izs > ibndz) 
                            src1z[isrc] = tzz[ixs*n1+izs-1];
                        if (izs < mod.nz+ibndz-1) 
                            src2z[isrc] = tzz[ixs*n1+izs+1];
                    }
                    else if (src.orient==5) { /* dipole + - */
                        src1z[isrc] = tzz[ixs*n1+izs];
                        src2z[isrc] = tzz[(ixs+1)*n1+izs];
                    }
                }
                else if (src.type==6) {
                    src1x[isrc] = vx[ixs*n1+izs];
                }
                else if (src.type==7) {
                    src1z[isrc] = vz[ixs*n1+izs];
                }
                
            }
        }
        else { /* Elastic scheme */

          	if (store) {
                if (verbose>=5) vmess("source at x=%d z=%d stored before free surface", ixs, izs);

              	if (src.type==1) {
                    if (src.orient==1) { /* monopole */
                        src1x[isrc] = txx[ixs*n1+izs];
                        src1z[isrc] = tzz[ixs*n1+izs];
                    }
                    else if (src.orient==2) { /* dipole +/- */
                        src1x[isrc] = txx[ixs*n1+izs];
                        src1z[isrc] = tzz[ixs*n1+izs];
                        src2x[isrc] = txx[ixs*n1+izs+1];
                        src2z[isrc] = tzz[ixs*n1+izs+1];
                    }
                    else if (src.orient==3) { /* dipole - + */
                        src1x[isrc] = txx[ixs*n1+izs];
                        src1z[isrc] = tzz[ixs*n1+izs];
                        src2x[isrc] = txx[(ixs-1)*n1+izs];
                        src2z[isrc] = tzz[(ixs-1)*n1+izs];
                    }
                    else if (src.orient==4) { /* dipole +/0/- */
                        if (izs > ibndz) {
                            src1x[isrc] = txx[ixs*n1+izs-1];
                            src1z[isrc] = tzz[ixs*n1+izs-1];
                        }
                        if (izs < mod.nz+ibndz-1) {
                            src1x[isrc] = txx[ixs*n1+izs+1];
                            src1z[isrc] = tzz[ixs*n1+izs+1];
                        }
                    }
                    else if (src.orient==5) { /* dipole + - */
                        src1x[isrc] = txx[ixs*n1+izs];
                        src1z[isrc] = tzz[ixs*n1+izs];
                        src2x[isrc] = txx[(ixs+1)*n1+izs];
                        src2z[isrc] = tzz[(ixs+1)*n1+izs];
                    }
              	}
              	else if (src.type==2) {
                    
                    /* Txz source */
                    if ((izs == ibndz) && bnd.top==1) {
                        src1x[isrc] = txz[(ixs-1)*n1+izs-1];
                        src2x[isrc] = txz[ixs*n1+izs-1];
                    }
                    else {
                        src1x[isrc] = txz[ixs*n1+izs];
                    }
                    /* possible dipole orientations for a txz source */
                    if (src.orient == 2) { /* dipole +/- */
                        src2x[isrc] = txz[ixs*n1+izs+1];
                    }
                    else if (src.orient == 3) { /* dipole - + */
                        src2x[isrc] = txz[(ixs-1)*n1+izs];
                    }
                    else if (src.orient == 4) { /*  dipole +/O/- */
                        /* correction: subtrace previous value to prevent z-1 values. */
                        src1x[isrc] = txz[ixs*n1+izs];
                        src2x[isrc] = txz[ixs*n1+izs+1];
                    }
                    else if (src.orient == 5) { /* dipole + - */
                        src2x[isrc] = txz[(ixs+1)*n1+izs];
                    }

              	}
               	else if (src.type==3) {
                   	src1z[isrc] = tzz[ixs*n1+izs];
               	}
               	else if (src.type==4) {
                   	src1x[isrc] = txx[ixs*n1+izs];
               	}
               	else if (src.type==5) {
                                        
                    src1x[isrc] = vx[ixs*n1+izs];
                    src1z[isrc] = vz[ixs*n1+izs];
                    src2x[isrc] = vx[ixs*n1+izs-1];
                    src2z[isrc] = vz[(ixs-1)*n1+izs];

                    /* determine second position of dipole */
                    if (src.orient == 2) { /* dipole +/- vertical */
                        izs += 1;
                        src1x[isrc+src.n] = vx[ixs*n1+izs];
                        src1z[isrc+src.n] = vz[ixs*n1+izs];
                        src2x[isrc+src.n] = vx[ixs*n1+izs-1];
                        src2z[isrc+src.n] = vz[(ixs-1)*n1+izs];
                    }
                    else if (src.orient == 3) { /* dipole - + horizontal */
                        ixs += 1;
                        src1x[isrc+src.n] = vx[ixs*n1+izs];
                        src1z[isrc+src.n] = vz[ixs*n1+izs];
                        src2x[isrc+src.n] = vx[ixs*n1+izs-1];
                        src2z[isrc+src.n] = vz[(ixs-1)*n1+izs];
                    }

				}
               	else if (src.type==6) {
                   	src1x[isrc] = vx[ixs*n1+izs];
               	}
               	else if (src.type==7) {
                   	src1z[isrc] = vz[ixs*n1+izs];
               	}
               	else if (src.type==8) {
                    
                    src1x[isrc] = vx[(ixs+1)*n1+izs];
                    src1z[isrc] = vz[ixs*n1+izs+1];
                    src2x[isrc] = vx[ixs*n1+izs];
                    src2z[isrc] = vz[ixs*n1+izs];
                    
                    /* determine second position of dipole */
                    if (src.orient == 2) { /* dipole +/- vertical */
                        izs += 1;
                        src1x[isrc+src.n] = vx[(ixs+1)*n1+izs];
                        src1z[isrc+src.n] = vz[ixs*n1+izs+1];
                        src2x[isrc+src.n] = vx[ixs*n1+izs];
                        src2z[isrc+src.n] = vz[ixs*n1+izs];
                    }
                    else if (src.orient == 3) { /* dipole - + horizontal */
                        ixs += 1;
                        src1x[isrc+src.n] = vx[(ixs+1)*n1+izs];
                        src1z[isrc+src.n] = vz[ixs*n1+izs+1];
                        src2x[isrc+src.n] = vx[ixs*n1+izs];
                        src2z[isrc+src.n] = vz[ixs*n1+izs];
                    }

               	} /* end of source.type */
           	}
		}
    }
    
    return 0;
}

    
    
int reStoreSourceOnSurface(modPar mod, srcPar src, bndPar bnd, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, int verbose)
{
    /**********************************************************************
     
     AUTHOR:
     Jan Thorbecke (janth@xs4all.nl)
     The Netherlands 
     
     ***********************************************************************/
    
	int   ixs, izs, isrc, is0;
    int   ibndz, ibndx, store;
	int   nx, nz, n1;
    
	nx  = mod.nx;
	nz  = mod.nz;
	n1  = mod.naz;
    
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
	}
	else {	
    	ibndz = mod.ioPz;
    	ibndx = mod.ioPx;
    	if (bnd.lef==4 || bnd.lef==2) ibndx += bnd.ntap;
	}

	/* restore source positions on the edge */
	is0 = -1*floor((src.n-1)/2);
#pragma omp	for private (isrc, ixs, izs, store) 
	for (isrc=0; isrc<src.n; isrc++) {
		/* calculate the source position */
		if (src.random || src.multiwav) {
			ixs = src.x[isrc] + ibndx;
			izs = src.z[isrc] + ibndz;
		}
		else { /* plane wave and point sources */
			ixs = ixsrc + ibndx + is0 + isrc;
			izs = izsrc + ibndz;
		}
        
        store=0;
		if ( (ixs <= ibndx+1)  && ISODD(bnd.lef)) store=1;
		if ( (ixs >= nx+ibndx) && ISODD(bnd.rig)) store=1;
		if ( (izs <= ibndz+1)  && ISODD(bnd.top)) store=1;
		if ( (izs >= nz+ibndz) && ISODD(bnd.bot)) store=1;
        
		if (mod.ischeme <= 2) { /* Acoustic scheme */
            
            if (store) {
                if (verbose>=5) vmess("source at x=%d z=%d restored at free surface", ixs, izs);

                /* Compressional source */
                if (src.type == 1) {
                    
                    if (src.orient==1) { /* monopole */
                        tzz[ixs*n1+izs]= src1z[isrc];
                    }
                    else if (src.orient==2) { /* dipole +/- */
                        tzz[ixs*n1+izs] = src1z[isrc];
                        tzz[ixs*n1+izs+1] = src2z[isrc];
                    }
                    else if (src.orient==3) { /* dipole - + */
                        tzz[ixs*n1+izs] = src1z[isrc];
                        tzz[(ixs-1)*n1+izs] = src2z[isrc];
                    }
                    else if (src.orient==4) { /* dipole +/0/- */
                        if (izs > ibndz) 
                            tzz[ixs*n1+izs-1] = src1z[isrc];
                        if (izs < mod.nz+ibndz-1) 
                            tzz[ixs*n1+izs+1] = src2z[isrc];
                    }
                    else if (src.orient==5) { /* dipole + - */
                        tzz[ixs*n1+izs] = src1z[isrc];
                        tzz[(ixs+1)*n1+izs] = src2z[isrc];
                    }
                }
                else if (src.type==6) {
                    vx[ixs*n1+izs] = src1x[isrc];
                }
                else if (src.type==7) {
                    vz[ixs*n1+izs] = src1z[isrc];
                }
                
            }
            
        }
        else { /* Elastic scheme */
            
          	if (store) {
                if (verbose>=5) vmess("source at x=%d z=%d restored at free surface", ixs, izs);

              	if (src.type==1) {
                    if (src.orient==1) { /* monopole */
                        txx[ixs*n1+izs] = src1x[isrc];
                        tzz[ixs*n1+izs] = src1z[isrc];
                    }
                    else if (src.orient==2) { /* dipole +/- */
                        txx[ixs*n1+izs] = src1x[isrc];
                        tzz[ixs*n1+izs] = src1z[isrc];
                        txx[ixs*n1+izs+1] = src2x[isrc];
                        tzz[ixs*n1+izs+1] = src2z[isrc];
                    }
                    else if (src.orient==3) { /* dipole - + */
                        txx[ixs*n1+izs] = src1x[isrc];
                        tzz[ixs*n1+izs] = src1z[isrc];
                        txx[(ixs-1)*n1+izs] = src2x[isrc];
                        tzz[(ixs-1)*n1+izs] = src2z[isrc];
                    }
                    else if (src.orient==4) { /* dipole +/0/- */
                        if (izs > ibndz) {
                            txx[ixs*n1+izs-1] = src1x[isrc];
                            tzz[ixs*n1+izs-1] = src1z[isrc];
                        }
                        if (izs < mod.nz+ibndz-1) {
                            txx[ixs*n1+izs+1] = src1x[isrc];
                            tzz[ixs*n1+izs+1] = src1z[isrc];
                        }
                    }
                    else if (src.orient==5) { /* dipole + - */
                        txx[ixs*n1+izs] = src1x[isrc];
                        tzz[ixs*n1+izs] = src1z[isrc];
                        txx[(ixs+1)*n1+izs] = src2x[isrc];
                        tzz[(ixs+1)*n1+izs] = src2z[isrc];
                    }
              	}
              	else if (src.type==2) {
                    
                    /* Txz source */
                    if ((izs == ibndz) && bnd.top==1) {
                        txz[(ixs-1)*n1+izs-1] = src1x[isrc];
                        txz[ixs*n1+izs-1] = src2x[isrc];
                    }
                    else {
                        txz[ixs*n1+izs] = src1x[isrc];
                    }
                    /* possible dipole orientations for a txz source */
                    if (src.orient == 2) { /* dipole +/- */
                        txz[ixs*n1+izs+1] = src2x[isrc];
                    }
                    else if (src.orient == 3) { /* dipole - + */
                        txz[(ixs-1)*n1+izs] = src2x[isrc];
                    }
                    else if (src.orient == 4) { /*  dipole +/O/- */
                        /* correction: subtrace previous value to prevent z-1 values. */
                        txz[ixs*n1+izs] = src1x[isrc];
                        txz[ixs*n1+izs+1] = src2x[isrc];
                    }
                    else if (src.orient == 5) { /* dipole + - */
                        txz[(ixs+1)*n1+izs] = src2x[isrc];
                    }
                    
              	}
               	else if (src.type==3) {
                   	tzz[ixs*n1+izs] = src1z[isrc];
               	}
               	else if (src.type==4) {
                   	txx[ixs*n1+izs] = src1x[isrc];
               	}
               	else if (src.type==5) {
                    
                    vx[ixs*n1+izs]= src1x[isrc];
                    vz[ixs*n1+izs] = src1z[isrc];
                    vx[ixs*n1+izs-1] = src2x[isrc];
                    vz[(ixs-1)*n1+izs] = src2z[isrc];
                    
                    /* determine second position of dipole */
                    if (src.orient == 2) { /* dipole +/- vertical */
                        izs += 1;
                        vx[ixs*n1+izs] = src1x[isrc+src.n];
                        vz[ixs*n1+izs] = src1z[isrc+src.n];
                        vx[ixs*n1+izs-1] = src2x[isrc+src.n];
                        vz[(ixs-1)*n1+izs] = src2z[isrc+src.n];
                    }
                    else if (src.orient == 3) { /* dipole - + horizontal */
                        ixs += 1;
                        vx[ixs*n1+izs] = src1x[isrc+src.n];
                        vz[ixs*n1+izs] = src1z[isrc+src.n];
                        vx[ixs*n1+izs-1] = src2x[isrc+src.n];
                        vz[(ixs-1)*n1+izs] = src2z[isrc+src.n];
                    }
                    
				}
               	else if (src.type==6) {
                   	vx[ixs*n1+izs] = src1x[isrc];
               	}
               	else if (src.type==7) {
                   	vz[ixs*n1+izs] = src1z[isrc];
               	}
               	else if (src.type==8) {
                    
                    vx[(ixs+1)*n1+izs] = src1x[isrc];
                    vz[ixs*n1+izs+1] = src1z[isrc];
                    vx[ixs*n1+izs] = src2x[isrc];
                    vz[ixs*n1+izs] = src2z[isrc];
                    
                    /* determine second position of dipole */
                    if (src.orient == 2) { /* dipole +/- vertical */
                        izs += 1;
                        vx[(ixs+1)*n1+izs] = src1x[isrc+src.n];
                        vz[ixs*n1+izs+1] = src1z[isrc+src.n];
                        vx[ixs*n1+izs] = src2x[isrc+src.n];
                        vz[ixs*n1+izs] = src2z[isrc+src.n];
                    }
                    else if (src.orient == 3) { /* dipole - + horizontal */
                        ixs += 1;
                        vx[(ixs+1)*n1+izs] = src1x[isrc+src.n];
                        vz[ixs*n1+izs+1] = src1z[isrc+src.n];
                        vx[ixs*n1+izs] = src2x[isrc+src.n];
                        vz[ixs*n1+izs] = src2z[isrc+src.n];
                    }
                    
               	}
           	}
		}
    }

      return 0;
}
