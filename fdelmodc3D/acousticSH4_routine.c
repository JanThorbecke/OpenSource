#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>

int taperEdges(modPar mod, bndPar bnd, float *vx, float *vz, int verbose);

int applySource(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *l2m, float **src_nwav, int verbose);

int acousticSH4_routine_(int *nxf,  int *nzf, int *ldz, int *it0, int *it1, int *src_type, wavPar wav, bndPar bnd, int *ixsrc, int *izsrc, float **src_nwav, float *tx, float *tz, float *vz, float *ro, float *mul, int verbose);


/*********************************************************************
 COMPUTATIONAL OVERVIEW OF THE 4th ORDER STAGGERED GRID: 
 
 The captial symbols T (=Txx,Tzz) Txz,Vx,Vz represent the actual grid
 The indices ix,iz are related to the T grid, so the capital T 
 symbols represent the actual modelled grid.
 
 one cel (iz,ix)
 |
 V                              extra column of vx,txz
 |
 -------                                           V
 | txz vz| txz vz  txz vz  txz vz  txz vz  txz vz txz
 |       |      
 | vx  t | vx  t   vx  t   vx  t   vx  t   vx  t  vx
 -------
 txz vz  txz vz  txz vz  txz vz  txz vz  txz vz  txz
 
 vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx
 |   |   |   |   |   |   | 
 txz vz  txz Vz--Txz-Vz--Txz-Vz  Txz-Vz  txz vz  txz
 |   |   |   |   |   |   |
 vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx
 |   |   |   |   |   |   |
 txz vz  txz Vz  Txz-Vz  Txz-Vz  Txz-Vz  txz vz  txz
 |   |   |   |   |   |   |
 vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx
 |   |   |   |   |   |   |
 txz vz  txz Vz  Txz-Vz  Txz-Vz  Txz-Vz  txz vz  txz
 |   |   |   |   |   |   |
 vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx
 
 txz vz  txz vz  txz vz  txz vz  txz vz  txz vz  txz
 
 vx  t   vx  t   vx  t   vx  t   vx  t   vx  t  vx
 
 txz vz  txz vz  txz vz  txz vz  txz vz  txz vz  txz  <--| 
 |
 extra row of txz/vz |
 
 AUTHOR:
 Jan Thorbecke (janth@xs4all.nl)
 The Netherlands 
 
 ***********************************************************************/

int main(int argc, char **argv)
{

}

int acousticSH4_routine_(int *nxf,  int *nzf, int *ldz, int *it0, int *it1, int *src_type, wavPar wav, bndPar bnd, int *ixsrc, int *izsrc, float **src_nwav, float *tx, float *tz, float *vz, float *ro, float *mul, int verbose)
{

	float c1, c2;
	float *tmps;
	int   ix, iz, ixs, izs, ibnd, store;
	int   nx, nz, n1;
	int   is0, isrc, ioXx, ioXz, ioZz, ioZx, ioPx, ioPz;


	c1 = 9.0/8.0; 
	c2 = -1.0/24.0;
	nx  = *nxf;
	nz  = *nzf;
	n1  = *ldz;

	ibnd = 1;

	ioXx=2;
	ioXz=ioXx-1;
	ioZz=2;
	ioZx=ioZz-1;
	ioPx=1;
	ioPz=ioPx;

#pragma omp parallel default (shared) \
shared (ro, mul, tx, tz, vz) \
shared (*it0, *it1, c1, c2) \
shared (shot, bnd, mod, src, wav, rec, ixsrc, izsrc, it, src_nwav, verbose)
{
    /* Main loop over the number of time steps */
    for (it=*it0; it<*it1; it++) {


        
	/* calculate vx for all grid points except on the virtual boundary*/
#pragma omp for private (ix, iz) nowait
	for (ix=ioXx; ix<nx+1; ix++) {
#pragma ivdep
		for (iz=ioXz; iz<nz+1; iz++) {
			tx[ix*n1+iz] -= mul[ix*n1+iz]*(
						c1*(vz[ix*n1+iz]     - vz[(ix-1)*n1+iz]) +
						c2*(vz[(ix+1)*n1+iz] - vz[(ix-2)*n1+iz]));
		}
	}

	/* calculate vz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iz) 
	for (ix=ioZx; ix<nx+1; ix++) {
#pragma ivdep
		for (iz=ioZz; iz<nz+1; iz++) {
			tz[ix*n1+iz] -= mul[ix*n1+iz]*(
						c1*(vz[ix*n1+iz]   - vz[ix*n1+iz-1]) +
						c2*(vz[ix*n1+iz+1] - vz[ix*n1+iz-2]));
		}
	}

	/* Add force source */
	if (src.type > 5) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, tx, tz, vz, NULL, NULL, ro, mul, src_nwav, verbose);
	}


	/* calculate p/tzz for all grid points except on the virtual boundary */
#pragma omp	for private (ix, iz)
#pragma ivdep
	for (ix=ioPx; ix<nx+1; ix++) {
#pragma ivdep
		for (iz=ioPz; iz<nz+1; iz++) {
			vz[ix*n1+iz] -= ro[ix*n1+iz]*(
						c1*(tx[(ix+1)*n1+iz] - tx[ix*n1+iz]) +
						c2*(tx[(ix+2)*n1+iz] - tx[(ix-1)*n1+iz]) +
						c1*(tz[ix*n1+iz+1]   - tz[ix*n1+iz]) +
						c2*(tz[ix*n1+iz+2]   - tz[ix*n1+iz-1]));
		}
	}

	/* Add stress source */
	if (src.type < 6) {
		 applySource(mod, src, wav, bnd, itime, *ixsrc, *izsrc, tx, tz, vz, NULL, NULL, ro, mul, src_nwav, verbose);
	}
    
/* Free surface: calculate free surface conditions for stresses */

	/* check if there are sources placed on the free surface */
	store=0;
	if (src.type==1 || src.type==6) {
		tmps = (float *)calloc(1, sizeof(float));
		is0 = -1*floor((src.n-1)/2);
		isrc=0;
        /* calculate the source position */
        ixs = *ixsrc + ibnd + is0 + isrc;
        izs = *izsrc + ibnd;
        if (ixs == 1) store=1;
        if (ixs == nx) store=1;
        if (izs == 1) store=1;
        if (izs == nz) store=1;
        if (store) {
            if (src.type==1) tmps[isrc] = vz[ixs*n1+izs];
            else tmps[isrc] = tx[ixs*n1+izs];
        }
		
	}

	if (bnd.free[0]) { /* free surface at top */
#pragma omp	for private (ix) nowait
		for (ix=1; ix<=nx; ix++) {
			iz = bnd.surface[ix-1];
			vz[ix*n1+iz] = 0.0;
		}
	}
	if (bnd.free[1]) { /* free surface at right */
#pragma omp	for private (iz) nowait
		for (iz=1; iz<=nz; iz++) {
			vz[nx*n1+iz] = 0.0;
		}
	}
	if (bnd.free[2]) { /* free surface at bottom */
#pragma omp	for private (ix) nowait
		for (ix=1; ix<=nx; ix++) {
			vz[ix*n1+nz] = 0.0;
		}
	}
	if (bnd.free[3]) { /* free surface at left */
#pragma omp	for private (iz) nowait
		for (iz=1; iz<=nz; iz++) {
			vz[n1+iz] = 0.0;
		}
	}

	/* restore source positions on the edge */
	if (src.type==1 || src.type==6) {
		if (store) {
#pragma omp	for private (isrc)
			for (isrc=0; isrc<src.n; isrc++) {
				/* calculate the source position */
                ixs = *ixsrc + ibnd + is0 + isrc;
                izs = *izsrc + ibnd;
				if (src.type==1) vz[ixs*n1+izs] += tmps[isrc];
				else tx[ixs*n1+izs] += tmps[isrc];
			}
		}
		free(tmps);
	}
        
    /* taper the edges of the model */
    taperEdges(mod, bnd, vx, vz, verbose);

    } /*end of time loop */
} /* end of OMP parallel */

	return 0;
}



int applySource(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, float *ro, float *l2m, float **src_nwav, int verbose)
{
    int is0, ibndz, ibndx;
    int isrc, ix, iz, n1;
    int id1, id2;
    float src_ampl, time, scl, dt;
    static int first=1;
    
    else if (src_type==7) {
        ibndz = mod.iorder/2;
        ibndx = mod.iorder/2-1;
    }
    else {  
        ibndz = mod.iorder/2-1;
        ibndx = mod.iorder/2-1;
    }
    
    n1   = mod.naz;
    dt   = mod.dt;
            
#pragma omp     for private (isrc, src_ampl, ix, iz, time, id1, id2, scl) 
        src_ampl=0.0;
        ix = ixsrc + ibndx;
        iz = izsrc + ibndz;
        time = itime*dt - src.tbeg[isrc];
        id1 = floor(time/dt);
        id2 = id1+1;
    
        /* delay not reached or no samples left in source wavelet? */
        if ( (time < 0.0) || ( (itime*dt) >= src.tend[isrc]) ) continue;
                
        src_ampl = src_nwav[0][id1]*(id2-time/dt) + src_nwav[0][id2]*(time/dt-id1);
        
        if (src_ampl==0.0) continue;
    
        if ( ((ix-ibndx)<0) || ((ix-ibndx)>mod.nx) ) continue; /* source outside grid */
                        
        /* source scaling factor to compensate for discretisation */
                
        src_ampl *= (1.0/mod.dx)*l2m[ix*n1+iz];
        
        /* Force source */
        
        if (src.type == 7) {
            vz[ix*n1+iz] += src_ampl*ro[ix*n1+iz]/(l2m[ix*n1+iz]);
        }
        else if (src.type == 2) {
            txz[ix*n1+iz] += src_ampl;
        }
        /* Tzz source */
        else if(src.type == 3) {
            tzz[ix*n1+iz] += src_ampl;
        } 
    
    return 0;
}



int taperEdges(modPar mod, bndPar bnd, float *vx, float *vz, int verbose)
{
    int   ix, iz, ibnd, ib, ntaper;
    int   nx, nz, n1;
    
    nx  = mod.nx;
    nz  = mod.nz;
    n1  = mod.naz;
    ibnd = mod.iorder/2-1;
    
    /* top */
    if (bnd.tap[0] > 0) {
        ntaper = bnd.tap[0];
        ib = (ntaper+ibnd-1);
#pragma omp for private(ix,iz)
        for (ix=ibnd; ix<nx+ibnd; ix++) {
#pragma ivdep
            for (iz=ibnd; iz<ibnd+ntaper; iz++) {
                vx[ix*n1+iz] *= bnd.tapx[ib-iz];
                vz[ix*n1+iz+1] *= bnd.tapz[ib-iz];
            }
        }
    }
    /* right */
    if (bnd.tap[1] > 0) {
        ntaper = bnd.tap[1];
        ib = (nx+ibnd-ntaper);
#pragma omp for private(ix,iz)
        for (ix=nx+ibnd-ntaper; ix<nx+ibnd; ix++) {
#pragma ivdep
            for (iz=ibnd; iz<nz+ibnd; iz++) {
                vx[ix*n1+iz] *= bnd.tapx[ix-ib];
                vz[ix*n1+iz] *= bnd.tapz[ix-ib];
            }
        }
    }
    /* bottom */
    if (bnd.tap[2] > 0) {
        ntaper = bnd.tap[2];
        ib = (nz+ibnd-ntaper);
#pragma omp for private(ix,iz)
        for (ix=ibnd; ix<nx+ibnd; ix++) {
#pragma ivdep
            for (iz=nz+ibnd-ntaper; iz<nz+ibnd; iz++) {
                vx[ix*n1+iz]   *= bnd.tapx[iz-ib];
                vz[ix*n1+iz+1] *= bnd.tapz[iz-ib];
            }
        }
    }
    /* left */
    if (bnd.tap[3] > 0) {
        ntaper = bnd.tap[3];
        ib = (ntaper+ibnd-1);
#pragma omp for private(ix,iz)
        for (ix=ibnd; ix<ntaper+ibnd; ix++) {
#pragma ivdep
            for (iz=ibnd; iz<nz+ibnd; iz++) {
                vx[ix*n1+iz] *= bnd.tapx[ib-ix];
                vz[ix*n1+iz] *= bnd.tapz[ib-ix];
            }
        }
    }
    
    return 0;
}




