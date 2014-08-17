#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc.h"

static float *pml_nxVx, *pml_nxVz, *pml_nxP;
static float *pml_nzVx, *pml_nzVz, *pml_nzP;

int applySource(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float **src_nwav, int verbose);

int acoustic4pml(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float **src_nwav, float *vx, float *vz, float *p, float *rox, float *roz, float *l2m, int verbose)
{
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

	float c1, c2;
	float *tmps;
	int   ix, iz, ixs, izs, ibnd, store;
	int   nx, nz, n1;
	int   is0, isrc, ioXx, ioXz, ioZz, ioZx, ioPx, ioPz;
    size_t sizem;


	c1 = 9.0/8.0; 
	c2 = -1.0/24.0;
	nx  = mod.nx;
	nz  = mod.nz;
	n1  = mod.naz;
	sizem=mod.naz*mod.nax;


	ibnd = mod.iorder/2-1;

	ioXx=mod.iorder/2;
	ioXz=ioXx-1;
	ioZz=mod.iorder/2;
	ioZx=ioZz-1;
	ioPx=mod.iorder/2-1;
	ioPz=ioPx;

        
    if (itime==0) {/* allocate arrays to store boundary values */
        pml_nxVx = (float *)calloc(mod.nt*(2*mod.nax),sizeof(float));
        pml_nxVz = (float *)calloc(mod.nt*(2*mod.nax),sizeof(float));
        pml_nxP  = (float *)calloc(mod.nt*(2*mod.nax),sizeof(float));
        pml_nzVx = (float *)calloc(mod.nt*(2*mod.naz),sizeof(float));
        pml_nzVz = (float *)calloc(mod.nt*(2*mod.naz),sizeof(float));
        pml_nzP  = (float *)calloc(mod.nt*(2*mod.naz),sizeof(float));
    }
    
    /* First pass though the modeling */

    /* inject the edges of the model */
    for (ix=ioXx; ix<nx+1; ix++) {
        iz = ioXz;
        vx[ix*n1+iz] -= pml_nxVx[ix];
        iz = nz-50;
        vx[ix*n1+iz] -= pml_nxVx[ix+mod.nax];
    }
    for (iz=ioXz; iz<nz+1; iz++) {
        ix = ioXx;
        vx[ix*n1+iz] -= pml_nzVx[iz];
        ix = nx;
        vx[ix*n1+iz] -= pml_nzVx[iz+mod.naz];
    }
    
    /* inject the edges of the model */
    for (ix=ioZx; ix<nx+1; ix++) {
        iz = ioZz;
        vz[ix*n1+iz] -= pml_nxVz[ix];
        iz = nz-50;
        vz[ix*n1+iz] -= pml_nxVz[ix+mod.nax];
    }
    for (iz=ioZz; iz<nz+1; iz++) {
        ix = ioZx;
        vz[ix*n1+iz] -= pml_nzVz[iz];
        ix = nx;
        vz[ix*n1+iz] -= pml_nzVz[iz+mod.naz];
    }

    /* inject the edges of the model */
    for (ix=ioPx; ix<nx+1; ix++) {
        iz = ioPz;
        p[ix*n1+iz] -= pml_nxP[ix];
        iz = nz-50;
        p[ix*n1+iz] -= pml_nxP[ix+mod.nax];
    }
    for (iz=ioPz; iz<nz+1; iz++) {
        ix = ioPx;
        p[ix*n1+iz] -= pml_nzP[iz];
        ix = nx;
        p[ix*n1+iz] -= pml_nzP[iz+mod.naz];
    }

	/* calculate vx for all grid points except on the virtual boundary*/
#pragma omp for private (ix, iz) nowait
	for (ix=ioXx; ix<nx+1; ix++) {
#pragma ivdep
		for (iz=ioXz; iz<nz+1; iz++) {
			vx[sizem+ix*n1+iz] -= rox[ix*n1+iz]*(
						c1*(p[sizem+ix*n1+iz]     - p[sizem+(ix-1)*n1+iz]) +
						c2*(p[sizem+(ix+1)*n1+iz] - p[sizem+(ix-2)*n1+iz]));
		}
	}


	/* calculate vz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iz) 
	for (ix=ioZx; ix<nx+1; ix++) {
#pragma ivdep
		for (iz=ioZz; iz<nz+1; iz++) {
			vz[sizem+ix*n1+iz] -= roz[ix*n1+iz]*(
						c1*(p[sizem+ix*n1+iz]   - p[sizem+ix*n1+iz-1]) +
						c2*(p[sizem+ix*n1+iz+1] - p[sizem+ix*n1+iz-2]));
		}
	}

	/* Add force source */
	if (src.type > 5) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, &vx[sizem], &vz[sizem], &p[sizem], NULL, NULL, rox, roz, l2m, src_nwav, verbose);
	}


    
	/* rigid boundary condition clears velocities on boundaries */
	if (bnd.rig[0]) { /* rigid surface at top */
#pragma omp	for private (ix, iz) nowait
#pragma ivdep
		for (ix=1; ix<=nx; ix++) {
/* ToDo			iz = bnd.surface[ix-ibnd]; */
			vx[sizem+ix*n1+1] = 0.0;
			vz[sizem+ix*n1+1] = -vz[sizem+ix*n1+2];
			vz[sizem+ix*n1+0] = -vz[sizem+ix*n1+3];
		}
	}
	if (bnd.rig[1]) { /* rigid surface at right */
#pragma omp	for private (ix, iz) nowait
#pragma ivdep
		for (iz=1; iz<=nz; iz++) {
            vz[sizem+nx*n1+iz]     = 0.0;
			vx[sizem+(nx+1)*n1+iz] = -vx[sizem+nx*n1+iz];
			vx[sizem+(nx+2)*n1+iz] = -vx[sizem+(nx-1)*n1+iz];
		}
	}
	if (bnd.rig[2]) { /* rigid surface at bottom */
#pragma omp	for private (ix, iz) nowait
#pragma ivdep
		for (ix=1; ix<=nx; ix++) {
			vx[sizem+ix*n1+nz]   = 0.0;
			vz[sizem+ix*n1+nz+1] = -vz[sizem+ix*n1+nz];
			vz[sizem+ix*n1+nz+2] = -vz[sizem+ix*n1+nz-1];
		}
	}
	if (bnd.rig[3]) { /* rigid surface at left */
#pragma omp	for private (ix, iz) nowait
#pragma ivdep
		for (iz=1; iz<=nz; iz++) {
            vz[sizem+n1+iz] = 0.0;
			vx[sizem+n1+iz] = -vx[sizem+2*n1+iz];
			vx[sizem+n1+iz] = -vx[sizem+3*n1+iz];
		}
	}


	/* calculate p/tzz for all grid points except on the virtual boundary */
#pragma omp	for private (ix, iz)
#pragma ivdep
	for (ix=ioPx; ix<nx+1; ix++) {
#pragma ivdep
		for (iz=ioPz; iz<nz+1; iz++) {
			p[sizem+ix*n1+iz] -= l2m[ix*n1+iz]*(
						c1*(vx[sizem+(ix+1)*n1+iz] - vx[sizem+ix*n1+iz]) +
						c2*(vx[sizem+(ix+2)*n1+iz] - vx[sizem+(ix-1)*n1+iz]) +
						c1*(vz[sizem+ix*n1+iz+1]   - vz[sizem+ix*n1+iz]) +
						c2*(vz[sizem+ix*n1+iz+2]   - vz[sizem+ix*n1+iz-1]));
		}
	}

	/* Add stress source */
	if (src.type < 6) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, &vx[sizem], &vz[sizem], &p[sizem], NULL, NULL, rox, roz, l2m, src_nwav, verbose);
	}
    

    
/* Free surface: calculate free surface conditions for stresses */

	/* check if there are sources placed on the free surface */
	store=0;
	if (src.type==1 || src.type==6) {
		tmps = (float *)calloc(src.n, sizeof(float));
		is0 = -1*floor((src.n-1)/2);
#pragma omp	for private (isrc) 
		for (isrc=0; isrc<src.n; isrc++) {
			/* calculate the source position */
			if (src.random) {
				ixs = src.x[isrc] + ibnd;
				izs = src.z[isrc] + ibnd;
			}
			else { /* plane wave and point sources */
				ixs = ixsrc + ibnd + is0 + isrc;
				izs = izsrc + ibnd;
			}
			if (ixs == 1) store=1;
			if (ixs == nx) store=1;
			if (izs == 1) store=1;
			if (izs == nz) store=1;
			if (store) {
				if (src.type==1) tmps[isrc] = p[sizem+ixs*n1+izs];
				else tmps[isrc] = vx[sizem+ixs*n1+izs];
			}
		}
	}

	if (bnd.free[0]) { /* free surface at top */
#pragma omp	for private (ix) nowait
		for (ix=1; ix<=nx; ix++) {
			iz = bnd.surface[ix-1];
			p[sizem+ix*n1+iz] = 0.0;
		}
	}
	if (bnd.free[1]) { /* free surface at right */
#pragma omp	for private (iz) nowait
		for (iz=1; iz<=nz; iz++) {
			p[sizem+nx*n1+iz] = 0.0;
		}
	}
	if (bnd.free[2]) { /* free surface at bottom */
#pragma omp	for private (ix) nowait
		for (ix=1; ix<=nx; ix++) {
			p[sizem+ix*n1+nz] = 0.0;
		}
	}
	if (bnd.free[3]) { /* free surface at left */
#pragma omp	for private (iz) nowait
		for (iz=1; iz<=nz; iz++) {
			p[sizem+n1+iz] = 0.0;
		}
	}

	/* restore source positions on the edge */
	if (src.type==1 || src.type==6) {
		if (store) {
#pragma omp	for private (isrc)
			for (isrc=0; isrc<src.n; isrc++) {
				/* calculate the source position */
				if (src.random) {
					ixs = src.x[isrc] + ibnd;
					izs = src.z[isrc] + ibnd;
				}
				else { /* plane wave sources */
					ixs = ixsrc + ibnd + is0 + isrc;
					izs = izsrc + ibnd;
				}
				if (src.type==1) p[sizem+ixs*n1+izs] += tmps[isrc];
				else vx[sizem+ixs*n1+izs] += tmps[isrc];
			}
		}
		free(tmps);
	}
    
    /* store the edges of the model */
    for (ix=ioXx; ix<nx+1; ix++) {
        iz = ioXz;
        pml_nxVx[ix]         = vx[sizem+ix*n1+iz];
        iz = nz-50;
        pml_nxVx[ix+mod.nax] = vx[sizem+ix*n1+iz];
    }
    for (iz=ioXz; iz<nz+1; iz++) {
        ix = ioXx;
        pml_nzVx[iz]         = vx[sizem+ix*n1+iz];
        ix = nx;
        pml_nzVx[iz+mod.naz] = vx[sizem+ix*n1+iz];
    }
    
    /* store the edges of the model */
	for (ix=ioZx; ix<nx+1; ix++) {
        iz = ioZz;
        pml_nxVz[ix]         = vz[sizem+ix*n1+iz];
        iz = nz-50;
        pml_nxVz[ix+mod.nax] = vz[sizem+ix*n1+iz];
    }
    for (iz=ioZz; iz<nz+1; iz++) {
        ix = ioZx;
        pml_nzVz[iz]         = vz[sizem+ix*n1+iz];
        ix = nx;
        pml_nzVz[iz+mod.naz] = vz[sizem+ix*n1+iz];
    }

    /* store the edges of the model */
	for (ix=ioPx; ix<nx+1; ix++) {
        iz = ioPz;
        pml_nxP[ix]         = p[sizem+ix*n1+iz];
        iz = nz-50;
        pml_nxP[ix+mod.nax] = p[sizem+ix*n1+iz];
    }
    for (iz=ioPz; iz<nz+1; iz++) {
        ix = ioPx;
        pml_nzP[iz]         = p[sizem+ix*n1+iz];
        ix = nx;
        pml_nzP[iz+mod.naz] = p[sizem+ix*n1+iz];
    }

    
    /****************************/
    
    if (itime!=0) { /* calculate the field with (hope)fully absorbing boundaries */
        
        /* Second pass through the modeling */
        
        /* calculate vx for all grid points except on the virtual boundary*/
#pragma omp for private (ix, iz) nowait
        for (ix=ioXx; ix<nx+1; ix++) {
#pragma ivdep
            for (iz=ioXz; iz<nz+1; iz++) {
                vx[ix*n1+iz] -= rox[ix*n1+iz]*(
                                               c1*(p[ix*n1+iz]     - p[(ix-1)*n1+iz]) +
                                               c2*(p[(ix+1)*n1+iz] - p[(ix-2)*n1+iz]));
            }
        }
        
        /* calculate vz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iz) 
        for (ix=ioZx; ix<nx+1; ix++) {
#pragma ivdep
            for (iz=ioZz; iz<nz+1; iz++) {
                vz[ix*n1+iz] -= roz[ix*n1+iz]*(
                                               c1*(p[ix*n1+iz]   - p[ix*n1+iz-1]) +
                                               c2*(p[ix*n1+iz+1] - p[ix*n1+iz-2]));
            }
        }
        
        /* Add force source */
        if (src.type > 5) {
            applySource(mod, src, wav, bnd, itime-1, ixsrc, izsrc, vx, vz, p, NULL, NULL, rox, roz, l2m, src_nwav, verbose);
        }
        
        
        
        /* rigid boundary condition clears velocities on boundaries */
        if (bnd.rig[0]) { /* rigid surface at top */
#pragma omp	for private (ix, iz) nowait
#pragma ivdep
            for (ix=1; ix<=nx; ix++) {
                /* ToDo			iz = bnd.surface[ix-ibnd]; */
                vx[ix*n1+1] = 0.0;
                vz[ix*n1+1] = -vz[ix*n1+2];
                vz[ix*n1+0] = -vz[ix*n1+3];
            }
        }
        if (bnd.rig[1]) { /* rigid surface at right */
#pragma omp	for private (ix, iz) nowait
#pragma ivdep
            for (iz=1; iz<=nz; iz++) {
                vz[nx*n1+iz]     = 0.0;
                vx[(nx+1)*n1+iz] = -vx[nx*n1+iz];
                vx[(nx+2)*n1+iz] = -vx[(nx-1)*n1+iz];
            }
        }
        if (bnd.rig[2]) { /* rigid surface at bottom */
#pragma omp	for private (ix, iz) nowait
#pragma ivdep
            for (ix=1; ix<=nx; ix++) {
                vx[ix*n1+nz]   = 0.0;
                vz[ix*n1+nz+1] = -vz[ix*n1+nz];
                vz[ix*n1+nz+2] = -vz[ix*n1+nz-1];
            }
        }
        if (bnd.rig[3]) { /* rigid surface at left */
#pragma omp	for private (ix, iz) nowait
#pragma ivdep
            for (iz=1; iz<=nz; iz++) {
                vz[n1+iz] = 0.0;
                vx[n1+iz] = -vx[2*n1+iz];
                vx[n1+iz] = -vx[3*n1+iz];
            }
        }
        
        
        /* calculate p/tzz for all grid points except on the virtual boundary */
#pragma omp	for private (ix, iz)
#pragma ivdep
        for (ix=ioPx; ix<nx+1; ix++) {
#pragma ivdep
            for (iz=ioPz; iz<nz+1; iz++) {
                p[ix*n1+iz] -= l2m[ix*n1+iz]*(
                                              c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
                                              c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]) +
                                              c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
                                              c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]));
            }
        }
        
        /* Add stress source */
        if (src.type < 6) {
            applySource(mod, src, wav, bnd, itime-1, ixsrc, izsrc, vx, vz, p, NULL, NULL, rox, roz, l2m, src_nwav, verbose);
        }
        
        
        
        /* Free surface: calculate free surface conditions for stresses */
        
        /* check if there are sources placed on the free surface */
        store=0;
        if (src.type==1 || src.type==6) {
            tmps = (float *)calloc(src.n, sizeof(float));
            is0 = -1*floor((src.n-1)/2);
#pragma omp	for private (isrc) 
            for (isrc=0; isrc<src.n; isrc++) {
                /* calculate the source position */
                if (src.random) {
                    ixs = src.x[isrc] + ibnd;
                    izs = src.z[isrc] + ibnd;
                }
                else { /* plane wave and point sources */
                    ixs = ixsrc + ibnd + is0 + isrc;
                    izs = izsrc + ibnd;
                }
                if (ixs == 1) store=1;
                if (ixs == nx) store=1;
                if (izs == 1) store=1;
                if (izs == nz) store=1;
                if (store) {
                    if (src.type==1) tmps[isrc] = p[ixs*n1+izs];
                    else tmps[isrc] = vx[ixs*n1+izs];
                }
            }
        }
        
        if (bnd.free[0]) { /* free surface at top */
#pragma omp	for private (ix) nowait
            for (ix=1; ix<=nx; ix++) {
                iz = bnd.surface[ix-1];
                p[ix*n1+iz] = 0.0;
            }
        }
        if (bnd.free[1]) { /* free surface at right */
#pragma omp	for private (iz) nowait
            for (iz=1; iz<=nz; iz++) {
                p[nx*n1+iz] = 0.0;
            }
        }
        if (bnd.free[2]) { /* free surface at bottom */
#pragma omp	for private (ix) nowait
            for (ix=1; ix<=nx; ix++) {
                p[ix*n1+nz] = 0.0;
            }
        }
        if (bnd.free[3]) { /* free surface at left */
#pragma omp	for private (iz) nowait
            for (iz=1; iz<=nz; iz++) {
                p[n1+iz] = 0.0;
            }
        }
        
        /* restore source positions on the edge */
        if (src.type==1 || src.type==6) {
            if (store) {
#pragma omp	for private (isrc)
                for (isrc=0; isrc<src.n; isrc++) {
                    /* calculate the source position */
                    if (src.random) {
                        ixs = src.x[isrc] + ibnd;
                        izs = src.z[isrc] + ibnd;
                    }
                    else { /* plane wave sources */
                        ixs = ixsrc + ibnd + is0 + isrc;
                        izs = izsrc + ibnd;
                    }
                    if (src.type==1) p[ixs*n1+izs] += tmps[isrc];
                    else vx[ixs*n1+izs] += tmps[isrc];
                }
            }
            free(tmps);
        }

    }
    
    /* inject the edges of the model */
    for (ix=ioXx; ix<nx+1; ix++) {
        iz = ioXz;
        vx[ix*n1+iz] -= pml_nxVx[ix];
        iz = nz-50;
        vx[ix*n1+iz] -= pml_nxVx[ix+mod.nax];
    }
    for (iz=ioXz; iz<nz+1; iz++) {
        ix = ioXx;
        vx[ix*n1+iz] -= pml_nzVx[iz];
        ix = nx;
        vx[ix*n1+iz] -= pml_nzVx[iz+mod.naz];
    }
    
    /* inject the edges of the model */
    for (ix=ioZx; ix<nx+1; ix++) {
        iz = ioZz;
        vz[ix*n1+iz] -= pml_nxVz[ix];
        iz = nz-50;
        vz[ix*n1+iz] -= pml_nxVz[ix+mod.nax];
    }
    for (iz=ioZz; iz<nz+1; iz++) {
        ix = ioZx;
        vz[ix*n1+iz] -= pml_nzVz[iz];
        ix = nx;
        vz[ix*n1+iz] -= pml_nzVz[iz+mod.naz];
    }

    /* inject the edges of the model */
    for (ix=ioPx; ix<nx+1; ix++) {
        iz = ioPz;
        p[ix*n1+iz] -= pml_nxP[ix];
        iz = nz-50;
        p[ix*n1+iz] -= pml_nxP[ix+mod.nax];
    }
    for (iz=ioPz; iz<nz+1; iz++) {
        ix = ioPx;
        p[ix*n1+iz] -= pml_nzP[iz];
        ix = nx;
        p[ix*n1+iz] -= pml_nzP[iz+mod.naz];
    }

    

	return 0;
}
