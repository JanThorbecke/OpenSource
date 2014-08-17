#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc.h"

/**
*  Taper the edges of the wavefield to suppress unwanted reflections from the sides of the model.
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/


int taperEdges(modPar mod, bndPar bnd, float *vx, float *vz, int verbose)
{
	int   ix, iz, ibnd, ib, ntaper;
	int   nx, nz, n1;

	nx  = mod.nx;
	nz  = mod.nz;
	n1  = mod.naz;
	ibnd = mod.iorder/2-1;

	/* top */
	if (bnd.ntap > 0) {
		ntaper = bnd.ntap;
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
	if (bnd.ntap > 0) {
		ntaper = bnd.ntap;
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
	if (bnd.ntap > 0) {
		ntaper = bnd.ntap;
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
	if (bnd.ntap > 0) {
		ntaper = bnd.ntap;
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
