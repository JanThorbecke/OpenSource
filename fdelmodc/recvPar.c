#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "fdelmodc.h"
#include "par.h"

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

/**
*  Calculates the receiver positions based on the input parameters
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/


void name_ext(char *filename, char *extension);

int recvPar(recPar *rec, float sub_x0, float sub_z0, float dx, float dz, int nx, int nz)
{
	float *xrcv1, *xrcv2, *zrcv1, *zrcv2;
	int ix0, ix1, iz0, iz1, i, ix, iz, ir, isign, verbose;
	float dxrcv, dzrcv, *dxr, *dzr, r, rr;
	float rrcv, dphi, oxrcv, ozrcv;
	float xrange, zrange;
	int Nx1, Nx2, Nz1, Nz2, Ndx, Ndz, iarray, nskip, nrec;
	int nxrcv, nzrcv, ncrcv, nrcv, max_nrec;
	float *xrcva, *zrcva;

	if(!getparint("verbose", &verbose)) verbose = 0;
	if (!getparint("max_nrec",&max_nrec)) max_nrec=15000;


	nrec=0;
	/* check if receiver positions on a circle are defined */
	if (getparfloat("rrcv", &rrcv)) {
		if (!getparfloat("dphi",&dphi)) dphi=2.0;
		if (!getparfloat("oxrcv",&oxrcv)) oxrcv=0.0;
		if (!getparfloat("ozrcv",&ozrcv)) ozrcv=0.0;
		ncrcv = NINT(360.0/dphi);
		
		for (ix=0; ix<ncrcv; ix++) {
			rec->xr[ix] = oxrcv-sub_x0+rrcv*cos(((ix*dphi)/360.0)*(2.0*M_PI));
			rec->zr[ix] = ozrcv-sub_z0+rrcv*sin(((ix*dphi)/360.0)*(2.0*M_PI));
			rec->x[ix] = NINT((oxrcv-sub_x0+rrcv*cos(((ix*dphi)/360.0)*(2.0*M_PI)))/dx);
			rec->z[ix] = NINT((ozrcv-sub_z0+rrcv*sin(((ix*dphi)/360.0)*(2.0*M_PI)))/dz);
			if (verbose>4) fprintf(stderr,"Receiver Circle: xrcv[%d]=%f zrcv=%f\n", ix, rec->xr[ix]+sub_x0, rec->zr[ix]+sub_z0);
		}
		nrec += ncrcv;
	}

	/* check if receiver array is defined */
	nxrcv = countparval("xrcva");
	nzrcv = countparval("zrcva");
	if (nxrcv != nzrcv) {
		verr("Number of receivers in array xrcva (%d), zrcva(%d) are not equal",nxrcv, nzrcv);
	}       
	
	if (nxrcv != 0) {
		/* receiver array is defined */
		xrcva = (float *)malloc(nxrcv*sizeof(float));
		zrcva = (float *)malloc(nxrcv*sizeof(float));
		getparfloat("xrcva", xrcva);
		getparfloat("zrcva", zrcva);
		for (ix=0; ix<nxrcv; ix++) {
			rec->xr[nrec+ix] = xrcva[ix]-sub_x0;
			rec->zr[nrec+ix] = zrcva[ix]-sub_z0;
			rec->x[nrec+ix] = NINT((xrcva[ix]-sub_x0)/dx);
			rec->z[nrec+ix] = NINT((zrcva[ix]-sub_z0)/dz);
			if (verbose>4) fprintf(stderr,"Receiver Array: xrcv[%d]=%f zrcv=%f\n", ix, rec->xr[nrec+ix]+sub_x0, rec->zr[nrec+ix]+sub_z0);
		}
		nrec += nxrcv;
		free(xrcva);
		free(zrcva);
	}
	
	/* linear receiver arrays */
	
	Nx1 = countparval("xrcv1");
	Nx2 = countparval("xrcv2");
	Nz1 = countparval("zrcv1");
	Nz2 = countparval("zrcv2");
	assert(Nx1==Nx2);
	assert(Nz1==Nz2);
	assert(Nx1==Nz1);
		
	if (nrec==0 && Nx1==0) { /* no receivers are defined use default linear array of receivers on top of model */
		Nx1=1;
	}
	
	if (Nx1!=0) {
		xrcv1 = (float *)malloc(Nx1*sizeof(float));
		xrcv2 = (float *)malloc(Nx1*sizeof(float));
		zrcv1 = (float *)malloc(Nx1*sizeof(float));
		zrcv2 = (float *)malloc(Nx1*sizeof(float));
		
		if(!getparfloat("xrcv1", xrcv1)) xrcv1[0]=sub_x0;
		if(!getparfloat("xrcv2", xrcv2)) xrcv2[0]=(nx-1)*dx+sub_x0;
		if(!getparfloat("zrcv1", zrcv1)) zrcv1[0]=sub_z0;
		if(!getparfloat("zrcv2", zrcv2)) zrcv2[0]=zrcv1[0];		
		
		Ndx = countparval("dxrcv");
		Ndz = countparval("dzrcv");

		dxr = (float *)malloc(Nx1*sizeof(float));
		dzr = (float *)malloc(Nx1*sizeof(float));
		if(!getparfloat("dxrcv", dxr)) dxr[0]=dx;
		if(!getparfloat("dzrcv", dzr)) dzr[0]=0.0;
		if ( (Ndx<=1) && (Ndz==0) ){ /* default values are set */
			for (i=1; i<Nx1; i++) {
				dxr[i] = dxr[0];
				dzr[i] = dzr[0];
			}
			Ndx=1;
		}
		else if ( (Ndz==1) && (Ndx==0) ){ /* default values are set */
			for (i=1; i<Nx1; i++) {
				dxr[i] = dxr[0];
				dzr[i] = dzr[0];
			}
			Ndz=1;
		}
		else { /* make sure that each array has dzrcv or dxrcv defined for each line or receivers */
			if (Ndx>1) assert(Ndx==Nx1);
			if (Ndz>1) assert(Ndz==Nx1);
		}
		
/*
		if ( (Ndx!=0) && (Ndz!=0) ) {
			vwarn("Both dzrcv and dxrcv are set: dxrcv value is used");
			Ndz=0;
			for (i=0; i<Nx1; i++) dzr[i] = 0.0;
		}
*/
		
		/* check if receiver arrays fit into model */
		for (iarray=0; iarray<Nx1; iarray++) {
			xrcv1[iarray] = MAX(sub_x0,      xrcv1[iarray]);
			xrcv1[iarray] = MIN(sub_x0+nx*dx,xrcv1[iarray]);
			xrcv2[iarray] = MAX(sub_x0,      xrcv2[iarray]);
			xrcv2[iarray] = MIN(sub_x0+nx*dx,xrcv2[iarray]);
			
			zrcv1[iarray] = MAX(sub_z0,      zrcv1[iarray]);
			zrcv1[iarray] = MIN(sub_z0+nz*dz,zrcv1[iarray]);
			zrcv2[iarray] = MAX(sub_z0,      zrcv2[iarray]);
			zrcv2[iarray] = MIN(sub_z0+nz*dz,zrcv2[iarray]);
		}


		/* calculate receiver array and store into rec->x,z */
		
		for (iarray=0; iarray<Nx1; iarray++) {
	
			xrange = (xrcv2[iarray]-xrcv1[iarray]); 
			zrange = (zrcv2[iarray]-zrcv1[iarray]); 
			if (dxr[iarray] != 0.0) {
				nrcv = NINT(abs(xrange/dxr[iarray]))+1;
				dxrcv=dxr[iarray];
				dzrcv = zrange/(nrcv-1);
				if (dzrcv != dzr[iarray]) {
					vwarn("For receiver array %d: calculated dzrcv=%f given=%f", iarray, dzrcv, dzr[iarray]);
					vwarn("The calculated receiver distance %f is used", dzrcv);
				}
			}
			else {
				if (dzr[iarray] == 0) {
					verr("For receiver array %d: receiver distance dzrcv is not given", iarray);
				}
				nrcv=NINT(abs(zrange/dzr[iarray]))+1;
				dxrcv = xrange/(nrcv-1);
				dzrcv = dzr[iarray];
				if (dxrcv != dxr[iarray]) {
					vwarn("For receiver array %d: calculated dxrcv=%f given=%f", iarray, dxrcv, dxr[iarray]);
					vwarn("The calculated receiver distance %f is used", dxrcv);
				}
			}

			// calculate coordinates
			for (ir=0; ir<nrcv; ir++) {
				rec->xr[nrec]=xrcv1[iarray]-sub_x0+ir*dxrcv;
				rec->zr[nrec]=zrcv1[iarray]-sub_z0+ir*dzrcv;

				rec->x[nrec]=NINT((rec->xr[nrec])/dx);
				rec->z[nrec]=NINT((rec->zr[nrec])/dz);
				nrec++;
			}
			if (nrec >= max_nrec) {
				verr("Number of receivers in arrays xrcv1,xrcv2,zrcv1,zrcv2 are larger than max: increas max_nrec %d",max_nrec);
			}
		}
	
		free(xrcv1);
		free(xrcv2);
		free(zrcv1);
		free(zrcv2);
		free(dxr);
		free(dzr);
	}
	rec->n = nrec;

	return 0;
}

