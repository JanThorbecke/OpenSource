#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "fdelmodc3D.h"
#include "par.h"

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

/**
*  Calculates the receiver positions based on the input parameters
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*
*   Ammendments:
*           Max Holicki changing the allocation receiver array (2-2016)
*           Joeri Brackenhoff adding the 3D extension
*           The Netherlands 
**/


void name_ext(char *filename, char *extension);

long recvPar3D(recPar *rec, float sub_x0, float sub_y0, float sub_z0, 
	float dx, float dy, float dz, long nx, long ny, long nz)
{
	float   *xrcv1, *xrcv2, *yrcv1, *yrcv2, *zrcv1, *zrcv2;
	long    i, ix, iy, ir, verbose;
	float   dxrcv, dyrcv, dzrcv, *dxr, *dyr, *dzr;
	float   rrcv, dphi, oxrcv, oyrcv, ozrcv, arcv;
	double  circ, h, a, b, e, s, xr, yr, zr, dr, srun, phase;
	float   xrange, yrange, zrange, sub_x1, sub_y1, sub_z1;
	long    Nx1, Nx2, Ny1, Ny2, Nz1, Nz2, Ndx, Ndy, Ndz, iarray, nrec, nh;
	long    nxrcv, nyrcv, nzrcv, ncrcv, nrcv, ntrcv, *nlxrcv, *nlyrcv;
	float   *xrcva, *yrcva, *zrcva;
	char*   rcv_txt;
	FILE    *fp;

	if (!getparlong("verbose", &verbose)) verbose = 0;

    /* Calculate Model Dimensions */
    sub_x1=sub_x0+(nx-1)*dx;
    sub_y1=sub_y0+(ny-1)*dy;
    sub_z1=sub_z0+(nz-1)*dz;

/* Compute how many receivers are defined and then allocate the receiver arrays */

    /* Receiver Array */
    nxrcv=countparval("xrcva");
    nyrcv=countparval("yrcva");
    nzrcv=countparval("zrcva");
    if (nxrcv!=nzrcv) verr("Number of receivers in array xrcva (%li), yrcva (%li), zrcva(%li) are not equal",nxrcv,nyrcv,nzrcv);
    if (verbose&&nxrcv) vmess("Total number of array receivers: %li",nxrcv);

    /* Linear Receiver Arrays */
	Nx1 = countparval("xrcv1");
	Nx2 = countparval("xrcv2");
	Ny1 = countparval("yrcv1");
	Ny2 = countparval("yrcv2");
	Nz1 = countparval("zrcv1");
	Nz2 = countparval("zrcv2");
    if (Nx1!=Nx2) verr("Number of receivers starting points in 'xrcv1' (%li) and number of endpoint in 'xrcv2' (%li) are not equal",Nx1,Nx2);
    if (Ny1!=Ny2) verr("Number of receivers starting points in 'yrcv1' (%li) and number of endpoint in 'yrcv2' (%li) are not equal",Ny1,Ny2);
    if (Nz1!=Nz2) verr("Number of receivers starting points in 'zrcv1' (%li) and number of endpoint in 'zrcv2' (%li) are not equal",Nz1,Nz2);
    if (Nx1!=Ny2) verr("Number of receivers starting points in 'xrcv1' (%li) and number of endpoint in 'yrcv2' (%li) are not equal",Nx1,Ny2);
    if (Nx1!=Nz2) verr("Number of receivers starting points in 'xrcv1' (%li) and number of endpoint in 'zrcv2' (%li) are not equal",Nx1,Nz2);

    rec->max_nrec=nyrcv*nxrcv;

	/* no receivers are defined use default linear array of receivers on top of model */
    if (!rec->max_nrec && Nx1==0) Nx1=1; // Default is to use top of model to record data
    if (!rec->max_nrec && Ny1==0) Ny1=1;
    if (!rec->max_nrec && Nz1==0) Nz1=1;

    if (Nx1) {
        /* Allocate Start & End Points of Linear Arrays */
        xrcv1=(float *)malloc(Nx1*sizeof(float));
        xrcv2=(float *)malloc(Nx1*sizeof(float));
        yrcv1=(float *)malloc(Nx1*sizeof(float));
        yrcv2=(float *)malloc(Nx1*sizeof(float));
        zrcv1=(float *)malloc(Nx1*sizeof(float));
        zrcv2=(float *)malloc(Nx1*sizeof(float));
        if (!getparfloat("xrcv1",xrcv1)) xrcv1[0]=sub_x0;
        if (!getparfloat("xrcv2",xrcv2)) xrcv2[0]=sub_x1;
        if (!getparfloat("yrcv1",yrcv1)) yrcv1[0]=sub_y0;
        if (!getparfloat("yrcv2",yrcv2)) yrcv2[0]=sub_y1;
        if (!getparfloat("zrcv1",zrcv1)) zrcv1[0]=sub_z0;
        if (!getparfloat("zrcv2",zrcv2)) zrcv2[0]=zrcv1[0];

		/* check if receiver arrays fit into model */
		for (iarray=0; iarray<Nx1; iarray++) {
			xrcv1[iarray] = MAX(sub_x0,      xrcv1[iarray]);
			xrcv1[iarray] = MIN(sub_x0+nx*dx,xrcv1[iarray]);
			xrcv2[iarray] = MAX(sub_x0,      xrcv2[iarray]);
			xrcv2[iarray] = MIN(sub_x0+nx*dx,xrcv2[iarray]);
			
			yrcv1[iarray] = MAX(sub_y0,      yrcv1[iarray]);
			yrcv1[iarray] = MIN(sub_y0+ny*dy,yrcv1[iarray]);
			yrcv2[iarray] = MAX(sub_y0,      yrcv2[iarray]);
			yrcv2[iarray] = MIN(sub_y0+ny*dy,yrcv2[iarray]);

			zrcv1[iarray] = MAX(sub_z0,      zrcv1[iarray]);
			zrcv1[iarray] = MIN(sub_z0+nz*dz,zrcv1[iarray]);
			zrcv2[iarray] = MAX(sub_z0,      zrcv2[iarray]);
			zrcv2[iarray] = MIN(sub_z0+nz*dz,zrcv2[iarray]);
		}

        /* Crop to Fit Model */
/* Max's addtion still have to check if it has the same fucntionality */
        for (iarray=0;iarray<Nx1;iarray++) {
            if (xrcv1[iarray]<sub_x0) {
                if (xrcv2[iarray]<sub_x0) {
                    verr("Linear array %li outside model bounds",iarray);
                }
				else {
                    vwarn("Cropping element %li of 'xrcv1' (%f) to model bounds (%f)",iarray,xrcv1[iarray],sub_x0);
                    xrcv1[iarray]=sub_x0;
                }
            } 
			else if (xrcv1[iarray] > sub_x1) {
                verr("Linear array %li outside model bounds",iarray);
            }
            if ( (xrcv2[iarray] < xrcv1[iarray]) ) {
                verr("Ill defined linear array %li, 'xrcv1' (%f) greater than 'xrcv2' (%f)",iarray,xrcv1[iarray],xrcv2[iarray]);
            }
			else if (xrcv2[iarray]>sub_x1) {
                vwarn("Cropping element %li of 'xrcv2' (%f) to model bounds (%f)",iarray,xrcv2[iarray],sub_x1);
                xrcv2[iarray]=sub_x1;
            }

            if (yrcv1[iarray]<sub_y0) {
                if (yrcv2[iarray]<sub_y0) {
                    verr("Linear array %li outside model bounds",iarray);
                }
				else {
                    vwarn("Cropping element %li of 'yrcv1' (%f) to model bounds (%f)",iarray,yrcv1[iarray],sub_y0);
                    yrcv1[iarray]=sub_y0;
                }
            } 
			else if (yrcv1[iarray] > sub_y1) {
                verr("Linear array %li outside model bounds",iarray);
            }
            if ( (yrcv2[iarray] < yrcv1[iarray]) ) {
                verr("Ill defined linear array %li, 'yrcv1' (%f) greater than 'yrcv2' (%f)",iarray,yrcv1[iarray],yrcv2[iarray]);
            }
			else if (yrcv2[iarray]>sub_y1) {
                vwarn("Cropping element %li of 'yrcv2' (%f) to model bounds (%f)",iarray,yrcv2[iarray],sub_y1);
                yrcv2[iarray]=sub_y1;
            }

            if (zrcv1[iarray] < sub_z0) {
                if (zrcv2[iarray] < sub_z0) {
                    verr("Linear array %li outside model bounds",iarray);
                }
				else {
               		vwarn("Cropping element %li of 'zrcv1' (%f) to model bounds (%f)",iarray,zrcv1[iarray],sub_z0);
                	zrcv1[iarray]=sub_z0;
                }
            }
			else if (zrcv1[iarray] > sub_z1) {
                verr("Linear array %li outside model bounds",iarray);
            }
            if ( (zrcv2[iarray] < zrcv1[iarray]) ) {
                verr("Ill defined linear array %li, 'zrcv1' (%f) greater than 'zrcv2' (%f)",iarray,zrcv1[iarray],zrcv2[iarray]);
            }
			else if (zrcv2[iarray]>sub_z1) {
                vwarn("Cropping element %li of 'xrcv2' (%f) to model bounds (%f)",iarray,zrcv2[iarray],sub_z1);
                zrcv2[iarray]=sub_z1;
            }
        }

        /* Get Sampling Rates */
		Ndx = countparval("dxrcv");
		Ndy = countparval("dyrcv");
		Ndz = countparval("dzrcv");

		dxr = (float *)malloc(Nx1*sizeof(float));
		dyr = (float *)malloc(Nx1*sizeof(float));
		dzr = (float *)malloc(Nx1*sizeof(float));
		if(!getparfloat("dxrcv", dxr)) dxr[0]=dx;
		if(!getparfloat("dyrcv", dyr)) dyr[0]=dy;
		if(!getparfloat("dzrcv", dzr)) dzr[0]=0.0;
		if ( (Ndx<=1) && (Ndy<=1) && (Ndz==0) ){ /* default values are set */
			for (i=1; i<Nx1; i++) {
				dxr[i] = dxr[0];
				dyr[i] = dyr[0];
				dzr[i] = dzr[0];
			}
			Ndx=1;
            Ndy=1;
			Ndz=1;
		}
		else if ( (Ndz==1) && (Ndx==0) && (Ndy==0) ){ /* default values are set */
			for (i=1; i<Nx1; i++) {
				dxr[i] = dxr[0];
				dyr[i] = dyr[0];
				dzr[i] = dzr[0];
			}
			Ndz=1;
            Ndy=1;
			Ndx=1;
		}
		else { /* make sure that each array has dzrcv or dxrcv defined for each line or receivers */
			if (Ndx!=Ndz) {
				verr("Number of 'dxrcv' (%li) is not equal to number of 'dzrcv' (%li) or 1",Ndx,Ndz);
			}
			if (Ndx!=Ndy) {
				verr("Number of 'dxrcv' (%li) is not equal to number of 'dyrcv' (%li) or 1",Ndx,Ndy);
			}
			if (Ndx!=Nx1 && Ndx!=1) {
				verr("Number of 'dxrcv' (%li) is not equal to number of starting points in 'xrcv1' (%li) or 1",Ndx,Nx1);
			}
			if (Ndy!=Ny1 && Ndy!=1) {
				verr("Number of 'dyrcv' (%li) is not equal to number of starting points in 'yrcv1' (%li) or 1",Ndy,Ny1);
			}
		}

		/* check consistency of receiver steps */
        for (iarray=0; iarray<Ndx; iarray++) {
            if (dxr[iarray]<0) {
				dxr[i]=dx;
				vwarn("'dxrcv' element %li (%f) is less than zero, changing it to %f'",iarray,dxr[iarray],dx);
			}
        }
        for (iarray=0; iarray<Ndy; iarray++) {
            if (dyr[iarray]<0) {
				dyr[i]=dx;
				vwarn("'dyrcv' element %li (%f) is less than zero, changing it to %f'",iarray,dyr[iarray],dy);
			}
        }
        for (iarray=0;iarray<Ndz;iarray++) {
            if (dzr[iarray]<0) {
				dzr[iarray]=dz;
				vwarn("'dzrcv' element %li (%f) is less than zero, changing it to %f'",iarray,dzr[iarray],dz);
			}
        }
        for (iarray=0;iarray<Ndx;iarray++){
            if (dxr[iarray]==0 && dzr[iarray]==0) {
                xrcv2[iarray]=xrcv1[iarray];
				dxr[iarray]=1.;
                vwarn("'dxrcv' element %li & 'dzrcv' element 1 are both 0.",iarray+1);
                vmess("Placing 1 receiver at (%li,%li)",xrcv1[iarray],zrcv1[iarray]);
            }
        }
        for (iarray=0;iarray<Ndx;iarray++){
            if (xrcv1[iarray]==xrcv2[iarray] && dxr[iarray]!=0) {
                dxr[iarray]=0.;
                vwarn("Linear array %li: 'xrcv1'='xrcv2' and 'dxrcv' is not 0. Setting 'dxrcv'=0",iarray+1);
            }
        }
        for (iarray=0;iarray<Ndx;iarray++){
            if (yrcv1[iarray]==yrcv2[iarray] && dyr[iarray]!=0) {
                dyr[iarray]=0.;
                vwarn("Linear array %li: 'yrcv1'='yrcv2' and 'dyrcv' is not 0. Setting 'dyrcv'=0",iarray+1);
            }
        }
        for (iarray=0;iarray<Ndx;iarray++){
            if (zrcv1[iarray]==zrcv2[iarray] && dzr[iarray]!=0.){
                dzr[iarray]=0.;
                vwarn("Linear array %li: 'zrcv1'='zrcv2' and 'dzrcv' is not 0. Setting 'dzrcv'=0",iarray+1);
            }
        }

        /* Calculate Number of Receivers */
		nrcv = 0;
        nlxrcv=(long *)malloc(Nx1*sizeof(long));
        nlyrcv=(long *)malloc(Nx1*sizeof(long));
		for (iarray=0; iarray<Nx1; iarray++) {
			xrange = (xrcv2[iarray]-xrcv1[iarray]); 
			yrange = (yrcv2[iarray]-yrcv1[iarray]); 
			zrange = (zrcv2[iarray]-zrcv1[iarray]); 
			if (dxr[iarray] != 0.0 && dyr[iarray] != 0.0) {
				nlxrcv[iarray] = NINT(fabs(xrange/dxr[iarray]))+1;
				nlyrcv[iarray] = NINT(fabs(yrange/dyr[iarray]))+1;
			}
			else if (dxr[iarray] != 0.0) {
				nlxrcv[iarray] = NINT(fabs(xrange/dxr[iarray]))+1;
				nlyrcv[iarray] = 1;
			}
			else if (dyr[iarray] != 0.0) {
				nlxrcv[iarray] = 1;
				nlyrcv[iarray] = NINT(fabs(yrange/dyr[iarray]))+1;
			}
			else {
				if (dzr[iarray] == 0) {
					verr("For receiver array %li: receiver distance dzrcv is not given", iarray);
				}
				nlxrcv[iarray] = NINT(fabs(zrange/dzr[iarray]))+1;
				nlyrcv[iarray] = NINT(fabs(zrange/dzr[iarray]))+1;
			}
            nrcv+=nlyrcv[iarray]*nlxrcv[iarray];
		}

        /* Calculate Number of Receivers */
        if (verbose) vmess("Total number of linear array receivers: %li",nrcv);
        if (!nrcv) {
            free(xrcv1);
            free(xrcv2);
            free(yrcv1);
            free(yrcv2);
            free(zrcv1);
            free(zrcv2);
            free(dxr);
            free(dyr);
            free(dzr);
            free(nlxrcv);
            free(nlyrcv);
        }
        rec->max_nrec+=nrcv;
    } 
	else {
		nrcv=0;
	}

/* allocate the receiver arrays */

    /* Total Number of Receivers */
    if (verbose) vmess("Total number of receivers: %li",rec->max_nrec);

    /* Allocate Arrays */
    rec->x  = (long *)calloc(rec->max_nrec,sizeof(long));
    rec->y  = (long *)calloc(rec->max_nrec,sizeof(long));
    rec->z  = (long *)calloc(rec->max_nrec,sizeof(long));
    rec->xr = (float *)calloc(rec->max_nrec,sizeof(float));
    rec->yr = (float *)calloc(rec->max_nrec,sizeof(float));
    rec->zr = (float *)calloc(rec->max_nrec,sizeof(float));

/* read in the receiver postions */

	nrec=0;
    /* Receiver Array */
	if (nxrcv != 0 && nyrcv != 0) {
		/* receiver array is defined */
		xrcva = (float *)malloc(nxrcv*sizeof(float));
		yrcva = (float *)malloc(nxrcv*sizeof(float));
		zrcva = (float *)malloc(nxrcv*sizeof(float));
		getparfloat("xrcva", xrcva);
		getparfloat("yrcva", yrcva);
		getparfloat("zrcva", zrcva);
		for (iy=0; iy<nyrcv; iy++) {
            for (ix=0; ix<nxrcv; ix++) {
                rec->xr[nrec+iy*nxrcv+ix] = xrcva[ix]-sub_x0;
                rec->yr[nrec+iy*nxrcv+ix] = yrcva[iy]-sub_y0;
                rec->zr[nrec+iy*nxrcv+ix] = zrcva[ix]-sub_z0;
                rec->x[nrec+iy*nxrcv+ix] = NINT((xrcva[ix]-sub_x0)/dx);
                rec->y[nrec+iy*nxrcv+ix] = NINT((yrcva[iy]-sub_y0)/dy);
                rec->z[nrec+iy*nxrcv+ix] = NINT((zrcva[ix]-sub_z0)/dz);
                if (verbose>4) fprintf(stderr,"Receiver Array: xrcv[%li]=%f yrcv[%li]=%f zrcv=%f\n", ix, rec->xr[nrec+ix]+sub_x0, iy, rec->yr[nrec+ix]+sub_y0, rec->zr[nrec+ix]+sub_z0);
            }
        }
		free(xrcva);
		free(yrcva);
		free(zrcva);
		nrec += nyrcv*nxrcv;
	}

    /* Linear Receiver Arrays */
    if (nrcv!=0) {
		xrcv1 = (float *)malloc(Nx1*sizeof(float));
		xrcv2 = (float *)malloc(Nx1*sizeof(float));
		yrcv1 = (float *)malloc(Nx1*sizeof(float));
		yrcv2 = (float *)malloc(Nx1*sizeof(float));
		zrcv1 = (float *)malloc(Nx1*sizeof(float));
		zrcv2 = (float *)malloc(Nx1*sizeof(float));
		
		if(!getparfloat("xrcv1", xrcv1)) xrcv1[0]=sub_x0;
		if(!getparfloat("xrcv2", xrcv2)) xrcv2[0]=(nx-1)*dx+sub_x0;
		if(!getparfloat("yrcv1", yrcv1)) yrcv1[0]=sub_y0;
		if(!getparfloat("yrcv2", yrcv2)) yrcv2[0]=(ny-1)*dy+sub_y0;
		if(!getparfloat("zrcv1", zrcv1)) zrcv1[0]=sub_z0;
		if(!getparfloat("zrcv2", zrcv2)) zrcv2[0]=zrcv1[0];		
		
		Ndx = countparval("dxrcv");
		Ndy = countparval("dyrcv");
		Ndz = countparval("dzrcv");

		dxr = (float *)malloc(Nx1*sizeof(float));
		dyr = (float *)malloc(Nx1*sizeof(float));
		dzr = (float *)malloc(Nx1*sizeof(float));
		if(!getparfloat("dxrcv", dxr)) dxr[0]=dx;
		if(!getparfloat("dyrcv", dyr)) dyr[0]=dy;
		if(!getparfloat("dzrcv", dzr)) dzr[0]=0.0;
		if ( (Ndx<=1) && (Ndy<=1) && (Ndz==0) ){ /* default values are set */
			for (i=1; i<Nx1; i++) {
				dxr[i] = dxr[0];
				dyr[i] = dyr[0];
				dzr[i] = dzr[0];
			}
			Ndx=1;
            Ndy=1;
		}
        else if ( (Ndx<=1) && (Ndy==0) && (Ndz==0) ){ /* default values are set */
			for (i=1; i<Nx1; i++) {
				dxr[i] = dxr[0];
				dyr[i] = dyr[0];
				dzr[i] = dzr[0];
			}
			Ndx=1;
		}
        else if ( (Ndy<=1) && (Ndx==0) && (Ndz==0) ){ /* default values are set */
			for (i=1; i<Nx1; i++) {
				dxr[i] = dxr[0];
				dyr[i] = dyr[0];
				dzr[i] = dzr[0];
			}
			Ndy=1;
		}
		else if ( (Ndz==1) && (Ndy==0) && (Ndx==0) ){ /* default values are set */
			for (i=1; i<Nx1; i++) {
				dxr[i] = dxr[0];
				dyr[i] = dyr[0];
				dzr[i] = dzr[0];
			}
			Ndz=1;
		}
		else { /* make sure that each array has dzrcv or dxrcv defined for each line or receivers */
			if (Ndx>1) assert(Ndx==Nx1);
			if (Ndy>1) assert(Ndy==Ny1);
			if (Ndz>1) assert(Ndz==Nx1);
		}
		
		/* check if receiver arrays fit into model */
		for (iarray=0; iarray<Nx1; iarray++) {
			xrcv1[iarray] = MAX(sub_x0,      xrcv1[iarray]);
			xrcv1[iarray] = MIN(sub_x0+nx*dx,xrcv1[iarray]);
			xrcv2[iarray] = MAX(sub_x0,      xrcv2[iarray]);
			xrcv2[iarray] = MIN(sub_x0+nx*dx,xrcv2[iarray]);

			yrcv1[iarray] = MAX(sub_y0,      yrcv1[iarray]);
			yrcv1[iarray] = MIN(sub_y0+ny*dy,yrcv1[iarray]);
			yrcv2[iarray] = MAX(sub_y0,      yrcv2[iarray]);
			yrcv2[iarray] = MIN(sub_y0+ny*dy,yrcv2[iarray]);
			
			zrcv1[iarray] = MAX(sub_z0,      zrcv1[iarray]);
			zrcv1[iarray] = MIN(sub_z0+nz*dz,zrcv1[iarray]);
			zrcv2[iarray] = MAX(sub_z0,      zrcv2[iarray]);
			zrcv2[iarray] = MIN(sub_z0+nz*dz,zrcv2[iarray]);
		}

		/* calculate receiver array and store into rec->x,y,z */

		for (iarray=0; iarray<Nx1; iarray++) {
			xrange = (xrcv2[iarray]-xrcv1[iarray]); 
			yrange = (yrcv2[iarray]-yrcv1[iarray]); 
			zrange = (zrcv2[iarray]-zrcv1[iarray]); 
			if (dxr[iarray] != 0.0) {
				nrcv = nlyrcv[iarray]*nlxrcv[iarray];
				dxrcv = dxr[iarray];
				dyrcv = yrange/(nlyrcv[iarray]-1);
				dzrcv = zrange/(nlxrcv[iarray]-1);
        fprintf(stderr,"recvPar3D.c 0 dzrcv=%f zrange=%f %d %d %d\n", dzrcv, zrange, iarray, nlxrcv[iarray], isnan(dzrcv));
				if (dyrcv != dyr[iarray] && !isnan(dyrcv)) {
					vwarn("For receiver array %li: calculated dyrcv=%f given=%f", iarray, dyrcv, dyr[iarray]);
					vwarn("The calculated receiver distance %f is used", dyrcv);
				}
				if (dzrcv != dzr[iarray] && !isnan(dzrcv)) {
					vwarn("For receiver array %li: calculated dzrcv=%f given=%f", iarray, dzrcv, dzr[iarray]);
					vwarn("The calculated receiver distance %f is used", dzrcv);
				}
			}
            else if (dyr[iarray] != 0.0) {
				nrcv = nlyrcv[iarray]*nlxrcv[iarray];
				dxrcv = xrange/(nlxrcv[iarray]-1);
				dyrcv = dyr[iarray];
				dzrcv = zrange/(nlxrcv[iarray]-1);
				if (dxrcv != dxr[iarray] && !isnan(dxrcv)) {
					vwarn("For receiver array %li: calculated dxrcv=%f given=%f", iarray, dxrcv, dxr[iarray]);
					vwarn("The calculated receiver distance %f is used", dxrcv);
				}
				if (dzrcv != dzr[iarray] && !isnan(dzrcv)) {
					vwarn("For receiver array %li: calculated dzrcv=%f given=%f", iarray, dzrcv, dzr[iarray]);
					vwarn("The calculated receiver distance %f is used", dzrcv);
				}
			}
			else {
				if (dzr[iarray] == 0) {
					verr("For receiver array %li: receiver distance dzrcv is not given", iarray);
				}
				nrcv = nlyrcv[iarray]*nlxrcv[iarray];
				dxrcv = xrange/(nrcv-1);
				dyrcv = yrange/(nrcv-1);
				dzrcv = dzr[iarray];
				if (dxrcv != dxr[iarray]) {
					vwarn("For receiver array %li: calculated dxrcv=%f given=%f", iarray, dxrcv, dxr[iarray]);
					vwarn("The calculated receiver distance %f is used", dxrcv);
				}
				if (dyrcv != dyr[iarray]) {
					vwarn("For receiver array %li: calculated dyrcv=%f given=%f", iarray, dyrcv, dyr[iarray]);
					vwarn("The calculated receiver distance %f is used", dyrcv);
				}
			}

			// calculate coordinates
			for (iy=0; iy<nlyrcv[iarray]; iy++) {
                for (ix=0; ix<nlxrcv[iarray]; ix++) {
                    rec->xr[nrec]=xrcv1[iarray]-sub_x0+ix*dxrcv;
                    rec->yr[nrec]=yrcv1[iarray]-sub_y0+iy*dyrcv;
                    rec->zr[nrec]=zrcv1[iarray]-sub_z0+ix*dzrcv;

                    rec->x[nrec]=NINT((rec->xr[nrec])/dx);
                    rec->y[nrec]=NINT((rec->yr[nrec])/dy);
                    rec->z[nrec]=NINT((rec->zr[nrec])/dz);
                    nrec++;
                }
            }
		}
		free(xrcv1);
		free(xrcv2);
		free(yrcv1);
		free(yrcv2);
		free(zrcv1);
		free(zrcv2);
		free(dxr);
		free(dyr);
		free(dzr);
        free(nlxrcv);
        free(nlyrcv);
	}
        fprintf(stderr,"recvPar3D.c dxrcv=%f\n", dxrcv);
        fprintf(stderr,"recvPar3D.c dyrcv=%f\n", dyrcv);
        fprintf(stderr,"recvPar3D.c dzrcv=%f\n", dzrcv);


    rec->n=rec->max_nrec;
	return 0;
}
