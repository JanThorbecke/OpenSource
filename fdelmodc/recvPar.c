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
*
*   Ammendments:
*           Max Holicki changing the allocation receiver array (2-2016)
*           The Netherlands 
**/


void name_ext(char *filename, char *extension);

int recvPar(recPar *rec, float sub_x0, float sub_z0, float dx, float dz, int nx, int nz)
{
	float *xrcv1, *xrcv2, *zrcv1, *zrcv2;
	int   i, ix, ir, verbose;
	float dxrcv, dzrcv, *dxr, *dzr;
	float rrcv, dphi, oxrcv, ozrcv, arcv;
	double circ, h, a, b, e, s, xr, zr, dr, srun, phase;
	float xrange, zrange, sub_x1, sub_z1;
	int Nx1, Nx2, Nz1, Nz2, Ndx, Ndz, iarray, nrec, nh;
	int nxrcv, nzrcv, ncrcv, nrcv, ntrcv, *nlrcv;
	float *xrcva, *zrcva;
	char* rcv_txt;
	FILE *fp;

	if (!getparint("verbose", &verbose)) verbose = 0;

    /* Calculate Model Dimensions */
    sub_x1=sub_x0+(nx-1)*dx;
    sub_z1=sub_z0+(nz-1)*dz;

/* Compute how many receivers are defined and then allocate the receiver arrays */

    /* Receivers on a Circle */
    if (getparfloat("rrcv",&rrcv)) {
        if (!getparfloat("dphi",&dphi)) dphi=2.0;
        ncrcv=NINT(360.0/dphi);
        if (verbose) vmess("Total number of receivers on a circle: %d",ncrcv);
    } 
	else {
		ncrcv=0;
	}

    /* Receivers from a File */
    ntrcv=0;
    if (!getparstring("rcv_txt",&rcv_txt)) rcv_txt=NULL;
    if (rcv_txt!=NULL) {
        /* Open text file */
        fp=fopen(rcv_txt,"r");
        assert(fp!=NULL);
        /* Get number of lines */
        while (!feof(fp)) if (fgetc(fp)=='\n') ntrcv++;
        fseek(fp,-1,SEEK_CUR);
        if (fgetc(fp)!='\n') ntrcv++; /* Checks if last line terminated by /n */
        if (verbose) vmess("Number of receivers in rcv_txt file: %d",ntrcv);
        rewind(fp);
    }

    /* Receiver Array */
    nxrcv=countparval("xrcva");
    nzrcv=countparval("zrcva");
    if (nxrcv!=nzrcv) verr("Number of receivers in array xrcva (%d), zrcva(%d) are not equal",nxrcv,nzrcv);
    if (verbose&&nxrcv) vmess("Total number of array receivers: %d",nxrcv);

    /* Linear Receiver Arrays */
	Nx1 = countparval("xrcv1");
	Nx2 = countparval("xrcv2");
	Nz1 = countparval("zrcv1");
	Nz2 = countparval("zrcv2");
    if (Nx1!=Nx2) verr("Number of receivers starting points in 'xrcv1' (%d) and number of endpoint in 'xrcv2' (%d) are not equal",Nx1,Nx2);
    if (Nz1!=Nz2) verr("Number of receivers starting points in 'zrcv1' (%d) and number of endpoint in 'zrcv2' (%d) are not equal",Nz1,Nz2);
    if (Nx1!=Nz2) verr("Number of receivers starting points in 'xrcv1' (%d) and number of endpoint in 'zrcv2' (%d) are not equal",Nx1,Nz2);

    rec->max_nrec=ncrcv+ntrcv+nxrcv;

	/* no receivers are defined use default linear array of receivers on top of model */
    if (!rec->max_nrec && Nx1==0) Nx1=1; // Default is to use top of model to record data

    if (Nx1) {
        /* Allocate Start & End Points of Linear Arrays */
        xrcv1=(float *)malloc(Nx1*sizeof(float));
        xrcv2=(float *)malloc(Nx1*sizeof(float));
        zrcv1=(float *)malloc(Nx1*sizeof(float));
        zrcv2=(float *)malloc(Nx1*sizeof(float));
        if (!getparfloat("xrcv1",xrcv1)) xrcv1[0]=sub_x0;
        if (!getparfloat("xrcv2",xrcv2)) xrcv2[0]=sub_x1;
        if (!getparfloat("zrcv1",zrcv1)) zrcv1[0]=sub_z0;
        if (!getparfloat("zrcv2",zrcv2)) zrcv2[0]=zrcv1[0];

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

        /* Crop to Fit Model */
/* Max's addtion still have to check if it has the same fucntionality */
        for (iarray=0;iarray<Nx1;iarray++) {
            if (xrcv1[iarray]<sub_x0) {
                if (xrcv2[iarray]<sub_x0) {
                    verr("Linear array %d outside model bounds",iarray);
                }
				else {
                    vwarn("Cropping element %d of 'xrcv1' (%f) to model bounds (%f)",iarray,xrcv1[iarray],sub_x0);
                    xrcv1[iarray]=sub_x0;
                }
            } 
			else if (xrcv1[iarray] > sub_x1) {
                verr("Linear array %d outside model bounds",iarray);
            }
            if ( (xrcv2[iarray] < xrcv1[iarray]) ) {
                verr("Ill defined linear array %d, 'xrcv1' (%f) greater than 'xrcv2' (%f)",iarray,xrcv1[iarray],xrcv2[iarray]);
            }
			else if (xrcv2[iarray]>sub_x1) {
                vwarn("Cropping element %d of 'xrcv2' (%f) to model bounds (%f)",iarray,xrcv2[iarray],sub_x1);
                xrcv2[iarray]=sub_x1;
            }

            if (zrcv1[iarray] < sub_z0) {
                if (zrcv2[iarray] < sub_z0) {
                    verr("Linear array %d outside model bounds",iarray);
                }
				else {
               		vwarn("Cropping element %d of 'zrcv1' (%f) to model bounds (%f)",iarray,zrcv1[iarray],sub_z0);
                	zrcv1[iarray]=sub_z0;
                }
            }
			else if (zrcv1[iarray] > sub_z1) {
                verr("Linear array %d outside model bounds",iarray);
            }
            if ( (zrcv2[iarray] < zrcv1[iarray]) ) {
                verr("Ill defined linear array %d, 'zrcv1' (%f) greater than 'zrcv2' (%f)",iarray,zrcv1[iarray],zrcv2[iarray]);
            }
			else if (zrcv2[iarray]>sub_z1) {
                vwarn("Cropping element %d of 'xrcv2' (%f) to model bounds (%f)",iarray,zrcv2[iarray],sub_z1);
                zrcv2[iarray]=sub_z1;
            }
        }

        /* Get Sampling Rates */
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
			Ndz=1;
		}
		else if ( (Ndz==1) && (Ndx==0) ){ /* default values are set */
			for (i=1; i<Nx1; i++) {
				dxr[i] = dxr[0];
				dzr[i] = dzr[0];
			}
			Ndz=1;
			Ndx=1;
		}
		else { /* make sure that each array has dzrcv or dxrcv defined for each line or receivers */
			if (Ndx!=Ndz) {
				verr("Number of 'dxrcv' (%d) is not equal to number of 'dzrcv' (%d) or 1",Ndx,Ndz);
			}
			if (Ndx!=Nx1 && Ndx!=1) {
				verr("Number of 'dxrcv' (%d) is not equal to number of starting points in 'xrcv1' (%d) or 1",Ndx,Nx1);
			}
		}

		/* check consistency of receiver steps */
        for (iarray=0; iarray<Ndx; iarray++) {
            if (dxr[iarray]<0) {
				dxr[i]=dx;
				vwarn("'dxrcv' element %d (%f) is less than zero, changing it to %f'",iarray,dxr[iarray],dx);
			}
        }
        for (iarray=0;iarray<Ndz;iarray++) {
            if (dzr[iarray]<0) {
				dzr[iarray]=dz;
				vwarn("'dzrcv' element %d (%f) is less than zero, changing it to %f'",iarray,dzr[iarray],dz);
			}
        }
        for (iarray=0;iarray<Ndx;iarray++){
            if (dxr[iarray]==0 && dzr[iarray]==0) {
                xrcv2[iarray]=xrcv1[iarray];
				dxr[iarray]=1.;
                vwarn("'dxrcv' element %d & 'dzrcv' element 1 are both 0.",iarray+1);
                vmess("Placing 1 receiver at (%d,%d)",xrcv1[iarray],zrcv1[iarray]);
            }
        }
        for (iarray=0;iarray<Ndx;iarray++){
            if (xrcv1[iarray]==xrcv2[iarray] && dxr[iarray]!=0) {
                dxr[iarray]=0.;
                vwarn("Linear array %d: 'xrcv1'='xrcv2' and 'dxrcv' is not 0. Setting 'dxrcv'=0",iarray+1);
            }
        }
        for (iarray=0;iarray<Ndx;iarray++){
            if (zrcv1[iarray]==zrcv2[iarray] && dzr[iarray]!=0.){
                dzr[iarray]=0.;
                vwarn("Linear array %d: 'zrcv1'='zrcv2' and 'dzrcv' is not 0. Setting 'dzrcv'=0",iarray+1);
            }
        }

        /* Calculate Number of Receivers */
		nrcv = 0;
        nlrcv=(int *)malloc(Nx1*sizeof(int));
		for (iarray=0; iarray<Nx1; iarray++) {
			xrange = (xrcv2[iarray]-xrcv1[iarray]); 
			zrange = (zrcv2[iarray]-zrcv1[iarray]); 
			if (dxr[iarray] != 0.0) {
				nlrcv[iarray] = NINT(fabs(xrange/dxr[iarray]))+1;
			}
			else {
				if (dzr[iarray] == 0) {
					verr("For receiver array %d: receiver distance dzrcv is not given", iarray);
				}
				nlrcv[iarray] = NINT(fabs(zrange/dzr[iarray]))+1;
			}
            nrcv+=nlrcv[iarray];
		}

        /* Calculate Number of Receivers */
/*
        nlrcv=(int *)malloc(Nx1*sizeof(int));
        if (!isnan(*xrcv1)) *nlrcv=MIN(NINT((*xrcv2-*xrcv1)/(*dxr)),NINT((*zrcv2-*zrcv1)/(*dzr)))+1;
        else *nlrcv=0;
        nrcv=*nlrcv;
        if (verbose>4 && nlrcv[iarray]!=0) vmess("Linear receiver array 1 has final bounds: (X: %f -> %f,Z: %f ->
%f)",xrcv1[iarray],xrcv1[iarray]+nlrcv[iarray]*(*dxr),zrcv1[iarray],zrcv1[iarray]+nlrcv[iarray]*(*dzr));
        if (Ndx>1) {
            for (iarray=1;iarray<Nx1;iarray++) {
                if (!isnan(xrcv1[iarray])) {
					nlrcv[iarray]=MIN(NINT((xrcv2[iarray]-xrcv1[iarray])/dxr[iarray]),NINT((zrcv2[iarray]-zrcv1[iarray])/dzr[iarray]))+1;
				}
                else {
					nlrcv[iarray]=0;
				}
                nrcv+=nlrcv[iarray];
                if (verbose>4&&nlrcv[iarray]!=0) vmess("Linear receiver array %d has final bounds: (X: %f -> %f,Z: %f ->
%f)",iarray,xrcv1[iarray],xrcv1[iarray]+nlrcv[iarray]*dxr[iarray],zrcv1[iarray],zrcv1[iarray]+nlrcv[iarray]*dzr[iarray]);
            }
        }
		 else {
            for (iarray=1;iarray<Nx1;iarray++) {
                if (!isnan(xrcv1[iarray])) nlrcv[iarray]=MIN(NINT((xrcv2[iarray]-xrcv1[iarray])/(*dxr)),NINT((zrcv2[iarray]-zrcv1[iarray])/(*dzr)))+1;
                else nlrcv[iarray]=0;
                nrcv+=nlrcv[iarray];
                if (verbose>4&&nlrcv[iarray]!=0) vmess("Linear receiver array %d has final bounds: (X: %f -> %f,Z: %f ->
%f)",iarray,xrcv1[iarray],xrcv1[iarray]+nlrcv[iarray]**dxr,zrcv1[iarray],zrcv1[iarray]+nlrcv[iarray]**dzr);
            }
        }
*/
        if (verbose) vmess("Total number of linear array receivers: %d",nrcv);
        if (!nrcv) {
            free(xrcv1);
            free(xrcv2);
            free(zrcv1);
            free(zrcv2);
            free(dxr);
            free(dzr);
            free(nlrcv);
        }
        rec->max_nrec+=nrcv;
    } 
	else {
		nrcv=0;
	}

/* allocate the receiver arrays */

    /* Total Number of Receivers */
    if (verbose) vmess("Total number of receivers: %d",rec->max_nrec);

    /* Allocate Arrays */
    rec->x  = (int *)calloc(rec->max_nrec,sizeof(int));
    rec->z  = (int *)calloc(rec->max_nrec,sizeof(int));
    rec->xr = (float *)calloc(rec->max_nrec,sizeof(float));
    rec->zr = (float *)calloc(rec->max_nrec,sizeof(float));

/* read in the receiver postions */

	nrec=0;
    /* Receivers on a Circle */
    if (ncrcv) {
		if (!getparfloat("oxrcv",&oxrcv)) oxrcv=0.0;
		if (!getparfloat("ozrcv",&ozrcv)) ozrcv=0.0;
		if (!getparfloat("arcv",&arcv)) {
			arcv=rrcv; 
			for (ix=0; ix<ncrcv; ix++) {
				rec->xr[ix] = oxrcv-sub_x0+rrcv*cos(((ix*dphi)/360.0)*(2.0*M_PI));
				rec->zr[ix] = ozrcv-sub_z0+arcv*sin(((ix*dphi)/360.0)*(2.0*M_PI));
				rec->x[ix] = NINT(rec->xr[ix]/dx);
				rec->z[ix] = NINT(rec->zr[ix]/dz);
				//rec->x[ix] = NINT((oxrcv-sub_x0+rrcv*cos(((ix*dphi)/360.0)*(2.0*M_PI)))/dx);
				//rec->z[ix] = NINT((ozrcv-sub_z0+arcv*sin(((ix*dphi)/360.0)*(2.0*M_PI)))/dz);
				if (verbose>4) fprintf(stderr,"Receiver Circle: xrcv[%d]=%f zrcv=%f\n", ix, rec->xr[ix]+sub_x0, rec->zr[ix]+sub_z0);
			}
		}
		else { /* an ellipse */
			/* simple numerical solution to find equidistant points on an ellipse */
			nh  = (ncrcv)*1000; /* should be fine enough for most configurations */
			h = 2.0*M_PI/nh;
			a = MAX(rrcv, arcv);
			b = MIN(rrcv, arcv);
			e = sqrt(a*a-b*b)/a;
			//fprintf(stderr,"a=%f b=%f e=%f\n", a, b, e);
			circ = 0.0;
			for (ir=0; ir<nh; ir++) {
				s = sin(ir*h);
				circ += sqrt(1.0-e*e*s*s);
			}
			circ = a*h*circ;
			//fprintf(stderr,"circ = %f circle=%f\n", circ, 2.0*M_PI*rrcv);
			/* define distance between receivers on ellipse */
			dr = circ/ncrcv;
			ix = 0;
			srun = 0.0;
			if (arcv >= rrcv) phase=0.0;
			else phase=0.5*M_PI;
			for (ir=0; ir<nh; ir++) {
				s = sin(ir*h);
				srun += sqrt(1.0-e*e*s*s);
				if (a*h*srun >= ix*dr ) {
					xr = rrcv*cos(ir*h+phase);
					zr = arcv*sin(ir*h+phase);
					rec->xr[ix] = oxrcv-sub_x0+xr;
					rec->zr[ix] = ozrcv-sub_z0+zr;
					rec->x[ix] = NINT(rec->xr[ix]/dx);
					rec->z[ix] = NINT(rec->zr[ix]/dz);
					if (verbose>4) fprintf(stderr,"Receiver Ellipse: xrcv[%d]=%f zrcv=%f\n", ix, rec->xr[ix]+sub_x0, rec->zr[ix]+sub_z0);
					ix++;
				}
				if (ix == ncrcv) break;
			}
		}

		/* check if receivers fit into the model otherwise clip to edges */
		for (ix=0; ix<ncrcv; ix++) {
			rec->x[ix] = MIN(nx-1, MAX(rec->x[ix], 0));
			rec->z[ix] = MIN(nz-1, MAX(rec->z[ix], 0));
		}
		nrec += ncrcv;
	}

    /* Receiver Text File */

    if (ntrcv) {
		/* Allocate arrays */
		xrcva = (float *)malloc(ntrcv*sizeof(float));
		zrcva = (float *)malloc(ntrcv*sizeof(float));
		/* Read in receiver coordinates */
		for (i=0;i<ntrcv;i++) {
			if (fscanf(fp,"%e %e\n",&xrcva[i],&zrcva[i])!=2) vmess("Receiver Text File: Can not parse coordinates on line %d.",i);
		}
		/* Close file */
		fclose(fp);
		/* Process coordinates */
		for (ix=0; ix<ntrcv; ix++) {
			rec->xr[nrec+ix] = xrcva[ix]-sub_x0;
			rec->zr[nrec+ix] = zrcva[ix]-sub_z0;
			rec->x[nrec+ix] = NINT((xrcva[ix]-sub_x0)/dx);
			rec->z[nrec+ix] = NINT((zrcva[ix]-sub_z0)/dz);
			if (verbose>4) vmess("Receiver Text Array: xrcv[%d]=%f zrcv=%f", ix, rec->xr[nrec+ix]+sub_x0, rec->zr[nrec+ix]+sub_z0);
		}
		free(xrcva);
		free(zrcva);
		nrec += ntrcv;
	}

    /* Receiver Array */
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
		free(xrcva);
		free(zrcva);
		nrec += nxrcv;
	}

    /* Linear Receiver Arrays */
    if (nrcv!=0) {
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
				nrcv = nlrcv[iarray];
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
				nrcv = nlrcv[iarray];
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
		}
		free(xrcv1);
		free(xrcv2);
		free(zrcv1);
		free(zrcv2);
		free(dxr);
		free(dzr);
        free(nlrcv);
	}

    rec->n=rec->max_nrec;
	return 0;
}
