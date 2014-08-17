#include <assert.h>
#include "segy.h"
#include "optim.h"

/* To Do add zi function of depth as in getrecvsp */

int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);

void getrecextr(int *xi, int *zi, int *nrec, int nx, int nz, float dx, float dz, float ox, float oz,  int *id0, int *id1, int *dstep, int verbose)
{
	FILE    *fp;
	size_t  nread;
	int     error, n1, n2;
	int		nrx, i, j, ndeltx, np, lint;
	int     nrz, ndeltz, nxm, ngath;
	int     boundary, ix, size;
	float   f1, f2, d1, d2, *tmpdata;
	float	xprev, deltx, zprev, deltz;
	float   zrecmin, zrecmax, zstart;
	float   xrcv1, xrcv2, zrcv1, zrcv2;
	float	*xrcv, *zrcv, dxrcv, dzrcv;
	float   xmin, xmax, scl;
	char    *file_int;
	segy	hdr;

/* Determine the lateral receiver positions */

	if(!getparint("lint", &lint)) lint = 1;
	if(!getparfloat("dxrcv",&dxrcv)) dxrcv = dx;
	if(!getparfloat("dzrcv",&dzrcv)) dzrcv = 0.0;
	if (fabs(dxrcv/dx - NINT(dxrcv/dx)) > 0.001) 
		fprintf(stderr,"    getrecextr: dxrcv does not fit on grid; rounded to nearest position\n");
	
	nrx = countparval("xrcv");
	*nrec = 0;
	if(nrx == 0) {
    	if(!getparfloat("xrcv1", &xrcv1)) xrcv1=ox;
    	if(!getparfloat("xrcv2", &xrcv2)) xrcv2=ox+(nx-1)*dx;
		if (xrcv1<ox) {
			fprintf(stderr,"xrcv1 defined outside the model !!\n");
			xrcv1 = ox;
		}

		if (xrcv1 < xrcv2) {
			if (dxrcv != 0.0) *nrec = NINT((xrcv2-xrcv1)/dxrcv)+1;
			else *nrec = 0;
		}
		for (i = 0; i < *nrec; i++) {
			xi[i] = NINT((xrcv1-ox+i*dxrcv)/dx);
			if (xi[i] >= nx) { 
				*nrec = i-1;
				fprintf(stderr,"xrcv2 defined outside the model !!\n");
				break;
			}
		}
	}
	else if (nrx == 1) {
		xrcv = (float *)malloc(nrx*sizeof(float));
		getparfloat("xrcv", xrcv);
		if (dzrcv == 0.0 ) *nrec = 1;
		else *nrec = NINT((nz-1)*dz/dzrcv)+1;
		for (i = 0; i < *nrec; i++) xi[i] = NINT((xrcv[0]-ox)/dx);
		free(xrcv);
	}
	else {
		xrcv = (float *)malloc(nrx*sizeof(float));
		getparfloat("xrcv",xrcv);
		for (i = 0; i < nrx; i++) xrcv[i] -= ox;
	
		if (lint == 1) {
			if(verbose) fprintf(stderr,"    getrecextr: interpolating receivers with interval %f between positions defined by xrcv\n",dxrcv);
			xprev = xrcv[0];
			np = 0;
			for (i = 1; i < nrx; i++) {
				deltx = xrcv[i] - xprev;
				ndeltx = NINT(ABS(deltx/dxrcv));
	
				for (j = 0; j <= ndeltx; j++) {
					xi[np++] = NINT((xprev + (j*dxrcv*deltx)/deltx)/dx);
				}
				xprev = xrcv[i];
				np--;
			}
			xi[np] = NINT(xrcv[nrx-1]/dx);
			*nrec = np+1;
		}
		else {
			if(verbose) fprintf(stderr,"    getrecextr: receivers at positions defined by xrcv\n");
			for (i = 0; i < nrx; i++) xi[i] = NINT(xrcv[i]/dx);
			*nrec = nrx;
		}
		free(xrcv);
	}

	
/* Determine the depth receiver positions */

	if(!getparstring("file_int", &file_int)) file_int=NULL;

	if (file_int == NULL) {
		nrz = countparval("zrcv");
		if(nrz == 0) {
    		if(!getparfloat("zrcv1", &zrcv1)) zrcv1=oz+(*id0)*dz;
    		if(!getparfloat("zrcv2", &zrcv2)) zrcv2=zrcv1;
			if (zrcv1<oz) {
				fprintf(stderr,"zrcv1 defined outside the model !!\n");
				zrcv1 = oz;
			}
			if (zrcv1 < zrcv2) {
				if (dzrcv != 0.0) *nrec = NINT((zrcv2-zrcv1)/dzrcv)+1;
				else dzrcv = (zrcv2-zrcv1)/(*nrec-1);
				for (i = 0; i < *nrec; i++) {
					zi[i] = NINT((zrcv1-oz+i*dzrcv)/dz);
					if (zi[i] > nz) { 
						*nrec = i-1;
						fprintf(stderr,"zrcv2 defined outside the model !!\n");
						break;
					}
				}
			}
			else {
				if(verbose) fprintf(stderr,"    getrecextr: no zrcv specified: receivers at z = %.2f\n", oz+*id0*dz);
				for (i = 0; i < *nrec; i++) zi[i] = *id0;
			}
		}
		else if (nrz == 1) {
			zrcv = (float *)malloc(nrz*sizeof(float));
			getparfloat("zrcv", zrcv);
			for (i = 0; i < *nrec; i++) zi[i] = NINT((zrcv[0]-oz)/dz);
			free(zrcv);
		}
		else {
			zrcv = (float *)malloc(nrz*sizeof(float));
			getparfloat("zrcv",zrcv);
			for (i = 0; i < nrz; i++) zrcv[i] -= oz;

			if (lint == 1) {
				if(verbose) fprintf(stderr,"    getrecextr: interpolating receivers with interval %f between positions defined by zrcv\n",dzrcv);
				zprev = zrcv[0];
				np = 0;
				for (i = 1; i < nrz; i++) {
					deltz = zrcv[i] - zprev;
					if (dzrcv == 0) ndeltz = *nrec;
					else ndeltz = NINT(ABS(deltz/dzrcv));
		
					for (j = 0; j <= ndeltz; j++) {
						if (NINT(deltz*1e4) != 0) zi[np++] = NINT((zprev + (j*dzrcv*deltz)/deltz)/dz);
						else zi[np++] = NINT(zprev/dz);
					}
					zprev = zrcv[i];
					np--;
				}
				zi[np] = NINT(zrcv[nrz-1]/dz);
			}
			else {
				if(verbose) fprintf(stderr,"    getrecextr: receivers at positions defined by zrcv\n");
				for (i = 0; i < nrz; i++) zi[i] = NINT(zrcv[i]/dz);
			}
			free(zrcv);
		}
	}
	else {

		if(!getparint("boundary", &boundary)) boundary=1;
        if (verbose) fprintf(stderr,"    getrecextr: receiver definition on boundary\n");

		boundary -= 1;
		error = getFileInfo(file_int, &n1, &n2, &ngath, &d1, &d2, &f1, &f2, &xmin, &xmax, &scl, &nxm);
		assert (error >= 0 );
		size = n1*n2;
		tmpdata = (float *)malloc(size*sizeof(float));
		
		fp = fopen( file_int, "r" );
		for (i=0; i<n2; i++) {
			nread = fread( &hdr, 1, TRCBYTES, fp );
			assert(nread == TRCBYTES);
			nread = fread( &tmpdata[i*n1], sizeof(float), hdr.ns, fp );
			assert(nread == hdr.ns);
		}
		fclose(fp);

		fprintf(stderr,"    getrecextr: assuming samples in interface file represent the x-axis\n");

		assert (n1 == nx);
		for(ix=0; ix<*nrec; ix++) 
			zi[ix] = NINT(tmpdata[boundary*nx+xi[ix]]/dz);
		free(tmpdata);
	}

	if (verbose>=3) {
		fprintf(stderr,"receivers are positioned at: \n");
		for (i=0; i<*nrec; i++) 
			fprintf(stderr," receiver no. %d at gridpoint: x=%d z=%d \n", i, xi[i], zi[i]);
	}

/* Check if zi fits in model */

	for(ix=0; ix<*nrec; ix++) {
		if (zi[ix] > nz) {
			fprintf(stderr,"z[%d]=%d defined outside the model !!\n", ix, zi[ix]);
			fprintf(stderr,"position reset to maximim depth in model.\n");
			zi[ix] = nz-1;
		}
	}

	if(!getparfloat("zstart",&zstart)) zstart = 0.0;
	*id0 = NINT(zstart/dz);
	zrecmin = zrecmax = zi[0];
	for (i=1; i<*nrec; i++) {
		zrecmin = MIN(zrecmin ,zi[i]);
		zrecmax = MAX(zrecmax ,zi[i]);
	}
	if (*id0 <= zrecmin) {
		*id1 = NINT(zrecmax);
		*dstep = 1;
	}
	else {
		*id1 = NINT(zrecmin);
		*dstep = -1;
	}

	return;
}
