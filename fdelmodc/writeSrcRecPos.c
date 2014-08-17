#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"par.h"
#include"fdelmodc.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

/**
* Writes the source and receiver positions into a gridded file,
* which has the same size as the input gridded model files. 
* Source positions have a value +1 and receivers -1.
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

int writesufile(char *filename, float *data, int n1, int n2, float f1, float f2, float d1, float d2);

int writeSrcRecPos(modPar *mod, recPar *rec, srcPar *src, shotPar *shot)
{
	FILE *fp;
	float *dum, sub_x0, sub_z0, dx, dz;
	int is, nx, nz, is0, ish, ix, iz, ndot, idx, idz;
	char tmpname[1024];

 	ndot = 2;
	nx = mod->nx;
	nz = mod->nz;
	dx = mod->dx;
	dz = mod->dz;
	sub_x0 = mod->x0;
	sub_z0 = mod->z0;

	/* write velocity field with positions of the sources */
	dum = (float *)calloc(nx*nz, sizeof(float));
	vmess("Positions: shot=%d src=%d rec=%d", shot->n, src->n, rec->n);
	/* source positions for random shots */
	if (src->random) {
		sprintf(tmpname,"SrcPositions%d.txt",src->n);
		fp = fopen(tmpname, "w+");
		for (is=0; is<src->n; is++) {
			for (idx=0; idx<=ndot; idx++) {
				for (idz=0; idz<=ndot; idz++) {
					dum[(MAX(0,src->x[is]-idx))*nz+MAX(0,src->z[is]-idz)] = 1.0;
					dum[(MAX(0,src->x[is]-idx))*nz+MIN(nz-1,src->z[is]+idz)] = 1.0;
					dum[(MIN(nx-1,src->x[is]+idx))*nz+MIN(nz-1,src->z[is]+idz)] = 1.0;
					dum[(MIN(nx-1,src->x[is]+idx))*nz+MAX(0,src->z[is]-idz)] = 1.0;
				}
			}
			fprintf(fp, "%f %f\n", src->z[is]*dz+sub_z0, src->x[is]*dx+sub_x0);
		}
		fclose(fp);
	}
	/* source positions for single shot sources with plane waves */
	else if (src->plane) {
    	is0 = -1*floor((src->n-1)/2);
		sprintf(tmpname,"SrcPositions%d.txt",shot->n);
		fp = fopen(tmpname, "w+");
		for (ish=0; ish<shot->n; ish++) {
			for (is=0; is<src->n; is++) {
				ix = shot->x[ish] + 1 + is0 + is;
				iz = shot->z[ish] + 1;
				dum[ix*nz+iz] = 1.0;
				dum[(MAX(0,ix-1))*nz+iz] = 1.0;
				dum[(MIN(nx-1,ix+1))*nz+iz] = 1.0;
				dum[ix*nz+MAX(0,iz-1)] = 1.0;
				dum[ix*nz+MIN(nz-1,iz+1)] = 1.0;
				fprintf(fp, "(%f, %f)\n", ix*dx+sub_x0, iz*dz+sub_z0);
			}
		}
		fclose(fp);
	}
	else {
		sprintf(tmpname,"SrcPositions%d.txt",shot->n);
		fp = fopen(tmpname, "w+");
		for (is=0; is<shot->n; is++) {
			for (idx=0; idx<=ndot; idx++) {
				for (idz=0; idz<=ndot; idz++) {
					dum[(MAX(0,shot->x[is]-idx))*nz+MAX(0,shot->z[is]-idz)] = 1.0;
					dum[(MAX(0,shot->x[is]-idx))*nz+MIN(nz-1,shot->z[is]+idz)] = 1.0;
					dum[(MIN(nx-1,shot->x[is]+idx))*nz+MIN(nz-1,shot->z[is]+idz)] = 1.0;
					dum[(MIN(nx-1,shot->x[is]+idx))*nz+MAX(0,shot->z[is]-idz)] = 1.0;
				}
			}
			fprintf(fp, "%f %f\n", shot->z[is]*dz+sub_z0, shot->x[is]*dx+sub_x0);
		}
		fclose(fp);
	}

	/* receiver positions */
	sprintf(tmpname,"RcvPositions%d.txt",rec->n);
	fp = fopen(tmpname, "w+");
	for (is=0; is<rec->n; is++) {
		dum[rec->x[is]*nz+rec->z[is]] = -1.0;
		dum[(MAX(0,rec->x[is]-1))*nz+rec->z[is]] = -1.0;
		dum[(MIN(nx-1,rec->x[is]+1))*nz+rec->z[is]] = -1.0;
		dum[rec->x[is]*nz+MAX(0,rec->z[is]-1)] = -1.0;
		dum[rec->x[is]*nz+MIN(nz-1,rec->z[is]+1)] = -1.0;
		if (rec->int_vx==3) {
			fprintf(fp, "(%f, %f)\n", rec->xr[is]*dx+sub_x0, rec->zr[is]*dz+sub_z0);
		}
		else {
			fprintf(fp, "(%f, %f)\n", rec->x[is]*dx+sub_x0, rec->z[is]*dz+sub_z0);
		}
	}
	fclose(fp);
	writesufile("SrcRecPositions.su", dum, nz, nx, sub_z0, sub_x0, dz, dx);
	free(dum);

	return 0;
}
