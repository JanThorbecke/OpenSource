#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include <fcntl.h>

// #define t0(x,y,z)   time0[nxz*(y) + nz*(x) + (z)]
// #define s0(x,y,z)   slow0[nxz*(y) + nz*(x) + (z)]
#define t0(x,y,z)   time0[nxy*(z) + nx*(y) + (x)]
#define s0(x,y,z)   slow0[nxy*(z) + nx*(y) + (x)]
#define	SQR(x)	((x) * (x))
#define	DIST(x,y,z,x1,y1,z1)	sqrt(SQR(x-(x1))+SQR(y-(y1)) + SQR(z-(z1)))

/* definitions from verbose.c */
extern void verr(char *fmt, ...);
extern void vwarn(char *fmt, ...);
extern void vmess(char *fmt, ...);

void src3D(float *time0, float *slow0, long nz, long nx, long ny, float h, float fxs, float fys, float fzs, long xs, long ys, long zs, long *cube)
{
	long
		srctype=1,	/* if 1, source is a point;
						2, source is on the walls of the data volume;
						3, source on wall, time field known; */
		srcwall,	/* if 1, source on x=0 wall, if 2, on x=nx-1 wall
						if 3, source on y=0 wall, if 4, on y=ny-1 wall
						if 5, source on z=0 wall, if 6, on z=nz-1 wall */
		xs1,			/* shot x position (in grid points) */
		ys1,			/* shot y position */
		zs1,			/* shot depth */
		xx, yy, zz,	/* Used to loop around xs, ys, zs coordinates	*/
		ii, i, j, k, 
		wfint, ofint,
		nxy, nyz, nxz, nxyz, nwall,
		NCUBE=2;
	float
		fxs1,	/* shot position in X (in real units)*/
		fys1,	/* shot position in Y (in real units)*/
		fzs1,	/* shot position in Z (in real units)*/
		*wall,
		/* maximum offset (real units) to compute */
		/* used in linear velocity gradient cube source */
		rx, ry, rz, dvz, dv, v0,
		rzc, rxyc, rz1, rxy1, rho, theta1, theta2,
		xsrc1, ysrc1, zsrc1;
	char
		*oldtfile,	/* file through which old travel times are input */
		*wallfile;   /* file containing input wall values of traveltimes */


	if(!getparlong("NCUBE",&NCUBE)) NCUBE=2;
	*cube = NCUBE;

	if(!getparlong("srctype",&srctype)) srctype=1;

	nxy = nx * ny;
	nyz = ny * nz;
	nxz = nx * nz;
	nxyz = nx * ny * nz;

	/* SET TIMES TO DUMMY VALUE */
	for(i=0;i<nxyz;i++) time0[i] = 1.0e10;

	if (srctype == 1) {			/*  VIDALE'S POINT SOURCE */
		/* FILL IN CUBE AROUND SOURCE POINT */
		/* HOLE'S NEW LINEAR VELOCITY GRADIENT CUBE (APRIL 1991)*/
		v0 = h/s0(xs,ys,zs);
		for (xx = xs-NCUBE; xx <= xs+NCUBE; xx++) {
			if (xx < 0 || xx >= nx)	continue; 
			for (yy = ys-NCUBE; yy <= ys+NCUBE; yy++) {
				if (yy < 0 || yy >= ny)	continue; 
				for (zz = zs-NCUBE; zz <= zs+NCUBE; zz++) {
					if (zz < 0 || zz >= nz)	continue; 
					if (zz == zs)
					  dvz = 1/s0(xx,yy,zz+1)-1/s0(xs,ys,zs);
					else
					  dvz = (1/s0(xx,yy,zz)-1/s0(xs,ys,zs))/(zz-zs);
					dv = fabs(dvz);
					if (dv == 0.)  {
					  t0(xx,yy,zz) = s0(xs,ys,zs)*DIST(xs,ys,zs,xx,yy,zz);
					  continue;
					}
					rzc = -v0/dv;
					rx = h*(xx - xs);
					ry = h*(yy - ys);
					rz = h*(zz - zs);
					rz1 = rz*dvz/dv;
					rxy1 = sqrt(rx*rx+ry*ry+rz*rz-rz1*rz1);
					if (rxy1<=h/1.e6)
					  t0(xx,yy,zz) = fabs(log((v0+dv*rz1)/v0)/dv);
					else {
					  rxyc = (rz1*rz1+rxy1*rxy1-2*rz1*rzc)/(2*rxy1);
					  rho = sqrt(rzc*rzc+rxyc*rxyc);
					  theta1 = asin(-rzc/rho);
					  /* can't handle asin(1.) ! */
					  if (fabs(rz1-rzc)>=rho)  rho=1.0000001*fabs(rz1-rzc);
					  theta2 = asin((rz1-rzc)/rho);
					  if (rxyc<0) theta1=M_PI-theta1;
					  if (rxyc<rxy1) theta2=M_PI-theta2;
					  t0(xx,yy,zz) = log(tan(theta2/2)/tan(theta1/2)) / dv;
				        }
				}
			}
		}
	}
	else if (srctype == 2) {		/*  HOLE'S EXTERNAL SOURCE */

		/* FILL IN WALLS' TIMES FROM EXTERNAL DATAFILE */
		read (wfint,wall,4*nwall);	/* READ X=0 WALL */
		if (wall[0]>-1.e-20) {
			ii = 0;
			for (k=0; k<nz; k++) {
				for (j=0; j<ny; j++) {
					t0(0,j,k) = wall[ii];
					ii++;
				}
			}
		}
		read (wfint,wall,4*nwall);	/* READ X=NX-1 WALL */
		if (wall[0]>-1.e-20) {
			ii = 0;
			for (k=0; k<nz; k++) {
				for (j=0; j<ny; j++) {
					t0(nx-1,j,k) = wall[ii];
					ii++;
				}
			}
		}
		read (wfint,wall,4*nwall);	/* READ Y=0 WALL */
		if (wall[0]>-1.e-20) {
			ii = 0;
			for (k=0; k<nz; k++) {
				for (i=0; i<nx; i++) {
					t0(i,0,k) = wall[ii];
					ii++;
				}
			}
		}
		read (wfint,wall,4*nwall);	/* READ Y=NY-1 WALL */
		if (wall[0]>-1.e-20) {
			ii = 0;
			for (k=0; k<nz; k++) {
				for (i=0; i<nx; i++) {
					t0(i,ny-1,k) = wall[ii];
					ii++;
				}
			}
		}
		read (wfint,wall,4*nwall);	/* READ Z=0 WALL */
		if (wall[0]>-1.e-20) {
			ii = 0;
			for (j=0; j<ny; j++) {
				for (i=0; i<nx; i++) {
					t0(i,j,0) = wall[ii];
					ii++;
				}
			}
		}
		read (wfint,wall,4*nwall);	/* READ Z=NZ-1 WALL */
		if (wall[0]>-1.e-20) {
			ii = 0;
			for (j=0; j<ny; j++) {
				for (i=0; i<nx; i++) {
					t0(i,j,nz-1) = wall[ii];
					ii++;
				}
			}
		}
	}
	else if (srctype == 3) {                /*  HOLE'S REDO OLD TIMES */
	        /* READ IN OLD TIME FILE */
	        if (srctype == 3)  read(ofint,time0,nxyz*4);
	}

	return;
}
