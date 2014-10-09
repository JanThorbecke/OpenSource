#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc.h"

/**
*  Stores the wavefield at the receiver positions.
*
*  On a staggered grid the fields are all on different positions, 
*  to compensate for that the rec.int_vx and rec.int_vz options
*  can be set.
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

int getRecTimes(modPar mod, recPar rec, bndPar bnd, int itime, int isam, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rec_vx, float *rec_vz, float *rec_txx, float *rec_tzz, float *rec_txz, float *rec_p, float *rec_pp, float *rec_ss, float *rec_udp, float *rec_udvz, int verbose)
{
	int n1, ibndx, ibndz;
	int irec, ix, iz, ix2, iz2, ix1, iz1;
	float rdz, rdx, C00, C10, C01, C11;
	float *vz_t, c1, c2, roz;

    ibndx = mod.ioPx;
    ibndz = mod.ioPz;
    if (bnd.lef==4 || bnd.lef==2) ibndx += bnd.ntap;
    if (bnd.top==4 || bnd.top==2) ibndz += bnd.ntap;
	n1    = mod.naz;

	if (!rec.n) return 0;

/***********************************************************************
* velocity or txz or potential registrations issues:
* rec_x and rec_z are related to actual txx/tzz/p positions.
* offsets from virtual boundaries must be taken into account.
*
* vx velocities have one sample less in x-direction
* vz velocities have one sample less in z-direction
* txz stresses have one sample less in z-direction and x-direction
*
* Note, in the acoustic scheme P is stored in the Tzz array.
***********************************************************************/

	for (irec=0; irec<rec.n; irec++) {
		iz = rec.z[irec]+ibndz;
		ix = rec.x[irec]+ibndx;
		iz1 = iz-1;
		ix1 = ix-1;
		iz2 = iz+1;
		ix2 = ix+1;
		/* interpolation to precise (not necessary on a grid point) position */
		if ((rec.int_vx==3) || (rec.int_vz==3)) {

			iz = (int)floorf(rec.zr[irec]/mod.dz)+ibndz;
			ix = (int)floorf(rec.xr[irec]/mod.dx)+ibndx;
			rdz = (rec.zr[irec] - (iz-ibndz)*mod.dz)/mod.dz;
			rdx = (rec.xr[irec] - (ix-ibndx)*mod.dx)/mod.dx;
			iz1 = iz-1;
			ix1 = ix-1;
			iz2 = iz+1;
			ix2 = ix+1;
			
			/*
			 // Interpolate according to Dirk Kraaijpool's scheme 
			 // Reference:  "Seismic ray fields and ray field maps : theory and algorithms" , 
			 // PhD thesis Utrecht University,Faculty of Geosciences, 2003) 
			 
			 C00 = tzz[ix*n1+iz]      + 0.5*((tzz[(ix+1)*n1+iz]   +tzz[(ix-1)*n1+iz]+ 
			 tzz[(ix  )*n1+iz+1] +tzz[(ix  )*n1+iz-1])/(2.0*mod.dx));
			 C10 = tzz[(ix+1)*n1+iz]  + 0.5*((tzz[(ix+2)*n1+iz]   +tzz[(ix  )*n1+iz]+
			 tzz[(ix+1)*n1+iz+1] +tzz[(ix+1)*n1+iz-1])/(2.0*mod.dz));
			 C01 = tzz[ix*n1+iz+1]    + 0.5*((tzz[(ix+1)*n1+iz+1] +tzz[(ix-1)*n1+iz+1]+
			 tzz[(ix)*n1+iz+2]   +tzz[(ix  )*n1+iz])/(2.0*mod.dx));
			 C11 = tzz[(ix+1)*n1+iz+1]+ 0.5*((tzz[(ix+2)*n1+iz+1] +tzz[(ix  )*n1+iz+1]+
			 tzz[(ix+1)*n1+iz+2] +tzz[(ix+1)*n1+iz])/(2.0*mod.dz));
			 */
			
			if (rec.type.p){
				/* bi-linear interpolation */
				C00 = tzz[ix*n1+iz];
				C10 = tzz[(ix+1)*n1+iz];
				C01 = tzz[ix*n1+iz+1];
				C11 = tzz[(ix+1)*n1+iz+1];
				rec_p[irec*rec.nt+isam] = C00*(1.0-rdx)*(1.0-rdz) + C10*rdx*(1.0-rdz) +
										  C01*(1.0-rdx)*rdz       + C11*rdx*rdz;
			}
			if (rec.type.txx) {
				C00 = txx[ix*n1+iz];
				C10 = txx[(ix+1)*n1+iz];
				C01 = txx[ix*n1+iz+1];
				C11 = txx[(ix+1)*n1+iz+1];
				rec_txx[irec*rec.nt+isam] = C00*(1.0-rdx)*(1.0-rdz) + C10*rdx*(1.0-rdz) +
											C01*(1.0-rdx)*rdz       + C11*rdx*rdz;
			}
			if (rec.type.tzz) {
				C00 = tzz[ix*n1+iz];
				C10 = tzz[(ix+1)*n1+iz];
				C01 = tzz[ix*n1+iz+1];
				C11 = tzz[(ix+1)*n1+iz+1];
				rec_tzz[irec*rec.nt+isam] = C00*(1.0-rdx)*(1.0-rdz) + C10*rdx*(1.0-rdz) +
											C01*(1.0-rdx)*rdz       + C11*rdx*rdz;
			}
			if (rec.type.txz) {
				C00 = txz[ix2*n1+iz2];
				C10 = txz[(ix2+1)*n1+iz2];
				C01 = txz[ix2*n1+iz2+1];
				C11 = txz[(ix2+1)*n1+iz2+1];
				rec_txz[irec*rec.nt+isam] = C00*(1.0-rdx)*(1.0-rdz) + C10*rdx*(1.0-rdz) +
											C01*(1.0-rdx)*rdz       + C11*rdx*rdz;
			}
			if (rec.type.pp) {
				C00 = (vx[ix2*n1+iz]-vx[ix*n1+iz] +
					   vz[ix*n1+iz2]-vz[ix*n1+iz])/mod.dx;
				C10 = (vx[(ix2+1)*n1+iz]-vx[(ix+1)*n1+iz] +
					   vz[(ix+1)*n1+iz2]-vz[(ix+1)*n1+iz])/mod.dx;
				C01 = (vx[ix2*n1+iz+1]-vx[ix*n1+iz+1] +
					   vz[ix*n1+iz2+1]-vz[ix*n1+iz+1])/mod.dx;
				C11 = (vx[(ix2+1)*n1+iz+1]-vx[(ix+1)*n1+iz+1] +
					   vz[(ix+1)*n1+iz2+1]-vz[(ix+1)*n1+iz+1])/mod.dx;
				rec_pp[irec*rec.nt+isam] = C00*(1.0-rdx)*(1.0-rdz) + C10*rdx*(1.0-rdz) +
										   C01*(1.0-rdx)*rdz       + C11*rdx*rdz;
			}
			if (rec.type.ss) {
				C00 = (vx[ix2*n1+iz2]-vx[ix2*n1+iz] -
					   (vz[ix2*n1+iz2]-vz[ix*n1+iz2]))/mod.dx;
				C10 = (vx[(ix2+1)*n1+iz2]-vx[(ix2+1)*n1+iz] -
						(vz[(ix2+1)*n1+iz2]-vz[(ix+1)*n1+iz2]))/mod.dx;
				C01 = (vx[ix2*n1+iz2+1]-vx[ix2*n1+iz+1] -
						(vz[ix2*n1+iz2+1]-vz[ix*n1+iz2+1]))/mod.dx;;
				C11 = (vx[(ix2+1)*n1+iz2+1]-vx[(ix2+1)*n1+iz+1] -
						(vz[(ix2+1)*n1+iz2+1]-vz[(ix+1)*n1+iz2+1]))/mod.dx;
				rec_ss[irec*rec.nt+isam] = C00*(1.0-rdx)*(1.0-rdz) + C10*rdx*(1.0-rdz) +
										   C01*(1.0-rdx)*rdz       + C11*rdx*rdz;
			}
			if (rec.type.vz) {
				C00 = vz[ix*n1+iz2];
				C10 = vz[(ix+1)*n1+iz2];
				C01 = vz[ix*n1+iz2+1];
				C11 = vz[(ix+1)*n1+iz2+1];
				rec_vz[irec*rec.nt+isam] = C00*(1.0-rdx)*(1.0-rdz) + C10*rdx*(1.0-rdz) +
										   C01*(1.0-rdx)*rdz       + C11*rdx*rdz;
			}
			if (rec.type.vx) {
				C00 = vx[ix2*n1+iz];
				C10 = vx[(ix2+1)*n1+iz];
				C01 = vx[ix2*n1+iz+1];
				C11 = vx[(ix2+1)*n1+iz+1];
				rec_vx[irec*rec.nt+isam] = C00*(1.0-rdx)*(1.0-rdz) + C10*rdx*(1.0-rdz) +
										   C01*(1.0-rdx)*rdz       + C11*rdx*rdz;
			}
			
		}
		else { /* read values directly from the grid points */
			if (rec.type.p)   rec_p[irec*rec.nt+isam] = tzz[ix*n1+iz];
			if (rec.type.txx) rec_txx[irec*rec.nt+isam] = txx[ix*n1+iz];
			if (rec.type.tzz) rec_tzz[irec*rec.nt+isam] = tzz[ix*n1+iz];
			if (rec.type.txz) {
				if (rec.int_vz == 2 || rec.int_vx == 2) {
					rec_txz[irec*rec.nt+isam] = 0.25*(
							txz[ix*n1+iz2]+txz[ix2*n1+iz2]+
							txz[ix*n1+iz]+txz[ix2*n1+iz]);
				}
				else {
					rec_txz[irec*rec.nt+isam] = txz[ix2*n1+iz2];
				}
			}
			if (rec.type.pp) {
				rec_pp[irec*rec.nt+isam] = (vx[ix2*n1+iz]-vx[ix*n1+iz] +
											vz[ix*n1+iz2]-vz[ix*n1+iz])/mod.dx;
			}
			if (rec.type.ss) {
				rec_ss[irec*rec.nt+isam] = (vx[ix2*n1+iz2]-vx[ix2*n1+iz] -
										   (vz[ix2*n1+iz2]-vz[ix*n1+iz2]))/mod.dx;
			}
			if (rec.type.vz) {
/* interpolate vz to vx position to the right and above of vz */
				if (rec.int_vz == 1) {
					rec_vz[irec*rec.nt+isam] = 0.25*(
							vz[ix*n1+iz2]+vz[ix2*n1+iz2]+
							vz[ix*n1+iz]+vz[ix2*n1+iz]);
				}
/* interpolate vz to Txx/Tzz position by taking the mean of 2 values */
				else if (rec.int_vz == 2) {
					rec_vz[irec*rec.nt+isam] = 0.5*(vz[ix*n1+iz2]+vz[ix*n1+iz]);
				}
				else {
//					fprintf(stderr,"getting Vz at x=%d z=%d value=%e\n", ix, iz, vz[ix*n1+iz2]);
					rec_vz[irec*rec.nt+isam] = vz[ix*n1+iz2];
				}
			}
			if (rec.type.vx) {
/* interpolate vx to vz position to the left and below of vx */
				if (rec.int_vx == 1) {
					rec_vx[irec*rec.nt+isam] = 0.25*(
							vx[ix2*n1+iz]+vx[ix2*n1+iz2]+
							vx[ix*n1+iz]+vx[ix*n1+iz2]);
				}
/* interpolate vx to Txx/Tzz position by taking the mean of 2 values */
				else if (rec.int_vx == 2) {
					rec_vx[irec*rec.nt+isam] = 0.5*(vx[ix2*n1+iz]+vx[ix*n1+iz]);
				}
				else {
					rec_vx[irec*rec.nt+isam] = vx[ix2*n1+iz];
				}
			}
		}

	} /* end of irec loop */

	/* store all x-values on z-level for P Vz for up-down decomposition */
	if (rec.type.ud) {
		iz = rec.z[0]+ibndz;
		iz2 = iz+1;
		vz_t = (float *)calloc(2*mod.nax,sizeof(float));
		/* P and Vz are staggered in time and need to correct for this */
		/* -1- compute Vz at next time step and average with current time step */
		c1 = 9.0/8.0;
		c2 = -1.0/24.0;
		roz = mod.dt/(mod.dx*rec.rho);
    	for (ix=mod.ioZx; ix<mod.ieZx; ix++) {
           	vz_t[ix] = vz[ix*n1+iz] - roz*(
                       	c1*(tzz[ix*n1+iz]   - tzz[ix*n1+iz-1]) +
                       	c2*(tzz[ix*n1+iz+1] - tzz[ix*n1+iz-2]));
           	vz_t[mod.nax+ix] = vz[ix*n1+iz2] - roz*(
                       	c1*(tzz[ix*n1+iz2]   - tzz[ix*n1+iz2-1]) +
                       	c2*(tzz[ix*n1+iz2+1] - tzz[ix*n1+iz2-2]));
       	}
		for (ix=0; ix<mod.nax; ix++) {
			/* -2- compute average in time and depth to get Vz at same depth and time as P */
			rec_udvz[ix*rec.nt+isam] = 0.25*(vz[ix*n1+iz2]+vz[ix*n1+iz]+vz_t[mod.nax+ix]+vz_t[ix]);
			rec_udp[ix*rec.nt+isam]  = tzz[ix*n1+iz];
		}
		free(vz_t);
	}

	return 0;
}
