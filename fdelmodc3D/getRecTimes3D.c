#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc3D.h"
#include"par.h"

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

long getRecTimes3D(modPar mod, recPar rec, bndPar bnd, long itime, long isam, float *vx, float *vy, float *vz, float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz, float *l2m, float *rox, float *roy, float *roz, float *rec_vx, float *rec_vy, float *rec_vz, float *rec_txx, float *rec_tyy, float *rec_tzz, float *rec_txz, float *rec_txy, float *rec_tyz, float *rec_p, float *rec_pp, float *rec_ss, float *rec_udp, float *rec_udvz, long verbose)
{
	long n1, n2, ibndx, ibndy, ibndz;
	long irec, ix, iy, iz, ix2, iy2, iz2, ix1, iy1, iz1;
	float dvx, dvy, dvz, rdz, rdy, rdx;
    float C000, C100, C010, C001, C110, C101, C011, C111;
	float *vz_t, c1, c2, lroz, field;

    ibndx = mod.ioPx;
    ibndy = mod.ioPy;
    ibndz = mod.ioPz;
    if (bnd.lef==4 || bnd.lef==2) ibndx += bnd.ntap;
    if (bnd.bac==4 || bnd.fro==2) ibndy += bnd.ntap;
    if (bnd.top==4 || bnd.top==2) ibndz += bnd.ntap;
	n1    = mod.naz;
    n2    = mod.nax;
	c1 = 9.0/8.0;
	c2 = -1.0/24.0;

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
		iy = rec.y[irec]+ibndy;
		ix = rec.x[irec]+ibndx;
		iz1 = iz-1;
		iy1 = iy-1;
		ix1 = ix-1;
		iz2 = iz+1;
		iy2 = iy+1;
		ix2 = ix+1;
		/* interpolation to precise (not necessary on a grid point) position */
		if ( rec.int_p==3 ) {

			iz = (long)floorf(rec.zr[irec]/mod.dz)+ibndz;
			iy = (long)floorf(rec.yr[irec]/mod.dy)+ibndy;
			ix = (long)floorf(rec.xr[irec]/mod.dx)+ibndx;
			rdz = (rec.zr[irec] - (iz-ibndz)*mod.dz)/mod.dz;
			rdy = (rec.yr[irec] - (iy-ibndy)*mod.dy)/mod.dy;
			rdx = (rec.xr[irec] - (ix-ibndx)*mod.dx)/mod.dx;
			iz1 = iz-1;
			iy1 = iy-1;
			ix1 = ix-1;
			iz2 = iz+1;
			iy2 = iy+1;
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
				C000 = tzz[iy*n1*n2+ix*n1+iz];
				C100 = tzz[iy*n1*n2+(ix+1)*n1+iz];
				C010 = tzz[iy*n1*n2+ix*n1+iz+1];
				C001 = tzz[(iy+1)*n1*n2+ix*n1+iz];
				C110 = tzz[iy*n1*n2+(ix+1)*n1+iz+1];
				C101 = tzz[(iy+1)*n1*n2+(ix+1)*n1+iz];
				C011 = tzz[(iy+1)*n1*n2+ix*n1+iz+1];
				C111 = tzz[(iy+1)*n1*n2+(ix+1)*n1+iz+1];
				rec_p[irec*rec.nt+isam] =   C000*(1.0-rdx)*(1.0-rdz)*(1.0-rdy) +
                                            C100*rdx*(1.0-rdz)*(1.0-rdy) +
										    C010*(1.0-rdx)*rdz*(1.0-rdy) +
										    C001*(1.0-rdx)*(1.0-rdz)*rdy +
                                            C110*rdx*rdz*(1.0-rdy) +
											C101*rdx*(1.0-rdz)*rdy +
											C011*(1.0-rdx)*rdz*rdy +
											C111*rdx*rdz*rdy;
			}
			if (rec.type.txx) {
				C000 = txx[iy*n1*n2+ix*n1+iz];
				C100 = txx[iy*n1*n2+(ix+1)*n1+iz];
				C010 = txx[iy*n1*n2+ix*n1+iz+1];
				C001 = txx[(iy+1)*n1*n2+ix*n1+iz];
				C110 = txx[iy*n1*n2+(ix+1)*n1+iz+1];
				C101 = txx[(iy+1)*n1*n2+(ix+1)*n1+iz];
				C011 = txx[(iy+1)*n1*n2+ix*n1+iz+1];
				C111 = txx[(iy+1)*n1*n2+(ix+1)*n1+iz+1];
				rec_txx[irec*rec.nt+isam] = C000*(1.0-rdx)*(1.0-rdz)*(1.0-rdy) +
                                            C100*rdx*(1.0-rdz)*(1.0-rdy) +
										    C010*(1.0-rdx)*rdz*(1.0-rdy) +
										    C001*(1.0-rdx)*(1.0-rdz)*rdy +
                                            C110*rdx*rdz*(1.0-rdy) +
											C101*rdx*(1.0-rdz)*rdy +
											C011*(1.0-rdx)*rdz*rdy +
											C111*rdx*rdz*rdy;
			}
			if (rec.type.tyy) {
				C000 = tyy[iy*n1*n2+ix*n1+iz];
				C100 = tyy[iy*n1*n2+(ix+1)*n1+iz];
				C010 = tyy[iy*n1*n2+ix*n1+iz+1];
				C001 = tyy[(iy+1)*n1*n2+ix*n1+iz];
				C110 = tyy[iy*n1*n2+(ix+1)*n1+iz+1];
				C101 = tyy[(iy+1)*n1*n2+(ix+1)*n1+iz];
				C011 = tyy[(iy+1)*n1*n2+ix*n1+iz+1];
				C111 = tyy[(iy+1)*n1*n2+(ix+1)*n1+iz+1];
				rec_tyy[irec*rec.nt+isam] = C000*(1.0-rdx)*(1.0-rdz)*(1.0-rdy) +
                                            C100*rdx*(1.0-rdz)*(1.0-rdy) +
										    C010*(1.0-rdx)*rdz*(1.0-rdy) +
										    C001*(1.0-rdx)*(1.0-rdz)*rdy +
                                            C110*rdx*rdz*(1.0-rdy) +
											C101*rdx*(1.0-rdz)*rdy +
											C011*(1.0-rdx)*rdz*rdy +
											C111*rdx*rdz*rdy;
			}
			if (rec.type.tzz) {
				C000 = tzz[iy*n1*n2+ix*n1+iz];
				C100 = tzz[iy*n1*n2+(ix+1)*n1+iz];
				C010 = tzz[iy*n1*n2+ix*n1+iz+1];
				C001 = tzz[(iy+1)*n1*n2+ix*n1+iz];
				C110 = tzz[iy*n1*n2+(ix+1)*n1+iz+1];
				C101 = tzz[(iy+1)*n1*n2+(ix+1)*n1+iz];
				C011 = tzz[(iy+1)*n1*n2+ix*n1+iz+1];
				C111 = tzz[(iy+1)*n1*n2+(ix+1)*n1+iz+1];
				rec_tzz[irec*rec.nt+isam] = C000*(1.0-rdx)*(1.0-rdz)*(1.0-rdy) +
                                            C100*rdx*(1.0-rdz)*(1.0-rdy) +
										    C010*(1.0-rdx)*rdz*(1.0-rdy) +
										    C001*(1.0-rdx)*(1.0-rdz)*rdy +
                                            C110*rdx*rdz*(1.0-rdy) +
											C101*rdx*(1.0-rdz)*rdy +
											C011*(1.0-rdx)*rdz*rdy +
											C111*rdx*rdz*rdy;
			}
			if (rec.type.txz) {
				C000 = txz[iy*n1*n2+ix*n1+iz];
				C100 = txz[iy*n1*n2+(ix+1)*n1+iz];
				C010 = txz[iy*n1*n2+ix*n1+iz+1];
				C001 = txz[(iy+1)*n1*n2+ix*n1+iz];
				C110 = txz[iy*n1*n2+(ix+1)*n1+iz+1];
				C101 = txz[(iy+1)*n1*n2+(ix+1)*n1+iz];
				C011 = txz[(iy+1)*n1*n2+ix*n1+iz+1];
				C111 = txz[(iy+1)*n1*n2+(ix+1)*n1+iz+1];
				rec_txz[irec*rec.nt+isam] = C000*(1.0-rdx)*(1.0-rdz)*(1.0-rdy) +
                                            C100*rdx*(1.0-rdz)*(1.0-rdy) +
										    C010*(1.0-rdx)*rdz*(1.0-rdy) +
										    C001*(1.0-rdx)*(1.0-rdz)*rdy +
                                            C110*rdx*rdz*(1.0-rdy) +
											C101*rdx*(1.0-rdz)*rdy +
											C011*(1.0-rdx)*rdz*rdy +
											C111*rdx*rdz*rdy;
			}
			if (rec.type.txy) {
				C000 = txy[iy*n1*n2+ix*n1+iz];
				C100 = txy[iy*n1*n2+(ix+1)*n1+iz];
				C010 = txy[iy*n1*n2+ix*n1+iz+1];
				C001 = txy[(iy+1)*n1*n2+ix*n1+iz];
				C110 = txy[iy*n1*n2+(ix+1)*n1+iz+1];
				C101 = txy[(iy+1)*n1*n2+(ix+1)*n1+iz];
				C011 = txy[(iy+1)*n1*n2+ix*n1+iz+1];
				C111 = txy[(iy+1)*n1*n2+(ix+1)*n1+iz+1];
				rec_txy[irec*rec.nt+isam] = C000*(1.0-rdx)*(1.0-rdz)*(1.0-rdy) +
                                            C100*rdx*(1.0-rdz)*(1.0-rdy) +
										    C010*(1.0-rdx)*rdz*(1.0-rdy) +
										    C001*(1.0-rdx)*(1.0-rdz)*rdy +
                                            C110*rdx*rdz*(1.0-rdy) +
											C101*rdx*(1.0-rdz)*rdy +
											C011*(1.0-rdx)*rdz*rdy +
											C111*rdx*rdz*rdy;
			}
			if (rec.type.tyz) {
				C000 = tyz[iy*n1*n2+ix*n1+iz];
				C100 = tyz[iy*n1*n2+(ix+1)*n1+iz];
				C010 = tyz[iy*n1*n2+ix*n1+iz+1];
				C001 = tyz[(iy+1)*n1*n2+ix*n1+iz];
				C110 = tyz[iy*n1*n2+(ix+1)*n1+iz+1];
				C101 = tyz[(iy+1)*n1*n2+(ix+1)*n1+iz];
				C011 = tyz[(iy+1)*n1*n2+ix*n1+iz+1];
				C111 = tyz[(iy+1)*n1*n2+(ix+1)*n1+iz+1];
				rec_tyz[irec*rec.nt+isam] = C000*(1.0-rdx)*(1.0-rdz)*(1.0-rdy) +
                                            C100*rdx*(1.0-rdz)*(1.0-rdy) +
										    C010*(1.0-rdx)*rdz*(1.0-rdy) +
										    C001*(1.0-rdx)*(1.0-rdz)*rdy +
                                            C110*rdx*rdz*(1.0-rdy) +
											C101*rdx*(1.0-rdz)*rdy +
											C011*(1.0-rdx)*rdz*rdy +
											C111*rdx*rdz*rdy;
			}
			if (rec.type.pp) {
				C000 = (vx[iy*n1*n2+ix2*n1+iz]-vx[iy*n1*n2+ix*n1+iz] +
						vy[iy2*n1*n2+ix*n1+iz]-vy[iy*n1*n2+ix*n1+iz] +
						vz[iy*n1*n2+ix*n1+iz2]-vz[iy*n1*n2+ix*n1+iz])/mod.dx;
				C100 = (vx[iy*n1*n2+(ix2+1)*n1+iz]-vx[iy*n1*n2+(ix+1)*n1+iz] +
						vy[iy2*n1*n2+(ix+1)*n1+iz]-vy[iy*n1*n2+(ix+1)*n1+iz] +
						vz[iy*n1*n2+(ix+1)*n1+iz2]-vz[iy*n1*n2+(ix+1)*n1+iz])/mod.dx;
				C010 = (vx[iy*n1*n2+ix2*n1+iz+1]-vx[iy*n1*n2+ix*n1+iz+1] +
						vy[iy2*n1*n2+ix*n1+iz+1]-vy[iy*n1*n2+ix*n1+iz+1] +
						vz[iy*n1*n2+ix*n1+iz2+1]-vz[iy*n1*n2+ix*n1+iz+1])/mod.dx;
				C001 = (vx[(iy+1)*n1*n2+ix2*n1+iz]-vx[(iy+1)*n1*n2+ix*n1+iz] +
						vy[(iy2+1)*n1*n2+ix*n1+iz]-vy[(iy+1)*n1*n2+ix*n1+iz] +
						vz[(iy+1)*n1*n2+ix*n1+iz2]-vz[(iy+1)*n1*n2+ix*n1+iz])/mod.dx;
				C110 = (vx[iy*n1*n2+(ix2+1)*n1+iz+1]-vx[iy*n1*n2+(ix+1)*n1+iz+1] +
						vy[iy2*n1*n2+(ix+1)*n1+iz+1]-vy[iy*n1*n2+(ix+1)*n1+iz+1] +
						vz[iy*n1*n2+(ix+1)*n1+iz2+1]-vz[iy*n1*n2+(ix+1)*n1+iz+1])/mod.dx;
				C101 = (vx[(iy+1)*n1*n2+(ix2+1)*n1+iz]-vx[(iy+1)*n1*n2+(ix+1)*n1+iz] +
						vy[(iy2+1)*n1*n2+(ix+1)*n1+iz]-vy[(iy+1)*n1*n2+(ix+1)*n1+iz] +
						vz[(iy+1)*n1*n2+(ix+1)*n1+iz2]-vz[(iy+1)*n1*n2+(ix+1)*n1+iz])/mod.dx;
				C011 = (vx[(iy+1)*n1*n2+ix2*n1+iz+1]-vx[(iy+1)*n1*n2+ix*n1+iz+1] +
						vy[(iy2+1)*n1*n2+ix*n1+iz+1]-vy[(iy+1)*n1*n2+ix*n1+iz+1] +
						vz[(iy+1)*n1*n2+ix*n1+iz2+1]-vz[(iy+1)*n1*n2+ix*n1+iz+1])/mod.dx;
				C111 = (vx[(iy+1)*n1*n2+(ix2+1)*n1+iz+1]-vx[(iy+1)*n1*n2+(ix+1)*n1+iz+1] +
						vy[(iy2+1)*n1*n2+(ix+1)*n1+iz+1]-vy[(iy+1)*n1*n2+(ix+1)*n1+iz+1] +
						vz[(iy+1)*n1*n2+(ix+1)*n1+iz2+1]-vz[(iy+1)*n1*n2+(ix+1)*n1+iz+1])/mod.dx;
				rec_pp[irec*rec.nt+isam] =  C000*(1.0-rdx)*(1.0-rdz)*(1.0-rdy) +
                                            C100*rdx*(1.0-rdz)*(1.0-rdy) +
										    C010*(1.0-rdx)*rdz*(1.0-rdy) +
										    C001*(1.0-rdx)*(1.0-rdz)*rdy +
                                            C110*rdx*rdz*(1.0-rdy) +
											C101*rdx*(1.0-rdz)*rdy +
											C011*(1.0-rdx)*rdz*rdy +
											C111*rdx*rdz*rdy;
			}
			if (rec.type.ss) {
				C000 = 	(vx[iy2*n1*n2+ix2*n1+iz2]-vx[iy*n1*n2+ix2*n1+iz] -
						(vy[iy2*n1*n2+ix2*n1+iz2]-vy[iy2*n1*n2+ix*n1+iz]) -
						(vz[iy2*n1*n2+ix2*n1+iz2]-vz[iy*n1*n2+ix*n1+iz2]))/mod.dx;
				C100 =	(vx[iy2*n1*n2+(ix2+1)*n1+iz2]-vx[iy*n1*n2+(ix2+1)*n1+iz] -
						(vy[iy2*n1*n2+(ix2+1)*n1+iz2]-vy[iy2*n1*n2+(ix+1)*n1+iz]) -
						(vz[iy2*n1*n2+(ix2+1)*n1+iz2]-vz[iy*n1*n2+(ix+1)*n1+iz2]))/mod.dx;
				C010 =	(vx[iy2*n1*n2+ix2*n1+iz2+1]-vx[iy*n1*n2+ix2*n1+iz+1] -
						(vy[iy2*n1*n2+ix2*n1+iz2+1]-vy[iy2*n1*n2+ix*n1+iz+1]) -
						(vz[iy2*n1*n2+ix2*n1+iz2+1]-vz[iy*n1*n2+ix*n1+iz2+1]))/mod.dx;
				C001 =	(vx[(iy2+1)*n1*n2+ix2*n1+iz2]-vx[(iy+1)*n1*n2+ix2*n1+iz] -
						(vy[(iy2+1)*n1*n2+ix2*n1+iz2]-vy[(iy2+1)*n1*n2+ix*n1+iz]) -
						(vz[(iy2+1)*n1*n2+ix2*n1+iz2]-vz[(iy+1)*n1*n2+ix*n1+iz2]))/mod.dx;
				C110 =	(vx[iy2*n1*n2+(ix2+1)*n1+iz2+1]-vx[iy*n1*n2+(ix2+1)*n1+iz+1] -
						(vy[iy2*n1*n2+(ix2+1)*n1+iz2+1]-vy[iy2*n1*n2+(ix+1)*n1+iz+1]) -
						(vz[iy2*n1*n2+(ix2+1)*n1+iz2+1]-vz[iy*n1*n2+(ix+1)*n1+iz2+1]))/mod.dx;
				C101 =	(vx[(iy2+1)*n1*n2+(ix2+1)*n1+iz2]-vx[(iy+1)*n1*n2+(ix2+1)*n1+iz] -
						(vy[(iy2+1)*n1*n2+(ix2+1)*n1+iz2]-vy[(iy2+1)*n1*n2+(ix+1)*n1+iz]) -
						(vz[(iy2+1)*n1*n2+(ix2+1)*n1+iz2]-vz[(iy+1)*n1*n2+(ix+1)*n1+iz2]))/mod.dx;
				C011 =	(vx[(iy2+1)*n1*n2+ix2*n1+iz2+1]-vx[(iy+1)*n1*n2+ix2*n1+iz+1] -
						(vy[(iy2+1)*n1*n2+ix2*n1+iz2+1]-vy[(iy2+1)*n1*n2+ix*n1+iz+1]) -
						(vz[(iy2+1)*n1*n2+ix2*n1+iz2+1]-vz[(iy+1)*n1*n2+ix*n1+iz2+1]))/mod.dx;
				C111 =	(vx[(iy2+1)*n1*n2+(ix2+1)*n1+iz2+1]-vx[(iy+1)*n1*n2+(ix2+1)*n1+iz+1] -
						(vy[(iy2+1)*n1*n2+(ix2+1)*n1+iz2+1]-vy[(iy2+1)*n1*n2+(ix+1)*n1+iz+1]) -
						(vz[(iy2+1)*n1*n2+(ix2+1)*n1+iz2+1]-vz[(iy+1)*n1*n2+(ix+1)*n1+iz2+1]))/mod.dx;
				rec_ss[irec*rec.nt+isam] =  C000*(1.0-rdx)*(1.0-rdz)*(1.0-rdy) +
                                            C100*rdx*(1.0-rdz)*(1.0-rdy) +
										    C010*(1.0-rdx)*rdz*(1.0-rdy) +
										    C001*(1.0-rdx)*(1.0-rdz)*rdy +
                                            C110*rdx*rdz*(1.0-rdy) +
											C101*rdx*(1.0-rdz)*rdy +
											C011*(1.0-rdx)*rdz*rdy +
											C111*rdx*rdz*rdy;
			}
			if (rec.type.vz) {
				C000 = vz[iy*n1*n2+ix*n1+iz2];
				C100 = vz[iy*n1*n2+(ix+1)*n1+iz2];
				C010 = vz[iy*n1*n2+ix*n1+iz2+1];
				C001 = vz[(iy+1)*n1*n2+ix*n1+iz2];
				C110 = vz[iy*n1*n2+(ix+1)*n1+iz2+1];
				C101 = vz[(iy+1)*n1*n2+(ix+1)*n1+iz2];
				C011 = vz[(iy+1)*n1*n2+ix*n1+iz2+1];
				C111 = vz[(iy+1)*n1*n2+(ix+1)*n1+iz2+1];
				rec_vz[irec*rec.nt+isam] =  C000*(1.0-rdx)*(1.0-rdz)*(1.0-rdy) +
                                            C100*rdx*(1.0-rdz)*(1.0-rdy) +
										    C010*(1.0-rdx)*rdz*(1.0-rdy) +
										    C001*(1.0-rdx)*(1.0-rdz)*rdy +
                                            C110*rdx*rdz*(1.0-rdy) +
											C101*rdx*(1.0-rdz)*rdy +
											C011*(1.0-rdx)*rdz*rdy +
											C111*rdx*rdz*rdy;
			}
			if (rec.type.vy) {
				C000 = vy[iy2*n1*n2+ix*n1+iz];
				C100 = vy[iy2*n1*n2+(ix+1)*n1+iz];
				C010 = vy[iy2*n1*n2+ix*n1+iz+1];
				C001 = vy[(iy2+1)*n1*n2+ix*n1+iz];
				C110 = vy[iy2*n1*n2+(ix+1)*n1+iz+1];
				C101 = vy[(iy2+1)*n1*n2+(ix+1)*n1+iz];
				C011 = vy[(iy2+1)*n1*n2+ix*n1+iz+1];
				C111 = vy[(iy2+1)*n1*n2+(ix+1)*n1+iz+1];
				rec_vy[irec*rec.nt+isam] =  C000*(1.0-rdx)*(1.0-rdz)*(1.0-rdy) +
                                            C100*rdx*(1.0-rdz)*(1.0-rdy) +
										    C010*(1.0-rdx)*rdz*(1.0-rdy) +
										    C001*(1.0-rdx)*(1.0-rdz)*rdy +
                                            C110*rdx*rdz*(1.0-rdy) +
											C101*rdx*(1.0-rdz)*rdy +
											C011*(1.0-rdx)*rdz*rdy +
											C111*rdx*rdz*rdy;
			}
			if (rec.type.vx) {
				C000 = vy[iy*n1*n2+ix2*n1+iz];
				C100 = vy[iy*n1*n2+(ix2+1)*n1+iz];
				C010 = vy[iy*n1*n2+ix2*n1+iz+1];
				C001 = vy[(iy+1)*n1*n2+ix2*n1+iz];
				C110 = vy[iy*n1*n2+(ix2+1)*n1+iz+1];
				C101 = vy[(iy+1)*n1*n2+(ix2+1)*n1+iz];
				C011 = vy[(iy+1)*n1*n2+ix2*n1+iz+1];
				C111 = vy[(iy+1)*n1*n2+(ix2+1)*n1+iz+1];
				rec_vx[irec*rec.nt+isam] =  C000*(1.0-rdx)*(1.0-rdz)*(1.0-rdy) +
                                            C100*rdx*(1.0-rdz)*(1.0-rdy) +
										    C010*(1.0-rdx)*rdz*(1.0-rdy) +
										    C001*(1.0-rdx)*(1.0-rdz)*rdy +
                                            C110*rdx*rdz*(1.0-rdy) +
											C101*rdx*(1.0-rdz)*rdy +
											C011*(1.0-rdx)*rdz*rdy +
											C111*rdx*rdz*rdy;
			}
			
		}
		else { /* read values directly from the grid points */
			if (verbose>=4 && isam==0) {
				vmess("Receiver %li read at gridpoint ix=%li iy=%li iz=%li",irec, ix, iy, iz);
			}
			/* interpolation of receivers to same time step is only done for acoustic scheme */
			if (rec.type.p) {
				if (rec.int_p == 1) {
					if (mod.ischeme == 1) { /* interpolate Tzz times -1/2 Dt backward to Vz times */
                        dvx = c1*(vx[iy*n1*n2+(ix+1)*n1+iz] - vx[iy*n1*n2+ix*n1+iz]) +
                              c2*(vx[iy*n1*n2+(ix+2)*n1+iz] - vx[iy*n1*n2+(ix-1)*n1+iz]);
                        dvy = c1*(vy[(iy+1)*n1*n2+ix*n1+iz] - vy[iy*n1*n2+ix*n1+iz]) +
                              c2*(vy[(iy+2)*n1*n2+ix*n1+iz] - vy[(iy-1)*n1*n2+ix*n1+iz]);
                        dvz = c1*(vz[iy*n1*n2+ix*n1+iz+1]   - vz[iy*n1*n2+ix*n1+iz]) +
                              c2*(vz[iy*n1*n2+ix*n1+iz+2]   - vz[iy*n1*n2+ix*n1+iz-1]);
                        field = tzz[iy*n1*n2+ix*n1+iz] + (1.0/2.0)*l2m[iy*n1*n2+ix*n1+iz]*(dvx+dvy+dvz);
                        dvx = c1*(vx[iy*n1*n2+(ix+1)*n1+iz1] - vx[iy*n1*n2+ix*n1+iz1]) +
                              c2*(vx[iy*n1*n2+(ix+2)*n1+iz1] - vx[iy*n1*n2+(ix-1)*n1+iz1]);
                        dvy = c1*(vy[(iy+1)*n1*n2+ix*n1+iz1] - vy[iy*n1*n2+ix*n1+iz1]) +
                              c2*(vy[(iy+2)*n1*n2+ix*n1+iz1] - vy[(iy-1)*n1*n2+ix*n1+iz1]);
                        dvz = c1*(vz[iy*n1*n2+ix*n1+iz1+1]   - vz[iy*n1*n2+ix*n1+iz1]) +
                              c2*(vz[iy*n1*n2+ix*n1+iz1+2]   - vz[iy*n1*n2+ix*n1+iz1-1]);
                        field += tzz[iy*n1*n2+ix*n1+iz1] + (1.0/2.0)*l2m[iy*n1*n2+ix*n1+iz1]*(dvx+dvy+dvz);
						rec_p[irec*rec.nt+isam] = 0.5*field;
					}
					else {
						rec_p[irec*rec.nt+isam] = 0.5*(tzz[iy*n1*n2+ix*n1+iz1]+tzz[iy*n1*n2+ix*n1+iz]);
					}
				}
				else if (rec.int_p == 2) {
					if (mod.ischeme == 1) { /* interpolate Tzz times -1/2 Dt backward to Vx times */
                        dvx = c1*(vx[iy*n1*n2+(ix+1)*n1+iz] - vx[iy*n1*n2+ix*n1+iz]) +
                              c2*(vx[iy*n1*n2+(ix+2)*n1+iz] - vx[iy*n1*n2+(ix-1)*n1+iz]);
                        dvy = c1*(vy[(iy+1)*n1*n2+ix*n1+iz] - vy[iy*n1*n2+ix*n1+iz]) +
                              c2*(vy[(iy+2)*n1*n2+ix*n1+iz] - vy[(iy-1)*n1*n2+ix*n1+iz]);
                        dvz = c1*(vz[iy*n1*n2+ix*n1+iz+1]   - vz[iy*n1*n2+ix*n1+iz]) +
                              c2*(vz[iy*n1*n2+ix*n1+iz+2]   - vz[iy*n1*n2+ix*n1+iz-1]);
                        field = tzz[iy*n1*n2+ix*n1+iz] + (1.0/2.0)*l2m[iy*n1*n2+ix*n1+iz]*(dvx+dvy+dvz);
                        dvx = c1*(vx[iy*n1*n2+(ix1+1)*n1+iz] - vx[iy*n1*n2+ix1*n1+iz]) +
                              c2*(vx[iy*n1*n2+(ix1+2)*n1+iz] - vx[iy*n1*n2+(ix1-1)*n1+iz]);
                        dvy = c1*(vy[(iy+1)*n1*n2+ix1*n1+iz] - vy[iy*n1*n2+ix1*n1+iz]) +
                              c2*(vy[(iy+2)*n1*n2+ix1*n1+iz] - vy[(iy-1)*n1*n2+ix1*n1+iz]);
                        dvz = c1*(vz[iy*n1*n2+ix1*n1+iz+1]   - vz[iy*n1*n2+ix1*n1+iz]) +
                              c2*(vz[iy*n1*n2+ix1*n1+iz+2]   - vz[iy*n1*n2+ix1*n1+iz-1]);
                        field += tzz[iy*n1*n2+ix1*n1+iz] + (1.0/2.0)*l2m[iy*n1*n2+ix1*n1+iz]*(dvx+dvy+dvz);
						rec_p[irec*rec.nt+isam] = 0.5*field;
					}
					else {
						rec_p[irec*rec.nt+isam] = 0.5*(tzz[iy*n1*n2+ix1*n1+iz]+tzz[iy*n1*n2+ix*n1+iz]);
					}
				}
				else {
					rec_p[irec*rec.nt+isam] = tzz[iy*n1*n2+ix*n1+iz];
				}
			}
			if (rec.type.txx) rec_txx[irec*rec.nt+isam] = txx[iy*n1*n2+ix*n1+iz];
			if (rec.type.tzz) rec_tzz[irec*rec.nt+isam] = tzz[iy*n1*n2+ix*n1+iz];
			if (rec.type.txz) { /* time interpolation to be done */
				if (rec.int_vz == 2 || rec.int_vx == 2) {
					rec_txz[irec*rec.nt+isam] = 0.125*(
							txz[iy*n1*n2+ix*n1+iz]  + txz[iy*n1*n2+ix*n1+iz2]+
							txz[iy2*n1*n2+ix*n1+iz] + txz[iy2*n1*n2+ix*n1+iz2]+
							txz[iy*n1*n2+ix2*n1+iz]  + txz[iy*n1*n2+ix2*n1+iz2]+
							txz[iy2*n1*n2+ix2*n1+iz] + txz[iy2*n1*n2+ix2*n1+iz2]);
				}
				else {
					rec_txz[irec*rec.nt+isam] = txz[iy2*n1*n2+ix2*n1+iz2];
				}
			}
			if (rec.type.pp) {
				rec_pp[irec*rec.nt+isam] = (vx[iy*n1*n2+ix2*n1+iz]-vx[iy*n1*n2+ix*n1+iz] +
											vy[iy2*n1*n2+ix*n1+iz]-vy[iy*n1*n2+ix*n1+iz] +
											vz[iy*n1*n2+ix*n1+iz2]-vz[iy*n1*n2+ix*n1+iz])/mod.dx;
			}
			if (rec.type.ss) {
				rec_ss[irec*rec.nt+isam] = (vx[iy2*n1*n2+ix2*n1+iz2]-vx[iy*n1*n2+ix2*n1+iz] -
										   (vy[iy2*n1*n2+ix2*n1+iz2]-vy[iy2*n1*n2+ix*n1+iz]) -
										   (vz[iy2*n1*n2+ix2*n1+iz2]-vz[iy*n1*n2+ix*n1+iz2]))/mod.dx;
			}
			if (rec.type.vz) {
/* interpolate vz to vx position to the right and above of vz */
				if (rec.int_vz == 1) {
					rec_vz[irec*rec.nt+isam] = 0.125*(
							vz[iy*n1*n2+ix*n1+iz2]+vz[iy*n1*n2+ix1*n1+iz2]+
							vz[iy*n1*n2+ix*n1+iz] +vz[iy*n1*n2+ix1*n1+iz]+
							vz[iy1*n1*n2+ix*n1+iz2]+vz[iy1*n1*n2+ix1*n1+iz2]+
							vz[iy1*n1*n2+ix*n1+iz] +vz[iy1*n1*n2+ix1*n1+iz]);
				}
/* interpolate vz to Txx/Tzz position by taking the mean of 2 values */
				else if (rec.int_vz == 2) {
					if (mod.ischeme == 1) { /* interpolate Vz times +1/2 Dt forward to P times */
                        field = vz[iy*n1*n2+ix*n1+iz] - 0.5*roz[iy*n1*n2+ix*n1+iz]*(
                        	c1*(tzz[iy*n1*n2+ix*n1+iz]   - tzz[iy*n1*n2+ix*n1+iz-1]) +
                        	c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
                        field += vz[iy*n1*n2+ix*n1+iz2] - 0.5*roz[iy*n1*n2+ix*n1+iz2]*(
                        	c1*(tzz[iy*n1*n2+ix*n1+iz2]   - tzz[iy*n1*n2+ix*n1+iz2-1]) +
                        	c2*(tzz[iy*n1*n2+ix*n1+iz2+1] - tzz[iy*n1*n2+ix*n1+iz2-2]));
						rec_vz[irec*rec.nt+isam] = 0.5*field;
					}
					else {
						rec_vz[irec*rec.nt+isam] = 0.5*(vz[iy*n1*n2+ix*n1+iz2]+vz[iy*n1*n2+ix*n1+iz]);
					}
				}
				else {
					rec_vz[irec*rec.nt+isam] = vz[iy*n1*n2+ix*n1+iz2];
					//rec_vz[irec*rec.nt+isam] = vz[ix*n1+iz];
					//fprintf(stderr,"isam=%li vz[%li]=%e vz[%li]=%e vz[%li]=%e \n",isam, iz-1,vz[ix*n1+iz-1],iz,vz[ix*n1+iz], iz+1, vz[ix*n1+iz+1]);
				}
			}
			if (rec.type.vx) {
/* interpolate vx to vz position to the left and below of vx */
				if (rec.int_vx == 1) {
					rec_vx[irec*rec.nt+isam] = 0.125*(
							vx[iy*n1*n2+ix2*n1+iz]+vx[iy*n1*n2+ix2*n1+iz1]+
							vx[iy*n1*n2+ix*n1+iz]+vx[iy*n1*n2+ix*n1+iz1]+
							vx[iy2*n1*n2+ix2*n1+iz]+vx[iy2*n1*n2+ix2*n1+iz1]+
							vx[iy2*n1*n2+ix*n1+iz]+vx[iy2*n1*n2+ix*n1+iz1]);
				}
/* interpolate vx to Txx/Tzz position by taking the mean of 2 values */
				else if (rec.int_vx == 2) {
					if (mod.ischeme == 1) { /* interpolate Vx times +1/2 Dt forward to P times */
            			field = vx[iy*n1*n2+ix*n1+iz] - 0.5*rox[iy*n1*n2+ix*n1+iz]*(
                			c1*(tzz[iy*n1*n2+ix*n1+iz]     - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
                			c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
            			field += vx[ix2*n1+iz] - 0.5*rox[ix2*n1+iz]*(
                			c1*(tzz[iy*n1*n2+ix2*n1+iz]     - tzz[iy*n1*n2+(ix2-1)*n1+iz]) +
                			c2*(tzz[iy*n1*n2+(ix2+1)*n1+iz] - tzz[iy*n1*n2+(ix2-2)*n1+iz]));
						rec_vx[irec*rec.nt+isam] = 0.5*field;
					}
					else {
						rec_vx[irec*rec.nt+isam] = 0.5*(vx[iy*n1*n2+ix2*n1+iz]+vx[iy*n1*n2+ix*n1+iz]);
					}
				}
				else {
					rec_vx[irec*rec.nt+isam] = vx[iy*n1*n2+ix2*n1+iz];
				}
			}
		}

	} /* end of irec loop */

	/* store all x-values on z-level for P Vz for up-down decomposition */
	if (rec.type.ud) {
		iz = rec.z[0]+ibndz;
		iz2 = iz+1;
		vz_t = (float *)calloc(2*mod.nax*mod.nay,sizeof(float));
		/* P and Vz are staggered in time and need to correct for this */
		/* -1- compute Vz at next time step and average with current time step */
		lroz = mod.dt/(mod.dx*rec.rho);
    	for (iy=mod.ioZy; iy<mod.ieZy; iy++) {
			for (ix=mod.ioZx; ix<mod.ieZx; ix++) {
				vz_t[iy*mod.nax+ix] = vz[iy*n1*n2+ix*n1+iz] - lroz*(
							c1*(tzz[iy*n1*n2+ix*n1+iz]   - tzz[iy*n1*n2+ix*n1+iz-1]) +
							c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
				vz_t[iy*mod.nax+mod.nay*mod.nax+ix] = vz[iy*n1*n2+ix*n1+iz2] - lroz*(
							c1*(tzz[iy*n1*n2+ix*n1+iz2]   - tzz[iy*n1*n2+ix*n1+iz2-1]) +
							c2*(tzz[iy*n1*n2+ix*n1+iz2+1] - tzz[iy*n1*n2+ix*n1+iz2-2]));
			}
		}
		for (iy=0; iy<mod.nay; iy++) {
			for (ix=0; ix<mod.nax; ix++) {
				/* -2- compute average in time and depth to get Vz at same depth and time as P */
				rec_udvz[iy*mod.nax*rec.nt+ix*rec.nt+isam] = 0.25*(vz[iy*n1*n2+ix*n1+iz2]+vz[iy*n1*n2+ix*n1+iz]+vz_t[iy*mod.nax+mod.nay*mod.nax+ix]+vz_t[iy*mod.nax+ix]);
				rec_udp[iy*mod.nax*rec.nt+ix*rec.nt+isam]  = tzz[iy*n1*n2+ix*n1+iz];
			}
		}
		free(vz_t);
	}

	return 0;
}
