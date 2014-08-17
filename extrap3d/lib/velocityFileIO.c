#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include "segy.h"
#include "Area.h"
#ifdef MPI
#include <mpi.h>
#endif

int openVelocityFile(char *file_vel, FILE **vel_file, Area *vel_area, int verbose)
{
	FILE *fp;
	segy hdr1, hdr2, hdr3;
	int nxv, nyv, nzv, iymin, iymax, ntraces, nxy, n1, pe=0, root_pe=0;
	int i,j,k,vx,vy,vz,iz, binary_file, ixmin,ixmax;
	int twoD;
	float xvmin, yvmin, zvmin, dxv, dyv, dzv;
	float scl, *dum, *zslice;
	char *tmp_dir, tmpname[1024];
	size_t  nread, nwrite, bytes, size, trace_sz;

#ifdef MPI
	MPI_Comm_rank( MPI_COMM_WORLD, &pe );
#endif

/* determine what kind of file is read in; binary or SU format */

	binary_file = (strstr(file_vel, ".su")==NULL);

	if(!getparstring("tmp_dir", &tmp_dir)) tmp_dir="/tmp";

	if(!getparint("twoD", &twoD)) twoD = 0;
	if (twoD) {
		fp = fopen( file_vel, "r" ); 
		assert( fp );

		if(!getparint("nyv", &nyv)) nyv = 1; assert( nyv > 0 );

		nread = fread( &hdr1, 1, TRCBYTES, fp );
		assert (nread == TRCBYTES);

		trace_sz = sizeof(float)*hdr1.ns+TRCBYTES;
		fseek ( fp, 0, SEEK_END );
		bytes = ftell(fp); 
		ntraces  = (int) (bytes/trace_sz);

		nzv = hdr1.ns;
		nxv = ntraces;
		if(!getparfloat("dxv", &dxv)) dxv = hdr1.d2; assert( dxv != 0 );
		if(!getparfloat("dyv", &dyv)) dyv = 1.; assert( dyv != 0 );
		if(!getparfloat("dzv", &dzv)) dzv = hdr1.d1; assert( dzv != 0 );
		if(!getparfloat("xvmin", &xvmin)) xvmin = hdr1.f2; 
		if(!getparfloat("yvmin", &yvmin)) yvmin = 0; 
		if(!getparfloat("zvmin", &zvmin)) zvmin = hdr1.f1; 

		fseek(fp, 0, SEEK_SET);
		j=0;
		dum = (float *)malloc(nxv*nyv*nzv*sizeof(float));
		while (j<ntraces) {
			nread = fread( &hdr1, 1, TRCBYTES, fp );
			if (nread == 0) break;
			nread = fread( &dum[j*hdr1.ns], sizeof(float), hdr1.ns, fp );
			assert (nread == hdr1.ns);
			j++;
		}
		fclose(fp);
		vz = 1;
		vx = 2;
		vy = 3;
		for (j=1; j<nyv; j++) {
			for (i=0; i<nxv*nzv; i++) {
				dum[j*nxv*nzv+i] = dum[i];
			}
		}
	}
	else {
/* Open velocity file and determine the size of the file */

	fp = fopen( file_vel, "r" ); 
	assert( fp );

	if(!getparint("vx", &vx)) vx = 0;
	if(!getparint("vy", &vy)) vy = 0;
	if(!getparint("vz", &vz)) vz = 0;
	if (vx>3 || vy>3 || vz>3 || vx<1 || vy<1 || vz<1 ) { 
		fprintf(stderr,"Error in setting the vx,vy and vz parameters,\n");
		fprintf(stderr,"each parameter must have a value of 1,2 or 3\n");
		fprintf(stderr,"They indicate the way the velocity file is stored.\n");
		fprintf(stderr,"For example vz=1 vx=2 vy=3 means that:\n");
		fprintf(stderr," - the sample direction (1) is depth\n");
		fprintf(stderr," - the trace direction (2) is the x-coordinate\n");
		fprintf(stderr," - the gather direction (3) is the y-coordinate\n");
		exit(1);
	}

	if (verbose) {
		fprintf(stderr,"\n    VELOCITY INFORMATION\n");
	}

	if ( binary_file ) {
		if (verbose) {
			fprintf(stderr," reading binary file\n");
		}
		if(!getparint("nxv", &nxv)) nxv = 1; assert( nxv > 0 );
		if(!getparint("nyv", &nyv)) nyv = 1; assert( nyv > 0 );
		if(!getparfloat("dxv", &dxv)) dxv = 1; assert( dxv != 0 );
		if(!getparfloat("dyv", &dyv)) dyv = 1; assert( dyv != 0 );
		if(!getparfloat("dzv", &dzv)) dzv = 1; assert( dzv != 0 );
		if(!getparfloat("xvmin", &xvmin)) xvmin = 0; 
		if(!getparfloat("yvmin", &yvmin)) yvmin = 0; 
		if(!getparfloat("zvmin", &zvmin)) zvmin = 0; 
		if(!getparint("nzv", &nzv)) {
			fseek(fp, 0, SEEK_END);
			bytes    = ftell(fp); 
			assert ( !(bytes % sizeof(float)) );
			size     = bytes/sizeof(float);
			nzv      = (int) (size/(size_t)(nxv*nyv));
			assert (size%(nxv*nyv) == 0);
			fseek(fp, 0, SEEK_SET);
		}
		/* read file */

		if (vx == 1) n1=nxv;
		if (vy == 1) n1=nyv;
		if (vz == 1) n1=nzv;

		fseek ( fp, 0, SEEK_END );
		bytes = ftell(fp); 
		ntraces  = (int) (bytes/(n1*sizeof(float)));

		if (bytes != nxv*nyv*nzv*sizeof(float)) {
			fprintf(stderr,"ERROR: Number of bytes(%d) in file not equal to nx*ny*nz*4: %d*%d*%d*4=%d bytes\n", bytes, nxv,nyv,nzv,nxv*nyv*nzv*4);
			return -1;
		}
		fseek(fp, 0, SEEK_SET);
		j=0;
		dum = (float *)malloc(nxv*nyv*nzv*sizeof(float));
		while (j<ntraces) {
			nread = fread( &dum[j*n1], sizeof(float), n1, fp );
			assert (nread == n1);
			j++;
		}
		fclose(fp);

	}
	else { /* SU-file */
		/* read first header */
		nread = fread( &hdr1, 1, TRCBYTES, fp );
		assert (nread == TRCBYTES);
		trace_sz = sizeof(float)*hdr1.ns+TRCBYTES;

		if (hdr1.scalco < 0) scl = 1.0/fabs(hdr1.scalco);
		else if (hdr1.scalco == 0) scl = 1.0;
		else scl = hdr1.scalco;

		/* read second header */
		dum = (float *)malloc(hdr1.ns*sizeof(float));
		nread = fread( &hdr2, 1, TRCBYTES, fp );
		assert (nread == TRCBYTES);
		free(dum);

		fseek ( fp, 0, SEEK_END );
		bytes = ftell(fp); 
		ntraces  = (int) (bytes/trace_sz);

		/* read last header */
		fseek ( fp, (long)(bytes-trace_sz), SEEK_SET );
		nread = fread( &hdr3, 1, TRCBYTES, fp );
		assert (nread == TRCBYTES);

		/* depending on vx,vy, vz the header info is intepreted */
		if (vx == 1) {
			if(!getparint("nxv", &nxv)) nxv = hdr1.ns; assert( nxv > 0 );
			if(!getparfloat("xvmin", &xvmin)) xvmin = hdr1.f1; 
			if(!getparfloat("dxv", &dxv)) dxv = hdr1.d1; assert( dxv != 0 );
			if (vy == 2) {
				if(!getparfloat("yvmin", &yvmin)) yvmin = hdr1.gy*scl; 
				if(!getparint("nyv", &nyv)) nyv = hdr1.trwf; assert( nyv > 0 );
				iymin = hdr1.gy;
				iymax = hdr3.gy;
				if(!getparfloat("dyv", &dyv)) dyv = scl*(iymax-iymin)/(nyv-1); 

				if(!getparfloat("dzv", &dzv)) {
					fprintf(stderr,"dzv should be given\n"); 
					exit(1);
				}
				nzv = ntraces/nyv;
				assert (ntraces % nyv == 0);
				if(!getparfloat("zvmin", &zvmin)) zvmin = 0; 
			}
			else {
				if(!getparfloat("dzv", &dzv)) dzv = hdr1.d2;
				if(!getparint("nzv", &nzv)) nzv = hdr1.trwf; assert( nzv > 0 );
				if(!getparfloat("zvmin", &zvmin)) zvmin = 0; 

				nyv = ntraces/nzv;
				assert (ntraces % nzv == 0);
				if(!getparfloat("yvmin", &yvmin)) yvmin = hdr1.gy*scl; 
				iymin = hdr1.gy;
				iymax = hdr3.gy;
				if(!getparfloat("dyv", &dyv)) dyv = scl*(iymax-iymin)/(nyv-1); 
			}
		}
		if (vx == 2) {
			if(!getparfloat("xvmin", &xvmin)) xvmin = hdr1.gx*scl; 
			if(!getparint("nxv", &nxv)) nxv = hdr1.trwf; assert( nxv > 0 );
			ixmin = hdr1.gx;
			ixmax = hdr3.gx;
			if(!getparfloat("dxv", &dxv)) dxv = scl*(ixmax-ixmin)/(nxv-1); 
			if (vy == 1) {
				if(!getparint("nyv", &nyv)) nyv = hdr1.ns; assert( nyv > 0 );
				if(!getparfloat("yvmin", &yvmin)) yvmin = hdr1.f1; 
				if(!getparfloat("dyv", &dyv)) dyv = hdr1.d1; assert( dyv != 0 );

				if(!getparfloat("dzv", &dzv)) {
					fprintf(stderr,"dzv should be given\n"); 
					exit(1);
				}
				nzv = ntraces/nxv;
				assert (ntraces % nxv == 0);
				if(!getparfloat("zvmin", &zvmin)) zvmin = 0; 
			}
			else {
				if(!getparfloat("dzv", &dzv)) dzv = hdr1.d1;
				if(!getparint("nzv", &nzv)) nzv = hdr1.ns;
				if(!getparfloat("zvmin", &zvmin)) zvmin = hdr1.f1; 

				nyv = ntraces/nxv;
				assert (ntraces % nxv == 0);
				if(!getparfloat("yvmin", &yvmin)) yvmin = hdr1.gy*scl; 
				iymin = hdr1.gy;
				iymax = hdr3.gy;
				if(!getparfloat("dyv", &dyv)) dyv = scl*(iymax-iymin)/(nyv-1); 
			}
		}
		if (vx == 3) {
			if (vy == 2) {
				if(!getparfloat("dzv", &dzv)) dzv = hdr1.d1;
				if(!getparint("nzv", &nzv)) nzv = hdr1.ns;
				if(!getparfloat("zvmin", &zvmin)) zvmin = hdr1.f1; 

				if(!getparfloat("yvmin", &yvmin)) yvmin = hdr1.gy*scl; 
				if(!getparint("nyv", &nyv)) nyv = hdr1.trwf; assert( nyv > 0 );
				iymin = hdr1.gy;
				iymax = hdr3.gy;
				if(!getparfloat("dyv", &dyv)) dyv = scl*(iymax-iymin)/(nyv-1); 

				nxv = ntraces/nyv;
				assert (ntraces % nyv == 0);
			}
			else {
				if(!getparint("nyv", &nyv)) nyv = hdr1.ns; assert( nyv > 0 );
				if(!getparfloat("yvmin", &yvmin)) yvmin = hdr1.f1; 
				if(!getparfloat("dyv", &dyv)) dyv = hdr1.d1; assert( dyv != 0 );

				if(!getparfloat("dzv", &dzv)) dzv = hdr1.d2;
				if(!getparint("nzv", &nzv)) nzv = hdr1.trwf; assert( nzv > 0 );
				if(!getparfloat("zvmin", &zvmin)) zvmin = 0; 

				nxv = ntraces/nzv;
				assert (ntraces % nzv == 0);
			}

			if(!getparfloat("xvmin", &xvmin)) xvmin = hdr1.gx*scl; 
			ixmin = hdr1.gx;
			ixmax = hdr3.gx;
			if(!getparfloat("dxv", &dxv)) dxv = scl*(ixmax-ixmin)/(nxv-1); 
		}

		/* read SU file */

		fseek(fp, 0, SEEK_SET);
		j=0;
		dum = (float *)malloc(nxv*nyv*nzv*sizeof(float));
		while (j<ntraces) {
			nread = fread( &hdr1, 1, TRCBYTES, fp );
			if (nread == 0) break;
			nread = fread( &dum[j*hdr1.ns], sizeof(float), hdr1.ns, fp );
			assert (nread == hdr1.ns);
			j++;
		}
		fclose(fp);
	}
	}

	if (verbose) {
		fprintf(stderr," xvmin = %.2f nxv = %d dxv = %.2f xvmax = %.2f\n", 
			xvmin, nxv, dxv, xvmin + (nxv-1)*dxv );
		fprintf(stderr," yvmin = %.2f nyv = %d dyv = %.2f yvmax = %.2f\n", 
			yvmin, nyv, dyv, yvmin + (nyv-1)*dyv );
		fprintf(stderr," zvmin = %.2f nzv = %d dzv = %.2f zvmax = %.2f\n", 
			zvmin, nzv, dzv, zvmin + (nzv-1)*dzv );
	}

	/* write temporary binary file which is stored per depth slice */

	sprintf(tmpname ,"%s/velocity%d.bin",tmp_dir,getpid());
	fprintf(stderr,"temporary file name at pe %d - %s\n", pe, tmpname);
	fp = fopen( tmpname, "w+" ); 

	zslice = (float *)malloc(nxv*nyv*sizeof(float));
	nxy = nxv*nyv;

	for (iz=0; iz<nzv; iz++) {
		if (vx == 1) {
			if (vy == 2) {
				for (j=0; j<nyv; j++) {
					for (i=0; i<nxv; i++) {
						zslice[j*nxv+i] = dum[iz*nxv*nyv+j*nxv+i];
					}
				}
			}
			else {
				for (j=0; j<nyv; j++) {
					for (i=0; i<nxv; i++) {
						zslice[j*nxv+i] = dum[j*nxv*nzv+iz*nxv+i];
					}
				}
			}
		}
		if (vx == 2) {
			if (vy == 1) {
				for (j=0; j<nyv; j++) {
					for (i=0; i<nxv; i++) {
						zslice[j*nxv+i] = dum[iz*nyv*nxv+i*nyv+j];
					}
				}
			}
			else {
				for (j=0; j<nyv; j++) {
					for (i=0; i<nxv; i++) {
						zslice[j*nxv+i] = dum[j*nzv*nxv+i*nzv+iz];
					}
				}
			}
		}
		if (vx == 3) {
			if (vy == 2) {
				for (j=0; j<nyv; j++) {
					for (i=0; i<nxv; i++) {
						zslice[j*nxv+i] = dum[i*nzv*nyv+j*nzv+iz];
					}
				}
			}
			else {
				for (j=0; j<nyv; j++) {
					for (i=0; i<nxv; i++) {
						zslice[j*nxv+i] = dum[i*nyv*nzv+iz*nyv+j];
					}
				}
			}
		}

		nwrite = fwrite( zslice, sizeof(float), nxy, fp );
		assert (nwrite == nxy);
	}

	free(zslice);
	free(dum);
	fflush(fp);

	*vel_file = fp;

	vel_area->xmin = xvmin;
	vel_area->ymin = yvmin;
	vel_area->zmin = zvmin;
	vel_area->ixmin = 0;
	vel_area->ixmax = nxv-1;
	vel_area->iymin = 0;
	vel_area->iymax = nyv-1;
	vel_area->nx    = nxv;
	vel_area->ny    = nyv;
	vel_area->nz    = nzv;
	vel_area->dx    = dxv;
	vel_area->dy    = dyv;
	vel_area->dz    = dzv;
	vel_area->sxy   = nxv*nyv;

	return 0;
}

void readVelocitySlice(FILE *fp, float *velocity, int iz, int nyv, int nxv)
{
	size_t offset, nread, pos;
	int nxy, j, i, err;
	int bytes;

	nxy = nxv*nyv;

	/* Read depth slices */
	offset = iz*nxy*sizeof(float);
	err = fseek(fp, offset, SEEK_SET);
	nread = fread( velocity, sizeof(float), nxy, fp );
	assert (nread == nxy);

	return;
}
