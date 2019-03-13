#include<mpi.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<assert.h>
#include<sys/time.h>
#include"par.h"
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"fdelmodc3D.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))


#define STRIPE_COUNT "4" /* must be an ascii string */
#define STRIPE_SIZE "1048576" /* 1 MB must be an ascii string */
//#define STRIPE_SIZE "268435456" /* 256 MB must be an ascii string */
#define C1 (9.0/8.0)
#define C2 (1.0/24.0)
#define Dx(f,ix,iy,iz,nz) C1*(f[iy+ix*nz+iz] - f[iy+(ix-1)*nz+iz]) - C2*(f[iy+(ix+1)*nz+iz] - f[iy+(ix-2)*nz+iz])
#define Dy(f,ix,iy,iz,nxz) C1*(f[iy*nxz+ix+iz] - f[(iy-1)*nxz+ix+iz]) - C2*(f[(iy+1)*nxz+ix+iz] - f[(iy-2)*nxz+ix+iz])
#define Dz(f,ix,iy,iz) C1*(f[iy+ix+iz] - f[iy+ix+iz-1]) - C2*(f[iy+ix+iz+1] - f[iy+ix+iz-2])

#define Dv(vx,vz,ix,iy,iz,nz,nxz) C1*((vx[iy*nxz+(ix+1)*nz+iz] - vx[iy*nxz+ix*nz+iz]) + \
				      (vy[(iy+1)*nxz+ix*nz+iz] - vy[iy*nxz+ix*nz+iz]) + \
				      (vz[iy*nxz+ix*nz+iz+1]   - vz[iy*nxz+ix*nz+iz])) - \
				  C2*((vx[iy*nxz+(ix+2)*nz+iz] - vx[iy*nxz+(ix-1)*nz+iz]) + \
				      (vy[(iy+2)*nxz+ix*nz+iz] - vy[(iy-1)*nxz+ix*nz+iz]) + \
				      (vz[iy*nxz+ix*nz+iz+2]   - vz[iy*nxz+ix*nz+iz-1]))

int getParameters3D(modPar *mod, recPar *rec, snaPar *sna, wavPar *wav, srcPar *src, shotPar *shot, bndPar *bnd, int verbose);
int readModel3D(modPar mod, bndPar bnd, float *rox, float *roz, float *l2m);

void vinit();
int updateVelocitiesHalo(float *vx, float *vy, float *vz, float *p, float *ro, int halo, int npx, int npy, int npz);
int updateVelocities(float *vx, float *vy, float *vz, float *p, float *ro, int halo, int npx, int npy, int npz);
int updatePressureHalo(float *vx, float *vy, float *vz, float *p, float *l2m, int halo, int npx, int npy, int npz);
int updatePressure(float *vx, float *vy, float *vz, float *p, float *l2m, int halo, int npx, int npy, int npz);
int exchangeHalo(float *leftRecv, float *leftSend, float *rightRecv, float *rightSend, int size, int leftrank, int rightrank, MPI_Request *reqRecv, MPI_Request *reqSend, int tag);
int newHaloVxVz(float *vx, float *vz, int npx, int npy, int npz, int halo, float *leftRecv, float *rightRecv, float *frontRecv, float *backRecv, float *topRecv, float *bottomRecv);
int newHaloP(float *p, int npx, int npy, int npz, int halo, float *leftRecv, float *rightRecv, float *frontRecv, float *backRecv, float *topRecv, float *bottomRecv);
int copyHaloVxVz(float *vx, float *vz, int npx, int npy, int npz, int halo, float *leftSend, float *rightSend, float *frontSend, float *backSend, float *topSend, float *bottomSend);
int copyHaloP(float *p, int npx, int npy, int npz, int halo, float *leftSend, float *rightSend, float *frontSend, float *backSend, float *topSend, float *bottomSend);
int waitForHalo(MPI_Request *reqRecv, MPI_Request *reqSend);
float gauss2time(float t, float f, float t0);
double wallclock_time(void);
void name_ext(char *filename, char *extension);


/* Self documentation */
char *sdoc[] = {
" ",
"   file_rcv=recv.su .. base name for receiver files",
"   file_snap=snap.su . base name for snapshot files",
"   nx=256 ............ number of samples in x-direction",
"   ny=nx ............. number of samples in y-direction",
"   nz=nx ............. number of samples in z-direction",
"   dx=5 .............. spatial sampling in x-direction",
"   dy=5 .............. spatial sampling in y-direction",
"   dz=5 .............. spatial sampling in z-direction",
"" ,
"   verbose=0 ......... silent mode; =1: display info",
" ",
"      Jan Thorbecke 2016",
"      Cray / TU Delft",
"      E-mail: janth@xs4all.nl ",
"",
NULL};

int main (int argc, char *argv[])
{
	modPar mod;
	recPar rec;
	snaPar sna;
	wavPar wav;
	srcPar src;
	bndPar bnd;
	shotPar shot;
	float *wavelet;
	int nx, ny, nz, dims[3], period[3], reorder, coord[3], ndims=3;
	int npx, npy, npz, halo, nt;
	int my_rank, size, source, dest, snapwritten;
	int left, right, front, back, top, bottom;
	int direction, displ, halosizex, halosizey, halosizez;
	int ioXx, ioXz, ioYx, ioYz, ioZz, ioZx, ioPx, ioPz;
	int it, ix, iy, iz, iyp, ixp, izp, isrc, ixsrc, iysrc, izsrc, c1, c2;
	int sizes[3], subsizes[3], starts[3];
	int gsizes[3], gsubsizes[3], gstarts[3];
	int error, rc, verbose;
	float fx, fy, fz, dx, dy, dz, flx, fly, flz;
	float *p, *vx, *vy, *vz, *rox, *roz, *roy, *l2m, hcp, hro, fac;
	float *leftRecv, *leftSend, *rightRecv, *rightSend;
	float *frontRecv, *frontSend, *backRecv, *backSend;
	float *topRecv, *topSend, *bottomRecv, *bottomSend;
	float dt, src_ampl, fmax, fpeaksrc, t0src, time, snaptime;
	double t00, t0, t1, t2, tcomm, tcomp, thcomp, tcopy, ttot, tio;
	char err_buffer[MPI_MAX_ERROR_STRING];
	int resultlen;
	MPI_Comm COMM_CART;
	MPI_Request reqSend[6], reqRecv[6];
	MPI_Status status[12];
	MPI_Datatype local_array, global_array;
	MPI_Offset disp;
	MPI_Info fileinfo;
	MPI_File fh;
	char filename[1000], *file_snap, *file_rcv;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	vinit();

	t0= wallclock_time();
	initargs(argc,argv);
	requestdoc(0);

	if (!getparint("verbose",&verbose)) verbose=0;
	if (!getparstring("file_snap",&file_snap)) file_snap="snap.su";
	if (!getparstring("file_rcv",&file_rcv)) file_rcv="recv.su";

	getParameters3D(&mod, &rec, &sna, &wav, &src, &shot, &bnd, verbose);

	/* divide 3D cube among available processors */
	dims[0]=0; 
	dims[1]=0;
	dims[2]=0;
	MPI_Dims_create(size, ndims, dims);

	/* dims contains the number of MPI-tasks in each direction */
	/* set number of grid points based on number of procs in dims */
	if (!getparint("nx",&nx)) nx=256;
	if (!getparint("ny",&ny)) ny=nx;
	if (!getparint("nz",&nz)) nz=nx;

	if (!getparfloat("dx",&dx)) dx=5.;
	if (!getparfloat("dy",&dy)) dy=5.;
	if (!getparfloat("dz",&dz)) dz=5.;

	halo = 2;

	/* new larger dimensions to fit with the domain-decomposition */
	nz=dims[0]*ceil(mod.naz/dims[0]);
	nx=dims[1]*ceil(mod.nax/dims[1]); 
	ny=dims[2]*ceil(mod.nay/dims[2]); 

//	dt=0.001;
	nt=4096;
	t0src=0.50;
	hcp=2000.0;
	hro=1800.0;
	tcomm=tcomp=thcomp=tcopy=tio=0.0;

	/* for stability 10 points per wavelenght */
	fmax=hcp/(mod.dx*8);
	dt=0.295*mod.dx/hcp;
	fpeaksrc=0.2*fmax; /* Ricker wavelet has peak at 1/3 of fmax */
	fac = mod.dt/mod.dx;

	fx=-mod.dx*nx/2; fy=-mod.dy*ny/2; fz=0;
	npz = 2*halo+nz/dims[0];
	npx = 2*halo+nx/dims[1];
	npy = 2*halo+ny/dims[2];
	wavelet = (float *)calloc(nt,sizeof(float));

	/* find which MPI-task has the source position */


	snaptime = t0src+1.80*npx*dx*0.5/hcp;
	snapwritten=0;
	nt = (int) 1.1*(t0src+snaptime)/dt;
	nt = (int) (t0src+1.5)/dt;

	if (verbose && my_rank==0) {
		fprintf(stderr,"fmax=%f fpeak=%f dt=%e\n", fmax, fpeaksrc, dt);
		fprintf(stderr,"nx=%d nprocx=%d ny=%d nprocy=%d nz=%d nprocz=%d\n", nx, dims[1], ny, dims[2], nz, dims[0]);
		fprintf(stderr,"npx=%d npy=%d npz=%d nt=%d time=%f\n", npx, npy, npz, nt, nt*dt);
		fprintf(stderr,"source expected at local boundary at %f seconds\n", npx*dx*0.5/hcp);
		fflush(stderr);
	}

	if (my_rank==0) fprintf(stderr,"memory per MPI task is %ld MB\n", (6*npx*npy*npz*4/(1024*1024)));

	/* allocate wavefields and medium properties for local grid domains */
	p   = (float *)calloc(npx*npy*npz,sizeof(float));
	vx  = (float *)calloc(npx*npy*npz,sizeof(float));
	vy  = (float *)calloc(npx*npy*npz,sizeof(float));
	vz  = (float *)calloc(npx*npy*npz,sizeof(float));

    /* read velocity and density files */

    readModel3D(mod, bnd, rox, roz, l2m);

/* for 2.5 D model npy=1 */
	rox = (float *)calloc(npx*npy*npz,sizeof(float));
	roy= (float *)calloc(npx*npy*npz,sizeof(float));
	roz= (float *)calloc(npx*npy*npz,sizeof(float));
	l2m = (float *)calloc(npx*npy*npz,sizeof(float));

	/* define homogeneus model */
	for (ix=0; ix<npx*1*npz; ix++) {
		rox[ix]  = fac/hro;
		roy[ix]  = fac/hro;
		roz[ix]  = fac/hro;
		l2m[ix] = fac*hcp*hcp*hro;
	}

	/* create cartesian domain decomposition */
	period[0]=0; 
	period[1]=0;
	period[2]=0;
	reorder=0;
	MPI_Cart_create(MPI_COMM_WORLD, 3, dims, period, reorder, &COMM_CART);

	/* find out coordinates of the rank */
	MPI_Cart_coords(COMM_CART, my_rank, 3, coord);
	flz = fz+(dz*nz/dims[0])*coord[0];
	flx = fx+(dx*nx/dims[1])*coord[1];
	fly = fy+(dy*ny/dims[2])*coord[2];
	if (verbose>=2) fprintf(stderr,"Rank %d coordinates are %d %d %d orig=(%5.2F, %5.2f, %5.2f) \n", my_rank, coord[0], coord[1], coord[2], flx, fly, flz);
	fflush(stderr);

	/* find out neighbours of the rank, MPI_PROC_NULL is a hard boundary of the model */ 
	displ=1;
	MPI_Cart_shift(COMM_CART, 1, 1, &left, &right);
	MPI_Cart_shift(COMM_CART, 2, 1, &top, &bottom);
	MPI_Cart_shift(COMM_CART, 0, 1, &front, &back);
	if (verbose>=2) fprintf(stderr, "Rank %d in direction 0 has LR neighbours %d %d FB %d %d TB %d %d\n", my_rank, left, right, front, back, top, bottom);
	fflush(stderr);

	/* allocate of halo areas */
	halosizex = npy*npz*halo;
	leftRecv  = (float *)calloc(3*halosizex,sizeof(float));
	rightRecv = (float *)calloc(3*halosizex,sizeof(float));
	leftSend  = (float *)calloc(3*halosizex,sizeof(float));
	rightSend = (float *)calloc(3*halosizex,sizeof(float));

	halosizey = npx*npz*halo;
	frontRecv = (float *)calloc(3*halosizey,sizeof(float));
	backRecv  = (float *)calloc(3*halosizey,sizeof(float));
	frontSend = (float *)calloc(3*halosizey,sizeof(float));
	backSend  = (float *)calloc(3*halosizey,sizeof(float));

	halosizez = npy*npx*halo;
	topRecv    = (float *)calloc(3*halosizez,sizeof(float));
	bottomRecv = (float *)calloc(3*halosizez,sizeof(float));
	topSend    = (float *)calloc(3*halosizez,sizeof(float));
	bottomSend = (float *)calloc(3*halosizez,sizeof(float));

	if (my_rank==0) fprintf(stderr,"memory per MPI task for halo exchange is %ld MB\n", ((12*(halosizex+halosizey+halosizez))*4/(1024*1024)));

	/* create subarrays(excluding halo areas) to write to file with MPI-IO */
	/* data in the local array */
	sizes[0]=npz; 
	sizes[1]=npx; 
	sizes[2]=npy;
	subsizes[0]=sizes[0]-2*halo; 
	subsizes[1]=sizes[1]-2*halo;  
	subsizes[2]=sizes[2]-2*halo;
	starts[0]=halo; 
	starts[1]=halo; 
	starts[2]=halo; 
	MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, 
                            MPI_FLOAT, &local_array); 
	MPI_Type_commit(&local_array);

	/* data in the global array */
	gsizes[0]=nz; 
	gsizes[1]=nx; 
	gsizes[2]=ny;
	gsubsizes[0]=subsizes[0]; 
	gsubsizes[1]=subsizes[1];  
	gsubsizes[2]=subsizes[2];
	gstarts[0]=subsizes[0]*coord[0]; 
	gstarts[1]=subsizes[1]*coord[1]; 
	gstarts[2]=subsizes[2]*coord[2]; 
	MPI_Type_create_subarray(3, gsizes, gsubsizes, gstarts, MPI_ORDER_C, 
                            MPI_FLOAT, &global_array); 
	MPI_Type_commit(&global_array);


	/* compute field of the inner grid excluding halo areas */
	ioXx=2;
	ioXz=ioXx-1;
	ioYx=2;
	ioYz=ioYx-1;
	ioZz=2;
	ioZx=ioZz-1;
	ioPx=1;
	ioPz=ioPx;

	t00 = wallclock_time();
	for (it=0; it<nt; it++) {
		time = it*dt;
		wavelet[it] = gauss2time(time,fpeaksrc,t0src);
	}
	if (my_rank==0) {
		FILE *fp;
		fp = fopen("src.bin", "w+");
		fwrite( wavelet, sizeof(float), nt, fp);
		fflush(fp);
		fclose(fp);
	}

/*
	nt =1;
			sprintf(filename,"snap_nz%d_nx%d_ny%d.bin",nz, nx, ny);

	for (ix=0; ix<npx*npy*npz; ix++) {
		p[ix]  = my_rank;
	}
			MPI_Info_create(&fileinfo);
			MPI_Info_set(fileinfo, "striping_factor", STRIPE_COUNT);
			MPI_Info_set(fileinfo, "striping_unit", STRIPE_SIZE);
			MPI_File_delete(filename, MPI_INFO_NULL);
			rc = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR|MPI_MODE_CREATE, fileinfo, &fh);
			if (rc != MPI_SUCCESS) {
  				fprintf(stderr, "could not open input file\n");
  				MPI_Abort(MPI_COMM_WORLD, 2);
			}
			disp = 0;
			rc = MPI_File_set_view(fh, disp, MPI_FLOAT, global_array, "native", fileinfo);
			if (rc != MPI_SUCCESS) {
				fprintf(stderr, "error setting view on results file\n");
				MPI_Abort(MPI_COMM_WORLD, 4);
			}
			rc = MPI_File_write_all(fh, p, 1, local_array, status);
			if (rc != MPI_SUCCESS) {
				MPI_Error_string(rc,err_buffer,&resultlen);
				fprintf(stderr,err_buffer);
				MPI_Abort(MPI_COMM_WORLD, 5);
			}
			MPI_File_close(&fh);
*/


	/* Main loop over the number of time steps */
	for (it=0; it<nt; it++) {
		t0 = wallclock_time();
		time = it*dt;
		//fprintf(stderr,"modeling time step %d for time %f\n", it, time);

		/* define source wavelet */
		wavelet[it] = gauss2time(time,fpeaksrc,t0src);

		/* update of grid values on halo areas */
		updateVelocitiesHalo(vx, vy, vz, p, rox, halo, npx, npy, npz);
		t1 = wallclock_time();
		thcomp += t1-t0;

		/* copy of halo areas  */
		copyHaloVxVz(vx, vz, npx, npy, npz, halo, leftSend, rightSend, frontSend, backSend, topSend, bottomSend);
		t2 = wallclock_time();
		tcopy += t2-t1;

		/* start a-synchronous communication of halo areas to neighbours */
		/* this is done first for Vx,Vz fields only */
		exchangeHalo(leftRecv, leftSend, rightRecv, rightSend, 2*halosizex, left, right, &reqRecv[0], &reqSend[0], 0);
		exchangeHalo(frontRecv, frontSend, backRecv, backSend, 2*halosizey, front, back, &reqRecv[2], &reqSend[2], 4);
		exchangeHalo(topRecv, topSend, bottomRecv, bottomSend, 2*halosizez, top, bottom, &reqRecv[4], &reqSend[4], 8);
		t1 = wallclock_time();
		tcomm += t1-t2;

		/* compute update on grid values excluding halo areas */
		updateVelocities(vx, vy, vz, p, rox, halo, npx, npy, npz);
		t2 = wallclock_time();
		tcomp += t2-t1;

		/* wait for Vx.Vz halo exchange */
		waitForHalo(&reqRecv[0], &reqSend[0]);
		t1 = wallclock_time();
		tcomm += t1-t2;

		/* copy of halo areas  back to arrays */
		newHaloVxVz(vx, vz, npx, npy, npz, halo, leftRecv, rightRecv, frontRecv, backRecv, topRecv, bottomRecv);
		t2 = wallclock_time();
		tcopy += t2-t1;
	
		/* add Force source on the Vz grid */
		src_ampl = wavelet[it];
	
		/* check if source position is in local domain */
		/* for the moment place a source in the middle of each domain */
		ixsrc = npx/2;
		iysrc = npy/2;
		izsrc = npz/2;
		isrc  = iysrc*npx*npz+ixsrc*npz+izsrc;
//		fprintf(stderr,"npz=%d npx=%d npy=%d isrc=%d\n", npz, npx, npy, isrc);
	
		/* source scaling factor to compensate for discretisation */
		src_ampl *= rox[isrc]*l2m[isrc]/(dt);
	
		/* Force source */
		//if (my_rank == 0) vz[isrc] += 0.25*src_ampl*ro[isrc]*dz;
		vz[isrc] += 0.25*src_ampl*rox[isrc]*dz;
	
		/* compute field on the grid of the halo areas */
		updatePressureHalo(vx, vy, vz, p, l2m, halo, npx, npy, npz);
		t1 = wallclock_time();
		thcomp += t1-t2;

		/* copy p-field and sent to neighbours */
		copyHaloP(p, npx, npy, npz, halo, leftSend, rightSend, frontSend, backSend, topSend, bottomSend);
		exchangeHalo(leftRecv, leftSend, rightRecv, rightSend, halosizex, left, right, &reqRecv[0], &reqSend[0], 0);
		exchangeHalo(frontRecv, frontSend, backRecv, backSend, halosizey, front, back, &reqRecv[2], &reqSend[2], 4);
		exchangeHalo(topRecv, topSend, bottomRecv, bottomSend, halosizez, top, bottom, &reqRecv[4], &reqSend[4], 8);
		t2 = wallclock_time();
		tcomm += t2-t1;

		/* compute update on grid values excluding halo areas */
		updatePressure(vx, vy, vz, p, l2m, halo, npx, npy, npz);
		t1 = wallclock_time();
		tcomp += t1-t2;

		/* wait for P halo exchange */
		waitForHalo(&reqRecv[0], &reqSend[0]);
		t2 = wallclock_time();
		tcomm += t2-t1;

		newHaloP(p, npx, npy, npz, halo, leftRecv, rightRecv, frontRecv, backRecv, topRecv, bottomRecv);
		t1 = wallclock_time();
		tcopy += t1-t2;
	
//		fprintf(stderr,"rank %d did time step %d in %f seconds\n", my_rank, it, t1-t0);
//		fflush(stderr);

		/* write snapshots to file */
//		if (time >= snaptime && !snapwritten) {
		if ((it+1)%100==0 ) {

			t1 = wallclock_time();
			sprintf(filename,"snap_nz%d_nx%d_ny%d_it%4d.bin",nz, nx, ny, it);

			MPI_Info_create(&fileinfo);
			MPI_Info_set(fileinfo, "striping_factor", STRIPE_COUNT);
			MPI_Info_set(fileinfo, "striping_unit", STRIPE_SIZE);
			MPI_File_delete(filename, MPI_INFO_NULL);
			rc = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR|MPI_MODE_CREATE, fileinfo, &fh);
			if (rc != MPI_SUCCESS) {
  				fprintf(stderr, "could not open input file\n");
  				MPI_Abort(MPI_COMM_WORLD, 2);
			}
			disp = 0;
			rc = MPI_File_set_view(fh, disp, MPI_FLOAT, global_array, "native", fileinfo);
			if (rc != MPI_SUCCESS) {
				fprintf(stderr, "error setting view on results file\n");
				MPI_Abort(MPI_COMM_WORLD, 4);
			}
			rc = MPI_File_write_all(fh, p, 1, local_array, status);
			if (rc != MPI_SUCCESS) {
				MPI_Error_string(rc,err_buffer,&resultlen);
				fprintf(stderr,err_buffer);
				MPI_Abort(MPI_COMM_WORLD, 5);
			}
			MPI_File_close(&fh);


/*			MPI_Info_create(&fileinfo);
			MPI_File_delete(filename, MPI_INFO_NULL);
			MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
			MPI_File_set_view(fh, 0, MPI_FLOAT, global_array, "native", MPI_INFO_NULL);
			MPI_File_write_all(fh, p, npz*npx*npy, local_array, status);
			MPI_File_close(&fh);
*/
			snapwritten+=1;
			t2 = wallclock_time();
			tio += t2-t1;
		}

	}
	ttot = wallclock_time() - t00;

	if (my_rank == 0) {
		fprintf(stderr,"rank %d total time in %f seconds\n", my_rank, ttot);
		fprintf(stderr,"rank %d comm  time in %f seconds\n", my_rank, tcomm);
		fprintf(stderr,"rank %d comp  time in %f seconds\n", my_rank, tcomp);
		fprintf(stderr,"rank %d hcomp time in %f seconds\n", my_rank, thcomp);
		fprintf(stderr,"rank %d copy  time in %f seconds\n", my_rank, tcopy);
		fprintf(stderr,"rank %d io    time in %f seconds\n", my_rank, tio);
		fprintf(stderr,"rank %d snaphsots written to file\n", snapwritten);
	}


	MPI_Finalize();
	return 0;
}



int updateVelocities(float *vx, float *vy, float *vz, float *p, float *ro, int halo, int npx, int npy, int npz)
{
	int ix, iy, iz, iyp, ixp, izp, c1, c2, nxz;
	int ixs, ixe, iys, iye, izs, ize;
	float DpDx, DpDy, DpDz;

	nxz=npx*npz;
	c1 = 9.0/8.0;
	c2 = -1.0/24.0;

	ixs=2*halo; ixe=npx-2*halo;
	iys=2*halo; iye=npy-2*halo;
	izs=2*halo; ize=npz-2*halo;

	/* calculate vx,vy,vz for all grid points except on the virtual boundary*/
#pragma omp for private (iy, ix, iz) nowait
#pragma ivdep
	for (iy=iys; iy<iye; iy++) {
		iyp=iy*nxz;
		for (ix=ixs; ix<ixe; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=izs; iz<ize; iz++) {
				DpDx = Dx(p,ix,iyp,iz,npz);
				DpDy = Dy(p,ixp,iy,iz,nxz);
				DpDz = Dz(p,ixp,iyp,iz);

				vz[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDz;
				vx[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDx;
				vy[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDy;
			}
		}
	}

	return 0;
}

int updateVelocitiesHalo(float *vx, float *vy, float *vz, float *p, float *ro, int halo, int npx, int npy, int npz)
{
	int ix, iy, iz, iyp, ixp, izp, c1, c2, nxz;
	float DpDx, DpDy, DpDz;

	nxz=npx*npz;
	c1 = 9.0/8.0;
	c2 = -1.0/24.0;

	/* calculate vx,vy,vz for all halo grid points */

	/* compute halo areas at left side */
#pragma omp for private (iy, ix, iz) nowait
	for (iy=halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		for (ix=halo; ix<2*halo; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				DpDx = Dx(p,ix,iyp,iz,npz);
				DpDy = Dy(p,ixp,iy,iz,nxz);
				DpDz = Dz(p,ixp,iyp,iz);

				vz[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDz;
				vx[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDx;
				vy[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDy;
			}
		}
	}

	/* compute halo areas at right side */
	for (iy=halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		for (ix=npx-2*halo; ix<npx-halo; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				DpDx = Dx(p,ix,iyp,iz,npz);
				DpDy = Dy(p,ixp,iy,iz,nxz);
				DpDz = Dz(p,ixp,iyp,iz);

				vz[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDz;
				vx[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDx;
				vy[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDy;
			}
		}
	}


	/* compute halo areas at front side */
	for (iy=halo; iy<2*halo; iy++) {
		iyp=iy*nxz;
		for (ix=2*halo; ix<npx-2*halo; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				DpDx = Dx(p,ix,iyp,iz,npz);
				DpDy = Dy(p,ixp,iy,iz,nxz);
				DpDz = Dz(p,ixp,iyp,iz);

				vz[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDz;
				vx[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDx;
				vy[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDy;
			}
		}
	}

	/* compute halo areas at back side */
	for (iy=npy-2*halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		for (ix=2*halo; ix<npx-2*halo; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				DpDx = Dx(p,ix,iyp,iz,npz);
				DpDy = Dy(p,ixp,iy,iz,nxz);
				DpDz = Dz(p,ixp,iyp,iz);

				vz[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDz;
				vx[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDx;
				vy[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDy;
			}
		}
	}

	/* compute halo areas at top side */
	for (iy=2*halo; iy<npy-2*halo; iy++) {
		iyp=iy*nxz;
		for (ix=2*halo; ix<npx-2*halo; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<2*halo; iz++) {
				DpDx = Dx(p,ix,iyp,iz,npz);
				DpDy = Dy(p,ixp,iy,iz,nxz);
				DpDz = Dz(p,ixp,iyp,iz);

				vz[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDz;
				vx[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDx;
				vy[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDy;
			}
		}
	}

	/* compute halo areas at bottom side */
	for (iy=2*halo; iy<npy-2*halo; iy++) {
		iyp=iy*nxz;
		for (ix=2*halo; ix<npx-2*halo; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=npz-2*halo; iz<npz-halo; iz++) {
				DpDx = Dx(p,ix,iyp,iz,npz);
				DpDy = Dy(p,ixp,iy,iz,nxz);
				DpDz = Dz(p,ixp,iyp,iz);

				vz[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDz;
				vx[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDx;
				vy[iyp+ixp+iz] += ro[iyp+ixp+iz]*DpDy;
			}
		}
	}

	return 0;
}

int updatePressure(float *vx, float *vy, float *vz, float *p, float *l2m, int halo, int npx, int npy, int npz)
{
	int ix, iy, iz, iyp, ixp, izp, c1, c2, nxz;
	int ixs, ixe, iys, iye, izs, ize;

	nxz=npx*npz;
	c1 = 9.0/8.0;
	c2 = -1.0/24.0;

	ixs=2*halo; ixe=npx-2*halo;
	iys=2*halo; iye=npy-2*halo;
	izs=2*halo; ize=npz-2*halo;

/* calculate p/tzz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iz)
#pragma ivdep
	for (iy=iys; iy<iye; iy++) {
		iyp=iy*nxz;
		for (ix=ixs; ix<ixe; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=izs; iz<ize; iz++) {
				p[iyp+ixp+iz] += l2m[iyp+ixp+iz]*(Dv(vx,vz,ix,iy,iz,npz,nxz));
			}
		}
	}

	return 0;
}

int updatePressureHalo(float *vx, float *vy, float *vz, float *p, float *l2m, int halo, int npx, int npy, int npz)
{
	int ix, iy, iz, iyp, ixp, izp, c1, c2, nxz;

	nxz=npx*npz;
	c1 = 9.0/8.0;
	c2 = -1.0/24.0;

	/* calculate p/tzz for all grid points except on the virtual boundary */

	/* compute halo areas at left side */
#pragma omp for private (iy, ix, iz) nowait
	for (iy=halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		for (ix=halo; ix<2*halo; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				p[iyp+ixp+iz] += l2m[iyp+ixp+iz]*(Dv(vx,vz,ix,iy,iz,npz,nxz));
			}
		}
	}

	/* compute halo areas at right side */
	for (iy=halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		for (ix=npx-2*halo; ix<npx-halo; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				p[iyp+ixp+iz] += l2m[iyp+ixp+iz]*(Dv(vx,vz,ix,iy,iz,npz,nxz));
			}
		}
	}


	/* compute halo areas at front side */
	for (iy=halo; iy<2*halo; iy++) {
		iyp=iy*nxz;
		for (ix=2*halo; ix<npx-2*halo; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				p[iyp+ixp+iz] += l2m[iyp+ixp+iz]*(Dv(vx,vz,ix,iy,iz,npz,nxz));
			}
		}
	}

	/* compute halo areas at back side */
	for (iy=npy-2*halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		for (ix=2*halo; ix<npx-2*halo; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				p[iyp+ixp+iz] += l2m[iyp+ixp+iz]*(Dv(vx,vz,ix,iy,iz,npz,nxz));
			}
		}
	}

	/* compute halo areas at top side */
	for (iy=2*halo; iy<npy-2*halo; iy++) {
		iyp=iy*nxz;
		for (ix=2*halo; ix<npx-2*halo; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<2*halo; iz++) {
				p[iyp+ixp+iz] += l2m[iyp+ixp+iz]*(Dv(vx,vz,ix,iy,iz,npz,nxz));
			}
		}
	}

	/* compute halo areas at bottom side */
	for (iy=2*halo; iy<npy-2*halo; iy++) {
		iyp=iy*nxz;
		for (ix=2*halo; ix<npx-2*halo; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=npz-2*halo; iz<npz-halo; iz++) {
				p[iyp+ixp+iz] += l2m[iyp+ixp+iz]*(Dv(vx,vz,ix,iy,iz,npz,nxz));
			}
		}
	}

	return 0;
}

int exchangeHalo(float *leftRecv, float *leftSend, float *rightRecv, float *rightSend, int size, int leftrank, int rightrank, MPI_Request *reqRecv, MPI_Request *reqSend, int tag)
{
	int error, my_rank, ltag;
	MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if (leftrank != MPI_PROC_NULL) {
		ltag = tag;
		error = MPI_Irecv(leftRecv, size, MPI_FLOAT, leftrank, ltag, MPI_COMM_WORLD, &reqRecv[0]);
		assert (error == MPI_SUCCESS);
//		fprintf(stderr,"rank %d recv data from %d left\n", my_rank, leftrank);
		ltag = tag+1;
		error = MPI_Isend(leftSend, size, MPI_FLOAT, leftrank, ltag, MPI_COMM_WORLD, &reqSend[0]);
		assert (error == MPI_SUCCESS);
//		fprintf(stderr,"rank %d send data to %d left\n", my_rank, leftrank);
	}
	else {
		reqRecv[0] = MPI_REQUEST_NULL;
		reqSend[0] = MPI_REQUEST_NULL;
	}

	if (rightrank != MPI_PROC_NULL) {
		ltag = tag+1;
		error = MPI_Irecv(rightRecv, size, MPI_FLOAT, rightrank, ltag, MPI_COMM_WORLD, &reqRecv[1]);
//		fprintf(stderr,"rank %d recv data from %d right\n", my_rank, rightrank);
		assert (error == MPI_SUCCESS);
		ltag = tag;
		error = MPI_Isend(rightSend, size, MPI_FLOAT, rightrank, ltag, MPI_COMM_WORLD, &reqSend[1]);
		assert (error == MPI_SUCCESS);
//		fprintf(stderr,"rank %d send data to %d right\n", my_rank, rightrank);
	}
	else {
		reqRecv[1] = MPI_REQUEST_NULL;
		reqSend[1] = MPI_REQUEST_NULL;
	}

	return 0;
}

int waitForHalo(MPI_Request *reqRecv, MPI_Request *reqSend)
{
	int i;
	MPI_Status status;
	int error;

	for (i=0; i<6; i++) {
		error = MPI_Wait(&reqSend[i], &status);
		assert (error == MPI_SUCCESS);
	}

//	MPI_Barrier(MPI_COMM_WORLD);

	for (i=0; i<6; i++) {
		error = MPI_Wait(&reqRecv[i], &status);
		assert (error == MPI_SUCCESS);
	}

	return 0;
}

int copyHaloVxVz(float *vx, float *vz, int npx, int npy, int npz, int halo, float *leftSend, float *rightSend, float *frontSend, float *backSend, float *topSend, float *bottomSend)
{
	int ix, iy, iz, ih, iyp, ixp, halosizex, halosizey, halosizez, nxz;

	nxz = npx*npz;

	/* copy halo areas at left side */
	halosizex = npy*npz*halo;
	for (iy=halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		for (ix=halo; ix<2*halo; ix++) {
			ixp=ix*npz;
			ih=(ix-halo)*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				leftSend[iy*npz*halo+ih+iz]             = vx[iyp+ixp+iz];
				leftSend[halosizex+iy*npz*halo+ih+iz]   = vz[iyp+ixp+iz];
			}
		}
	}

	/* copy halo areas at right side */
	halosizex = npy*npz*halo;
	for (iy=halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		for (ix=npx-2*halo; ix<npx-halo; ix++) {
			ixp=ix*npz;
			ih=(ix-(npx-2*halo))*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				rightSend[iy*npz*halo+ih+iz]             = vx[iyp+ixp+iz];
				rightSend[halosizex+iy*npz*halo+ih+iz]   = vz[iyp+ixp+iz];
			}
		}
	}


	/* copy halo areas at front side */
	halosizey = npx*npz*halo;
	for (iy=halo; iy<2*halo; iy++) {
		iyp=iy*nxz;
		ih=(iy-halo)*nxz;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				frontSend[ih+ixp+iz]             = vx[iyp+ixp+iz];
				frontSend[halosizey+ih+ixp+iz]   = vz[iyp+ixp+iz];
			}
		}
	}

	/* copy halo areas at back side */
	for (iy=npy-2*halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		ih=(iy-(npy-2*halo))*nxz;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				backSend[ih+ixp+iz]             = vx[iyp+ixp+iz];
				backSend[halosizey+ih+ixp+iz]   = vz[iyp+ixp+iz];
			}
		}
	}

	/* copy halo areas at top side */
	halosizez = npy*npx*halo;
	for (iy=0; iy<npy; iy++) {
		iyp=iy*nxz;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<2*halo; iz++) {
				ih=iz-halo;
				topSend[iy*npx*halo+ix*halo+ih]             = vx[iyp+ixp+iz];
				topSend[halosizez+iy*npx*halo+ix*halo+ih]   = vz[iyp+ixp+iz];
			}
		}
	}

	/* copy halo areas at bottom side */
	for (iy=0; iy<npy; iy++) {
		iyp=iy*nxz;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=npz-2*halo; iz<npz-halo; iz++) {
				ih=(iz-(npz-2*halo));
				bottomSend[iy*npx*halo+ix*halo+ih]             = vx[iyp+ixp+iz];
				bottomSend[halosizez+iy*npx*halo+ix*halo+ih]   = vz[iyp+ixp+iz];
			}
		}
	}

	return 0;
}

int copyHaloP(float *p, int npx, int npy, int npz, int halo, float *leftSend, float *rightSend, float *frontSend, float *backSend, float *topSend, float *bottomSend)
{
	int ix, iy, iz, ih, iyp, ixp, halosizex, halosizey, halosizez, nxz;

	nxz = npx*npz;

	/* copy halo areas at left side */
	halosizex = npy*npz*halo;
	for (iy=halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		for (ix=halo; ix<2*halo; ix++) {
			ixp=ix*npz;
			ih=(ix-halo)*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				leftSend[iy*npz*halo+ih+iz]             = p[iyp+ixp+iz];
			}
		}
	}

	/* copy halo areas at right side */
	halosizex = npy*npz*halo;
	for (iy=halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		for (ix=npx-2*halo; ix<npx-halo; ix++) {
			ixp=ix*npz;
			ih=(ix-(npx-2*halo))*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				rightSend[iy*npz*halo+ih+iz]             = p[iyp+ixp+iz];
			}
		}
	}


	/* copy halo areas at front side */
	halosizey = npx*npz*halo;
	for (iy=halo; iy<2*halo; iy++) {
		iyp=iy*nxz;
		ih=(iy-halo)*nxz;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				frontSend[ih+ixp+iz]             = p[iyp+ixp+iz];
			}
		}
	}

	/* copy halo areas at back side */
	for (iy=npy-2*halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		ih=(iy-(npy-2*halo))*nxz;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				backSend[ih+ixp+iz]             = p[iyp+ixp+iz];
			}
		}
	}

	/* copy halo areas at top side */
	halosizez = npy*npx*halo;
	for (iy=0; iy<npy; iy++) {
		iyp=iy*nxz;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<2*halo; iz++) {
				ih=iz-halo;
				topSend[iy*npx*halo+ix*halo+ih]             = p[iyp+ixp+iz];
			}
		}
	}

	/* copy halo areas at bottom side */
	for (iy=0; iy<npy; iy++) {
		iyp=iy*nxz;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=npz-2*halo; iz<npz-halo; iz++) {
				ih=(iz-(npz-2*halo));
				bottomSend[iy*npx*halo+ix*halo+ih]             = p[iyp+ixp+iz];
			}
		}
	}

	return 0;
}

/* copy communicated halo areas back to compute grids */
int newHaloVxVz(float *vx, float *vz, int npx, int npy, int npz, int halo, float *leftRecv, float *rightRecv, float *frontRecv, float *backRecv, float *topRecv, float *bottomRecv)
{
	int ix, iy, iz, ih, iyp, ixp, halosizex, halosizey, halosizez, nxz;

	nxz = npx*npz;

	/* copy halo areas at left side */
	halosizex = npy*npz*halo;
	for (iy=halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		for (ix=0; ix<halo; ix++) {
			ixp=ix*npz;
			ih=ixp;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				vx[iyp+ixp+iz] = leftRecv[iy*npz*halo+ih+iz];
				vz[iyp+ixp+iz] = leftRecv[halosizex+iy*npz*halo+ih+iz];
			}
		}
	}

	/* copy halo areas at right side */
	halosizex = npy*npz*halo;
	for (iy=halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		for (ix=npx-halo; ix<npx; ix++) {
			ixp=ix*npz;
			ih=(ix-(npx-halo))*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				vx[iyp+ixp+iz] = rightRecv[iy*npz*halo+ih+iz];
				vz[iyp+ixp+iz] = rightRecv[halosizex+iy*npz*halo+ih+iz];
			}
		}
	}


	/* copy halo areas at front side */
	halosizey = npx*npz*halo;
	for (iy=0; iy<halo; iy++) {
		iyp=iy*nxz;
		ih=iyp;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				vx[iyp+ixp+iz] = frontRecv[ih+ixp+iz];
				vz[iyp+ixp+iz] = frontRecv[halosizey+ih+ixp+iz];
			}
		}
	}

	/* copy halo areas at back side */
	for (iy=npy-halo; iy<npy; iy++) {
		iyp=iy*nxz;
		ih=(iy-(npy-halo))*nxz;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				vx[iyp+ixp+iz] = backRecv[ih+ixp+iz];
				vz[iyp+ixp+iz] = backRecv[halosizey+ih+ixp+iz];
			}
		}
	}

	/* copy halo areas at top side */
	halosizez = npy*npx*halo;
	for (iy=0; iy<npy; iy++) {
		iyp=iy*nxz;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=0; iz<halo; iz++) {
				ih=iz;
				vx[iyp+ixp+iz] = topRecv[iy*npx*halo+ix*halo+ih];
				vz[iyp+ixp+iz] = topRecv[halosizez+iy*npx*halo+ix*halo+ih];
			}
		}
	}

	/* copy halo areas at bottom side */
	for (iy=0; iy<npy; iy++) {
		iyp=iy*nxz;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=npz-halo; iz<npz; iz++) {
				ih=(iz-(npz-halo));
				vx[iyp+ixp+iz] = bottomRecv[iy*npx*halo+ix*halo+ih];
				vz[iyp+ixp+iz] = bottomRecv[halosizez+iy*npx*halo+ix*halo+ih];
			}
		}
	}

	return 0;
}

/* copy communicated halo areas back to compute grids */
int newHaloP(float *p, int npx, int npy, int npz, int halo, float *leftRecv, float *rightRecv, float *frontRecv, float *backRecv, float *topRecv, float *bottomRecv)
{
	int ix, iy, iz, ih, iyp, ixp, halosizex, halosizey, halosizez, nxz;

	nxz = npx*npz;

	/* copy halo areas at left side */
	halosizex = npy*npz*halo;
	for (iy=halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		for (ix=0; ix<halo; ix++) {
			ixp=ix*npz;
			ih=ixp;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				p[iyp+ixp+iz] = leftRecv[iy*npz*halo+ih+iz];
			}
		}
	}

	/* copy halo areas at right side */
	halosizex = npy*npz*halo;
	for (iy=halo; iy<npy-halo; iy++) {
		iyp=iy*nxz;
		for (ix=npx-halo; ix<npx; ix++) {
			ixp=ix*npz;
			ih=(ix-(npx-halo))*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				p[iyp+ixp+iz] = rightRecv[iy*npz*halo+ih+iz];
			}
		}
	}


	/* copy halo areas at front side */
	halosizey = npx*npz*halo;
	for (iy=0; iy<halo; iy++) {
		iyp=iy*nxz;
		ih=iyp;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				p[iyp+ixp+iz] = frontRecv[ih+ixp+iz];
			}
		}
	}

	/* copy halo areas at back side */
	for (iy=npy-halo; iy<npy; iy++) {
		iyp=iy*nxz;
		ih=(iy-(npy-halo))*nxz;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=halo; iz<npz-halo; iz++) {
				p[iyp+ixp+iz] = backRecv[ih+ixp+iz];
			}
		}
	}

	/* copy halo areas at top side */
	halosizez = npy*npx*halo;
	for (iy=0; iy<npy; iy++) {
		iyp=iy*nxz;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=0; iz<halo; iz++) {
				ih=iz;
				p[iyp+ixp+iz] = topRecv[iy*npx*halo+ix*halo+ih];
			}
		}
	}

	/* copy halo areas at bottom side */
	for (iy=0; iy<npy; iy++) {
		iyp=iy*nxz;
		for (ix=0; ix<npx; ix++) {
			ixp=ix*npz;
#pragma ivdep
			for (iz=npz-halo; iz<npz; iz++) {
				ih=(iz-(npz-halo));
				p[iyp+ixp+iz] = bottomRecv[iy*npx*halo+ix*halo+ih];
			}
		}
	}

	return 0;
}
float gauss2time(float t, float f, float t0)
{
    float value, time;

	time = t-t0;
    value = ((1.0-2.0*M_PI*M_PI*f*f*time*time)*exp(-M_PI*M_PI*f*f*time*time));
    return value;
}

