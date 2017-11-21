#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include"par.h"
#include"raytime.h"
#include "segy.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

typedef struct _icoord { /* 3D coordinate integer */
    int z;
    int x;
    int y;
} icoord;

typedef struct _fcoord { /* 3D coordinate float */
    float z;
    float x;
    float y;
} fcoord;

double wallclock_time(void);

void name_ext(char *filename, char *extension);

void threadAffinity(void);

int getParameters(modPar *mod, recPar *rec, srcPar *src, shotPar *shot, rayPar *ray, int verbose);

int getWaveParameter(float *slowness, icoord size, float dgrid, fcoord s, fcoord r, rayPar ray, fcoord *T, float *Jr);

void applyMovingAverageFilter(float *slowness, icoord size, int window, int dim, float *averageModel);

int readModel(modPar mod, float *velocity, float *slowness, int nw);

int defineSource(wavPar wav, srcPar src, modPar mod, float **src_nwav, int reverse, int verbose);

int writeSrcRecPos(modPar *mod, recPar *rec, srcPar *src, shotPar *shot);


int raytime(float *time, float *ampl, int *xnx, float *xrcv, float *xsrc, float *zsrc)
{
	modPar mod;
	recPar rec;
	srcPar src;
	shotPar shot;
	rayPar ray;
    float *velocity, *slowness, *smooth;
	double t0, t1, t2, tinit, tray, tio;
	size_t size;
	int nw, n1, ix, iz, ir, ixshot, izshot;
	int irec;
    fcoord coordsx, coordgx, Time;
    icoord grid;
    float Jr;
    segy hdr;
    char filetime[1024], fileamp[1024];
    size_t  nwrite;
	int verbose;
    FILE *fpt, *fpa;

	t0= wallclock_time();
	//initargs(argc,argv);
	requestdoc(0);

	if (!getparint("verbose",&verbose)) verbose=0;
	getParameters(&mod, &rec, &src, &shot, &ray, verbose);

	/* allocate arrays for model parameters: the different schemes use different arrays */

	n1 = mod.nz;
    nw = ray.smoothwindow;

	velocity = (float *)calloc(mod.nx*mod.nz,sizeof(float));
	slowness = (float *)calloc((mod.nx+2*nw)*(mod.nz+2*nw),sizeof(float));
//    slowness = (float *)calloc(mod.nx*mod.nz,sizeof(float));

	/* read velocity and density files */

	readModel(mod, velocity, slowness, nw);

	/* allocate arrays for wavefield and receiver arrays */

	size = shot.n*rec.n;
    //time = (float *)calloc(size,sizeof(float));
    //ampl = (float *)calloc(size,sizeof(float));

	/* Sinking source and receiver arrays: 
	   If P-velocity==0 the source and receiver 
	   postions are placed deeper until the P-velocity changes. 
	   Setting the option rec.sinkvel only sinks the receiver position 
       (not the source) and uses the velocity 
	   of the first receiver to sink through to the next layer. */

/* sink receivers to value different than sinkvel */
	for (ir=0; ir<rec.n; ir++) {
		iz = rec.z[ir];
		ix = rec.x[ir];
		while(velocity[(ix)*n1+iz] == rec.sinkvel) iz++;
		rec.z[ir]=iz+rec.sinkdepth;
		rec.zr[ir]=rec.zr[ir]+(rec.z[ir]-iz)*mod.dz;
//		rec.zr[ir]=rec.z[ir]*mod.dz;
		if (verbose>3) vmess("receiver position %d at grid[ix=%d, iz=%d] = (x=%f z=%f)", ir, ix, rec.z[ir], rec.xr[ir]+mod.x0, rec.zr[ir]+mod.z0);
	}
/*
*/

/* sink sources to value different than zero */
	for (izshot=0; izshot<shot.nz; izshot++) {
		for (ixshot=0; ixshot<shot.nx; ixshot++) {
			iz = shot.z[izshot];
			ix = shot.x[ixshot];
			while(velocity[(ix)*n1+iz] == 0.0) iz++;
			shot.z[izshot]=iz+src.sinkdepth; 
		}
	}

	if (verbose>3) writeSrcRecPos(&mod, &rec, &src, &shot);

    /* smooth slowness grid */
    grid.x = mod.nx;
    grid.z = mod.nz;
    grid.y = 1;
    if ( (ray.smoothwindow) != 0 ) { /* smooth slowness */ 
        smooth = (float *)calloc(grid.x*grid.z,sizeof(float));
        applyMovingAverageFilter(slowness, grid, ray.smoothwindow, 2, smooth);
        memcpy(slowness,smooth,grid.x*grid.z*sizeof(float));
        free(smooth);
	}

    /* prepare output file and headers */
    /*strcpy(filetime, rec.file_rcv);
    name_ext(filetime, "_time");
    fpt = fopen(filetime, "w");
    assert(fpt != NULL);

	if (ray.geomspread) {
        strcpy(fileamp, rec.file_rcv);
        name_ext(fileamp, "_amp");
        fpa = fopen(fileamp, "w");
        assert(fpa != NULL);
	}

    hdr.dt     = (unsigned short)1;
    hdr.scalco = -1000;
    hdr.scalel = -1000;
    hdr.trid   = 1;
    hdr.trwf   = shot.n;
    hdr.ns     = rec.n;*/

	t1=wallclock_time();
	tinit = t1-t0;
    tray=0.0;
    tio=0.0;

	/* Outer loop over number of shots */
	for (izshot=0; izshot<shot.nz; izshot++) {
		for (ixshot=0; ixshot<shot.nx; ixshot++) {

	        t2=wallclock_time();
        	if (verbose) {
            	vmess("Modeling source %d at gridpoints ix=%d iz=%d", (izshot*shot.nx)+ixshot, shot.x[ixshot], shot.z[izshot]);
            	vmess(" which are actual positions x=%.2f z=%.2f", mod.x0+mod.dx*shot.x[ixshot], mod.z0+mod.dz*shot.z[izshot]);
            	vmess("Receivers at gridpoint x-range ix=%d - %d", rec.x[0], rec.x[rec.n-1]);
            	vmess(" which are actual positions x=%.2f - %.2f", mod.x0+rec.xr[0], mod.x0+rec.xr[rec.n-1]);
            	vmess("Receivers at gridpoint z-range iz=%d - %d", rec.z[0], rec.z[rec.n-1]);
            	vmess(" which are actual positions z=%.2f - %.2f", mod.z0+rec.zr[0], mod.z0+rec.zr[rec.n-1]);
        	}

        	coordsx.x = mod.x0+shot.x[ixshot]*mod.dx;
        	coordsx.z = mod.z0+shot.z[izshot]*mod.dz;
        	coordsx.y = 0;

			xnx[ (izshot*shot.nx)+ixshot] = rec.n;
            xsrc[(izshot*shot.nx)+ixshot] = mod.x0+mod.dx*shot.x[ixshot];
            zsrc[(izshot*shot.nx)+ixshot] = mod.z0+mod.dz*shot.z[izshot];

	        t1=wallclock_time();
            tio += t1-t2;
#pragma omp parallel for default(shared) \
private (coordgx,irec,Time,Jr) 
        	for (irec=0; irec<rec.n; irec++) {
            	coordgx.x=mod.x0+rec.xr[irec];
            	coordgx.z=mod.z0+rec.zr[irec];
            	coordgx.y = 0;
	
            	getWaveParameter(slowness, grid, mod.dx, coordsx, coordgx, ray, &Time, &Jr);

				xrcv[((izshot*shot.nx)+ixshot)*rec.n + irec] = (mod.x0+rec.x[0]*mod.dx) + ((rec.x[1]-rec.x[0])*mod.dx*((float)irec));
            	time[((izshot*shot.nx)+ixshot)*rec.n + irec] = Time.x + Time.y + Time.z;
            	ampl[((izshot*shot.nx)+ixshot)*rec.n + irec] = Jr;
            	if (verbose>4) vmess("shot=%f,%f receiver at %f,%f T0=%f T1=%f T2=%f Jr=%f",coordsx.x, coordsx.z, coordgx.x, coordgx.z, Time.x, Time.y, Time.z, Jr); 
        	}
	        t2=wallclock_time();
            tray += t2-t1;

        	/*hdr.sx     = 1000*(mod.x0+mod.dx*shot.x[ixshot]);
        	hdr.sdepth = 1000*(mod.z0+mod.dz*shot.z[izshot]);
        	hdr.selev  = (int)(-1000.0*(mod.z0+mod.dz*shot.z[izshot]));
        	hdr.fldr   = ((izshot*shot.nx)+ixshot)+1;
        	hdr.tracl  = ((izshot*shot.nx)+ixshot)+1;
        	hdr.tracf  = ((izshot*shot.nx)+ixshot)+1;
        	hdr.ntr    = shot.n;
        	hdr.d1     = (rec.x[1]-rec.x[0])*mod.dx;
        	hdr.f1     = mod.x0+rec.x[0]*mod.dx;
        	hdr.d2     = (shot.x[1]-shot.x[0])*mod.dx;
        	hdr.f2     = mod.x0+shot.x[0]*mod.dx;
    
        	nwrite = fwrite( &hdr, 1, TRCBYTES, fpt);
        	assert(nwrite == TRCBYTES);
        	nwrite = fwrite( &time[((izshot*shot.nx)+ixshot)*rec.n], sizeof(float), rec.n, fpt);
        	assert(nwrite == rec.n);
	    	fflush(fpt);
	    	if (ray.geomspread) {
            	nwrite = fwrite( &hdr, 1, TRCBYTES, fpa);
            	assert(nwrite == TRCBYTES);
            	nwrite = fwrite( &ampl[((izshot*shot.nx)+ixshot)*rec.n], sizeof(float), rec.n, fpa);
            	assert(nwrite == rec.n);
	        	fflush(fpa);
        	}*/
	        t1=wallclock_time();
            tio += t1-t2;
    	}
	} /* end of loop over number of shots */
	//fclose(fpt);
	//if (ray.geomspread) fclose(fpa);

	t1= wallclock_time();
	if (verbose) {
		vmess("*******************************************");
		vmess("************* runtime info ****************");
		vmess("*******************************************");
		vmess("Total compute time ray-tracing    = %.2f s.", t1-t0);
		vmess("   - intializing arrays and model = %.3f", tinit);
		vmess("   - ray tracing                  = %.3f", tray);
		vmess("   - writing data to file         = %.3f", tio);
	}

	/* free arrays */

	//initargs(argc,argv); /* this will free the arg arrays declared */
	free(velocity);
	free(slowness);
	
	return 0;
}
