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

int getParameters(modPar *mod, recPar *rec, snaPar *sna, wavPar *wav, srcPar *src, shotPar *shot, bndPar *bnd, rayPar *ray, int verbose);

int getWaveParameter(float *slowness, icoord size, float dgrid, fcoord s, fcoord r, rayPar ray, fcoord *T, float *Jr);

int readModel(modPar mod, bndPar bnd, float *velocity, float *slowness);

int defineSource(wavPar wav, srcPar src, modPar mod, float **src_nwav, int reverse, int verbose);

int writeSrcRecPos(modPar *mod, recPar *rec, srcPar *src, shotPar *shot);

int raytime(float *time, float *ampl, int *xnx, float *xrcv, float *xsrc, float *zsrc)
{
	modPar mod;
	recPar rec;
	snaPar sna;
	wavPar wav;
	srcPar src;
	bndPar bnd;
	shotPar shot;
	rayPar ray;
	float **src_nwav;
	float *rox, *roz, *l2m, *lam, *mul;
	float *tss, *tes, *tep, *p, *q, *r;
	float *vx, *vz, *tzz, *txz, *txx;
	float *rec_vx, *rec_vz, *rec_p;
    float *velocity, *slowness;
	float *rec_txx, *rec_tzz, *rec_txz;
	float *rec_pp, *rec_ss;
	float *rec_udp, *rec_udvz;
	float *beam_vx, *beam_vz, *beam_p;
	float *beam_txx, *beam_tzz, *beam_txz;
	float *beam_pp, *beam_ss;	
	float sinkvel;
	double t0, t1, t2, t3, tt, tinit;
	size_t size, sizem, nsamp;
	int n1, ix, iz, ir, ixshot, izshot, i;
	int ioPx, ioPz;
	int it0, it1, its, it, fileno, isam;
	int ixsrc, izsrc, irec;
    int nRayStep;
    fcoord coordsx, coordgx, Time;
    icoord grid;
    float Jr;
    segy hdr;
    char filetime[1024], fileamp[1024];
    size_t  nwrite;
	int verbose;
    FILE *fpt, *fpa;
    double ddt;

	t0= wallclock_time();
	//initargs(argc,argv);
	//requestdoc(0);

	//if (!getparint("verbose",&verbose)) verbose=0;
	getParameters(&mod, &rec, &sna, &wav, &src, &shot, &bnd, &ray, verbose);

	/* allocate arrays for model parameters: the different schemes use different arrays */

	n1 = mod.nz;
	sizem=mod.nx*mod.nz;

	velocity = (float *)calloc(mod.nx*mod.nz,sizeof(float));
    slowness = (float *)calloc(mod.nx*mod.nz,sizeof(float));

	/* read velocity and density files */

	readModel(mod, bnd, velocity, slowness);

	/* read and/or define source wavelet(s) */

//	defineSource(wav, src, mod, src_nwav, mod.grid_dir, verbose);

	/* allocate arrays for wavefield and receiver arrays */

	size = shot.n*rec.n;
    //time = (float *)calloc(size,sizeof(float));
    //ampl = (float *)calloc(size,sizeof(float));

	t1= wallclock_time();
	if (verbose) {
		tinit = t1-t0;
		vmess("*******************************************");
		vmess("************* runtime info ****************");
		vmess("*******************************************");
		vmess("CPU time for intializing arrays and model = %f", tinit);
	}

	/* Sinking source and receiver arrays: 
	   If P-velocity==0 the source and receiver 
	   postions are placed deeper until the P-velocity changes. 
	   The free-surface position is stored in bnd.surface[ix].
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
	}*/

    /*ddt        = (double)mod.dt;
    hdr.dt     = (unsigned short)lround((((double)1.0e6*ddt*rec.skipdt)));
    hdr.scalco = -1000;
    hdr.scalel = -1000;
    hdr.trid   = 1;
    hdr.trwf   = shot.n;
    hdr.ns     = rec.n;*/

	/* Outer loop over number of shots */
	for (izshot=0; izshot<shot.nz; izshot++) {
		for (ixshot=0; ixshot<shot.nx; ixshot++) {

        	if (verbose) {
				vmess("Modeling source %d at gridpoints ix=%d iz=%d", (izshot*shot.n)+ixshot, shot.x[ixshot], shot.z[izshot]);
                vmess(" which are actual positions x=%.2f z=%.2f", mod.x0+mod.dx*shot.x[ixshot], mod.z0+mod.dz*shot.z[izshot]);
            	vmess("Receivers at gridpoint x-range ix=%d - %d", rec.x[0], rec.x[rec.n-1]);
            	vmess(" which are actual positions x=%.2f - %.2f", mod.x0+rec.xr[0], mod.x0+rec.xr[rec.n-1]);
            	vmess("Receivers at gridpoint z-range iz=%d - %d", rec.z[0], rec.z[rec.n-1]);
            	vmess(" which are actual positions z=%.2f - %.2f", mod.z0+rec.zr[0], mod.z0+rec.zr[rec.n-1]);
        	}

        	coordsx.x = mod.x0+shot.x[ixshot]*mod.dx;
        	coordsx.z = mod.z0+shot.z[izshot]*mod.dz;
        	coordsx.y = 0;
        	grid.x = mod.nx;
        	grid.z = mod.nz;
        	grid.y = 1;

			xnx[(izshot*shot.nx)+ixshot]  = rec.n;
        	xsrc[(izshot*shot.nx)+ixshot] = mod.x0+mod.dx*shot.x[ixshot];
        	zsrc[(izshot*shot.nx)+ixshot] = mod.z0+mod.dz*shot.z[izshot];

        	for (irec=0; irec<rec.n; irec++) {
            	coordgx.x=mod.x0+rec.xr[irec];
            	coordgx.z=mod.z0+rec.zr[irec];
            	coordgx.y = 0;

            	getWaveParameter(slowness, grid, mod.dx, coordsx, coordgx, ray, &Time, &Jr);

				xrcv[((izshot*shot.nx)+ixshot)*rec.n + irec] = (mod.x0+rec.x[0]*mod.dx) + ((rec.x[1]-rec.x[0])*mod.dx*((float)irec));
            	time[((izshot*shot.nx)+ixshot)*rec.n + irec] = Time.x + Time.y + Time.z;
            	ampl[((izshot*shot.nx)+ixshot)*rec.n + irec] = Jr;
            	fprintf(stderr,"shot=%f,%f receiver at %f,%f T0=%f T1=%f T2=%f Jr=%f\n",coordsx.x, coordsx.z, coordgx.x, coordgx.z, Time.x, Time.y, Time.z, Jr); 
        	}

        /*hdr.sx     = 1000*(mod.x0+mod.dx*shot.x[ishot]);
        hdr.sdepth = 1000*(mod.z0+mod.dz*shot.z[ishot]);
        hdr.selev  = (int)(-1000.0*(mod.z0+mod.dz*shot.z[ishot]));
        hdr.fldr   = ishot+1;
        hdr.tracl  = ishot+1;
        hdr.tracf  = ishot+1;
        hdr.ntr    = shot.n;
        hdr.d1     = (rec.x[1]-rec.x[0])*mod.dx;
        hdr.f1     = mod.x0+rec.x[0]*mod.dx;
        hdr.d2     = (shot.x[1]-shot.x[0])*mod.dx;
        hdr.f2     = mod.x0+shot.x[0]*mod.dx;
    
        nwrite = fwrite( &hdr, 1, TRCBYTES, fpt);
        assert(nwrite == TRCBYTES);
        nwrite = fwrite( &time[ishot*rec.n], sizeof(float), rec.n, fpt);
        assert(nwrite == rec.n);
	    fflush(fpt);
	    if (ray.geomspread) {
            nwrite = fwrite( &hdr, 1, TRCBYTES, fpa);
            assert(nwrite == TRCBYTES);
            nwrite = fwrite( &ampl[ishot*rec.n], sizeof(float), rec.n, fpa);
            assert(nwrite == rec.n);
	        fflush(fpa);
        }*/
    	}
	} /* end of loop over number of shots */
	//fclose(fpt);
	//if (ray.geomspread) fclose(fpa);

	t1= wallclock_time();
	if (verbose) {
		vmess("Total compute time ray-tracing = %.2f s.", t1-t0);
	}

	/* free arrays */

	//initargs(argc,argv); /* this will free the arg arrays declared */
	free(velocity);
	free(slowness);
	
	return 0;
}
