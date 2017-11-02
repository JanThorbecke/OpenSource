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

int readModel(modPar mod, float *velocity, float *slowness);

int defineSource(wavPar wav, srcPar src, modPar mod, float **src_nwav, int reverse, int verbose);

int writeSrcRecPos(modPar *mod, recPar *rec, srcPar *src, shotPar *shot);

/* Self documentation */
char *sdoc[] = {
" ",
" raytime - Jesper Spetzler ray-trace modeling ",
" ",
" IO PARAMETERS:",
"   file_cp= .......... P (cp) velocity file",
"   file_src= ......... file with source signature",
"   file_rcv=recv.su .. base name for receiver files",
"   dx= ............... read from model file: if dx==0 then dx= can be used to set it",
"   dz= ............... read from model file: if dz==0 then dz= can be used to set it",
"" ,
" RAY TRACING PARAMETERS:",
"   smoothwindow=0 .... if set lenght of 2/3D smoothing window on slowness",
"   useT2=0 ........... 1: compute more accurate T2 pertubation correction",
"   geomspread=1 ...... 1: compute Geometrical Spreading Factor",
"   nraystep=5 ........ number of points on ray",
" OPTIONAL PARAMETERS:",
"   ischeme=3 ......... 1=acoustic, 2=visco-acoustic 3=elastic, 4=visco-elastic",
"   sinkdepth=0 ....... receiver grid points below topography (defined bij cp=0.0)",
"   sinkdepth_src=0 ... source grid points below topography (defined bij cp=0.0)",
"   sinkvel=0 ......... use velocity of first receiver to sink through to next layer",
"   verbose=0 ......... silent mode; =1: display info",
" ",
" SHOT AND GENERAL SOURCE DEFINITION:",
"   xsrc=middle ....... x-position of (first) shot ",
"   zsrc=zmin ......... z-position of (first) shot ",
"   nshot=1 ........... number of shots to model",
"   dxshot=dx ......... if nshot > 1: x-shift in shot locations",
"   dzshot=0 .......... if nshot > 1: z-shift in shot locations",
"   xsrca= ............ defines source array x-positions",
"   zsrca= ............ defines source array z-positions",
"   wav_random=1 ...... 1 generates (band limited by fmax) noise signatures ",
"   src_at_rcv=1 ...... inject wavefield at receiver coordinates (1), inject at source (0)",
"" ,
" PLANE WAVE SOURCE DEFINITION:",
"   plane_wave=0 ...... model plane wave with nsrc= sources",
"   nsrc=1 ............ number of sources per (plane-wave) shot ",
"",
" RANDOM SOURCE DEFINITION FOR SEISMIC INTERFEROMTERY:",
"   src_random=0 ...... 1 enables nsrc random sources positions in one modeling",
"   nsrc=1 ............ number of sources to use for one shot",
"   xsrc1=0 ........... left bound for x-position of sources",
"   xsrc2=0 ........... right bound for x-position of sources",
"   zsrc1=0 ........... left bound for z-position of sources",
"   zsrc2=0 ........... right bound for z-position of sources",
"   amplitude=0 ....... distribution of source amplitudes",
"   distribution=0 .... random function for amplitude and tlength 0=flat 1=Gaussian ",
"   seed=10 ........... seed for start of random sequence ",
"" ,
" RECEIVER SELECTION:",
"   xrcv1=xmin ........ first x-position of linear receiver array(s)",
"   xrcv2=xmax ........ last x-position of linear receiver array(s)",
"   dxrcv=dx .......... x-position increment of receivers in linear array(s)",
"   zrcv1=zmin ........ first z-position of linear receiver array(s)",
"   zrcv2=zrcv1 ....... last z-position of linear receiver array(s)",
"   dzrcv=0.0 ......... z-position increment of receivers in linear array(s)",
"   xrcva= ............ defines receiver array x-positions",
"   zrcva= ............ defines receiver array z-positions",
"   rrcv= ............. radius for receivers on a circle ",
"   arcv= ............. vertical arc-lenght for receivers on a ellipse (rrcv=horizontal)",
"   oxrcv=0.0 ......... x-center position of circle",
"   ozrcv=0.0 ......... z-center position of circle",
"   dphi=2 ............ angle between receivers on circle ",
"   rcv_txt=........... text file with receiver coordinates. Col 1: x, Col. 2: z",
"   rec_ntsam=nt ...... maximum number of time samples in file_rcv files",
"",
"      Jan Thorbecke 2017",
"      TU Delft",
"      E-mail: janth@xs4all.nl ",
"",
NULL};


int main(int argc, char **argv)
{
	modPar mod;
	recPar rec;
	srcPar src;
	shotPar shot;
	rayPar ray;
    float *velocity, *slowness;
	double t0, t1, tinit;
	size_t size;
	int n1, ix, iz, ir, ixshot, izshot;
	int irec;
    fcoord coordsx, coordgx, Time;
    icoord grid;
    float Jr, *ampl, *time;
    segy hdr;
    char filetime[1024], fileamp[1024];
    size_t  nwrite;
	int verbose;
    FILE *fpt, *fpa;

	t0= wallclock_time();
	initargs(argc,argv);
	requestdoc(0);

	if (!getparint("verbose",&verbose)) verbose=0;
	getParameters(&mod, &rec, &src, &shot, &ray, verbose);

	/* allocate arrays for model parameters: the different schemes use different arrays */

	n1 = mod.nz;

	velocity = (float *)calloc(mod.nx*mod.nz,sizeof(float));
    slowness = (float *)calloc(mod.nx*mod.nz,sizeof(float));

	/* read velocity and density files */

	readModel(mod, velocity, slowness);

	/* read and/or define source wavelet(s) */

//	defineSource(wav, src, mod, src_nwav, mod.grid_dir, verbose);

	/* allocate arrays for wavefield and receiver arrays */

	size = shot.n*rec.n;
    time = (float *)calloc(size,sizeof(float));
    ampl = (float *)calloc(size,sizeof(float));

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
    strcpy(filetime, rec.file_rcv);
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
    hdr.ns     = rec.n;

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

#pragma omp parallel for default(shared) \
private (coordgx) \
private (irec,Time,Jr) 
        	for (irec=0; irec<rec.n; irec++) {
            	coordgx.x=mod.x0+rec.xr[irec];
            	coordgx.z=mod.z0+rec.zr[irec];
            	coordgx.y = 0;
	
            	getWaveParameter(slowness, grid, mod.dx, coordsx, coordgx, ray, &Time, &Jr);

            	time[((izshot*shot.nx)+ixshot)*rec.n + irec] = Time.x + Time.y + Time.z;
            	ampl[((izshot*shot.nx)+ixshot)*rec.n + irec] = Jr;
            	if (verbose>4) vmess("shot=%f,%f receiver at %f,%f T0=%f T1=%f T2=%f Jr=%f",coordsx.x, coordsx.z, coordgx.x, coordgx.z, Time.x, Time.y, Time.z, Jr); 
        	}

        	hdr.sx     = 1000*(mod.x0+mod.dx*shot.x[ixshot]);
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
        	}
    	}
	} /* end of loop over number of shots */
	fclose(fpt);
	if (ray.geomspread) fclose(fpa);

	t1= wallclock_time();
	if (verbose) {
		vmess("Total compute time ray-tracing = %.2f s.", t1-t0);
	}

	/* free arrays */

	initargs(argc,argv); /* this will free the arg arrays declared */
	free(velocity);
	free(slowness);
	
	return 0;
}
