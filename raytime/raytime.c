#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include <genfft.h>
#include"par.h"
#include"raytime.h"
#include "segy.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

double wallclock_time(void);

void name_ext(char *filename, char *extension);

void threadAffinity(void);

int getParameters(modPar *mod, recPar *rec, srcPar *src, shotPar *shot, rayPar *ray, int verbose);

int getWaveParameter(float *slowness, icoord size, float dgrid, fcoord s, fcoord r, rayPar ray, fcoord *T, float *Jr);

void applyMovingAverageFilter(float *slowness, icoord size, int window, int dim, float *averageModel);

int readModel(modPar mod, float *velocity, float *slowness, int nw);

int defineSource(wavPar wav, srcPar src, modPar mod, float **src_nwav, int reverse, int verbose);

int writeSrcRecPos(modPar *mod, recPar *rec, srcPar *src, shotPar *shot);

void vidale(float *ttime, float *slow, icoord *isrc, icoord grid, float dx, int order, int mzrcv);

/* Self documentation */
char *sdoc[] = {
" ",
" raytime - Jesper Spetzler ray-trace modeling ",
" ",
" IO PARAMETERS:",
"   file_cp= .......... P (cp) velocity file",
"   file_src= ......... file with source signature",
"   file_rcv=recv.su .. base name for receiver files",
"   file_rcvtime= ..... receiver file in x-t",
"   dx= ............... read from model file: if dx==0 then dx= can be used to set it",
"   dz= ............... read from model file: if dz==0 then dz= can be used to set it",
"   nt=1024 ........... number of time-samples in file_rcvtime",
"" ,
" RAY TRACING PARAMETERS:",
"   method=jesper ..... calculation method (jesper, fd) ",
"   smoothwindow=0 .... if set lenght of 2/3D smoothing window on slowness",
"   useT2=0 ........... 1: compute more accurate T2 pertubation correction",
"   geomspread=1 ...... 1: compute Geometrical Spreading Factor",
"   nraystep=5 ........ number of points on ray",
"   sbox=1 ............ radius of inner straight ray (fd method)",
" OPTIONAL PARAMETERS:",
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
    float *velocity, *slowness, *smooth, *trueslow, **inter;
	double t0, t1, t2, tinit, tray, tio;
	size_t size;
	int nw, n1, ix, iz, ir, ixshot, izshot;
	int nt, ntfft, nfreq, ig;
	int irec, sbox, ipos, nrx, nrz, nr;
    fcoord coordsx, coordgx, Time;
    icoord grid, isrc; 
    float Jr, *ampl, *time, *ttime, *ttime_p, cp_average, *wavelet, dw, dt;
	float dxrcv, dzrcv, rdelay, tr;
    segy hdr;
    char filetime[1024], fileamp[1024], *method, *file_rcvtime, *file_src;
    size_t  nwrite, nread;
	int verbose;
    complex *cmute, *cwav;
    FILE *fpt, *fpa, *fpwav, *fprcv;

	t0= wallclock_time();
	initargs(argc,argv);
	requestdoc(0);

	if(!getparint("verbose",&verbose)) verbose=0;
    if(!getparint("sbox", &sbox)) sbox = 1;
    if(!getparstring("method", &method)) method="jesper";
	if (!getparfloat("rec_delay",&rdelay)) rdelay=0.0;

	getParameters(&mod, &rec, &src, &shot, &ray, verbose);

    /* read file_src if file_rcvtime is defined */

    if (!getparstring("file_rcvtime",&file_rcvtime)) file_rcvtime=NULL;

	if (file_rcvtime != NULL) {
    	if (!getparstring("file_src",&file_src)) file_src=NULL;
		if (!getparfloat("dt",&dt)) dt=0.004;
		if (file_src != NULL ) {
        	fpwav = fopen( file_src, "r" );
        	assert( fpwav != NULL);
        	nread = fread( &hdr, 1, TRCBYTES, fpwav );
        	assert(nread == TRCBYTES);
			ntfft = optncr(MAX(hdr.ns, rec.nt));
 			wavelet = (float *)calloc(ntfft,sizeof(float));
			/* read first trace */
        	nread = fread(wavelet, sizeof(float), hdr.ns, fpwav);
        	assert (nread == hdr.ns);
        	fclose(fpwav);
		}
		else {
			ntfft = optncr(rec.nt);
 			wavelet = (float *)calloc(ntfft,sizeof(float));
			wavelet[0] = 1.0;
		}
    	nfreq = ntfft/2+1;
    	cwav    = (complex *)calloc(nfreq,sizeof(complex));
    	cmute   = (complex *)calloc(nfreq,sizeof(complex));
        rc1fft(wavelet,cwav,ntfft,-1);
    	dw      = 2*M_PI/(ntfft*dt);
	}

	/* allocate arrays for model parameters: the different schemes use different arrays */

	n1 = mod.nz;
    if(!strcmp(method,"fd")) nw = 0;
    else nw = ray.smoothwindow;

	velocity = (float *)calloc(mod.nx*mod.nz,sizeof(float));
	slowness = (float *)calloc((mod.nx+2*nw)*(mod.nz+2*nw),sizeof(float));
    trueslow = (float *)calloc(mod.nx*mod.nz,sizeof(float));

    if(!strcmp(method,"fd")) {
		ttime = (float *)calloc(mod.nx*mod.nz,sizeof(float));
	}

	/* read velocity and density files */
	readModel(mod, velocity, slowness, nw);

	/* allocate arrays for wavefield and receiver arrays */

	size = shot.n*rec.n;
    time = (float *)calloc(size,sizeof(float));
    ampl = (float *)calloc(size,sizeof(float));

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
		vmess("   - method for ray-tracing       = %s", method);
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
    if ( nw != 0 ) { /* smooth slowness */ 
        smooth = (float *)calloc(grid.x*grid.z,sizeof(float));
        applyMovingAverageFilter(slowness, grid, nw, 2, smooth);
        memcpy(slowness,smooth,grid.x*grid.z*sizeof(float));
        free(smooth);
	}

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
	if (file_rcvtime != NULL) {
        fprcv = fopen(file_rcvtime, "w");
        assert(fprcv != NULL);
	}

    memset(&hdr,0,sizeof(hdr));
    hdr.scalco = -1000;
    hdr.scalel = -1000;
    hdr.trid   = 1;

	t1=wallclock_time();
	tinit = t1-t0;
    tray=0.0;
    tio=0.0;

	/* Outer loop over number of shots */
	for (izshot=0; izshot<shot.nz; izshot++) {
		for (ixshot=0; ixshot<shot.nx; ixshot++) {

	        t2=wallclock_time();
        	if (verbose) {
            	vmess("Modeling source %d at gridpoints ix=%d iz=%d", (izshot*shot.n)+ixshot, shot.x[ixshot], shot.z[izshot]);
            	vmess(" which are actual positions x=%.2f z=%.2f", mod.x0+mod.dx*shot.x[ixshot], mod.z0+mod.dz*shot.z[izshot]);
            	vmess("Receivers at gridpoint x-range ix=%d - %d", rec.x[0], rec.x[rec.n-1]);
            	vmess(" which are actual positions x=%.2f - %.2f", mod.x0+rec.xr[0], mod.x0+rec.xr[rec.n-1]);
            	vmess("Receivers at gridpoint z-range iz=%d - %d", rec.z[0], rec.z[rec.n-1]);
            	vmess(" which are actual positions z=%.2f - %.2f", mod.z0+rec.zr[0], mod.z0+rec.zr[rec.n-1]);
        	}

        	coordsx.x = shot.x[ixshot]*mod.dx;
        	coordsx.z = shot.z[izshot]*mod.dz;
        	coordsx.y = 0;

	        t1=wallclock_time();
            tio += t1-t2;

            if (!strcmp(method,"jesper")) {
#pragma omp parallel for default(shared) \
private (coordgx,irec,Time,Jr) 
        		for (irec=0; irec<rec.n; irec++) {
            		coordgx.x=rec.xr[irec];
            		coordgx.z=rec.zr[irec];
            		coordgx.y = 0;
		
            		getWaveParameter(slowness, grid, mod.dx, coordsx, coordgx, ray, &Time, &Jr);
	
            		time[((izshot*shot.nx)+ixshot)*rec.n + irec] = Time.x + Time.y + 0.5*Time.z;
            		ampl[((izshot*shot.nx)+ixshot)*rec.n + irec] = Jr;
            		if (verbose>4) vmess("JS: shot=%f,%f receiver at %f,%f T0=%f T1=%f T2=%f Jr=%f",coordsx.x, coordsx.z, coordgx.x, coordgx.z, Time.x, Time.y, Time.z, Jr); 
        		}
			}
            else if(!strcmp(method,"fd")) {
	            int mzrcv;

                isrc.x = shot.x[ixshot];
                isrc.y = 0;
                isrc.z = shot.z[izshot];

                mzrcv = 0;
                for (irec = 0; irec < rec.n; irec++) mzrcv = MAX(rec.z[irec], mzrcv);

                vidale(ttime,slowness,&isrc,grid,mod.dx,sbox, mzrcv);
        		for (irec=0; irec<rec.n; irec++) {
            		coordgx.x=mod.x0+rec.xr[irec];
            		coordgx.z=mod.z0+rec.zr[irec];
            		coordgx.y = 0;
					ipos = ((izshot*shot.nx)+ixshot)*rec.n + irec;
	
            		time[ipos] = ttime[rec.z[irec]*mod.nx+rec.x[irec]];
					/* compute average velocity between source and receiver */
 					nrx = (rec.x[irec]-isrc.x);
 					nrz = (rec.z[irec]-isrc.z);
 					nr = abs(nrx) + abs(nrz);
					cp_average = 0.0;
        			for (ir=0; ir<nr; ir++) {
						ix = isrc.x + floor((ir*nrx)/nr);
						iz = isrc.z + floor((ir*nrz)/nr);
						//fprintf(stderr,"ir=%d ix=%d iz=%d velocity=%f\n", ir, ix, iz, velocity[ix*mod.nz+iz]);
						cp_average += velocity[ix*mod.nz+iz];
					}
					cp_average = cp_average/((float)nr);
            		ampl[ipos] = sqrt(time[ipos]*cp_average);
            		if (verbose>4) vmess("FD: shot=%f,%f receiver at %f(%d),%f(%d) T=%e V=%f Ampl=%f",coordsx.x, coordsx.z, coordgx.x, rec.x[irec], coordgx.z, rec.z[irec], time[ipos], cp_average, ampl[ipos]); 
        		}
            }
	        t2=wallclock_time();
            tray += t2-t1;

        	hdr.sx     = 1000*(mod.x0+mod.dx*shot.x[ixshot]);
        	hdr.sdepth = 1000*(mod.z0+mod.dz*shot.z[izshot]);
        	hdr.selev  = (int)(-1000.0*(mod.z0+mod.dz*shot.z[izshot]));
        	hdr.fldr   = ((izshot*shot.nx)+ixshot)+1;
        	hdr.tracl  = ((izshot*shot.nx)+ixshot)+1;
        	hdr.tracf  = ((izshot*shot.nx)+ixshot)+1;
        	hdr.ntr    = shot.n;
    		hdr.dt     = (unsigned short)1;
    		hdr.trwf   = shot.n;
    		hdr.ns     = rec.n;
        	//hdr.d1     = (rec.x[1]-rec.x[0])*mod.dx; // discrete
        	hdr.d1     = (rec.xr[1]-rec.xr[0]);
        	hdr.f1     = mod.x0+rec.x[0]*mod.dx;
        	hdr.d2     = (shot.x[MIN(shot.n-1,1)]-shot.x[0])*mod.dx;
        	hdr.f2     = mod.x0+shot.x[0]*mod.dx;
			dt_tmp = (fabs(hdr.d1*((float)hdr.scalco)));
			hdr.dt	   = (unsigned short)dt_tmp;
    
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
			if (file_rcvtime != NULL) {
    			hdr.ns     = rec.nt;
    			hdr.trwf   = rec.n;
    			hdr.ntr    = ((izshot*shot.nx)+ixshot+1)*rec.n;
    			hdr.dt     = dt*1000000;
    			hdr.d1     = dt;
        		hdr.f1     = 0.0;
        		hdr.d2     = (rec.xr[1]-rec.xr[0]);
    			hdr.f2     = mod.x0+rec.x[0]*mod.dx;

        		for (irec=0; irec<rec.n; irec++) {
					ipos = ((izshot*shot.nx)+ixshot)*rec.n + irec;
        			hdr.tracf  = irec+1;
        			hdr.tracl  = ((izshot*shot.nx)+ixshot*shot.nz)+irec+1;
        			hdr.gx     = 1000*(mod.x0+rec.xr[irec]);
        			hdr.offset = (rec.xr[irec]-shot.x[ixshot]*mod.dx);
        			hdr.gelev  = (int)(-1000*(mod.z0+rec.zr[irec]));

					tr = time[ipos]+rdelay;
                    for (ig=0; ig<nfreq; ig++) {
                        cmute[ig].r = (cwav[ig].r*cos(ig*dw*tr-M_PI/4.0)-cwav[ig].i*sin(ig*dw*tr-M_PI/4.0))/(ntfft*ampl[ipos]);
                        cmute[ig].i = (cwav[ig].i*cos(ig*dw*tr-M_PI/4.0)+cwav[ig].r*sin(ig*dw*tr-M_PI/4.0))/(ntfft*ampl[ipos]);
                    }
                	cr1fft(cmute,wavelet,ntfft,-1);
        			nwrite = fwrite( &hdr, 1, TRCBYTES, fprcv);
        			nwrite = fwrite( wavelet, sizeof(float), rec.nt, fprcv );
				}
			}
	        t1=wallclock_time();
            tio += t1-t2;
    	} /* end of ixshot loop */
	} /* end of loop over number of shots */
	fclose(fpt);
	if (file_rcvtime != NULL) fclose(fprcv);
	if (ray.geomspread) fclose(fpa);

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

	initargs(argc,argv); /* this will free the arg arrays declared */
	free(velocity);
	free(slowness);
	
	return 0;
}
