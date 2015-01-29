#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include"par.h"
#include"fdelmodc.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

double wallclock_time(void);

int getEmParameters(modPar *mod, recPar *rec, snaPar *sna, wavPar *wav, srcPar *src, shotPar *shot, bndPar *bnd, int verbose);

int readEmModel(modPar mod, bndPar bnd, float *eprs, float *ksigma, float *mu);

int defineSource(wavPar wav, srcPar src, float **src_nwav, int reverse, int verbose);

int writeSrcRecPos(modPar *mod, recPar *rec, srcPar *src, shotPar *shot);

int em4(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float **src_nwav, float *hz, float *hx, float *Ey, float *eprs, float *ksigma, float *mu, int verbose);

int getRecTimes(modPar mod, recPar rec, bndPar bnd, int itime, int isam, float *vx, float *vz, float *tzz, float *txx, 
	float *txz, float *rec_vx, float *rec_vz, float *rec_txx, float *rec_tzz, float *rec_txz, 
	float *rec_p, float *rec_pp, float *rec_ss, int verbose);

int writeEmRec(recPar rec, modPar mod, int ixsrc, int izsrc, int nsam, int ishot, int fileno, 
			 float *rec_vx, float *rec_vz, float *rec_txx, float *rec_tzz, float *rec_txz, 
			 float *rec_p, float *rec_pp, float *rec_ss, int verbose);

int writeEmSnapTimes(modPar mod, snaPar sna, int ixsrc, int izsrc, int itime, 
				   float *vx, float *vz, float *tzz, float *txx, float *txz, int verbose);

int getBeamTimes(modPar mod, snaPar sna, float *vx, float *vz, float *tzz, float *txx, float *txz, 
				 float *beam_vx, float *beam_vz, float *beam_txx, float *beam_tzz, float *beam_txz, 
				 float *beam_p, float *beam_pp, float *beam_ss, int verbose);

int writeBeams(modPar mod, snaPar sna, int ixsrc, int izsrc, int ishot, int fileno, 
			   float *beam_vx, float *beam_vz, float *beam_txx, float *beam_tzz, float *beam_txz, 
			   float *beam_p, float *beam_pp, float *beam_ss, int verbose);

int taperEdges(modPar mod, bndPar bnd, float *vx, float *vz, int verbose);

/* Self documentation */
char *sdoc[] = {
" ",
" fdemmodc - electro-magnetic GPR finite difference wavefield modeling ",
" ",
" IO PARAMETERS:",
"   file_er= .......... Er file: griddel model for relative permittivity",
"   file_ks= .......... Ksigma file: griddel model for conductivity",
"   file_src= ......... file with source signature",
"   file_rcv=recv.su .. base name for receiver files",
"   file_snap=snap.su . base name for snapshot files",
"   file_beam=beam.su . base name for beam fields ",
"   dx= ............... read from model file: if dx==0 then dx= can be used to set it",
"   dz= ............... read from model file: if dz==0 then dz= can be used to set it",
"   dt= ............... read from file_src: if dt==0 then dt= can be used to set it",
"" ,
" OPTIONAL PARAMETERS:",
//"   ischeme=3 ......... 1=acoustic, 2=visco-acoustic 3=elastic, 4=visco-elastic",
"   tmod=(nt-1)*dt .... total registration time (nt from file_src)",
"   ntaper=0 .......... length of taper at edges of model",
"   tapfact=0.30 ...... taper strength: larger value gets stronger taper",
//"   For the 4 boundaries the options are:  2=pml 3=rigid 4=taper",
//"   top=4 ............ type of boundary on top edge of model",
//"   left=4 ........... type of boundary on left edge of model",
//"   right=4 .......... type of boundary on right edge of model",
//"   bottom=4 ......... type of boundary on bottom edge of model",

//"   tapleft=0 ......... =1: taper left edge of model",
//"   tapright=0 ........ =1: taper right edge of model",
//"   taptop=0 .......... =1: taper top edge of model",
//"   tapbottom=0 ....... =1: taper bottom edge of model",
//"   grid_dir=0 ........ direction of time modeling (1=reverse time)",
"   sinkdepth=0 ....... receiver grid points below topography (defined bij cp=0.0)",
"   sinkdepth_src=0 ... source grid points below topography (defined bij cp=0.0)",
"   sinkvel=0 ......... use velocity of first receiver to sink through to next layer",
"   beam=0 ............ calculate energy beam of wavefield in model",
"   verbose=0 ......... silent mode; =1: display info",
" ",
" SHOT AND GENERAL SOURCE DEFINITION:",
"   src_type=1 ........ 1=Ey 6=Hz 7=Hx",
"   src_orient=1 ...... orientation of the source",
"                     - 1=monopole",
"                     - 2=dipole +/- vertical oriented",
"                     - 3=dipole - + horizontal oriented",
//"                     - 4=dipole +/0/-",
//"                     - 5=dipole + -",
"   xsrc=middle ....... x-position of (first) shot ",
"   zsrc=zmin ......... z-position of (first) shot ",
"   nshot=1 ........... number of shots to model",
"   dxshot=dx ......... if nshot > 1: x-shift in shot locations",
"   dzshot=0 .......... if nshot > 1: z-shift in shot locations",
"   xsrca= ............ defines source array x-positions",
"   zsrca= ............ defines source array z-positions",
"   wav_random=1 ...... 1 generates (band limited by fmax) noise signatures ",
"   fmax=from_src ..... maximum frequency in wavelet",
"   src_multiwav=0 .... use traces in file_src as areal source",
"   src_injectionrate=0 set to 1 to use injection rate source",

"" ,
" PLANE WAVE SOURCE DEFINITION:",
"   plane_wave=0 ...... model plane wave with nsrc= sources",
"   nsrc=1 ............ number of sources per (plane-wave) shot ",
"   src_angle=0 ....... angle of plane source array",
"   src_velo=1500 ..... velocity to use in src_angle definition",
"   src_window=0 ...... length of taper at edges of source array",
"",
" RANDOM SOURCE DEFINITION FOR SEISMIC INTERFEROMTERY:",
"   src_random=0 ...... 1 enables nsrc random sources positions in one modeling",
"   nsrc=1 ............ number of sources to use for one shot",
"   xsrc1=0 ........... left bound for x-position of sources",
"   xsrc2=0 ........... right bound for x-position of sources",
"   zsrc1=0 ........... left bound for z-position of sources",
"   zsrc2=0 ........... right bound for z-position of sources",
"   tsrc1=0.0 ......... begin time interval for random sources being triggered",
"   tsrc2=tmod ........ end time interval for random sources being triggered",
"   tactive=tsrc2 ..... end time for random sources being active",
"   tlength=tsrc2-tsrc1 average duration of random source signal",
"   length_random=1 ... duration of source is rand*tlength",
"   amplitude=0 ....... distribution of source amplitudes",
"   distribution=0 .... random function for amplitude and tlength 0=flat 1=Gaussian ",
"   seed=10 ........... seed for start of random sequence ",
"" ,
" SNAP SHOT SELECTION:",
"   tsnap1=0.1 ........ first snapshot time (s)",
"   tsnap2=0.0 ........ last snapshot time (s)",
"   dtsnap=0.1 ........ snapshot time interval (s)",
"   dxsnap=dx ......... sampling in snapshot in x-direction",
"   xsnap1=0 .......... first x-position for snapshots area",
"   xsnap2=0 .......... last x-position for snapshot area",
"   dzsnap=dz ......... sampling in snapshot in z-direction",
"   zsnap1=0 .......... first z-position for snapshots area",
"   zsnap2=0 .......... last z-position for snapshot area",
"   sna_type_ey=1 ..... Ey registration _sey",
"   sna_type_hx=1 ..... Hx registration _shx",
"   sna_type_hz=0 ..... Hz registration _shz",
"   sna_hzhxtime=0 .... registration of hz/hx times",
"                       The fd scheme is also staggered in time.",
"                       Time at which hz/hx snapshots are written:",
"                     - 0=previous hz/hx relative to ey at time t",
"                     - 1=next     hz/hx relative to ey at time t",
"" ,
" RECEIVER SELECTION:",
"   xrcv1=xmin ........ first x-position of linear receiver array(s)",
"   xrcv2=xmax ........ last x-position of linear receiver array(s)",
"   dxrcv=dx .......... x-position increment of receivers in linear array(s)",
"   zrcv1=zmin ........ first z-position of linear receiver array(s)",
"   zrcv2=zrcv1 ....... last z-position of linear receiver array(s)",
"   dzrcv=0.0 ......... z-position increment of receivers in linear array(s)",
"   dtrcv=.004 ........ desired sampling in receiver data (seconds)",
"   max_nrec=10000 .... maximum number of receivers",
"   xrcva= ............ defines receiver array x-positions",
"   zrcva= ............ defines receiver array z-positions",
"   rrcv= ............. radius for receivers on a circle ",
"   oxrcv=0.0 ......... x-center position of circle",
"   ozrcv=0.0 ......... z-center position of circle",
"   dphi=2 ............ angle between receivers on circle ",
"   rec_ntsam=nt ...... maximum number of time samples in file_rcv files",
"   rec_delay=0 ....... time in seconds to start recording",
"   rec_type_ey=1 ..... Ey registration _rey",
"   rec_type_hx=1 ..... Vz registration _rhx",
"   rec_type_hz=0 ..... Vx registration _rhz",
"   rec_int_hz=0  ..... interpolation of Hz receivers",
"                     - 0=Hz->Hz (no interpolation)",
"                     - 1=Hz->Hx",
"                     - 2=Hz->Ey",
"                     - 3=Hz->receiver position",
"   rec_int_hx=0 ...... interpolation of Hx receivers",
"                     - 0=Hx->Hx (no interpolation)",
"                     - 1=Hx->He",
"                     - 2=Hx->Ey",
"                     - 3=Hx->receiver position",
"" ,
"",
"      Jan Thorbecke 2014",
"      TU Delft",
"      E-mail: janth@xs4all.nl ",
"",
NULL};

/*
This code is largely based on fdelmodc where the acoustic4 algorithm has been adapted for GPR modeling 
To map the seismic parameter and names to EM data the following scheme is used:

Acoustic      => EM
P-velocity cp => epsilon_r
              => kappa_sigma attenuation
density rho   => mu_0 = 4.0*M_PI*10e-7 
kappa         => epsilon
Vx            => Hz
Vz            => -Hx
P             => Ey
2014 */


int main(int argc, char **argv)
{
	modPar mod;
	recPar rec;
	snaPar sna;
	wavPar wav;
	srcPar src;
	bndPar bnd;
	shotPar shot;
	float **src_nwav;
	float *mu, *ksigma, *eprs, *lam, *mul;
	float *tss, *tes, *tep, *p, *q, *r;
	float *vx, *vz, *tzz, *txz, *txx;
	float *rec_vx, *rec_vz, *rec_p;
	float *rec_txx, *rec_tzz, *rec_txz;
	float *rec_pp, *rec_ss;
	float *beam_vx, *beam_vz, *beam_p;
	float *beam_txx, *beam_tzz, *beam_txz;
	float *beam_pp, *beam_ss;	
	float sinkvel;
	double t0, t1, t2, t3, tt, tinit;
	size_t size, sizem, nsamp, memsize;
	int n1, ix, iz, ir, ishot, i;
	int ioPx, ioPz;
	int it0, it1, its, it, fileno, isam;
	int ixsrc, izsrc;
	int verbose;

	t0= wallclock_time();
	initargs(argc,argv);
	requestdoc(0);

	if (!getparint("verbose",&verbose)) verbose=0;
	getEmParameters(&mod, &rec, &sna, &wav, &src, &shot, &bnd, verbose);

	/* allocate arrays for model parameters: the different schemes use different arrays */

	n1 = mod.naz;
	sizem=mod.nax*mod.naz;

	ksigma = (float *)calloc(sizem,sizeof(float));
	eprs = (float *)calloc(sizem,sizeof(float));
	mu = (float *)calloc(sizem,sizeof(float));

	/* read velocity and density files */

	readEmModel(mod, bnd, eprs, ksigma, mu);

	/* read and/or define source wavelet(s) */

	/* Using a random source, which can have a random length 
	   for each source position, a pointer array with variable 
	   length (wav.nsamp[i]) is used.
	   The total length of all the source lengths together is wav.nst */
	
	if (wav.random) {
		src_nwav = (float **)calloc(wav.nx,sizeof(float *));
		src_nwav[0] = (float *)calloc(wav.nst,sizeof(float));
		assert(src_nwav[0] != NULL);
		nsamp = 0;
		for (i=0; i<wav.nx; i++) {
			src_nwav[i] = (float *)(src_nwav[0] + nsamp);
			nsamp += wav.nsamp[i];
		}
	}
	else {
		src_nwav = (float **)calloc(wav.nx,sizeof(float *));
		src_nwav[0] = (float *)calloc(wav.nt*wav.nx,sizeof(float));
		assert(src_nwav[0] != NULL);
		for (i=0; i<wav.nx; i++) {
			src_nwav[i] = (float *)(src_nwav[0] + wav.nt*i);
		}
	}

	defineSource(wav, src, src_nwav, mod.grid_dir, verbose);

	/* allocate arrays for wavefield and receiver arrays */

	vx  = (float *)calloc(sizem,sizeof(float));
	vz  = (float *)calloc(sizem,sizeof(float));
	tzz = (float *)calloc(sizem,sizeof(float)); /* =P field for acoustic */
	
	size = rec.n*rec.nt;
	if (rec.type.vz)  rec_vz  = (float *)calloc(size,sizeof(float));
	if (rec.type.vx)  rec_vx  = (float *)calloc(size,sizeof(float));
	if (rec.type.p)   rec_p   = (float *)calloc(size,sizeof(float));
	
	if(sna.beam) {
		size = sna.nz*sna.nx;
		if (sna.type.vz)  beam_vz  = (float *)calloc(size,sizeof(float));
		if (sna.type.vx)  beam_vx  = (float *)calloc(size,sizeof(float));
		if (sna.type.p)   beam_p   = (float *)calloc(size,sizeof(float));
	}

	t1= wallclock_time();
	if (verbose) {
		tinit = t1-t0;
		vmess("*******************************************");
		vmess("************* runtime info ****************");
		vmess("*******************************************");
		vmess("CPU time for intializing arrays and model = %f", tinit);
	}

	if (verbose>3) writeSrcRecPos(&mod, &rec, &src, &shot);

	/* Outer loop over number of shots */
	for (ishot=0; ishot<shot.n; ishot++) {

		izsrc = shot.z[ishot];
		ixsrc = shot.x[ishot];
		fileno= 0;

		memset(vx,0,sizem*sizeof(float));
		memset(vz,0,sizem*sizeof(float));
		memset(tzz,0,sizem*sizeof(float));
		if (mod.ischeme==2) {
			memset(q,0,sizem*sizeof(float));
		}
		if (mod.ischeme>2) {
			memset(txz,0,sizem*sizeof(float));
			memset(txx,0,sizem*sizeof(float));
		}
		if (mod.ischeme==4) {
			memset(r,0,sizem*sizeof(float));
			memset(p,0,sizem*sizeof(float));
			memset(q,0,sizem*sizeof(float));
		}
		if (verbose) {
			if (!src.random) {
				vmess("Modeling source %d at gridpoints ix=%d iz=%d", ishot, shot.x[ishot], shot.z[ishot]);
				vmess(" which are actual positions x=%.2f z=%.2f", mod.x0+mod.dx*shot.x[ishot], mod.z0+mod.dz*shot.z[ishot]);
			}
			vmess("Receivers at gridpoint x-range ix=%d - %d", rec.x[0], rec.x[rec.n-1]);
			vmess(" which are actual positions x=%.2f - %.2f", mod.x0+rec.xr[0], mod.x0+rec.xr[rec.n-1]);
			vmess("Receivers at gridpoint z-range iz=%d - %d", rec.z[0], rec.z[rec.n-1]);
			vmess(" which are actual positions z=%.2f - %.2f", mod.z0+rec.zr[0], mod.z0+rec.zr[rec.n-1]);
		}

		if (mod.grid_dir) { /* reverse time modeling */
			it0=-mod.nt+1;
			it1=0;
			its=-1;

			it0=0;
			it1=mod.nt;
			its=1;
		}
		else {
			it0=0;
			it1=mod.nt;
			its=1;
		}

		/* Main loop over the number of time steps */
		for (it=it0; it<it1; it++) {

#pragma omp parallel default (shared) \
shared (ksigma, eprs, lam, mu, txx, txz, tzz, vx, vz) \
shared (tss, tep, tes, r, q, p) \
shared (tinit, it0, it1, its) \
shared(beam_vx, beam_vz, beam_txx, beam_tzz, beam_txz, beam_p, beam_pp, beam_ss) \
shared(rec_vx, rec_vz, rec_txx, rec_tzz, rec_txz, rec_p, rec_pp, rec_ss) \
private (tt, t2, t3) \
shared (shot, bnd, mod, src, wav, rec, ixsrc, izsrc, it, src_nwav, verbose)
{
			em4(mod, src, wav, bnd, it, ixsrc, izsrc,src_nwav, 
				vx, vz, tzz, eprs, ksigma, mu, verbose);

			/* write samples to file if rec.nt samples are calculated */

#pragma omp master
{
			if ( (((it-rec.delay) % rec.skipdt)==0) && (it >= rec.delay) ) {
				int writeToFile, itwritten;

				writeToFile = ! ( (((it-rec.delay)/rec.skipdt)+1)%rec.nt );
				itwritten   = fileno*(rec.nt)*rec.skipdt;
				isam        = (it-rec.delay-itwritten)/rec.skipdt;

				/* store time at receiver positions */
				getRecTimes(mod, rec, bnd, it, isam, vx, vz, tzz, txx, txz, 
					rec_vx, rec_vz, rec_txx, rec_tzz, rec_txz, 
					rec_p, rec_pp, rec_ss, verbose);
			
				/* at the end of modeling a shot, write receiver array to output file(s) */
				if (writeToFile && (it+rec.skipdt <= it1-1) ) {
					fileno = ( ((it-rec.delay)/rec.skipdt)+1)/rec.nt;
					writeEmRec(rec, mod, ixsrc, izsrc, isam+1, ishot, fileno,
						rec_vx, rec_vz, rec_txx, rec_tzz, rec_txz, 
						rec_p, rec_pp, rec_ss, verbose);
				}
			}

			/* write snapshots to output file(s) */
			if (sna.nsnap) {
				writeEmSnapTimes(mod, sna, ixsrc, izsrc, it, vx, vz, tzz, txx, txz, verbose);
			}

			/* calculate beams */
			if(sna.beam) {
				getBeamTimes(mod, sna, vx, vz, tzz, txx,  txz, 
					beam_vx, beam_vz, beam_txx, beam_tzz, beam_txz, 
					beam_p, beam_pp, beam_ss, verbose);
			}
}
					
			/* taper the edges of the model */
//			taperEdges(mod, bnd, vx, vz, verbose);

#pragma omp master
{
			if (verbose) {
				if (it==it0+100*its) t2=wallclock_time();
				if (it==(it0+500*its)) {
					t3=wallclock_time();
					tt=(t3-t2)*(((it1-it0)*its)/400.0);
					vmess("Estimated compute time = %.2f s. per shot.",tt);
					vmess("Estimated total compute time = %.2f s.",tinit+shot.n*tt);
				}
			}
}
} /* end of OpenMP parallel section */

		} /* end of loop over time steps it */

		/* write output files: receivers and or beams */
		if (fileno) fileno++;
		
		if (rec.scale==1) { /* scale receiver with distance src-rcv */
			float xsrc, zsrc, Rrec, rdx, rdz;
			int irec;
			xsrc=mod.x0+mod.dx*ixsrc;
			zsrc=mod.z0+mod.dz*izsrc;
			for (irec=0; irec<rec.n; irec++) {
				rdx=mod.x0+rec.xr[irec]-xsrc;
				rdz=mod.z0+rec.zr[irec]-zsrc;
				Rrec = sqrt(rdx*rdx+rdz*rdz);
				fprintf(stderr,"Rec %d is scaled with distance %f R=%.2f,%.2f S=%.2f,%.2f\n", irec, Rrec,rdx,rdz,xsrc,zsrc);
				for (it=0; it<rec.nt; it++) {
					rec_p[irec*rec.nt+it] *= sqrt(Rrec);
				}
			}
		}
		writeEmRec(rec, mod, ixsrc, izsrc, isam+1, ishot, fileno,
			rec_vx, rec_vz, rec_txx, rec_tzz, rec_txz, 
			rec_p, rec_pp, rec_ss, verbose);
		
		writeBeams(mod, sna, ixsrc, izsrc, ishot, fileno, 
				   beam_vx, beam_vz, beam_txx, beam_tzz, beam_txz, 
				   beam_p, beam_pp, beam_ss, verbose);
		

	} /* end of loop over number of shots */


	t1= wallclock_time();
	if (verbose) {
		vmess("Total compute time FD modelling = %.2f s.", t1-t0);
	}

	/* free arrays */
	
	free(ksigma);
	free(eprs);
	free(src_nwav[0]);
	free(src_nwav);
	free(vx);
	free(vz);
	free(tzz);
	if (rec.type.vz)  free(rec_vz);
	if (rec.type.vx)  free(rec_vx);
	if (rec.type.p)   free(rec_p);
	if(sna.beam) {
		if (sna.type.vz)  free(beam_vz);
		if (sna.type.vx)  free(beam_vx);
		if (sna.type.p)   free(beam_p);
	}
	
	return 0;
}
