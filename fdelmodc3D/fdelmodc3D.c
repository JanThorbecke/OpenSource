#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include"par.h"
#include"fdelmodc3D.h"
#ifdef MPI
#include <mpi.h>
#endif

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

double wallclock_time(void);

void threadAffinity(void);

float ***alloc3float(modPar mod);
void free3float(float ***p);

long getParameters3D(modPar *mod, recPar *rec, snaPar *sna, wavPar *wav, srcPar *src,
	shotPar *shot, bndPar *bnd, long verbose);

long readModel3D(modPar mod, bndPar bnd, float ***rox, float ***roy, float ***roz,
	float ***l2m, float ***lam, float ***muu, float *tss, float *tes, float *tep);

long defineSource3D(wavPar wav, srcPar src, modPar mod, recPar rec, float **src_nwav,
	long reverse, long verbose);

long writeSrcRecPos3D(modPar *mod, recPar *rec, srcPar *src, shotPar *shot);

long acoustic6_3D(modPar mod, srcPar src, wavPar wav, bndPar bnd, long itime,
    long ixsrc, long iysrc, long izsrc, float **src_nwav, float *vx, float *vy, float *vz, 
	float *p, float ***rox, float ***roy, float ***roz, float ***l2m, long verbose);

long acoustic4_3D(modPar mod, srcPar src, wavPar wav, bndPar bnd, long itime,
	long ixsrc, long iysrc, long izsrc, float **src_nwav, float *vx, float *vy, float *vz,
	float *p, float ***rox, float ***roy, float ***roz, float ***l2m, long verbose);

long acoustic2_3D(modPar mod, srcPar src, wavPar wav, bndPar bnd, long itime,
    long ixsrc, long iysrc, long izsrc, float **src_nwav, float *vx, float *vy, float *vz,
    float *p, float ***rox, float ***roy, float ***roz, float ***l2m, long verbose);

long acoustic4_qr_3D(modPar mod, srcPar src, wavPar wav, bndPar bnd, long itime, 
    long ixsrc, long iysrc, long izsrc, float **src_nwav, float *vx, float *vy, float *vz,
    float *p, float ***rox, float ***roy, float ***roz, float ***l2m, long verbose);

long acousticSH4_3D(modPar mod, srcPar src, wavPar wav, bndPar bnd, long itime,
	long ixsrc, long iysrc, long izsrc, float **src_nwav, float *tx, float *ty, float *tz,
	float *vz, float ***rox, float ***roy, float ***roz, float ***mul, long verbose);

long elastic4dc_3D(modPar mod, srcPar src, wavPar wav, bndPar bnd, long itime, 
	long ixsrc, long iysrc, long izsrc, float **src_nwav, float *vx, float *vy, float *vz,
	float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz,
    float ***rox, float ***roy, float ***roz, float ***l2m, float ***lam, float ***mul, long verbose);

long getRecTimes3D(modPar mod, recPar rec, bndPar bnd, long itime, long isam,
	float *vx, float *vy, float *vz, float *tzz, float *tyy, float *txx,
	float *txz, float *txy, float *tyz, float ***l2m, float ***rox, float ***roy, float ***roz,
	float *rec_vx, float *rec_vy, float *rec_vz, float *rec_txx, float *rec_tyy, float *rec_tzz,
	float *rec_txz, float *rec_txy, float *rec_tyz, float *rec_p, float *rec_pp, float *rec_ss,
	float *rec_udp, float *rec_udvz, long verbose);

long writeRec3D(recPar rec, modPar mod, bndPar bnd, wavPar wav, long ixsrc, long iysrc, long izsrc,
	long nsam, long ishot, long fileno, float *rec_vx, float *rec_vy, float *rec_vz,
	float *rec_txx, float *rec_tyy, float *rec_tzz, float *rec_txz,  float *rec_tyz,
	float *rec_txy, float *rec_p, float *rec_pp, float *rec_ss,
	float *rec_udp, float *rec_udvz, long verbose);

long writeSnapTimes3D(modPar mod, snaPar sna, bndPar bnd, wavPar wav,
	long ixsrc, long iysrc, long izsrc, long itime, float *vx, float *vy, float *vz,
	float *tzz, float *tyy, float *txx, float *txz, float *tyz, float *txy, long verbose);

long getBeamTimes3D(modPar mod, snaPar sna, float *vx, float *vy, float *vz, float *tzz,
	float *tyy, float *txx, float *txz, float *tyz, float *txy, float *beam_vx, float *beam_vy,
	float *beam_vz, float *beam_txx, float *beam_tyy, float *beam_tzz, float *beam_txz,
	float *beam_tyz, float *beam_txy, float *beam_p, float *beam_pp, float *beam_ss, long verbose);

long writeBeams3D(modPar mod, snaPar sna, long ixsrc, long iysrc, long izsrc, long ishot,
	long fileno, float *beam_vx, float *beam_vy, float *beam_vz, float *beam_txx, float *beam_tyy,
	float *beam_tzz, float *beam_txz, float *beam_tyz, float *beam_txy, 
	float *beam_p, float *beam_pp, float *beam_ss, long verbose);

long allocStoreSourceOnSurface3D(srcPar src);

long freeStoreSourceOnSurface3D(void);

/* Self documentation */
char *sdoc[] = {
" ",
" fdelmodc - elastic acoustic finite difference wavefield modeling in 3D",
" ",
" IO PARAMETERS:",
"   file_cp= .......... P (cp) velocity file",
"   file_cs= .......... S (cs) velocity file",
"   file_den= ......... density (ro) file",
"   file_src= ......... file with source signature",
"   file_rcv=recv.su .. base name for receiver files",
"   file_snap=snap.su . base name for snapshot files",
"   file_beam=beam.su . base name for beam fields ",
"   dx= ............... read from model file: if dx==0 then dx= can be used to set it",
"   dy= ............... read from model file: if dy==0 then dy= can be used to set it",
"   dz= ............... read from model file: if dz==0 then dz= can be used to set it",
"   dt= ............... read from file_src: if dt is set it will interpolate file_src to dt sampling",
"" ,
" OPTIONAL PARAMETERS:",
"   ischeme=3 ......... 1=acoustic, 2=visco-acoustic 3=elastic, 4=visco-elastic, 5=double-couple",
"   tmod=(nt-1)*dt .... total modeling time (nt from file_src)",
"   ntaper=0 .......... length of taper in points at edges of model",
"   npml=35 ........... length of PML layer in points at edges of model",
"   R=1e-4 ............ the theoretical reflection coefficient at PML boundary",
"   m=2.0 ............. scaling order of the PML sigma function ",
"   tapfact=0.30 ...... taper strength: larger value gets stronger taper",
"   For the 4 boundaries the options are:  1=free 2=pml 3=rigid 4=taper",
"   top=1 ............. type of boundary on top edge of model",
"   left=4 ............ type of boundary on left edge of model",
"   right=4 ........... type of boundary on right edge of model",
"   bottom=4 .......... type of boundary on bottom edge of model",
"   front=4 ........... type of boundary on front edge of model",
"   back=4 ............ type of boundary on back edge of model",
"   grid_dir=0 ........ direction of time modeling (1=reverse time)",
"   Qp=15 ............. global Q-value for P-waves in visco-elastic (ischeme=2,4)",
"   file_qp= .......... model file Qp values as function of depth",
"   Qs=Qp ............. global Q-value for S-waves in visco-elastic (ischeme=4)",
"   file_qs= .......... model file Qs values as function of depth",
"   fw=0.5*fmax ....... central frequency for which the Q's are used",
"   sinkdepth=0 ....... receiver grid points below topography (defined bij cp=0.0)",
"   sinkdepth_src=0 ... source grid points below topography (defined bij cp=0.0)",
"   sinkvel=0 ......... use velocity of first receiver to sink through to next layer",
"   beam=0 ............ calculate energy beam of wavefield in model",
"   disable_check=0 ... disable stabilty and dispersion check and continue modeling",
"   verbose=0 ......... silent mode; =1: display info",
" ",
" SHOT AND GENERAL SOURCE DEFINITION:",
"   src_type=1 ........ 1=P 2=Txz 3=Tzz 4=Txx 5=S-pot 6=Fx 7=Fz 8=P-pot 9=double-couple",
"   src_orient=1 ...... orientation of the source",
"                     - 1=monopole",
"                     - 2=dipole +/- vertical oriented",
"                     - 3=dipole - + horizontal oriented",
"   dip=0 ............. dip for double-couple source",
"   strike=0 .......... strike for double-couple source",
"   rake=0 ............ rake for double-couple source",
"   xsrc=middle ....... x-position of (first) shot ",
"   ysrc=middle ....... y-position of (first) shot ",
"   zsrc=zmin ......... z-position of (first) shot ",
"   nshot=1 ........... number of shots to model",
"   dxshot=dx ......... if nshot > 1: x-shift in shot locations",
"   dyshot=0 .......... if nshot > 1: y-shift in shot locations",
"   dzshot=0 .......... if nshot > 1: z-shift in shot locations",
"   xsrca= ............ defines source array x-positions",
"   ysrca= ............ defines source array y-positions",
"   zsrca= ............ defines source array z-positions",
"   src_txt=........... text file with source coordinates. Col 1: x, Col. 2: z",
"   wav_random=1 ...... 1 generates (band limited by fmax) noise signatures ",
"   fmax=from_src ..... maximum frequency in wavelet",
"   src_multiwav=0 .... use traces in file_src as areal source",
"   src_at_rcv=1 ...... inject wavefield at receiver coordinates (1), inject at source (0)",
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
"   ysrc1=0 ........... left bound for y-position of sources",
"   ysrc2=0 ........... right bound for y-position of sources",
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
"   dysnap=dy ......... sampling in snapshot in y-direction",
"   ysnap1=0 .......... first y-position for snapshots area",
"   ysnap2=0 .......... last y-position for snapshot area",
"   dzsnap=dz ......... sampling in snapshot in z-direction",
"   zsnap1=0 .......... first z-position for snapshots area",
"   zsnap2=0 .......... last z-position for snapshot area",
"   snapwithbnd=0 ..... write snapshots with absorbing boundaries",
"   sna_type_p=1 ...... p registration _sp",
"   sna_type_vz=1 ..... Vz registration _svz",
"   sna_type_vy=0 ..... Vy registration _svy",
"   sna_type_vx=0 ..... Vx registration _svx",
"   sna_type_txx=0 .... Txx registration _stxx",
"   sna_type_tzz=0 .... Tzz registration _stzz",
"   sna_type_txz=0 .... Txz registration _stxz",
"   sna_type_pp=0 ..... P (divergence) registration _sP",
"   sna_type_ss=0 ..... S (curl) registration _sS",
"   sna_vxvztime=0 .... registration of vx/vx times",
"                       The fd scheme is also staggered in time.",
"                       Time at which vx/vz snapshots are written:",
"                     - 0=previous vx/vz relative to txx/tzz/txz at time t",
"                     - 1=next     vx/vz relative to txx/tzz/txz at time t",
"" ,
" RECEIVER SELECTION:",
"   xrcv1=xmin ........ first x-position of linear receiver array(s)",
"   xrcv2=xmax ........ last x-position of linear receiver array(s)",
"   dxrcv=dx .......... x-position increment of receivers in linear array(s)",
"   yrcv1=ymin ........ first y-position of linear receiver array(s)",
"   yrcv2=ymax ........ last y-position of linear receiver array(s)",
"   dyrcv=dy .......... y-position increment of receivers in linear array(s)",
"   zrcv1=zmin ........ first z-position of linear receiver array(s)",
"   zrcv2=zrcv1 ....... last z-position of linear receiver array(s)",
"   dzrcv=0.0 ......... z-position increment of receivers in linear array(s)",
"   dtrcv=.004 ........ desired sampling in receiver data (seconds)",
"   xrcva= ............ defines receiver array x-positions",
"   yrcva= ............ defines receiver array y-positions",
"   zrcva= ............ defines receiver array z-positions",
"   rrcv= ............. radius for receivers on a circle ",
"   arcv= ............. vertical arc-lenght for receivers on a ellipse (rrcv=horizontal)",
"   oxrcv=0.0 ......... x-center position of circle",
"   ozrcv=0.0 ......... z-center position of circle",
"   dphi=2 ............ angle between receivers on circle ",
"   rcv_txt=........... text file with receiver coordinates. Col 1: x, Col. 2: z",
"   rec_ntsam=nt ...... maximum number of time samples in file_rcv files",
"   rec_delay=0 ....... time in seconds to start recording: recorded time = tmod - rec_delay",
"   rec_type_p=1 ...... p registration _rp",
"   rec_type_vz=1 ..... Vz registration _rvz",
"   rec_type_vy=0 ..... Vy registration _rvy",
"   rec_type_vx=0 ..... Vx registration _rvx",
"   rec_type_txx=0 .... Txx registration _rtxx",
"   rec_type_tzz=0 .... Tzz registration _rtzz",
"   rec_type_txz=0 .... Txz registration _rtxz",
"   rec_type_pp=0 ..... P (divergence) registration _rP",
"   rec_type_ss=0 ..... S (curl) registration _rS",
"   rec_type_ud=0 ..... 1:pressure normalized decomposition in up and downgoing waves _ru, _rd",
"   ................... 2:particle velocity normalized decomposition in up and downgoing waves _ru, _rd",
"   kangle= ........... maximum wavenumber angle for decomposition",
"   rec_int_vx=0  ..... interpolation of Vx receivers",
"                     - 0=Vx->Vx (no interpolation)",
"                     - 1=Vx->Vz",
"                     - 2=Vx->Txx/Tzz(P)",
"                     - 3=Vx->receiver position",
"   rec_int_vy=0  ..... interpolation of Vy receivers",
"                     - 0=Vy->Vy (no interpolation)",
"                     - 1=Vy->Vz",
"                     - 2=Vy->Tyy/Tzz(P)",
"                     - 3=Vy->receiver position",
"   rec_int_vz=0 ...... interpolation of Vz receivers",
"                     - 0=Vz->Vz (no interpolation)",
"                     - 1=Vz->Vx",
"                     - 2=Vz->Txx/Tzz(P)",
"                     - 3=Vz->receiver position",
"   rec_int_p=0  ...... interpolation of P/Tzz receivers",
"                     - 0=P->P (no interpolation)",
"                     - 1=P->Vz",
"                     - 2=P->Vx",
"                     - 3=P->receiver position",
"" ,
" NOTES: For viscoelastic media dispersion and stability are not always",
" guaranteed by the calculated criteria, especially for Q values smaller than 13",
"",
"      Original fdelmodc by",
"      Jan Thorbecke 2011",
"      TU Delft",
"      E-mail: janth@xs4all.nl ",
"      3D extension by",
"      Joeri Brackenhoff 2019",
"      TU Delft",
"      E-mail: jabrackenhoff@tudelft.nl ",
"      2015  Contributions from Max Holicki",
"",
NULL};


int main(int argc, char **argv)
{
	modPar mod;
	recPar rec;
	snaPar sna;
	wavPar wav;
	srcPar src;
	bndPar bnd;
	shotPar shot;
	FILE 	*fp;
	float **src_nwav;
	float ***rox, ***roy, ***roz, ***l2m, ***lam, ***mul;
	float *tss, *tes, *tep, *p, *q, *r;
	float *vx, *vy, *vz, *tzz, *tyy, *txz, *txy, *tyz, *txx;
	float *rec_vx, *rec_vy, *rec_vz, *rec_p;
	float *rec_txx, *rec_tyy, *rec_tzz, *rec_txz, *rec_txy, *rec_tyz;
	float *rec_pp, *rec_ss;
	float *rec_udp, *rec_udvz;
	float *beam_vx, *beam_vy, *beam_vz, *beam_p;
	float *beam_txx, *beam_tyy, *beam_tzz, *beam_txz, *beam_txy, *beam_tyz;
	float *beam_pp, *beam_ss;	
	float sinkvel, npeshot;
	double t0, t1, t2, t3, tt, tinit;
	size_t size, sizem, nsamp;
	long n1, n2, ix, iy, iz, ir, ishot, i;
	long ioPx, ioPy, ioPz;
	long it0, it1, its, it, fileno, isam;
	long ixsrc, iysrc, izsrc, is0, is1;
	long verbose;
#ifdef MPI
	long     npes, pe;

	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &npes );
	MPI_Comm_rank( MPI_COMM_WORLD, &pe );
#else
	long     npes, pe;
	npes = 1;
	pe   = 0;
#endif


	t0= wallclock_time();
	initargs(argc,argv);
	requestdoc(0);

	if (!getparlong("verbose",&verbose)) verbose=0;
	getParameters3D(&mod, &rec, &sna, &wav, &src, &shot, &bnd, verbose);

	/* allocate arrays for model parameters: the different schemes use different arrays */ 

	n1 = mod.naz;
	n2 = mod.nax;
	sizem=mod.nax*mod.naz*mod.nay;

	rox = (float ***)alloc3float(mod);
	roy = (float ***)alloc3float(mod);
	roz = (float ***)alloc3float(mod);
	l2m = (float ***)alloc3float(mod);
	if (mod.ischeme==2) {
		tss = (float *)calloc(sizem,sizeof(float));
		tep = (float *)calloc(sizem,sizeof(float));
		q = (float *)calloc(sizem,sizeof(float));
	}
	if (mod.ischeme>2) {
		lam = (float ***)alloc3float(mod);
		mul = (float ***)alloc3float(mod);
	}
	if (mod.ischeme==4) {
		tss = (float *)calloc(sizem,sizeof(float));
		tes = (float *)calloc(sizem,sizeof(float));
		tep = (float *)calloc(sizem,sizeof(float));
		r = (float *)calloc(sizem,sizeof(float));
		p = (float *)calloc(sizem,sizeof(float));
		q = (float *)calloc(sizem,sizeof(float));
	}
	allocStoreSourceOnSurface3D(src);

	/* read velocity and density files */

	readModel3D(mod, bnd, rox, roy, roz, l2m, lam, mul, tss, tes, tep);

	// tss = (float *)calloc(sizem,sizeof(float));
	// for (iz=0; iz<mod.naz; iz++) {
	// 	for (iy=0; iy<mod.nay; iy++) {
	// 		for (ix=0; ix<mod.nax; ix++) {
	// 			tss[iy*n2*n1+ix*n1+iz] = roy[iy][ix][iz];
	// 		}
	// 	}
	// }
	// vmess("nax:%li nay:%li naz:%li",mod.nax,mod.nay,mod.naz);
	// fp = fopen("roy.bin","w+");
	// fwrite(tss,sizem,sizeof(float),fp);
	// fclose(fp);
	// for (iz=0; iz<mod.naz; iz++) {
	// 	for (iy=0; iy<mod.nay; iy++) {
	// 		for (ix=0; ix<mod.nax; ix++) {
	// 			tss[iy*n2*n1+ix*n1+iz] = rox[iy][ix][iz];
	// 		}
	// 	}
	// }
	// vmess("nax:%li nay:%li naz:%li",mod.nax,mod.nay,mod.naz);
	// fp = fopen("rox.bin","w+");
	// fwrite(tss,sizem,sizeof(float),fp);
	// fclose(fp);
	// for (iz=0; iz<mod.naz; iz++) {
	// 	for (iy=0; iy<mod.nay; iy++) {
	// 		for (ix=0; ix<mod.nax; ix++) {
	// 			tss[iy*n2*n1+ix*n1+iz] = roz[iy][ix][iz];
	// 		}
	// 	}
	// }
	// vmess("nax:%li nay:%li naz:%li",mod.nax,mod.nay,mod.naz);
	// fp = fopen("roz.bin","w+");
	// fwrite(tss,sizem,sizeof(float),fp);
	// fclose(fp);
	// free(tss);
	// return 0;

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

	defineSource3D(wav, src, mod, rec, src_nwav, mod.grid_dir, verbose);

	/* allocate arrays for wavefield and receiver arrays */

	vx  = (float *)calloc(sizem,sizeof(float));
	vy  = (float *)calloc(sizem,sizeof(float));
	vz  = (float *)calloc(sizem,sizeof(float));
	tzz = (float *)calloc(sizem,sizeof(float)); /* =P field for acoustic */
	if (mod.ischeme>2) {
		txz = (float *)calloc(sizem,sizeof(float));
		txy = (float *)calloc(sizem,sizeof(float));
		tyz = (float *)calloc(sizem,sizeof(float));
		txx = (float *)calloc(sizem,sizeof(float));
		tyy = (float *)calloc(sizem,sizeof(float));
	}
	
	size = rec.n*rec.nt;
	if (rec.type.vz)  rec_vz  = (float *)calloc(size,sizeof(float));
	if (rec.type.vy)  rec_vy  = (float *)calloc(size,sizeof(float));
	if (rec.type.vx)  rec_vx  = (float *)calloc(size,sizeof(float));
	if (rec.type.p)   rec_p   = (float *)calloc(size,sizeof(float));
	if (rec.type.txx) rec_txx = (float *)calloc(size,sizeof(float));
	if (rec.type.tyy) rec_tyy = (float *)calloc(size,sizeof(float));
	if (rec.type.tzz) rec_tzz = (float *)calloc(size,sizeof(float));
	if (rec.type.txz) rec_txz = (float *)calloc(size,sizeof(float));
	if (rec.type.txy) rec_txy = (float *)calloc(size,sizeof(float));
	if (rec.type.tyz) rec_tyz = (float *)calloc(size,sizeof(float));
	if (rec.type.pp)  rec_pp  = (float *)calloc(size,sizeof(float));
	if (rec.type.ss)  rec_ss  = (float *)calloc(size,sizeof(float));
    if (rec.type.ud) { 
		rec_udvz  = (float *)calloc(mod.nax*rec.nt,sizeof(float));
		rec_udp   = (float *)calloc(mod.nax*rec.nt,sizeof(float));
	}
	/* get velocity and density at first receiver location */
	//ir = mod.ioZz + rec.z[0]+(rec.x[0]+mod.ioZx)*n1+(rec.y[0]+mod.ioZy)*n1*n2;
	rec.rho = mod.dt/(mod.dx*roz[rec.y[0]+mod.ioZy][rec.x[0]+mod.ioZx][mod.ioZz+rec.z[0]]);
	rec.cp  = sqrt(l2m[rec.y[0]+mod.ioZy][rec.x[0]+mod.ioZx][mod.ioZz+rec.z[0]]*(roz[rec.y[0]+mod.ioZy][rec.x[0]+mod.ioZx][mod.ioZz+rec.z[0]]))*mod.dx/mod.dt;

	
	if(sna.beam) {
		size = sna.nz*sna.nx;
		if (sna.type.vz)  beam_vz  = (float *)calloc(size,sizeof(float));
		if (sna.type.vy)  beam_vy  = (float *)calloc(size,sizeof(float));
		if (sna.type.vx)  beam_vx  = (float *)calloc(size,sizeof(float));
		if (sna.type.p)   beam_p   = (float *)calloc(size,sizeof(float));
		if (sna.type.txx) beam_txx = (float *)calloc(size,sizeof(float));
		if (sna.type.tyy) beam_tyy = (float *)calloc(size,sizeof(float));
		if (sna.type.tzz) beam_tzz = (float *)calloc(size,sizeof(float));
		if (sna.type.txz) beam_txz = (float *)calloc(size,sizeof(float));
		if (sna.type.txy) beam_txy = (float *)calloc(size,sizeof(float));
		if (sna.type.tyz) beam_tyz = (float *)calloc(size,sizeof(float));
		if (sna.type.pp)  beam_pp  = (float *)calloc(size,sizeof(float));
		if (sna.type.ss)  beam_ss  = (float *)calloc(size,sizeof(float));
	}

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


    ioPx=mod.ioPx;
    ioPy=mod.ioPy;
    ioPz=mod.ioPz;
    if (bnd.lef==4 || bnd.lef==2) ioPx += bnd.ntap;
    if (bnd.fro==4 || bnd.fro==2) ioPy += bnd.ntap;
    if (bnd.top==4 || bnd.top==2) ioPz += bnd.ntap;
	if (rec.sinkvel) sinkvel=l2m[rec.y[0]+ioPy][rec.x[0]+ioPx][rec.z[0]+ioPz];
	else sinkvel = 0.0;

/* sink receivers to value different than sinkvel */
	for (ir=0; ir<rec.n; ir++) {
		iz = rec.z[ir];
		iy = rec.y[ir];
		ix = rec.x[ir];
		while(l2m[iy+ioPy][ix+ioPx][iz+ioPz] == sinkvel) iz++;
		rec.z[ir]=iz+rec.sinkdepth;
		rec.zr[ir]=rec.zr[ir]+(rec.z[ir]-iz)*mod.dz;
		if (verbose>3) vmess("receiver position %li at grid[ix=%li, iy=%li iz=%li] = (x=%f y=%f z=%f)", ir, ix+ioPx, iy+ioPy, rec.z[ir]+ioPz, rec.xr[ir]+mod.x0, rec.yr[ir]+mod.y0, rec.zr[ir]+mod.z0);
		//if (verbose>3) vmess("receiver position %li at grid[ix=%li, iy=%li iz=%li] = (x=%f y=%f z=%f)", ir, ix+ioPx, iy+ioPy, ioPz, rec.xr[ir]+mod.x0, rec.yr[ir]+mod.y0, rec.zr[ir]+mod.z0);
	}

/* sink sources to value different than zero */
	for (ishot=0; ishot<shot.n; ishot++) {
		iz = shot.z[ishot];
		iy = shot.y[ishot];
		ix = shot.x[ishot];
		while(l2m[iy+ioPy][ix+ioPx][iz+ioPz] == 0.0) iz++;
		shot.z[ishot]=iz+src.sinkdepth; 
	}

	/* scan for free surface boundary in case it has a topography */
	for (iy=0; iy<mod.ny; iy++) {
		for (ix=0; ix<mod.nx; ix++) {
			iz = ioPz;
			while(l2m[iy+ioPy][ix+ioPx][iz] == 0.0) iz++;
			bnd.surface[(iy+ioPy)*n2+ix+ioPx] = iz;
			if ((verbose>3) && (iz != ioPz)) vmess("Topography surface x=%.2f y=%.2f z=%.2f", mod.x0+mod.dx*ix, mod.y0+mod.dy*iy, mod.z0+mod.dz*(iz-ioPz));
		}
	}

	for (iy=0; iy<ioPy; iy++) {
		for (ix=0; ix<ioPx; ix++) {
			bnd.surface[iy*n2+ix] = bnd.surface[ioPy*n2+ioPx];
		}
		for (ix=ioPx+mod.nx; ix<mod.iePx; ix++) {
			bnd.surface[iy*n2+ix] = bnd.surface[ioPy*n2+mod.iePx-1];
		}
	}
	for (iy=ioPy+mod.ny; iy<mod.iePy; iy++) {
		for (ix=0; ix<ioPx; ix++) {
			bnd.surface[iy*n2+ix] = bnd.surface[(mod.iePy-1)*n2+ioPx];
		}
		for (ix=ioPx+mod.nx; ix<mod.iePx; ix++) {
			bnd.surface[iy*n2+ix] = bnd.surface[(mod.iePy-1)*n2+mod.iePx-1];
		}
	}
	if (verbose>3) writeSrcRecPos3D(&mod, &rec, &src, &shot);

	/* Outer loop over number of shots */
#ifdef MPI
    npeshot = MAX((((float)shot.n)/((float)npes)), 1.0);
    is0=ceil(pe*npeshot);
    is1=MIN(ceil((pe+1)*npeshot), shot.n);
    if (verbose>1) vmess("MPI: pe=%li does shots is0 %li - is1 %li\n", pe, is0, is1);
#else
	is0=0;
	is1=shot.n;
#endif


	for (ishot=is0; ishot<is1; ishot++) {

		izsrc = shot.z[ishot];
		iysrc = shot.y[ishot];
		ixsrc = shot.x[ishot];
		fileno= 0;

		memset(vx,0,sizem*sizeof(float));
		memset(vy,0,sizem*sizeof(float));
		memset(vz,0,sizem*sizeof(float));
		memset(tzz,0,sizem*sizeof(float));
		if (mod.ischeme==2) {
			memset(q,0,sizem*sizeof(float));
		}
		if (mod.ischeme>2) {
			memset(txz,0,sizem*sizeof(float));
			memset(txy,0,sizem*sizeof(float));
			memset(tyz,0,sizem*sizeof(float));
			memset(txx,0,sizem*sizeof(float));
			memset(tyy,0,sizem*sizeof(float));
		}
		if (mod.ischeme==4) {
			memset(r,0,sizem*sizeof(float));
			memset(p,0,sizem*sizeof(float));
			memset(q,0,sizem*sizeof(float));
		}
		if (verbose) {
			if (!src.random) {
				vmess("Modeling source %li at gridpoints ix=%li iy=%li iz=%li", ishot, shot.x[ishot], shot.y[ishot], shot.z[ishot]);
				vmess(" which are actual positions x=%.2f y=%.2f z=%.2f", mod.x0+mod.dx*shot.x[ishot], mod.y0+mod.dy*shot.y[ishot], mod.z0+mod.dz*shot.z[ishot]);
			}
			vmess("Receivers at gridpoint x-range ix=%li - %li", rec.x[0], rec.x[rec.n-1]);
			vmess(" which are actual positions x=%.2f - %.2f", mod.x0+rec.xr[0], mod.x0+rec.xr[rec.n-1]);
			vmess("Receivers at gridpoint y-range iy=%li - %li", rec.y[0], rec.y[rec.n-1]);
			vmess(" which are actual positions y=%.2f - %.2f", mod.y0+rec.yr[0], mod.y0+rec.yr[rec.n-1]);
			vmess("Receivers at gridpoint z-range iz=%li - %li", rec.z[0], rec.z[rec.n-1]);
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

					vmess("running time step time = %d.",it);

#pragma omp parallel default (shared) \
shared (rox, roy, roz, l2m, lam, mul, txx, tyy, tzz, txz, tyz, txy, vx, vy, vz) \
shared (tss, tep, tes, r, q, p) \
shared (tinit, it0, it1, its) \
shared(beam_vx, beam_vy, beam_vz, beam_txx, beam_tyy, beam_tzz, beam_txz, beam_tyz, beam_txy, beam_p, beam_pp, beam_ss) \
shared(rec_vx, rec_vy, rec_vz, rec_txx, rec_tyy, rec_tzz, rec_txz, rec_tyz, rec_txy, rec_p, rec_pp, rec_ss) \
shared (tt, t2, t3) \
shared (shot, bnd, mod, src, wav, rec, ixsrc, iysrc, izsrc, it, src_nwav, verbose)
{
			if (it==it0) {
				threadAffinity();
			}
			switch ( mod.ischeme ) {
				case -1 : /* Acoustic dissipative media FD kernel */
					acoustic4_qr_3D(mod, src, wav, bnd, it, ixsrc, iysrc, izsrc, src_nwav,
									vx, vy, vz, tzz, rox, roy, roz, l2m, verbose);
					break;
				case 1 : /* Acoustic FD kernel */
					if (mod.iorder==2) {
						acoustic2_3D(mod, src, wav, bnd, it, ixsrc, iysrc, izsrc, src_nwav,
								vx, vy, vz, tzz, rox, roy, roz, l2m, verbose);
					}
					else if (mod.iorder==4) {
                        if (mod.sh) {
                            acousticSH4_3D(mod, src, wav, bnd, it, ixsrc, iysrc, izsrc, src_nwav, 
                                    vx, vy, vz, tzz, rox, roy, roz, l2m, verbose);
                        }
                        else {
							acoustic4_3D(mod, src, wav, bnd, it, ixsrc, iysrc, izsrc, src_nwav,
									vx, vy, vz, tzz, rox, roy, roz, l2m, verbose);
                        }
					}
					else if (mod.iorder==6) {
							acoustic6_3D(mod, src, wav, bnd, it, ixsrc, iysrc, izsrc, src_nwav,
									vx, vy, vz, tzz, rox, roy, roz, l2m, verbose);
					}
					break;
				case 2 : /* Visco-Acoustic FD kernel */
						vmess("Visco-Acoustic order 4 not yet available");
					break;
				case 3 : /* Elastic FD kernel */
                    if (mod.iorder==4) {
						vmess("Elastic order 4 not yet available");
					}
					else if (mod.iorder==6) {
						vmess("Elastic order 6 not yet available");
                    }
					break;
				case 4 : /* Visco-Elastic FD kernel */
						vmess("Visco-Elastic order 4 not yet available");
					break;
				case 5 : /* Elastic FD kernel with S-velocity set to zero*/
						elastic4dc_3D(mod, src, wav, bnd, it, ixsrc, iysrc, izsrc, src_nwav,
							vx, vy, vz, tzz, tyy, txx, txz, txy, tyz, rox, roy, roz, l2m, lam, mul, verbose);
					break;
			}


			/* write samples to file if rec.nt samples are calculated */

#pragma omp master
{
			if ( (((it-rec.delay) % rec.skipdt)==0) && (it >= rec.delay) ) {
				long writeToFile, itwritten;

				writeToFile = ! ( (((it-rec.delay)/rec.skipdt)+1)%rec.nt );
				itwritten   = fileno*(rec.nt)*rec.skipdt;
                /* Note that time step it=0 (t=0 for t**-fields t=-1/2 dt for v*-field) is not recorded */
				isam        = (it-rec.delay-itwritten)/rec.skipdt+1;
				/* store time at receiver positions */
				getRecTimes3D(mod, rec, bnd, it, isam, vx, vy, vz,
					tzz, tyy, txx, txz, txy, tyz,
					l2m, rox, roy, roz,
					rec_vx, rec_vy, rec_vz,
					rec_txx, rec_tyy, rec_tzz, rec_txz, rec_txy, rec_tyz,
					rec_p, rec_pp, rec_ss, rec_udp, rec_udvz, verbose); 
			
				/* at the end of modeling a shot, write receiver array to output file(s) */
				if (writeToFile && (it+rec.skipdt <= it1-1) ) {
					fileno = ( ((it-rec.delay)/rec.skipdt)+1)/rec.nt;
					writeRec3D(rec, mod, bnd, wav, ixsrc, iysrc, izsrc, isam+1, ishot, fileno, 
             			rec_vx, rec_vy, rec_vz, rec_txx, rec_tyy, rec_tzz, rec_txz, rec_tyz, rec_txy,
             			rec_p, rec_pp, rec_ss, rec_udp, rec_udvz, verbose);
				}
			}

			/* write snapshots to output file(s) */
			if (sna.nsnap) {
				writeSnapTimes3D(mod, sna, bnd, wav, ixsrc, iysrc, izsrc, it,
								vx, vy, vz, tzz, tyy, txx, txz, tyz, txy, verbose);
			}

			/* calculate beams */
			if(sna.beam) {
				getBeamTimes3D(mod, sna, vx, vy, vz, tzz, tyy, txx, txz, tyz, txy,
				 	beam_vx, beam_vy, beam_vz, beam_txx, beam_tyy, beam_tzz, beam_txz, beam_tyz, beam_txy,
					beam_p, beam_pp, beam_ss, verbose);
			}
}
					
#pragma omp master
{
			if (verbose) {
				if (it==(it0+100*its)) t2=wallclock_time();
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
			float xsrc, ysrc, zsrc, Rrec, rdx, rdy, rdz;
			long irec;
			xsrc=mod.x0+mod.dx*ixsrc;
			ysrc=mod.y0+mod.dy*iysrc;
			zsrc=mod.z0+mod.dz*izsrc;
			for (irec=0; irec<rec.n; irec++) {
				rdx=mod.x0+rec.xr[irec]-xsrc;
				rdy=mod.y0+rec.yr[irec]-ysrc;
				rdz=mod.z0+rec.zr[irec]-zsrc;
				Rrec = sqrt(rdx*rdx+rdy*rdy+rdz*rdz);
				fprintf(stderr,"Rec %li is scaled with distance %f R=%.2f,%.2f,%.2f S=%.2f,%.2f,%.2f\n", irec, Rrec,rdx,rdy,rdz,xsrc,ysrc,zsrc);
				for (it=0; it<rec.nt; it++) {
					rec_p[irec*rec.nt+it] *= sqrt(Rrec);
				}
			}
		}
		writeRec3D(rec, mod, bnd, wav, ixsrc, iysrc, izsrc, isam+1, ishot, fileno, 
            rec_vx, rec_vy, rec_vz, rec_txx, rec_tyy, rec_tzz, rec_txz, rec_tyz, rec_txy,
            rec_p, rec_pp, rec_ss, rec_udp, rec_udvz, verbose);
		
		writeBeams3D(mod, sna, ixsrc, iysrc, izsrc, ishot, fileno, 
			   		beam_vx, beam_vy, beam_vz, beam_txx, beam_tyy, beam_tzz, beam_txz, beam_tyz, beam_txy, 
			   		beam_p, beam_pp, beam_ss, verbose);
		

	} /* end of loop over number of shots */


	t1= wallclock_time();
	if (verbose) {
		vmess("Total compute time FD modelling = %.2f s.", t1-t0);
	}

	/* free arrays */

	initargs(argc,argv); /* this will free the arg arrays declared */
	free3float(rox);
	free3float(roy);
	free3float(roz);
	free3float(l2m);
	free(src_nwav[0]);
	free(src_nwav);
	free(vx);
	free(vy);
	free(vz);
	free(tzz);
	freeStoreSourceOnSurface3D();
	if (rec.type.vz)  free(rec_vz);
	if (rec.type.vy)  free(rec_vy);
	if (rec.type.vx)  free(rec_vx);
	if (rec.type.p)   free(rec_p);
	if (rec.type.txx) free(rec_txx);
	if (rec.type.tyy) free(rec_tyy);
	if (rec.type.tzz) free(rec_tzz);
	if (rec.type.txz) free(rec_txz);
	if (rec.type.tyz) free(rec_tyz);
	if (rec.type.txy) free(rec_txy);
	if (rec.type.pp)  free(rec_pp);
	if (rec.type.ss)  free(rec_ss);
	if (rec.type.ud)  {
		free(rec_udvz);
		free(rec_udp);
	}
	if(sna.beam) {
		if (sna.type.vz)  free(beam_vz);
		if (sna.type.vy)  free(beam_vy);
		if (sna.type.vx)  free(beam_vx);
		if (sna.type.p)   free(beam_p);
		if (sna.type.txx) free(beam_txx);
		if (sna.type.tyy) free(beam_tyy);
		if (sna.type.tzz) free(beam_tzz);
		if (sna.type.txz) free(beam_txz);
		if (sna.type.tyz) free(beam_tyz);
		if (sna.type.txy) free(beam_txy);
		if (sna.type.pp)  free(beam_pp);
		if (sna.type.ss)  free(beam_ss);
	}
	
	if (mod.ischeme==2) {
		free(tss);
		free(tep);
		free(q);
	}
	if (mod.ischeme>2) {
		free3float(lam);
		free3float(mul);
		free(txz);
		free(tyz);
		free(txy);
		free(txx);
		free(tyy);
	}
	if (mod.ischeme==4) {
		free(tss);
		free(tes);
		free(tep);
		free(r);
		free(p);
		free(q);
	}

#ifdef MPI  
    MPI_Finalize();
#endif

	return 0;
}
