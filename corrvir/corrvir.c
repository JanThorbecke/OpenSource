#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "par.h"
#include "corr.h"
#include "segy.h" 

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

double wallclock_time(void);

int omp_get_max_threads(void);
int omp_get_num_threads(void);
int omp_get_thread_num(void);
void omp_set_num_threads(int num_threads);

typedef struct { /* complex number */
        float r,i;
} complex;

void read_sutrace_at_position(FILE *fp, int itrace, int ircvnum, int isrcnum, complex *tracedata, complex *tracebuffer, int *trace_in_buffer, size_t *nbuffer, int ntfft, int nt, int nfreq, size_t nbufmax, size_t trace_sz, double *tread, double *tfft, double *tmemcpy);

void get_sutrace_at_position(FILE *fp, size_t itrace, complex *tracedata, complex *tracebuffer, int ntfft, int nt, int nfreq, size_t nbufmax, size_t trace_sz, double *tread, double *tfft, double *tmemcpy);

int optncr(int n);
void cr1_fft(complex *cdata, float *data, int n, int sign);
void rc1_fft(float *rdata, complex *cdata, int n, int sign);

int getFileInfo(char *file_name, int *n1, int *n2, float *d1, float *d2, int verbose);

int computeFFT(float *datain, int ld1, int nStart, int ntfft, int nfreq, int nstation, complex *cA, complex *cB, double *tfft, int verbose);

int correlate(complex *cmaster, complex *cslaves, int nfreq, int ncor, double *tcorr, int verbose);
int coherence(complex *cmaster, complex *cslaves, int nfreq, int ncor, float reps, float epsmax, double *tcorr, int verbose);

/**************
* ntc output samples of correlation result
* note that nt (the number of samples read by the IO routine)
* should be 2*ntc and a number efficient for FFT's
*/

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" corrvir - correlation to compute virtual shot records",
"  ",
" corrvir file_shots= file_out= [optional parameters]",
"  ",
" Required parameters: ",
" ",
"   file_shots= ....... filename of shots used to calculate virtual shots",
"   file_out= ......... filename of the output file with the virtual shot records",
"   nsources= ......... maximum number of actual shots in a common-receiver gather",
"   nreceivers= ....... maximum number of receiver positions",
"   xsrcv= ............ defines virtual-source x-positions",
"   ysrcv= ............ defines virtual-source y-positions",
" ",
" Optional parameters: ",
" ",
"   dtolerr = 0.1 ....... tolerated error distance (in m) between above requested virtual-source position and actual receiver positions",
"   fmax = 70 ........... maximum frequency to use in correlation",
"   nbm = 4 ............. amount of memory (in GB) to buffer traces for reusage",
"   nsrc_ctr = 1 ...... . minumum number of sources to contribute to Fresnel stack",
"   normsrc = 0 ....... . normalize each correlated trace (before source summation)",
"   normalize = 1 ....... normalize the virtual trace with the number of contributing sources",
"   cohr = 0 ............ use cross-coherence:  {f1(w)*f2(w)/(|f1(w)|*|f2(w)|+eps)} ",
"   reps = 0.0 .......... relative stabilization factor for cohr; eps=reps*|f1(w)|*|f2(w)|",
"   epsmax = 0.1 ........ cut off range above which spectrum is flattened",
"   src_sel = 0 ......... use all sources that are contributing ",
"           = 1 ......... use only sources that are outside square area enclosing master and slave positions",
"           = 2 ......... use only sources that are in donut",	
"   bmin = 1 ............ define inner circle of donut",
"   bmax = 4 ............ define outer circle of donut",
" ",
" OUTPUT DEFINITION ",
"   causal=1 .......... output causal(1), non-causal(2), both(3), or summed(4)",
"   verbose=0 ......... silent option; >0 displays info",
" ",
" Notes: ",
"    ntc is the number of output samples of correlation result",
"    nt  is the number of samples per trace read by the IO routine, should be 2*ntc and a number efficient for FFT's",
"    For causal=3: t=0 is sample nt/2 and trace goes from [-nt/2*dt - 0 - nt/2*dt]",
" ",
" Authors  : Jan Thorbecke : 2013 (j.w.thorbecke@tudelft.nl)",
"            Boris Boullenger : 2015 (b.boullenger@tudelft.nl)",
" ",
NULL};
/**************** end self doc ***********************************/

int main (int argc, char **argv)
{
	int	i, j, k, ret, nshots, maxnfiles, nfiles;
	int	size, n1, n2, ntfft, nf, ifreq;
	int     verbose, causal;
	int     nt, nstation, ntc, nb, istation, jstation;
	int     pgsz, istep, nsrc_ctr, nsrc, nsrcMaster, nsrcSlave;
	size_t  outputSize, nwrite, cdatainSize, datainSize, cdataoutSize;
	size_t  jstep, lsz, offset, trace_sz, nread;
	int     ierr, itrace_out, done, cohr, normalize, normsrc;
	int	storedrcv, storedsrc, src_sel, read_trace;
	int 	nsources, nreceivers, jsource, ivirtual, ivrcv, ivsrc, isrcs, isrcm, ivs;
	int	sxv, syv, sx, sy, svelev, sxm, sym, sxs, sys;
	int	cx3, cy3, gxv, gyv, gvelev, gx, gy, gelev;
	double	cx1, cy1, cx2, cy2;
	int	distmc2, distsc2, distsc2min, distsc2max;
	int 	bmin, bmax;  
	int 	rpeg, speg, vrpeg, vspeg;
	int	xmin, xmax, ymin, ymax;
	int     *trace_in_buffer;
	size_t  nbuffer, nbufmax, nbufmemory;
	size_t ntrace, bufsize, nsamp, ipos, nfreq;
	int nbm;
	crgPos *rcvSrc;
	int ircv, *isrc, *nsrc_crg, *vsrc;
	float	dx, fmax, df, d1, d2, scl, fscalco, sclfft;
	float   dt, reps, epsmax;
	float 	dtolerr;
	float	*r, *vtrace;
	float 	diffx, diffy, dd, ddmin;
	float   *xsrcv, *ysrcv;
	float   rms2, sclrms;
	complex tmpcc; 
	int     nxsrcv, nysrcv, nvsrc;
	int	 	scalco;
	complex *c, *cmaster, *cslaves, *tracebuffer;
	double  t0, t1, t2, t3, ttotal, tinit, tread, tfft, tcorr, twrite, tlogic, trdloc, tftloc, tcorrloc, tmemcpy, tselect, tmemset, tloop, troutine0, troutine;
	double  tread_ivs, tfft_ivs, tcorr_ivs, twrite_ivs, tlogic_ivs, t1_ivs, t3_ivs, t1_ivs_ivr, t2_ivs_ivr, t3_ivs_ivr, t4_ivs_ivr, t5_ivs_ivr;
	char	*file_shots, *file_out, filename[1024], basename[1200], base[1200], *pbase;
	int     pe=0, root_pe=0, npes=1, ipe, size_s;
	float *trace;
	FILE *fpin, *fpout;
	segy hdr;
	segy *segy_hdrs;
	segy outputhdr;

	t0 = wallclock_time();
	tinit = twrite = tread = tcorr = tfft = tlogic = 0.0;
	initargs(argc, argv);
	requestdoc(1);
	
	if (!getparint("verbose", &verbose)) verbose = 0;
	if (!getparstring("file_shots", &file_shots)) file_shots=NULL;
	assert(file_shots != NULL);
	if (!getparstring("file_out", &file_out)) file_out=NULL;
	if (!getparint("nsources", &nsources)) nsources = 0;
	assert(nsources!=0);
	if (!getparint("nreceivers", &nreceivers)) nreceivers = 0;
	assert(nreceivers!=0);
	if (!getparint("nsrc_ctr", &nsrc_ctr)) nsrc_ctr = 1;
	if (!getparint("normsrc", &normsrc)) normsrc = 0;
	if (!getparint("causal", &causal)) causal = 1;
	if (!getparint("normalize", &normalize)) normalize = 1;
	if (!getparint("cohr", &cohr)) cohr = 0;
	if (!getparint("src_sel", &src_sel)) src_sel = 0;
	if (!getparint("bmin", &bmin)) bmin = 1;
	if (!getparint("nbm", &nbm)) nbm = 0;
	if (!getparint("bmax", &bmax)) bmax = 4;
	if (!getparfloat("reps", &reps)) reps = 0.0;
	if (!getparfloat("epsmax", &epsmax)) epsmax = 0.1;
	if (!getparfloat("fmax", &fmax)) fmax = 70;
	if (!getparfloat("dtolerr", &dtolerr)) dtolerr = 0.1;
	nxsrcv = countparval("xsrcv");
	nysrcv = countparval("ysrcv");
	if (nxsrcv != nysrcv) {
		verr("Number of sources in array xsrcv (%d), ysrcv(%d) are not equal",nxsrcv, nysrcv);
	}
	xsrcv = (float *)malloc(nxsrcv*sizeof(float));
	ysrcv = (float *)malloc(nxsrcv*sizeof(float));

	getparfloat("xsrcv", xsrcv);
	getparfloat("ysrcv", ysrcv);

	getFileInfo(file_shots, &n1, &n2, &d1, &d2, verbose);
	ntrace = n2;
	nt = n1;
	dt = d1;
	fpin = fopen(file_shots, "r");
	assert( fpin != NULL);
	nread = fread( &hdr, 1, TRCBYTES, fpin );
	assert(nread == TRCBYTES);
	fpout = fopen(file_out, "w+");

	fprintf(stderr,"nt=%d dt=%f ntrace=%ld\n",nt, dt, (long) ntrace);

    ntfft = optncr(2*nt);
	nf    = ntfft/2+1;
	df    = 1.0/(ntfft*dt);
	nfreq = MIN(nf,(int)(fmax/df)+1);
	sclfft= 1.0/((float)ntfft);
	trace_sz = sizeof(float)*(n1)+TRCBYTES;
	segy_hdrs = (segy *)calloc(ntrace,sizeof(segy));
	assert(segy_hdrs != NULL);
	ret = fseeko(fpin , 0, SEEK_SET);

	t1 = wallclock_time();
	tinit += t1-t0;

	/* check if data is fully buffered */
	if (nbm==0) { /* read only headers */
		nbufmax = 0;
		for (i=0; i<ntrace; i++) {
			if(i % 100000 == 0) fprintf(stderr,"i=%d out of %d\n", i, ntrace);
			offset = i*trace_sz;
			ret = fseeko(fpin , offset, SEEK_SET);
			if (ret<0) perror("fseeko");
   			nread = fread( &segy_hdrs[i], 1, TRCBYTES, fpin );
		}
	} 
	else { /* buffer is allocated to store all data traces */
        bufsize = (size_t)(ntrace*nfreq*sizeof(complex));
		nbufmax = ntrace;
		fprintf(stderr,"memory allocated to buffer all traces = %.2f GB\n", (float)(bufsize/(1024*1024*1024)));
		tracebuffer = (complex *)malloc(bufsize);
		assert (tracebuffer != NULL);

		/*================ Read geometry of all traces in file_shots ================*/

		r = (float *)calloc(ntfft,sizeof(float));
		c = (complex *)calloc((nf),sizeof(complex));
		for (i=0; i<ntrace; i++) {
			if(i % 100000 == 0) fprintf(stderr,"i=%d out of %d\n", i, ntrace);
			offset = i*trace_sz;
			ret = fseeko(fpin , offset, SEEK_SET);
			if (ret<0) perror("fseeko");
   			nread = fread( &segy_hdrs[i], 1, TRCBYTES, fpin );
			assert(nread == TRCBYTES);
			
			nsamp =  segy_hdrs[i].ns;
			memset(r,0,ntfft*sizeof(float));
       		nread = fread( r, sizeof(float), nt, fpin );
       		assert(nread == nt);
        	rc1_fft(r,c,ntfft,-1);
        	memcpy(&tracebuffer[i*nfreq],&c[0],nfreq*sizeof(complex));
		}
		free(r);
		free(c);
	}
	t2 = wallclock_time();
	tread += t2-t1;

	
	/*================ Sort traces in common receiver gathers ================*/
	
	/* nreceivers is total number of receiver positions which is (much) smaller than the number of traces in the file */
	/* Every receiver position has a number of sources: common receiver gather */

	rcvSrc = (crgPos *)calloc(nreceivers,sizeof(crgPos));
	for (i=0; i<nreceivers; i++) {
		rcvSrc[i].src = (traceCoord *)calloc(nsources,sizeof(traceCoord));
	}
		
	rcvSrc[0].gx = segy_hdrs[0].gx;
	rcvSrc[0].gy = segy_hdrs[0].gy;
	rcvSrc[0].rpeg = segy_hdrs[0].gdel;
	rcvSrc[0].gelev = segy_hdrs[0].gelev;
	rcvSrc[0].src[0].x = segy_hdrs[0].sx;
	rcvSrc[0].src[0].y = segy_hdrs[0].sy;
	rcvSrc[0].src[0].peg = segy_hdrs[0].sdel;
	scalco = segy_hdrs[0].scalco;
	rcvSrc[0].src[0].fpos = 0;
	rcvSrc[0].nsrc = 1;
	ircv=1;
	if (scalco < 0) fscalco = 1.0/(-1.0*scalco);
	else fscalco = (float)scalco;

	for (jstation=0; jstation<ntrace; jstation++) {
		gx = segy_hdrs[jstation].gx;
		gy = segy_hdrs[jstation].gy;
		rpeg = segy_hdrs[jstation].gdel;
		gelev = segy_hdrs[jstation].gelev;
		sx = segy_hdrs[jstation].sx;
		sy = segy_hdrs[jstation].sy;
		speg = segy_hdrs[jstation].sdel;
		
		/* bookkeeping: find out how many sources are contributing to each receiver position */
		storedrcv = 0;
		for (i=0; i<ircv; i++) {
			if ( (rcvSrc[i].gx == gx) && (rcvSrc[i].gy == gy) ) {
				storedrcv = 1;
				storedsrc = 0;
				for (j=0; j<rcvSrc[i].nsrc; j++) {
					if ( (rcvSrc[i].src[j].x == sx) && (rcvSrc[i].src[j].y == sy) ) {
						storedsrc=1;
						break;
					}
				}
				if ( !storedsrc ) {
					rcvSrc[i].src[rcvSrc[i].nsrc].x = sx;
					rcvSrc[i].src[rcvSrc[i].nsrc].y = sy;
					rcvSrc[i].src[rcvSrc[i].nsrc].peg = speg;
					rcvSrc[i].src[rcvSrc[i].nsrc].fpos = jstation;
					rcvSrc[i].nsrc += 1;
				}
				break;
			}
		}
		if ( !storedrcv ) {
			if (verbose>=2) fprintf(stderr,"for jstation = %d  number of receiver positions %d\n",jstation, ircv+1);
			rcvSrc[ircv].gx = gx;
			rcvSrc[ircv].gy = gy;
			rcvSrc[ircv].rpeg = rpeg;
			rcvSrc[ircv].gelev = gelev;
			rcvSrc[ircv].src[0].x = sx;
			rcvSrc[ircv].src[0].y = sy;
			rcvSrc[ircv].src[0].peg = speg;
			rcvSrc[ircv].src[0].fpos = jstation;
			rcvSrc[ircv].nsrc = 1;
			ircv++;
		}
	}
	fprintf(stderr,"number of receiver positions in file = %d \n",ircv);
	nreceivers=ircv;
	nsources = 0;
	nsrc_crg = (int *)calloc(nreceivers,sizeof(int));
	for (i=0; i<nreceivers; i++) {
		nsrc_crg[i] = rcvSrc[i].nsrc;
		nsources = MAX(nsources,rcvSrc[i].nsrc);
        	if(verbose>=2) fprintf(stderr,"number of sources in receiver position %d = %d \n",i+1, rcvSrc[i].nsrc);
	}

	
	
	/*================ Find receiver positions corresponding to requested virtual-sources ================*/
	if (nxsrcv) {
		nvsrc = nxsrcv;
		vsrc = (int *)malloc(nvsrc*sizeof(int));
		for (j=0; j<nxsrcv; j++) {
			diffx = xsrcv[j]- rcvSrc[0].gx*fscalco;
			diffy = ysrcv[j]- rcvSrc[0].gy*fscalco;
			ddmin = diffx*diffx+diffy*diffy;
			vsrc[j] = 0;
			for (i=1; i<nreceivers; i++) {
				diffx = xsrcv[j]- rcvSrc[i].gx*fscalco;
				diffy = ysrcv[j]- rcvSrc[i].gy*fscalco;
				dd = diffx*diffx+diffy*diffy;
				if (dd < ddmin ) {
					ddmin = dd; 
					vsrc[j] = i;
				}
			}
			fprintf(stderr,"ddmin for requested vsrc position (%.2f,%.2f) = %.2f\n",xsrcv[j],ysrcv[j],ddmin);
			if (ddmin <= dtolerr) {
				xsrcv[j]=rcvSrc[vsrc[j]].gx*fscalco;
				ysrcv[j]=rcvSrc[vsrc[j]].gy*fscalco;
				fprintf(stderr,"Requested vsrc position coincide with rcv position (%.2f,%.2f)\n",xsrcv[j],ysrcv[j]);		 	
			}
			else {
			fprintf(stderr,"Error: Requested vsrc position do not coincide with a rcv position\n");
			}
		}
	}
	else { /* all receivers are made virtual sources */
		nvsrc = nreceivers;
		vsrc = (int *)malloc(nvsrc*sizeof(int));
		for (i=0; i<nreceivers; i++) {
			vsrc[i] = i;
		}
	}

	/*================ initializations ================*/

	
	/*================ for nvirtual shots find suitable receivers ================*/
	t1 = wallclock_time();
	tinit += t1-t2;	

	#pragma omp parallel default(shared) \
	private(cmaster,cslaves,vtrace,r,c, t1_ivs, t3_ivs) \
	private(ivrcv, gxv, gyv, vrpeg, gvelev, nsrcSlave) \
	private(xmin,ymin,xmax,ymax, cx1, cy1, cx2, cy2, cx3, cy3, distmc2, distsc2min, distsc2max, distsc2) \
	private(nsrc, trdloc, tftloc, tcorrloc) \
	private(isrcm, isrcs, sxm, sym, sxs, sys, read_trace) \
	private(jsource, rms2, ifreq, tmpcc, sclrms, scl) \
	private(j, outputhdr, itrace_out, nwrite)
{

	cmaster = (complex *)malloc(nsources*nfreq*sizeof(complex));
	cslaves = (complex *)malloc(nsources*nfreq*sizeof(complex));
	vtrace = (float *)malloc(nt*sizeof(float));
	r = (float *)malloc(ntfft*sizeof(float));
	c = (complex *)malloc((ntfft/2+1)*sizeof(complex));
	//itrace_out=0;

        fprintf(stderr,"Number of OpenMP threads set  = %d number=%d\n", omp_get_max_threads(), omp_get_thread_num());

	ivs = 0;
	while (ivs<nvsrc) {/* loop over the number of virtual source positions to be created */
		
		t1_ivs = wallclock_time();
		//tread_ivs = tfft_ivs = tcorr_ivs = twrite_ivs = tlogic_ivs = tmemcpy = tselect = tmemset = tloop = troutine = 0.0;

		ivsrc=vsrc[ivs];
		sxv = rcvSrc[ivsrc].gx;
		syv = rcvSrc[ivsrc].gy;
		vspeg = rcvSrc[ivsrc].rpeg;
		svelev = rcvSrc[ivsrc].gelev;
		nsrcMaster = nsrc_crg[ivsrc];

		#pragma omp for
		for (ivrcv=0; ivrcv<nreceivers; ivrcv++) {
			
			//t1_ivs_ivr = wallclock_time();
			
			gxv = rcvSrc[ivrcv].gx;
			gyv = rcvSrc[ivrcv].gy;
			vrpeg = rcvSrc[ivrcv].rpeg;
			gvelev = rcvSrc[ivrcv].gelev;
			nsrcSlave = nsrc_crg[ivrcv];
			
			if (src_sel == 1) {
				xmin = MIN(gxv,sxv);
				ymin = MIN(gyv,syv);
				xmax = MAX(gxv,sxv);
				ymax = MAX(gyv,syv);
			}

			if (src_sel == 2) {
				cx1 = 0.5*(sxv+gxv);
				cy1 = 0.5*(syv+gyv);
				cx2 = ceil(cx1);
				cy2 = ceil(cy1);
				cx3 = cx2;
				cy3 = cy2;
				distmc2 = (sxv-cx3)*(sxv-cx3)+(syv-cy3)*(syv-cy3);
				distsc2min = distmc2+bmin*bmin*distmc2;
				distsc2max = distmc2+bmax*bmax*distmc2;
				distsc2 = (sxv-cx3)*(sxv-cx3)+(syv-cy3)*(syv-cy3);
			}

			//fprintf(stderr,"ivsrc=%d sxv=%d nsrcMaster=%d ivrcv=%d gxv=%d nsrcSlave=%d\n", ivsrc, sxv, nsrcMaster, ivrcv, gxv, nsrcSlave);
			
			nsrc = 0;
			memset(cslaves,0,nfreq*nsources*sizeof(complex));
			memset(cmaster,0,nfreq*nsources*sizeof(complex));
			//t0 = wallclock_time();
			//tmemset += t0-t1_ivs_ivr;

			//trdloc = tftloc = tcorrloc = 0.0;
			for (isrcm=0; isrcm<nsrcMaster; isrcm++) { /* for all traces find which traces are recorded by both virtual source and receiver position */
				sxm = rcvSrc[ivsrc].src[isrcm].x;
				sym = rcvSrc[ivsrc].src[isrcm].y;
				for (isrcs=0; isrcs<nsrcSlave; isrcs++) {
					sxs = rcvSrc[ivrcv].src[isrcs].x;
					sys = rcvSrc[ivrcv].src[isrcs].y;
					if ( (sxm == sxs) && (sym == sys) ) { /* master and slave have same contributing source */

						/* if the source is outside the inclosing area of the virtual-receiver and virtual-source coordinates 
				 		   then the source is accepted for contribution to the final summation */
						if (src_sel == 0) {
							read_trace=1;
						}
						else if (src_sel == 1) {
							read_trace = 0;
							if ( !( (sxs>xmin) && (sxs<xmax) ) ) {
								if ( !( (sys>ymin) && (sys<ymax) ) ) {
									read_trace=1;
								}
							}
						}
						else if (src_sel == 2) {
							read_trace = 0;
							if ( !( (distsc2 < distsc2min) && (distsc2 > distsc2max) ) ) {
								read_trace = 1;
							}
						}
						else if (src_sel == 3) {
							read_trace = 0;
							if ( (sxv < gxv) && (syv < gyv) ) {
								if ( (sxm < sxv) && (sym < syv) ) {
									read_trace = 1;
								}
							}
							else if ( (sxv < gxv) && (syv > gyv) ) {
								if ( (sxm < sxv) && (sym > syv) ) {
									read_trace = 1;
								}
							}
							else if ( (sxv > gxv) && (syv < gyv) ) {
								if ( (sxm > sxv) && (sym < syv) ) {
									read_trace = 1;
								}
							}
							else if ( (sxv > gxv) && (syv > gyv) ) {
								if ( (sxm > sxv) && (sym > syv) ) {
									read_trace = 1;
								}
							}
							
						}
						else if (src_sel == 4) {
							read_trace = 0;
							if ( (sxv < gxv) && (syv < gyv) ) {
								if ( (sxm < gxv) && (sym < gyv) ) {
									read_trace = 1;
								}
							}
							else if ( (sxv < gxv) && (syv > gyv) ) {
								if ( (sxm < gxv) && (sym > gyv) ) {
									read_trace = 1;
								}
							}
							else if ( (sxv > gxv) && (syv < gyv) ) {
								if ( (sxm > gxv) && (sym < gyv) ) {
									read_trace = 1;
								}
							}
							else if ( (sxv > gxv) && (syv > gyv) ) {
								if ( (sxm > gxv) && (sym > gyv) ) {
									read_trace = 1;
								}
							}
						}
						if (read_trace) { /* read data from file for sources that passed the tests above */
//#pragma omp critical
//{
							//troutine0 = wallclock_time();
							/* read master trace */
							//memcpy(&cmaster[nsrc*nfreq], &tracebuffer[rcvSrc[ivsrc].src[isrcm].fpos*nfreq],nfreq*sizeof(complex));
							get_sutrace_at_position(fpin, rcvSrc[ivsrc].src[isrcm].fpos, &cmaster[nsrc*nfreq], tracebuffer, ntfft, nt, nfreq, nbufmax, trace_sz, &tread, &tfft, &tmemcpy);
							//read_sutrace_at_position(fpin, rcvSrc[ivsrc].src[isrcm].fpos, ivsrc, isrcm, &cmaster[nsrc*nfreq], tracebuffer, trace_in_buffer, &nbuffer, ntfft, nt, nfreq, nbufmax, trace_sz, &trdloc, &tftloc, &tmemcpy);

							/* read slave trace */
							//ipos = (rcvSrc[ivrcv].src[isrcs].fpos)*nfreq;
							//memcpy(&cslaves[nsrc*nfreq], &tracebuffer[ipos],nfreq*sizeof(complex));
							get_sutrace_at_position(fpin, rcvSrc[ivrcv].src[isrcs].fpos, &cslaves[nsrc*nfreq], tracebuffer, ntfft, nt, nfreq, nbufmax, trace_sz, &tread, &tfft, &tmemcpy);
							//read_sutrace_at_position(fpin, rcvSrc[ivrcv].src[isrcs].fpos, ivrcv, isrcs, &cslaves[nsrc*nfreq], tracebuffer, trace_in_buffer, &nbuffer, ntfft, nt, nfreq, nbufmax, trace_sz, &trdloc, &tftloc, &tmemcpy);
							//troutine += wallclock_time()-troutine0;
//}
							nsrc++;
						}
					} /* end of if test that slave and master have same contributing source */
				} /* end of nsrcSlave loop */
			} /* end of nsrcMaster loop */
			//tread_ivs += trdloc;
			//tfft_ivs += tftloc;
			//tloop += wallclock_time()-t0;
			//t1 = wallclock_time();
			//fprintf(stderr,"ivrcv = %d out of %d\n",ivrcv, nreceivers);

			if (nsrc < nsrc_ctr) nsrc = -1; /* only compute virtual receiver when there are sufficient active sources contributing */

			/* correlation of virtual source ivsrc with virtual receiver ivcrv */
			if (cohr) {
				coherence(cmaster, cslaves, nfreq, nsrc, reps, epsmax, &tcorrloc, verbose);
			}
			else {
				correlate(cmaster, cslaves, nfreq, nsrc, &tcorrloc, verbose);
			}
			tcorr_ivs += wallclock_time()-t1;
			
			/* summation over the contributing sources */
			memset(c,0,nf*sizeof(complex));
			for (jsource=0; jsource<nsrc; jsource++) {
				if (normsrc == 1) {
					rms2 = 0.0;
					for (ifreq=0; ifreq<nfreq; ifreq++) {
						tmpcc = cslaves[jsource*nfreq+ifreq];
						rms2 += tmpcc.r*tmpcc.r+tmpcc.i*tmpcc.i;
					}
					if (rms2 > 0) sclrms = 1.0/sqrt(rms2);
					else sclrms = 1.0;
				}
				else (sclrms = 1.0);
				for (ifreq=0; ifreq<nfreq; ifreq++) {
					c[ifreq].r += cslaves[jsource*nfreq+ifreq].r*sclrms;
					c[ifreq].i += cslaves[jsource*nfreq+ifreq].i*sclrms;
				}
			}
		
			//t2_ivs_ivr = wallclock_time();

			cr1_fft(&c[0],r,ntfft,1); /* transform virtual trace back to time */
			
			//t3_ivs_ivr = wallclock_time();
			//tfft_ivs += t3_ivs_ivr-t2_ivs_ivr;
			
			if (normalize && nsrc>0) scl=sclfft/nsrc; /* normalize by the number of contributing sources */
			else scl = sclfft;
			
			if (causal==1) {
				for (j=0;j<nt; j++) {
					vtrace[j] = r[j]*scl;
				}
			}
			else if (causal==2) {
				vtrace[0] = r[0]*scl;
				for (j=1;j<nt; j++) {
					vtrace[j] = r[ntfft-j]*scl;
				}
			}
			else if (causal==3) {
				for (j=1;j<=(nt/2); j++) {
					vtrace[nt/2-j] = r[ntfft-j]*scl;
				}
				for (j=nt/2;j<nt; j++) {
					vtrace[j] = r[j-nt/2]*scl;
				}
			}
			else if (causal==4) {
				vtrace[0] = r[0]*scl;
				for (j=1;j<nt; j++) {
					vtrace[j] = 0.5*(r[ntfft-j] + r[j])*scl;
				}
			} /* one virtual trace (virtual source and virtualreceiver) is now calculated and is written to disk below */
			
			//t4_ivs_ivr = wallclock_time();
			//tselect += t4_ivs_ivr-t3_ivs_ivr;
			
			memcpy(&outputhdr,&segy_hdrs[(rcvSrc[ivrcv].src[0].fpos)],TRCBYTES);
			outputhdr.fldr = ivsrc+1;
			itrace_out=ivsrc*nreceivers+ivrcv+1;
			outputhdr.tracl = itrace_out; // wrong value in principle
			outputhdr.tracf = ivrcv+1;
			outputhdr.sx = sxv;
			outputhdr.sy = syv;
			outputhdr.selev = svelev;
			outputhdr.gx = gxv;
			outputhdr.gy = gyv;
			outputhdr.sdel = vspeg;
			outputhdr.gdel = vrpeg;
			outputhdr.gelev = gvelev;
			outputhdr.ns = nt;
			// outputhdr.offset = (gxv - sxv)/1000;
			#pragma omp critical
{
			nwrite = fwrite( &outputhdr, 1, TRCBYTES, fpout );
			assert (nwrite == TRCBYTES);
			nwrite = fwrite( &vtrace[0], sizeof(float), nt, fpout );
			assert (nwrite == nt);
			fflush(fpout);
}
		} /* end of virtual receiver loop */

		#pragma omp barrier /* to get timing when all threads are done */
		#pragma omp master
{
		if (verbose>=3) {
			t3_ivs = wallclock_time();
        		fprintf(stderr,"Number of OpenMP threads set  = %d number=%d\n", omp_get_max_threads(), omp_get_thread_num());
			//tlogic_ivs = ((t3_ivs-t1_ivs)-tread_ivs-tfft_ivs-tcorr_ivs-twrite_ivs-tmemcpy);
			fprintf(stderr,"************* Timings ************* vshot= %d (ivs = %d)\n", vspeg, ivs);
			//fprintf(stderr,"CPU-time read data         = %.3f\n", tread_ivs);
			//fprintf(stderr,"CPU-time FFT's             = %.3f\n", tfft_ivs);
			//fprintf(stderr,"CPU-time memcpy            = %.3f\n", tmemcpy);
			//fprintf(stderr,"CPU-time memset            = %.3f\n", tmemset);
			//fprintf(stderr,"CPU-time select            = %.3f\n", tselect);
			//fprintf(stderr,"CPU-time loop              = %.3f\n", tloop);
			//fprintf(stderr,"CPU-time routine           = %.3f\n", troutine);
			//fprintf(stderr,"CPU-time correlation       = %.3f\n", tcorr_ivs);
			//fprintf(stderr,"CPU-time write data        = %.3f\n", twrite_ivs);
			//fprintf(stderr,"CPU-time logic             = %.3f\n", tlogic_ivs);
			fprintf(stderr,"Total CPU-time this shot   = %.3f\n", t3_ivs-t1_ivs);
			fflush(stderr);
		}
		ivs++;
}
		#pragma omp barrier
	} /* end of virtual source loop */

	free(cmaster);
	free(cslaves);
	free(vtrace);
	free(r);
	free(c);

} /* end of parallel region */

	fclose(fpin);
	fclose(fpout);

	if (verbose) {
		t3 = wallclock_time();
		//tlogic = ((t3-t1)-tread-tfft-tcorr-twrite);
		ttotal = t3-t1;
		fprintf(stderr,"************* Total Timings ************* \n");
		//fprintf(stderr,"CPU-time initilization     = %.3f\n", tinit);
		//fprintf(stderr,"CPU-time read data         = %.3f\n", tread);
		//fprintf(stderr,"CPU-time FFT's             = %.3f\n", tfft);
		//fprintf(stderr,"CPU-time correlation       = %.3f\n", tcorr);
		//fprintf(stderr,"CPU-time write data        = %.3f\n", twrite);
		//fprintf(stderr,"CPU-time logic             = %.3f\n", tlogic);
		fprintf(stderr,"Total CPU-time             = %.3f\n", ttotal);
	}
/*================ end ================*/
	return 0;
}

void get_sutrace_at_position(FILE *fp, size_t itrace, complex *tracedata, complex *tracebuffer, int ntfft, int nt, int nfreq, size_t nbufmax, size_t trace_sz, double *tread, double *tfft, double *tmemcpy)
{
	static size_t freebuffer=0;
	size_t i, inbuffer, offset, nread, nwrite;
	int ret;
	float *trace;
	float *r;
	double t0, t1, t2, t3, t0_t, t1_t;
	complex *c;
	FILE *fpt;
	char filename[1024];

	if (nbufmax == 0) {
		t0_t = wallclock_time();	
		//trace = (float *)calloc(nt, sizeof(float));
		r = (float *)malloc(ntfft*sizeof(float));
		c = (complex *)malloc((ntfft/2+1)*sizeof(complex));

		/***** read trace from file *****/
       	offset = itrace*trace_sz+TRCBYTES;
       	ret = fseeko(fp, offset, SEEK_SET);
       	if (ret<0) perror("fseeko");
		memset(r,0,ntfft*sizeof(float));
       	nread = fread(r, sizeof(float), nt, fp);
       	assert(nread == nt);
		t1_t = wallclock_time();
		*tread += t1_t-t0_t;

		/***** transform to frequency domain *****/

		//memcpy(r,&trace[0],nt*sizeof(float)); /* might not be needed if segy_read_traces_c is nice */
		rc1_fft(r,c,ntfft,-1);
		memcpy(&tracedata[0].r,&c[0].r,nfreq*sizeof(complex));

		t2 = wallclock_time();
		*tfft += t2-t1_t;

		//free(trace);
		free(r);
		free(c);
	}
	else {  /**** point to trace from buffer *****/
		t0 = wallclock_time();
		//tracedata = (complex *)&tracebuffer[itrace*nfreq]; 
		memcpy(tracedata,&tracebuffer[itrace*nfreq],nfreq*sizeof(complex));
		t1 = wallclock_time();
		*tmemcpy += t1-t0;
	}
	return;
}

/* older implementation when part od traces is kept in memory buffer */

void read_sutrace_at_position_old(FILE *fp, int itrace, int ircvnum, int isrcnum, complex *tracedata, complex *tracebuffer, int *trace_in_buffer, size_t *nbuffer, int ntfft, int nt, int nfreq, size_t nbufmax, size_t trace_sz, double *tread, double *tfft, double *tmemcpy)
{
	static size_t freebuffer=0;
	size_t i, inbuffer, offset, nread, nwrite;
	int ret;
	float *trace;
	float *r;
	double t0, t1, t2, t3, t0_t, t1_t;
	complex *c;
	FILE *fpt;
	char filename[1024];

	inbuffer = -1;
	//fprintf(stderr,"reading trace %d\n", itrace);

	t0_t = wallclock_time();	
	for (i=0; i<nbufmax; i++) {
		if (itrace == trace_in_buffer[i]) {
			inbuffer = i;
			break;
		}
	}
	t1_t = wallclock_time();
	*tread += t1_t-t0_t;
	//fprintf(stderr,"itrace=%8d inbuffer=%8d nbufmax=%d ircvnum=%8d isrcnum=%8d time=%f\n", itrace, inbuffer, nbufmax, ircvnum, isrcnum, t1_t-t0_t);

	if (inbuffer == -1) {
		t0 = wallclock_time();
		trace = (float *)calloc(nt, sizeof(float));
		r = (float *)malloc(ntfft*sizeof(float));
		c = (complex *)malloc((ntfft/2+1)*sizeof(complex));

		/***** read trace from file *****/
        	offset = itrace*trace_sz+TRCBYTES;
        	ret = fseeko(fp, offset, SEEK_SET);
        	if (ret<0) perror("fseeko");
        	nread = fread(trace, sizeof(float), nt, fp);
        	assert(nread == nt);

		/***** transform to frequency domain *****/
		memset(r,0,ntfft*sizeof(float));
		memcpy(r,&trace[0],nt*sizeof(float)); /* might not be needed if segy_read_traces_c is nice */
		t1 = wallclock_time();
		*tread += t1-t0;

		rc1_fft(r,c,ntfft,-1);

		t2 = wallclock_time();
		*tfft += t2-t1;

		memcpy(&tracedata[0].r,&c[0].r,nfreq*sizeof(complex));

		/***** reset buffer counter when buffer is full *****/
		if (nbufmax != 0) {
			t1 = wallclock_time();
			if (freebuffer<nbufmax-1) {
				freebuffer++;
			}
			else freebuffer=0;
		//fprintf(stderr,"freebuffer = %d\n", freebuffer);
		//fprintf(stderr,"nbuffer = %d\n", *nbuffer);

			memcpy(&tracebuffer[freebuffer*nfreq].r,&c[0].r,nfreq*sizeof(complex));
			trace_in_buffer[freebuffer]=itrace;
			t2 = wallclock_time();
			*tmemcpy += t2-t1;
		}

		free(trace);
		free(r);
		free(c);

		t1 = wallclock_time();
		*tread += t1-t2;
	}
	else {  /**** read trace from buffer *****/
		t0 = wallclock_time();
		tracedata = (complex *)&tracebuffer[inbuffer*nfreq]; /* this is better, but need some more changes TODO */
		//fprintf(stderr,"read from inbuffer=%d nbufmax=%d\n", inbuffer,nbufmax);

		//memcpy(&tracedata[0].r,&tracebuffer[inbuffer*nfreq].r,nfreq*sizeof(complex));
		
		t1 = wallclock_time();
		*tread += t1-t0;
	}
	return;
}

