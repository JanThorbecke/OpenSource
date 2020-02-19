#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "segy.h"
#include "zfpmar.h"
#include <assert.h>
#include <zfp.h>

typedef struct { /* complex number */
        float r,i;
} complex;

#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

int optncr(int n);
long zfpcompress(float* data, long nx, long ny, long nz, double tolerance, zfpmar zfpm, FILE *file);
void cc1fft(complex *data, int n, int sign);
void rc1fft(float *rdata, complex *cdata, int n, int sign);
long writeData3D(FILE *fp, float *data, segy *hdrs, long n1, long n2);
long getFileInfo3D(char *filename, long *n1, long *n2, long *n3, long *ngath, float *d1, float *d2, float *d3, float *f1, float *f2, float *f3,
    float *sclsxgxsygy, long *nxm);

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" TWtransform - Transform data from uncompressed time domain to compressed frequency domain",
" ",
" TWtransform file_in= file_out= [optional parameters]",
" ",
" Required parameters: ",
" ",
"   file_in= ................. File containing the uncompressed time domain data",
"   file_out= ................ Output for the (compressed) frequency domain data",
" ",
" Optional parameters: ",
" ",
"   verbose=1 ................ silent option; >0 displays info",
"   fmin=0 ................... minimum frequency in the output",
"   fmax=70 .................. maximum frequency in the output",
"   mode=1 ................... sign of the frequency transform",
"   zfp=0 .................... (=1) compress the transformed data using zfp",
"   tolerance=1e-3 ........... accuracy of the zfp compression,",
"   smaller values give more accuracy to the compressed data but will decrease the compression rate",
"   weight=2.0 ............... scaling of the reflection data",
" ",
" ",
" author  : Joeri Brackenhoff : (j.a.brackenhoff@tudelft.nl)",
" author  : Jan Thorbecke     : (j.w.thorbecke@tudelft.nl)",
" ",
NULL};
/**************** end self doc ***********************************/

int main (int argc, char **argv)
{
	FILE    *fp, *fp_out;
	segy    hdr, *hdr_out;
	size_t  nread;
	long    fldr_shot, sx_shot, sy_shot, itrace, one_shot, igath, iw;
	long    end_of_file, nt, ntfft, nxy, nx, ny, nshots, ret, ntraces;
	long    nfreq, nw_low, nw_high, nw, mode, verbose, zfp;
	long	inx, iny, gy;
	float   scl, *trace, dt, dx, dy, ft, fx, fy, fmin, fmax, *cdata, scale;
	double	tolerance;
    char    *file_in, *file_out;
	complex *ctrace;
	zfptop	zfpt;
	zfpmar  zfpm;

    initargs(argc, argv);
    requestdoc(1);

    if (!getparstring("file_in", &file_in)) file_in = NULL;
        if (file_in==NULL) verr("No file_in is given");
    if (!getparstring("file_out", &file_out)) file_out = NULL;
        if (file_out==NULL) verr("No file_out is given");
    if (!getparfloat("fmin", &fmin)) fmin = 0;
    if (!getparfloat("fmax", &fmax)) fmax = 70;
    if (!getpardouble("tolerance", &tolerance)) tolerance = 1e-3;
    if (!getparfloat("weight", &scale)) scale = 2.0;
    if (!getparlong("mode", &mode)) mode = 1;
    if (!getparlong("verbose", &verbose)) verbose = 1;
    if (!getparlong("zfp", &zfp)) zfp = 0;

    nshots = 0;
    ret = getFileInfo3D(file_in, &nt, &nx, &ny, &nshots, &dt, &dx, &dy, &ft, &fx, &fy, &scl, &ntraces);

    ntfft = loptncr(nt); 
    nfreq = ntfft/2+1;
    nw_low = (long)MIN((fmin*ntfft*dt), nfreq-1);
    nw_low = MAX(nw_low, 1);
    nw_high = MIN((long)(fmax*ntfft*dt), nfreq-1);
    nw  = nw_high - nw_low + 1;

    nxy = nx*ny;

	/* Reading first header  */

	if (file_in == NULL) fp = stdin;
	else fp = fopen( file_in, "r" );
	if ( fp == NULL ) {
		fprintf(stderr,"input file %s has an error\n", file_in);
		perror("error in opening file: ");
		fflush(stderr);
		return -1;
	}

	if (file_out == NULL) fp_out = stdin;
	else fp_out = fopen( file_out, "w+" );
	if ( fp_out == NULL ) {
		fprintf(stderr,"input file %s has an error\n", file_out);
		perror("error in opening file: ");
		fflush(stderr);
		return -1;
	}

	fseek(fp, 0, SEEK_SET);
	nread = fread( &hdr, 1, TRCBYTES, fp );
	assert(nread == TRCBYTES);

	fseek(fp, 0, SEEK_SET);

	nt = hdr.ns;
	dt = hdr.dt/(1E6);

	trace   	= (float *)calloc(ntfft,sizeof(float));
	ctrace  	= (complex *)malloc(ntfft*sizeof(complex));
	cdata  		= (float *)calloc(nw*2*nxy,sizeof(float));
	hdr_out		= (segy *) calloc(nxy,sizeof(segy));

	end_of_file = 0;
	one_shot    = 1;
	igath       = 0;

	if (zfp) {
		vmess("zfp compression applied");
		zfpt.dx 		= dx;
		zfpt.dy 		= dy;
		zfpt.dz 		= dt;
		zfpt.ndim 		= 3;
		zfpt.ns 		= nshots;
		zfpt.scale 		= hdr.scalco;
		zfpt.nt			= ntfft;
		zfpt.fmin		= fmin;
		zfpt.fmax		= fmax;
		zfpt.nz			= 2*nw;
		zfpt.tolerance	= tolerance;
		zfpt.fz			= ft;
		zfpt.fx			= fx;
		zfpt.fy			= fy;

		nread = fwrite(&zfpt, 1, TOPBYTES, fp_out);
		assert(nread == TOPBYTES);
	}
	else {
		vmess("no zfp compression applied");
	}

	/* Read shots in file */

	while (!end_of_file) {

		/* start reading data (shot records) */
		itrace	= 0;
		iny 	= 1;
		nread = fread( &hdr, 1, TRCBYTES, fp );
		if (nread != TRCBYTES) { /* no more data in file */
			break;
		}

		sx_shot		= hdr.sx;
        sy_shot		= hdr.sy;
		gy			= hdr.gy;
		fldr_shot	= hdr.fldr;
		if (zfp) {
			zfpm.gx	= hdr.gx;
			zfpm.gy	= hdr.gy;
			zfpm.sx	= hdr.sx;
			zfpm.sy	= hdr.sy;
			zfpm.sz	= hdr.selev;
		}
		while (one_shot) {

			if (hdr.gy != gy) {
				iny++;
				gy = hdr.gy;
			}

			nread = fread( trace, sizeof(float), nt, fp );
			assert (nread == hdr.ns);

			/* transform to frequency domain */
			if (ntfft > hdr.ns) 
			memset( &trace[nt-1], 0, sizeof(float)*(ntfft-nt) );

			rc1fft(trace,ctrace,(int)ntfft,-1);
			for (iw=0; iw<nw; iw++) {
				cdata[itrace*nw*2+(iw*2)]	= scale*ctrace[nw_low+iw].r;
				cdata[itrace*nw*2+(iw*2)+1]	= scale*mode*ctrace[nw_low+iw].i;
			}
			itrace++;

            hdr.ep				= ntfft;
            hdr.ns      		= 2*nw;
            hdr.unscale 		= fmin;
			hdr.ungpow  		= fmax;
			hdr_out[itrace-1] 	= hdr;

			/* read next hdr of next trace */
			nread = fread( &hdr, 1, TRCBYTES, fp );
			if (nread != TRCBYTES) { 
				one_shot = 0;
				end_of_file = 1;
				break;
			}
			if ((sx_shot != hdr.sx) || (sy_shot != hdr.sy) || (fldr_shot != hdr.fldr)) break;
		}
		inx = itrace/iny;
		if (verbose) {
			vmess("finished reading shot x=%li y=%li (%li) with %li traces (nx=%li ny=%li) and weight=%.3f",sx_shot,sy_shot,igath,itrace,inx,iny,scale);
		}

		if (zfp) {
			zfpcompress(cdata,nx,ny,2*nw,tolerance,zfpm,fp_out);
		}
		else {
			ret = writeData3D(fp_out, (float *)&cdata[0], hdr_out, 2*nw, itrace);
            if (ret < 0 ) verr("error on writing output file.");
		}

		if (itrace != 0) { /* end of shot record */
			fseek( fp, -TRCBYTES, SEEK_CUR );
			igath++;
		}
		else {
			end_of_file = 1;
		}
	}
	fclose(fp_out);
	free(ctrace);
	free(trace);

	exit(0);
}

long zfpcompress(float* data, long nx, long ny, long nz, double tolerance, zfpmar zfpm, FILE *file)
{

	zfp_field*			field = NULL;
	zfp_stream* 		zfp = NULL;
	bitstream* 			stream = NULL;
	void* 				fi = NULL;
	void* 				fo = NULL;
	void* 				buffer = NULL;
	size_t 				rawsize = 0;
	size_t 				zfpsize = 0;
	size_t 				bufsize = 0;
	size_t				nwrite;
	zfp_exec_policy 	exec = zfp_exec_serial;

	zfp = zfp_stream_open(NULL);
	field = zfp_field_alloc();

	zfp_field_set_pointer(field, (void *)data);

	zfp_field_set_type(field, zfp_type_float);
	zfp_field_set_size_3d(field, (uint)nz, (uint)nx, (uint)ny);

	zfp_stream_set_accuracy(zfp, tolerance);

	if (!zfp_stream_set_execution(zfp, exec)) {
    	fprintf(stderr, "serial execution not available\n");
    	return EXIT_FAILURE;
    }

	bufsize = zfp_stream_maximum_size(zfp, field);
	if (!bufsize) {
      fprintf(stderr, "invalid compression parameters\n");
      return EXIT_FAILURE;
    }

	buffer = malloc(bufsize);
	if (!buffer) {
      fprintf(stderr, "cannot allocate memory\n");
      return EXIT_FAILURE;
    }

	stream = stream_open(buffer, bufsize);
    if (!stream) {
      fprintf(stderr, "cannot open compressed stream\n");
      return EXIT_FAILURE;
    }
    zfp_stream_set_bit_stream(zfp, stream);

	if (!zfp_stream_set_execution(zfp, exec)) {
        fprintf(stderr, "serial execution not available\n");
        return EXIT_FAILURE;
    }

    zfpsize = zfp_compress(zfp, field);
	if (zfpsize == 0) {
      fprintf(stderr, "compression failed\n");
      return EXIT_FAILURE;
    }

	zfpm.nx = nx;
	zfpm.ny = ny;
	zfpm.compsize = zfpsize;

	// file = fopen(zfppath, "wb");
	// if (file==NULL) {
	// 	fprintf(stderr,"input file %s has an error\n", zfppath);
	// 	perror("error in opening file: ");
	// 	fflush(stderr);
	// 	return -1;
	// }
	nwrite = fwrite(&zfpm, 1, MARBYTES, file);
	assert(nwrite == MARBYTES);
	if (fwrite(buffer, 1, zfpsize, file) != zfpsize) {
        fprintf(stderr, "cannot write compressed file\n");
        return EXIT_FAILURE;
    }
    // fclose(file);

	return 1;
}
