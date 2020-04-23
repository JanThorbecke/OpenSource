#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "segy.h"
#include <assert.h>
#include "zfpmar.h"
#include <zfp.h>

typedef struct { /* complex number */
        float r,i;
} complex;

#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

int optncr(int n);
void cc1fft(complex *data, int n, int sign);
void rc1fft(float *rdata, complex *cdata, int n, int sign);
long zfpdecompress(float* data, long nx, long ny, long nz, long comp, double tolerance, FILE *file);

long readShotData3Dzfp(char *filename, float *xrcv, float *yrcv, float *xsrc, float *ysrc, float *zsrc,
	long *xnx, complex *cdata, long nw, long nshots, long nx, long ny, float scale, long verbose)
{
	FILE    *fp;
	size_t  nread;
	long    sx_shot, sy_shot, iw, nxy, ishot, inx, iny, ix, iy;
	float   *trace, dx, dy, scl;
    double  tolerance;
    zfpmar  zfpm;
    zfptop  zfpt;

    nxy = nx*ny;

	/* Reading first header  */

	if (filename == NULL) fp = stdin;
	else fp = fopen( filename, "r" );
	if ( fp == NULL ) {
		fprintf(stderr,"input file %s has an error\n", filename);
		perror("error in opening file: ");
		fflush(stderr);
		return -1;
	}

	fseek(fp, 0, SEEK_SET);
	nread = fread( &zfpt, 1, TOPBYTES, fp );
	assert(nread == TOPBYTES);

	if (verbose) vmess("Reading in zfp frequency traces");

    tolerance   = zfpt.tolerance;
    dx          = zfpt.dx;
    dy          = zfpt.dy;

    if (zfpt.scale < 0.0) scl = 1.0/fabs((float)zfpt.scale);
	else if (zfpt.scale == 0.0) scl = 1.0;
	else scl = zfpt.scale;

    trace  = (float *)calloc(2*nw,sizeof(float));

	/* Read shots in file */

	for (ishot=0; ishot<nshots; ishot++) {

		/* start reading data (shot records) */
		nread = fread( &zfpm, 1, MARBYTES, fp );
		if (nread != MARBYTES) { /* no more data in file */
			break;
		}

        inx         = zfpm.nx;
        iny         = zfpm.ny;
		sx_shot     = zfpm.sx;
        sy_shot     = zfpm.sy;
		xsrc[ishot] = sx_shot*scl;
		ysrc[ishot] = sy_shot*scl;
		zsrc[ishot] = zfpm.sz*scl;
		xnx[ishot]  = inx*iny;

        trace  = (float *)realloc(trace, inx*iny*2*nw*sizeof(float));
        zfpdecompress(trace, inx, iny, 2*nw, zfpm.compsize , tolerance, fp);

		for (iy = 0; iy < iny; iy++) {
            for (ix = 0; ix < inx; ix++) {
                for (iw = 0; iw < nw; iw++) {
                    cdata[ishot*nxy*nw+iw*nxy+iy*inx+ix].r = trace[iy*inx*nw*2+ix*nw*2+(iw*2)];
                    cdata[ishot*nxy*nw+iw*nxy+iy*inx+ix].i = trace[iy*inx*nw*2+ix*nw*2+(iw*2)+1];
                }

                xrcv[ishot*nxy+iy*inx+ix] = (zfpm.gx*scl)+(ix*dx);
                yrcv[ishot*nxy+iy*inx+ix] = (zfpm.gy*scl)+(iy*dy);
            }
		}
		if (verbose>2) {
			vmess("finished reading shot x=%li y=%li (%li) with %li traces (nx=%li ny=%li)",sx_shot,sy_shot,ishot,inx*iny,inx,iny);
		}
	}

	fclose(fp);
	free(trace);

	return 0;
}

long zfpdecompress(float* data, long nx, long ny, long nz, long comp, double tolerance, FILE *file)
{
	zfp_field*			field = NULL;
	zfp_stream* 		zfp = NULL;
	bitstream* 			stream = NULL;
	zfp_exec_policy 	exec = zfp_exec_serial;
	size_t				nread, compsize;
	void 				*buffer;

	zfp = zfp_stream_open(NULL);
  	field = zfp_field_alloc();
	compsize = comp;

	buffer = malloc(compsize);
	if (!buffer) {
      fprintf(stderr, "cannot allocate memory\n");
      return EXIT_FAILURE;
    }
	nread = fread((uchar*)buffer, 1, compsize, file);
	assert(nread==compsize);

	stream = stream_open(buffer, compsize);
    if (!stream) {
      fprintf(stderr, "cannot open compressed stream\n");
      return EXIT_FAILURE;
    }
    zfp_stream_set_bit_stream(zfp, stream);

	zfp_field_set_type(field, zfp_type_float);
	zfp_field_set_size_3d(field, (uint)nz, (uint)nx, (uint)ny);

	zfp_stream_set_accuracy(zfp, tolerance);

	if (!zfp_stream_set_execution(zfp, exec)) {
    	fprintf(stderr, "serial execution not available\n");
    	return EXIT_FAILURE;
    }

	zfp_stream_rewind(zfp);

	if (!zfp_stream_set_execution(zfp, exec)) {
		fprintf(stderr, "serial execution not available\n");
		return EXIT_FAILURE;
	}

	zfp_field_set_pointer(field, (void *)data);

	while (!zfp_decompress(zfp, field)) {
      /* fall back on serial decompression if execution policy not supported */
      if (zfp_stream_execution(zfp) != zfp_exec_serial) {
        if (!zfp_stream_set_execution(zfp, zfp_exec_serial)) {
          fprintf(stderr, "cannot change execution policy\n");
          return EXIT_FAILURE;
        }
      }
      else {
        fprintf(stderr, "decompression failed\n");
        return EXIT_FAILURE;
      }
    }

	return 1;
}