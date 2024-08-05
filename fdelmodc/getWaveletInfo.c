#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "par.h"
#include "segy.h"

/**
*  reads file which contain the source wavelets and computes sampling interval
*  and tries to estimate the maximum frequency content.
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
typedef struct _dcomplexStruct { /* complex number */
    double r,i;
} dcomplex;
#endif/* complex */

int optncr(int n);
void rc1fft(float *rdata, complex *cdata, int n, int sign);

#define     MAX(x,y) ((x) > (y) ? (x) : (y))
#define     MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int getWaveletInfo(char *file_src, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *fmax, int *nxm, int verbose)
{
    FILE    *fp;
    size_t  nread, trace_sz;
    off_t   bytes;
    int     ret, one_shot, ntraces;
    int     optn, nfreq, i, iwmax;
    float   *trace, df;
	float   ampl, amplmax, tampl, tamplmax; 
    complex *ctrace;
    segy hdr;
    
    if (file_src == NULL) return 0; /* Input pipe can not be handled */
    else fp = fopen( file_src, "r" );
    assert( fp != NULL);
    if (strstr(file_src, ".su") == NULL) { /* assume a binary file is read in and read number of samples in file */
        ret = fseeko( fp, 0, SEEK_END );
        bytes = ftello( fp );
        *n1 = (int)(bytes/sizeof(float));
        *f1 = 0.0;
        *f2 = 0.0;
        ntraces = 1;
		if(!getparint("ntraces",&ntraces)) vwarn("For binary file_src ntraces is set to 1");
		if(!getparfloat("dt",d1)) verr("For binary file_src dt must be given");
		if(!getparfloat("fmax",fmax)) verr("For binary file_src fmax in wavelet must be given");
        *n2 = *nxm = ntraces;
        *n1 = *n1 / ntraces;
		*d2 = 1.0;
		return 0;
    }

    nread = fread( &hdr, 1, TRCBYTES, fp );
    assert(nread == TRCBYTES);
    ret = fseeko( fp, 0, SEEK_END );
	if (ret<0) perror("fseeko");
    bytes = ftello( fp );

    *n1 = hdr.ns;
    if (hdr.trid == 1 || hdr.dt != 0) {
        *d1 = ((float) hdr.dt)*1.e-6;
        if (hdr.delrt != 0) *f1 = ((float) hdr.delrt)/1000.;
        else *f1 = hdr.f1;
		if (*d1 == 0.0) *d1 = hdr.d1;
    }
    else {
        *d1 = hdr.d1;
        *f1 = hdr.f1;
    }
    *f2 = hdr.f2;

    trace_sz = (size_t)(sizeof(float)*(*n1)+TRCBYTES);
    ntraces  = (int) (bytes/trace_sz);
	*n2 = ntraces;

    /* check to find out number of traces in shot gather */

	optn  = optncr(*n1);
	nfreq = optn/2 + 1;
    df    = 1.0/(optn*(*d1));
	ctrace = (complex *)malloc(nfreq*sizeof(complex));
    one_shot = 1;
    trace = (float *)malloc(optn*sizeof(float));
    fseeko( fp, TRCBYTES, SEEK_SET );

    while (one_shot) {
		memset(trace,0,optn*sizeof(float));
        nread = fread( trace, sizeof(float), *n1, fp );
        assert (nread == *n1);
		tamplmax = 0.0;
		for (i=0;i<(*n1);i++) {
			tampl = fabsf(trace[i]);
			if (tampl > tamplmax) tamplmax = tampl;
		}
		if (trace[0]*1e-3 > tamplmax) {
			fprintf(stderr,"WARNING: file_src has a large amplitude %f at t=0\n", trace[0]);
			fprintf(stderr,"This will introduce high frequencies and can cause dispersion.\n");
		}

		/* estimate maximum frequency assuming amplitude spectrum is smooth */
		rc1fft(trace,ctrace,optn,1);

		/* find maximum amplitude */
		amplmax = 0.0;
		iwmax = 0;
		for (i=0;i<nfreq;i++) {
			ampl = sqrt(ctrace[i].r*ctrace[i].r+ctrace[i].i*ctrace[i].i);
			if (ampl > amplmax) {
				amplmax = ampl;
				iwmax = i;
			}
		}

		/* from the maximum amplitude position look for the largest frequency
         * which has an amplitude 400 times weaker than the maximum amplitude */
		for (i=iwmax;i<nfreq;i++) {
			ampl = sqrt(ctrace[i].r*ctrace[i].r+ctrace[i].i*ctrace[i].i);
			if (400*ampl < amplmax) {
				*fmax = (i-1)*df;
				break;
			}
		}

        nread = fread( &hdr, 1, TRCBYTES, fp );
        if (nread==0) break;
    }
	*nxm = (int)ntraces;

	if (verbose>2) {
		vmess("For file %s", file_src);
		vmess("nt=%d nx=%d", *n1, *n2);
		vmess("dt=%f dx=%f", *d1, *d2);
		vmess("fmax=%f", *fmax);
		vmess("tstart=%f", *f1);
	}

    fclose(fp);
    free(trace);
    free(ctrace);

    return 0;
}
