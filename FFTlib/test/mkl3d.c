#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <genfft.h>
#include "mkl_dfti.h"

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

int main (int argc, char **argv)
{
    long    nt, nz, ny, nx, ntr, ix, iy, it, is, iz, pos, file_det, nzs, dt;
    float   *indata, *rdata;
    complex *cdata, *cdataT;
    long    nf, nft, nkx, nky, xorig, yorig;

    DFTI_DESCRIPTOR_HANDLE my_desc_handle = NULL;
    MKL_LONG status;
    MKL_LONG dim_sizes[3];
    MKL_LONG strides_out[4];
    MKL_LONG strides_in[4];

    nt=8;
    nft=nt;
    nf=nft/2+1;
    nkx=3;
    nky=4;
    xorig=0;
    yorig=0;
    xorig=nx/2;
    yorig=ny/2;
    rdata        = (float *)calloc(nkx*nky*nft,sizeof(float));
    cdata        = (complex *)calloc(nkx*nky*nf,sizeof(complex));
    rdata[1] = 1.0;

    dim_sizes[0] = nky; dim_sizes[1]=nkx; dim_sizes[2]=nft;
    strides_in[0] = 0; strides_in[1]=nkx*nft; strides_in[2]=nft; strides_in[3]=1;
    strides_out[0] = 0; strides_out[1]=nkx*nf; strides_out[2]=nf; strides_out[3]=1;

    status = DftiCreateDescriptor(&my_desc_handle, DFTI_SINGLE, DFTI_REAL, 3, dim_sizes);
    status = DftiSetValue(my_desc_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    status = DftiSetValue(my_desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status = DftiSetValue(my_desc_handle, DFTI_OUTPUT_STRIDES, strides_out);
    status = DftiSetValue(my_desc_handle, DFTI_INPUT_STRIDES, strides_in);

    status = DftiCommitDescriptor(my_desc_handle);

    /* forward transform from libgenfft.a */
    yxt2wkykx(rdata, cdata, nft, nkx, nky, nft, nkx, nky, xorig, yorig);

    fprintf(stderr,"*** FFTlib after forward transform [nf][nky][nkx] ***\n");
    for (it=0; it<nf; it++) {
    for (iy=0; iy<nky; iy++) {
    for (ix=0; ix<nkx; ix++) {
	    fprintf(stderr,"out[%d][%d][%d]=%e %e\n", it, iy, ix, cdata[it*nkx*nky+iy*nkx+ix].r, cdata[it*nkx*nky+iy*nkx+ix].i );
    }
    }
    }

    /* backward transform from libgenfft.a */
    wkykx2yxt(cdata, rdata, nft, nkx, nky, nft, nkx, nky, xorig, yorig);

    fprintf(stderr,"*** FFTlib after backward transform [nky][nkx][nt] ***\n");
    for (iy=0; iy<nky; iy++) {
    for (ix=0; ix<nkx; ix++) {
    for (it=0; it<nft; it++) {
	    fprintf(stderr,"out[%d][%d][%d]=%e \n", iy, ix, it, rdata[iy*nkx*nft+ix*nft+it] );
    }
    }
    }

    memset(&cdata[0].r, 0, sizeof(complex)*nkx*nky*nf);
    //memset(&rdata[0], 0, sizeof(float)*nkx*nky*nft);
    //rdata[1] = 1.0;

    /* forward transform from MKL */
    status = DftiComputeForward(my_desc_handle, rdata, (MKL_Complex8 *)&cdata[0]);

    fprintf(stderr,"*** MKL after forward transform [nky][nkx][nf] ***\n");
    for (iy=0; iy<nky; iy++) {
    for (ix=0; ix<nkx; ix++) {
    for (it=0; it<nf; it++) {
	    fprintf(stderr,"out[%d][%d][%d]=%e %e\n", iy, ix, it, cdata[iy*nkx*nf+ix*nf+it].r, cdata[iy*nkx*nf+ix*nf+it].i );
    }
    }
    }

    /* tranpose to [nf][nky][nkx] */
    /* the signs of the FFT kernels in MKL will be different than in yxt2wkykx(time forward -1 x,y 1) */
    cdataT        = (complex *)calloc(nkx*nky*nf,sizeof(complex));
    fprintf(stderr,"*** MKL after forward transform tranposed to [nf][nky][nkx] ***\n");
    for (it=0; it<nf; it++) {
    for (iy=0; iy<nky; iy++) {
    for (ix=0; ix<nkx; ix++) {
	    cdataT[it*nkx*nky+iy*nkx+ix] = cdata[iy*nkx*nf+ix*nf+it];
	    fprintf(stderr,"out[%d][%d][%d]=%e %e\n", it, iy, ix, cdataT[it*nkx*nky+iy*nkx+ix].r, cdataT[it*nkx*nky+iy*nkx+ix].i );
    }
    }
    }

    /* backward transform from MKL */
    strides_out[0] = 0; strides_out[1]=nkx*nft; strides_out[2]=nft; strides_out[3]=1;
    strides_in[0] = 0; strides_in[1]=nkx*nf; strides_in[2]=nf; strides_in[3]=1;
    status = DftiSetValue(my_desc_handle, DFTI_INPUT_STRIDES, strides_in);
    status = DftiSetValue(my_desc_handle, DFTI_OUTPUT_STRIDES, strides_out);
    status = DftiCommitDescriptor(my_desc_handle);
    status = DftiComputeBackward(my_desc_handle, (MKL_Complex8 *)&cdata[0], rdata);

    fprintf(stderr,"*** MKL after backward transform [nky][nkx][nt] ***\n");
    for (iy=0; iy<nky; iy++) {
    for (ix=0; ix<nkx; ix++) {
    for (it=0; it<nft; it++) {
	    fprintf(stderr,"out[%d][%d][%d]=%e \n", iy, ix, it, rdata[iy*nkx*nft+ix*nft+it] );
    }
    }
    }
    
    status = DftiFreeDescriptor(&my_desc_handle);
    
    return 0;
}

