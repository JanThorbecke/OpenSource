#include "genfft.h"
#if defined(MKL)
#include "mkl_dfti.h"
#endif

/**
NAME:   wkykx2yxt

DESCRIPTION:
        3 Dimensional complex to real FFT (omega,kx -> x,t) with scaling

USAGE:
        void wkykx2yxt(complex *cdata, REAL *rdata, long nt, long nx, long ny, 
	        long ldt, long ldx, long ldy, long xorig, long yorig)
INPUT:
        - *cdata: complex 3D input array [nf][nky][nkx]
        -     nt: number of real (fast) output samples
        -     nx: number of complex (slow) samples to be transformed
        -     ny: number of complex (slow) samples to be transformed
        -    ldy: leading dimension (number of complex samples)
        -    ldx: leading dimension (number of complex samples)
        -    ldt: leading dimension (number of real samples)
        -  xorig: trace number of origin of x-axis first trace # 0 
        -  yorig: trace number of origin of y-axis first trace # 0 

OUTPUT: - *rdata: real 3D output array scaled [ny][nx][nt]

AUTHOR:
           Jan Thorbecke (janth@xs4all.nl)
           The Netherlands
----------------------------------------------------------------------*/

void wkykx2yxt(complex *cdata, REAL *rdata, long nt, long nx, long ny, long ldt, long ldx, long ldy, long xorig, long yorig)
{
    int     i, j, k, ld1, sign, nxorig, nyorig;
    REAL   scl, *rdum;
    complex *cdum;
#if defined(MKL)
    DFTI_DESCRIPTOR_HANDLE my_desc_handle = NULL;
    MKL_LONG status;
    MKL_LONG dim_sizes[3];
    MKL_LONG strides_out[4];
    MKL_LONG strides_in[4];
#endif

    assert ( (xorig >= 0) && (xorig < nx) );
    assert ( (yorig >= 0) && (yorig < ny) );
    scl = 1.0/(nt*nx*ny);
    ld1 = (nt+2)/2;
    nxorig = nx-xorig;
    nyorig = ny-yorig;

#if defined(MKL)
    //fprintf(stderr,"Using MKL for wkykx2yxt.c\n");
    cdum = (complex *)malloc(ld1*ldx*ldy*sizeof(complex));
    if (cdum == NULL) fprintf(stderr,"wkykx2yxt: memory allocation error\n");

    rdum = (float *)malloc(ldt*ldx*ldy*sizeof(float));
    if (rdum == NULL) fprintf(stderr,"wkykx2yxt: memory allocation error\n");
    
    /* tranpose cdata to [nky][nkx][nf] */
    /* the signs of the FFT kernels in MKL will be different than in yxt2wkykx(time forward -1 x,y 1) */
    for (j=0; j<ld1; j++) {
        for (k=0; k<ldy; k++) {
            for (i=0; i<ldx; i++) {
                cdum[k*ldx*ld1+i*ld1+j] = cdata[j*ldx*ldy+k*ldx+i];
            }
        }
    }

    dim_sizes[0]=ny;  dim_sizes[1]=nx;        dim_sizes[2]=nt;
    strides_in[0]=0;  strides_in[1]=ldx*ld1;  strides_in[2]=ld1; strides_in[3]=1;
    strides_out[0]=0; strides_out[1]=ldx*ldt; strides_out[2]=ldt; strides_out[3]=1;

    status = DftiCreateDescriptor(&my_desc_handle, DFTI_SINGLE, DFTI_REAL, 3, dim_sizes);
    status = DftiSetValue(my_desc_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    status = DftiSetValue(my_desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status = DftiSetValue(my_desc_handle, DFTI_OUTPUT_STRIDES, strides_out);
    status = DftiSetValue(my_desc_handle, DFTI_INPUT_STRIDES, strides_in);
    status = DftiCommitDescriptor(my_desc_handle);

    /* backward transform from MKL */
    status = DftiComputeBackward(my_desc_handle, (MKL_Complex8 *)&cdum[0], rdum);

    free(cdum);
    status = DftiFreeDescriptor(&my_desc_handle);

    /* set back to xorig and yorig */
    //rdum [ny][nx][nf] rearrange nkx and nky with respect to x(y)orig
    for (k = 0; k < nyorig; k++) {
        for (i = 0; i < nxorig; i++) {
            for (j = 0; j < nt; j++) {
                rdata[(k+yorig)*nx*ldt+(i+xorig)*ldt+j] = rdum[k*ldt*ldx+i*ldt+j]*scl;
            }
        }
        for(i = 0; i < xorig; i++) {
            for (j = 0; j < nt; j++) {
                rdata[(k+yorig)*nx*ldt+i*ldt+j] = rdum[k*ldt*ldx+(i+nxorig)*ldt+j]*scl;
            }
        }
    }
    for (k = 0; k < yorig; k++) {
        for (i = 0; i < nxorig; i++) {
            for (j = 0; j < nt; j++) {
                rdata[k*nx*ldt+(i+xorig)*ldt+j] = rdum[(k+nyorig)*ldt*ldx+i*ldt+j]*scl;
            }
        }
        for(i = 0; i < xorig; i++) {
            for (j = 0; j < nt; j++) {
                rdata[k*nx*ldt+i*ldt+j] = rdum[(k+nyorig)*ldt*ldx+(i+nxorig)*ldt+j]*scl;
            }
        }
    }
    free(rdum);

#else
    cdum = (complex *)malloc(ld1*ldx*ldy*sizeof(complex));    
    if (cdum == NULL) fprintf(stderr,"wkykx2yxt: memory allocation error\n");

    sign = -1;
    // cdum [nkx][nf][nky] Tranform over nky
    for (i = 0; i < ldx; i++) {
        for (j = 0; j < ld1; j++) {
            for (k = 0; k < ldy; k++) {
                cdum[i*ld1*ldy+j*ldy+k] = cdata[j*ldy*ldx+k*ldx+i];
            }
        }
        ccmfft(&cdum[i*ld1*ldy], ny, ld1, ldy, sign);
    }

    // cdata [nf][nky][nkx] Tranform over nkx
    for (k = 0; k < ldy; k++) {
        for (j = 0; j < ld1; j++) {
            for (i = 0; i < ldx; i++) {
                cdata[k*ld1*ldx+j*ldx+i] = cdum[i*ld1*ldy+j*ldy+k];
            }
        }
        ccmfft(&cdata[k*ld1*ldx], nx, ld1, ldx, sign);
    }
    free(cdum);

    cdum = (complex *)malloc(ld1*nx*ny*sizeof(complex));
    if (cdum == NULL) fprintf(stderr,"wkx2xt: memory allocation error\n");
    
    //cdum [ny][nx][nf] rearrange nkx and nky with respect to x(y)orig
    for (j = 0; j < ld1; j++) {
        for (k = 0; k < nyorig; k++) {
            for (i = 0; i < nxorig; i++) {
                cdum[(k+yorig)*nx*ld1+(i+xorig)*ld1+j] = cdata[k*ld1*ldx+j*ldx+i];
            }
            for(i = 0; i < xorig; i++) {
                cdum[(k+yorig)*nx*ld1+i*ld1+j] = cdata[k*ld1*ldx+j*ldx+i+nxorig];
            }
        }
        for (k = 0; k < yorig; k++) {
            for (i = 0; i < nxorig; i++) {
                cdum[k*nx*ld1+(i+xorig)*ld1+j] = cdata[(k+nyorig)*ld1*ldx+j*ldx+i];
            }
            for(i = 0; i < xorig; i++) {
                cdum[k*nx*ld1+i*ld1+j] = cdata[(k+nyorig)*ld1*ldx+j*ldx+i+nxorig];
            }
        }
    }

    sign = 1;
    crmfft(cdum, rdata, nt, nx*ny, ld1, ldt, sign);

    for(i = 0; i < nx*ny*ldt; i++) rdata[i] *= scl;

    free(cdum);
#endif


    return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define nwkykx2yxt    FNAME(WKYKX2YXTF)
#else
#define nwkykx2yxt    FNAME(wkykx2yxtf)
#endif

void nwkykx2yxt(complex *cdata, REAL *rdata, int *nt, int *nx, int *ny, int *ldt, int *ldx, int *ldy, int *xorig, int *yorig)
{

    wkykx2yxt(cdata, rdata, *nt, *nx, *ny, *ldt, *ldx, *ldy, *xorig, *yorig);

    return;
}
