#include "genfft.h"
#if defined(MKL)
#include "mkl_dfti.h"
#endif

/**
NAME:   yxt2wkykx

DESCRIPTION:
        3 Dimensional real to complex FFT (y, x,t -> omega,ky, kx) with scaling

USAGE:
        void yxt2wkykx(REAL *rdata, complex *cdata, long nt, long nx, long ny, 
	      long ldt, long ldx, long ldy, long xorig, long yorig)
INPUT:
        - *rdata: real 3D input array [ny][nx][nt]
        -     nt: number of real (fast) output samples
        -     nx: number of complex (slow) samples to be transformed
        -     ny: number of complex (slow) samples to be transformed
        -    ldt: leading dimension (number of real samples)
        -    ldx: leading dimension (number of complex samples)
        -    ldy: leading dimension (number of complex samples)
        -  xorig: trace number of origin of x-axis first trace # 0
        -  yorig: trace number of origin of y-axis first trace # 0

OUTPUT: - *cdata: complex 3D output array scaled [nf][nky][nkx]

AUTHOR:
           Jan Thorbecke (janth@xs4all.nl)
           The Netherlands
----------------------------------------------------------------------*/

void yxt2wkykx(REAL *rdata, complex *cdata, long nt, long nx, long ny, long ldt, long ldx, long ldy, long xorig, long yorig)
{
    long     i, j, k, ld1, sign;
    complex *cdum;
    REAL *rdum;
#if defined(MKL)
    DFTI_DESCRIPTOR_HANDLE my_desc_handle = NULL;
    MKL_LONG status;
    MKL_LONG dim_sizes[3];
    MKL_LONG strides_out[4];
    MKL_LONG strides_in[4];
#endif

    assert ( (xorig >= 0) && (xorig < nx) );
    assert ( (yorig >= 0) && (yorig < ny) );
    ld1 = (nt+2)/2;

#if defined(MKL)
    //fprintf(stderr,"Using MKL for yxt2wkykx.c \n");
    cdum = (complex *)malloc(ld1*ldx*ldy*sizeof(complex));
    if (cdum == NULL) fprintf(stderr,"yxt2wkykx: memory allocation error\n");

    if (xorig==0 && yorig==0) {
        rdum = rdata;
    }
    else {
        rdum = (float *)malloc(ldt*ldx*ldy*sizeof(float));
        if (rdum == NULL) fprintf(stderr,"yxt2wkykx: memory allocation error\n");
        for(k = yorig; k < ny; k++) {
            for(i = xorig; i < nx; i++) {
                for(j = 0; j < nt; j++) {
                    rdum[(k-yorig)*ldt*ldx+(i-xorig)*ldt+j] = rdata[k*nx*ldt+i*ldt+j];
                }
            }
            for(i = 0; i < xorig; i++) {
                for(j = 0; j < nt; j++) {
                    rdum[(k-yorig)*ldt*ldx+(nx-xorig+i)*ldt+j] = rdata[k*nx*ldt+i*ldt+j];
                }
            }
        }
        for(k = 0; k < yorig; k++) {
            for(i = xorig; i < nx; i++) {
                for(j = 0; j < nt; j++) {
                    rdum[(ny-yorig+k)*ldt*ldx+(i-xorig)*ldt+j] = rdata[k*nx*ldt+i*ldt+j];
                }
            }
            for(i = 0; i < xorig; i++) {
                for(j = 0; j < nt; j++) {
                    rdum[(ny-yorig+k)*ldt*ldx+(nx-xorig+i)*ldt+j] = rdata[k*nx*ldt+i*ldt+j];
                }
            }
        }
    }
    
    dim_sizes[0]=ny;  dim_sizes[1]=nx;       dim_sizes[2]=nt;
    strides_in[0]=0;  strides_in[1]=ldx*ldt; strides_in[2]=ldt;  strides_in[3]=1;
    strides_out[0]=0; strides_out[1]=ldx*ld1; strides_out[2]=ld1; strides_out[3]=1;

    status = DftiCreateDescriptor(&my_desc_handle, DFTI_SINGLE, DFTI_REAL, 3, dim_sizes);
    status = DftiSetValue(my_desc_handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    status = DftiSetValue(my_desc_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status = DftiSetValue(my_desc_handle, DFTI_OUTPUT_STRIDES, strides_out);
    status = DftiSetValue(my_desc_handle, DFTI_INPUT_STRIDES, strides_in);
    status = DftiCommitDescriptor(my_desc_handle);

    status = DftiComputeForward(my_desc_handle, rdum, (MKL_Complex8 *)&cdum[0]);

    /* tranpose to [nf][nky][nkx] */
    /* the signs of the FFT kernels in MKL will be different than in yxt2wkykx(time forward -1 x,y 1) */
    //fprintf(stderr,"*** MKL after forward transform tranposed to [nf][nky][nkx] ***\n");
    for (j=0; j<ld1; j++) {
        for (k=0; k<ny; k++) {
            for (i=0; i<nx; i++) {
                cdata[j*ldx*ldy+k*ldx+i] = cdum[k*ldx*ld1+i*ld1+j];
            }
        }
    }
    free(cdum);
    status = DftiFreeDescriptor(&my_desc_handle);
    if (xorig!=0 || yorig!=0) free(rdum);

#else
    cdum = (complex *)malloc(ld1*nx*ny*sizeof(complex));
    if (cdum == NULL) fprintf(stderr,"yxt2wkykx: memory allocation error\n");

    sign = -1;
    rcmfft(rdata, cdum, nt, nx*ny, ldt, ld1, sign);

    // cdata[nky][nf][nkx] = cdum[ny][nx][nf]
    for(j = 0; j < ld1; j++) {
        for(k = yorig; k < ny; k++) {
            for(i = xorig; i < nx; i++) {
                cdata[(k-yorig)*ld1*ldx+j*ldx+(i-xorig)] = cdum[k*nx*ld1+i*ld1+j];
            }
            for(i = 0; i < xorig; i++) {
                cdata[(k-yorig)*ld1*ldx+j*ldx+(nx-xorig+i)] = cdum[k*nx*ld1+i*ld1+j];
            }
        }
        for(k = 0; k < yorig; k++) {
            for(i = xorig; i < nx; i++) {
                cdata[(ny-yorig+k)*ld1*ldx+j*ldx+(i-xorig)] = cdum[k*nx*ld1+i*ld1+j];
            }
            for(i = 0; i < xorig; i++) {
                cdata[(ny-yorig+k)*ld1*ldx+j*ldx+(nx-xorig+i)] = cdum[k*nx*ld1+i*ld1+j];
            }
        }
    }
    
    free(cdum);
    cdum = (complex *)malloc(ld1*ldx*ldy*sizeof(complex));

    sign = 1;
    // cdum[nkx][nf][nky] = cdata[nky][nf][nkx]
    for (k = 0; k < ldy; k++) {
        ccmfft(&cdata[k*ld1*ldx], nx, ld1, ldx, sign);
        for (j = 0; j < ld1; j++) {
            for (i = 0; i < ldx; i++) {
                cdum[i*ld1*ldy+j*ldy+k] = cdata[k*ld1*ldx+j*ldx+i];
            }
        }
    }

    // cdata [nf][nky][nkx] = cdum[nkx][nf][nky]
    for (i = 0; i < ldx; i++) {
        ccmfft(&cdum[i*ld1*ldy], ny, ld1, ldy, sign);
        for (j = 0; j < ld1; j++) {
            for (k = 0; k < ldy; k++) {
                cdata[j*ldy*ldx+k*ldx+i] = cdum[i*ld1*ldy+j*ldy+k];
            }
        }
    }
    free(cdum);
#endif


    return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define nynxt2wkykx FNAME(YXT2WKYKXF)
#else
#define nynxt2wkykx FNAME(yxt2wkykxf)
#endif

void nynxt2wkykx(REAL *rdata, complex *cdata, int *nt, int *nx, int *ny, int *ldt, int *ldx, int *ldy, int *xorig, int *yorig)
{

    yxt2wkykx(rdata, cdata, *nt, *nx, *ny, *ldt, *ldx, *ldy, *xorig, *yorig);

    return;
}
