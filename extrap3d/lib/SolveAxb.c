#include <math.h>
#include <stdlib.h>
#include <assert.h>

#define SGELS sgels_

/***************************************************************
*
* Solving Ax = b for real numbers by QR decomposition using LAPACK 
* routines.
*
* Jan Thorbecke June-1994
*
*****************************************************************/

void SolveAxb(float *a, int nrow, int ncol, float *b, float *x)
{
	int 	i, lda, ldb, nrhs, lwork, info;
	float	*matrix, *work, *bl;
	char    *trans;
		
	lwork   = nrow+ncol;
	matrix  = (float *)malloc(nrow*ncol*sizeof(float));
	work  	= (float *)malloc(3*ncol*sizeof(float));
	bl     	= (float *)malloc(2*ncol*sizeof(float));

	for(i = 0; i < ncol*nrow; i++) matrix[i] = a[i];
	lda = nrow;

	for(i = 0; i < nrow; i++) bl[i] = b[i];
	ldb = nrow;

	trans = "N";
	nrhs  = 1;
	SGELS(trans,&nrow,&ncol,&nrhs,matrix,&lda,bl,&ldb,work,&lwork,&info);
	assert (info == 0);

	for(i = 0; i < ncol; i++) x[i] = bl[i];

	free(matrix);
	free(work);
	free(bl);

	return;
}
