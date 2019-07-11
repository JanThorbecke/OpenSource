#include <stdlib.h>
#include"fdelmodc3D.h"

/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* allocate a 3-d array of floats */
void ***alloc3float(modPar mod)
{
    size_t i3,i2,n1,n2,n3,size;
    void ***p;

	n1 = mod.naz;
	n2 = mod.nax;
	n3 = mod.nay;
	size = sizeof(float);

/* allocate pointer arrays of n3 and n2 */
    if ((p=(void***)malloc(n3*sizeof(void**)))==NULL)
        return NULL;
    if ((p[0]=(void**)malloc(n3*n2*sizeof(void*)))==NULL) {
        free(p);
        return NULL;
    }

/* allocate data size and assign pointers to n1, n2 and n3 */
	if (mod.nfx==1 && mod.nfy==1 ) { // 1D model 
		fprintf(stderr,"Allocating 1D model \n");
        if ((p[0][0]=(void*)malloc(n1*size))==NULL) {
            free(p);
            return NULL;
        }
        for (i3=0; i3<n3; i3++) {
            p[i3] = p[0]+n2*i3;
            for (i2=0; i2<n2; i2++)
                p[i3][i2] = (void*)p[0][0];
        }
	}
	else if (mod.nfy==1) { // 2D model
		fprintf(stderr,"Allocating 2D model ny=%d with nx=%d and nz=%d\n",n3,n2,n1);
        if ((p[0][0]=(void*)malloc(n2*n1*size))==NULL) {
            free(p);
            return NULL;
        }
        for (i3=0; i3<n3; i3++) {
            p[i3] = p[0]+n2*i3;
            for (i2=0; i2<n2; i2++)
                p[i3][i2] = (char*)p[0][0]+size*n1*i2;
        }
	}
	else { // 3D model 
		//fprintf(stderr,"Allocating 3D model \n");
        if ((p[0][0]=(void*)malloc(n3*n2*n1*size))==NULL) {
            free(p[0]);
            free(p);
            return NULL;
        }
        for (i3=0; i3<n3; i3++) {
            p[i3] = p[0]+n2*i3;
            for (i2=0; i2<n2; i2++)
                p[i3][i2] = (char*)p[0][0]+size*n1*(i2+n2*i3);
        }
	}
    return p;
}

void free3float(void ***p)
{
    free(p[0][0]);
    free(p[0]);
    free(p);
}

