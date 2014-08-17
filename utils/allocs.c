#include<stdlib.h>

/**
* Allocation routines for 2 Dimensional flaoting point arrays 
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/


/* allocate a 2-d array */
void **alloc2 (size_t n1, size_t n2, size_t size)
{
    size_t i2;
    void **p;

    if ((p=(void**)malloc(n2*sizeof(void*)))==NULL)
        return NULL;
    if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
        free(p);
        return NULL;
    }
    for (i2=0; i2<n2; i2++)
        p[i2] = (char*)p[0]+size*n1*i2;
    return p;
}

/* free a 2-d array */
void free2 (void **p)
{
    free(p[0]);
    free(p);
}


/* allocate a 2-d array of floats */
float **alloc2float(size_t n1, size_t n2)
{
    return (float**)alloc2(n1,n2,sizeof(float));
}

/* free a 2-d array of floats */
void free2float(float **p)
{
    free2((void**)p);
}

