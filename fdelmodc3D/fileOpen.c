#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include "segy.h"

/**
*  File handling routines 
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

void name_ext(char *filename, char *extension);

FILE *fileOpen(char *file, char *ext, int append)
{
	FILE *fp;
	char filename[1024];

	strcpy(filename, file);
	name_ext(filename, ext);
	if (append) fp = fopen(filename, "a");
   	else fp = fopen(filename, "w");
   	assert(fp != NULL);
	
	return fp;
}

int traceWrite(segy *hdr, float *data, int n, FILE *fp) 
{
    size_t  nwrite;

    nwrite = fwrite( hdr, 1, TRCBYTES, fp);
    assert(nwrite == TRCBYTES);
    nwrite = fwrite( data, sizeof(float), n, fp);
    assert(nwrite == n);

	return 0;
}

