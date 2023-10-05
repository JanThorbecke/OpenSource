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
#ifdef MPI
#include <mpi.h>
#endif

/**
*  File handling routines 
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

void name_ext(char *filename, char *extension);

#ifdef MPI
MPI_File fileOpen(char *file, char *ext, int append)
{
    MPI_File fh;
    MPI_Offset disp, offset, file_size;
    MPI_Datatype etype, ftype, buftype;
    MPI_Info info;
    MPI_Status status;
    char filename[1024];
    int result;

//    int pe;
//    MPI_Comm_rank( MPI_COMM_WORLD, &pe );
//    fprintf(stderr,"PE %d open file %s\n", pe, file);

    strcpy(filename, file);
    name_ext(filename, ext);
    result = MPI_File_open(MPI_COMM_WORLD, filename,
       MPI_MODE_APPEND | MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    assert(result == MPI_SUCCESS);

    return fh;
}

int traceWrite(segy *hdr, float *data, int n, long long offset, MPI_File fh) 
{
    MPI_Offset off;
    MPI_Status status;
    int result;

    off = (MPI_Offset) offset;
    MPI_File_write_at(fh, off, hdr, TRCBYTES, MPI_CHAR, &status);
    assert(result == MPI_SUCCESS);

    off = off + TRCBYTES;
    MPI_File_write_at(fh, off, data, n, MPI_FLOAT, &status);
    assert(result == MPI_SUCCESS);

    return 0;
}

void fileClose(MPI_File fh)
{
    int result;

    result = MPI_File_close(&fh);
    assert(result == MPI_SUCCESS);

    return;
}
#else
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

int traceWrite(segy *hdr, float *data, int n, long long offset, FILE *fp) 
{
    size_t  nwrite;

    nwrite = fwrite( hdr, 1, TRCBYTES, fp);
    assert(nwrite == TRCBYTES);
    nwrite = fwrite( data, sizeof(float), n, fp);
    assert(nwrite == n);
    fflush(fp);

    return 0;
}

void fileClose(FILE *fp)
{
    fclose(fp);
}
#endif

