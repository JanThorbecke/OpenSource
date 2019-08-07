#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "par.h"
#include "segy.h"

/**
*  Writes an 2D array to a SU file
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

#define TRCBYTES 240

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define ISODD(n) ((n) & 01)
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int writesufile(char *filename, float *data, size_t n1, size_t n2, float f1, float f2, float d1, float d2)
{
	FILE    *file_out;
	size_t  nwrite, itrace;
	int     ns;
	segy    *hdr;

/* Read in parameters */

	
	if (n1 > USHRT_MAX) {
		vwarn("Output file %s: number of samples is truncated from %d to USHRT_MAX.", filename, n1);
	}
	ns = MIN(n1,USHRT_MAX);

	file_out = fopen( filename, "w+" );
	assert( file_out );

	hdr = (segy *)calloc(1,TRCBYTES);
	hdr->ns = ns;
	hdr->dt = NINT(1000000*(d1));
	hdr->d1 = d1;
	hdr->d2 = d2;
	hdr->f1 = f1;
	hdr->f2 = f2;
	hdr->fldr = 1;
	hdr->trwf = n2;

	for (itrace=0; itrace<n2; itrace++) {
		hdr->tracl = itrace+1;
		nwrite = fwrite( hdr, 1, TRCBYTES, file_out );
		assert (nwrite == TRCBYTES);
		nwrite = fwrite( &data[itrace*n1], sizeof(float), ns, file_out );
		assert (nwrite == ns);
	} 
	fclose(file_out);
	free(hdr);

	return 0;
}

