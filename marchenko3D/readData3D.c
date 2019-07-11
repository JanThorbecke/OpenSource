#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "segy.h"

/**
* reads SU file and returns header and 2D array
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/


long readData3D(FILE *fp, float *data, segy *hdrs, long n1)
{
	size_t nread;
	long oneshot, itrace, sx, sy, fldr, gy, nx, ny;
	segy hdr;

	oneshot = 1;
	itrace  = 0;
    ny      = 1;

	nread = fread(&hdrs[0], 1, TRCBYTES, fp);
	if (nread == 0) return 0; /* end of file reached */

	if (n1==0) n1 = hdrs[0].ns;
	sx   = hdrs[0].sx;
    sy   = hdrs[0].sy;
	fldr = hdrs[0].fldr;
    gy   = hdrs[0].gy;
	while (oneshot) {
		nread = fread(&data[itrace*n1], sizeof(float), n1, fp);
		assert (nread == n1);
		itrace++;
		nread = fread(&hdr, 1, TRCBYTES, fp);
		if (nread == 0) break;
		assert(nread == TRCBYTES);
		if ( (sx != hdr.sx) || (sy != hdr.sy) || (fldr != hdr.fldr)) { /* end of shot record */
			fseek( fp, -TRCBYTES, SEEK_CUR );
			break;
		}
		if ( (gy != hdr.gy)) {
            ny++;
            gy = hdr.gy;
        }
		memcpy(&hdrs[itrace], &hdr, TRCBYTES);
	}
    nx = itrace/ny;

	return nx, ny;
}
