#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <genfft.h>

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

long getFileInfo3D(char *filename, long *n1, long *n2, long *n3, long *ngath, float *d1, float *d2, float *d3,
    float *f1, float *f2, float *f3, float *sclsxgxsygy, long *nxm);
double wallclock_time(void);
long writeData3D(FILE *fp, float *data, segy *hdrs, long n1, long n2);
long readSnapData3D(char *filename, float *data, segy *hdrs, long nsnaps, 
    long nx, long ny, long nz, long sx, long ex, long sy, long ey, long sz, long ez);

char *sdoc[] = {
" ",
" combine_induced - Combine induced seismicity results together ",
" ",
" authors  : Joeri Brackenhoff	: (J.A.Brackenhoff@tudelft.nl)",
"		   : Jan Thorbecke 		: (janth@xs4all.nl)",
" ",
" Required parameters: ",
"",
"   file_in= ................. File containing the first data",
" ",
" Optional parameters: ",
" ",
"   file_out= ................ Filename of the output",
"   numb= .................... integer number of first file",
"   dnumb= ................... integer number of increment in files",
"	nzmax= ................... Maximum number of files read",
NULL};

int main (int argc, char **argv)
{
	FILE    *fp_in, *fp_out;
	char    *fin, *fout, *ptr, fbegin[100], fend[100], fins[100], fin2[100], numb1[100];
	float   *indata, *outdata, fz, fy, fx, shift, dtshift, dt_time;
	float   dt, dy, dx, t0, x0, y0, sclsxgx, dt2, dy2, dx2, t02, x02, y02, sclsxgx2, dxrcv, dyrcv, dzrcv;
	long    nshots, nt, ny, nx, ntraces, nshots2, nt2, ny2, nx2, ntraces2, ix, iy, it, is, iz, pos, file_det, nxs, nys, nzs;
	long    numb, dnumb, ret, nzmax, verbose, nshot_out, ishift, nshift;
	segy    *hdr_in, *hdr_bin, *hdr_out;

	initargs(argc, argv);
	requestdoc(1);

    /*----------------------------------------------------------------------------*
    *   Get the parameters passed to the function 
    *----------------------------------------------------------------------------*/
	if (!getparstring("file_in", &fin)) fin = NULL;
    if (!getparstring("file_out", &fout)) fout = "out.su";
	if (!getparlong("numb", &numb)) numb=0;
	if (!getparlong("dnumb", &dnumb)) dnumb=0;
	if (!getparlong("nzmax", &nzmax)) nzmax=0;
	if (!getparlong("verbose", &verbose)) verbose=0;
	if (!getparlong("nshift", &nshift)) nshift=0;
	if (!getparfloat("shift", &shift)) shift=0.0;
	if (!getparfloat("dtshift", &dtshift)) dtshift=0.0;
	if (!getparfloat("dt_time", &dt_time)) dt_time=0.004;
	if (fin == NULL) verr("Incorrect downgoing input");

    /*----------------------------------------------------------------------------*
    *   Determine the position of the number in the string
    *   and split the file into beginning, middle and end
    *----------------------------------------------------------------------------*/
	if (dnumb < 1) dnumb = 1;
	sprintf(numb1,"%li",numb);
	ptr  = strstr(fin,numb1);
    pos = ptr - fin + 1;
    sprintf(fbegin,"%*.*s", pos-1, pos-1, fin);
   	sprintf(fend,"%s", fin+pos);

    /*----------------------------------------------------------------------------*
    *   Determine the amount of files that are present
    *----------------------------------------------------------------------------*/
	file_det = 1;
	nzs=0;
	while (file_det) { // Check for a file with the filename
        sprintf(fins,"%li",nzs*dnumb+numb);
        sprintf(fin,"%s%s%s",fbegin,fins,fend);
        fp_in = fopen(fin, "r");
        if (fp_in == NULL) { // If the filename does not exist
            if (nzs == 0) { // The filename is wrong to begin with
                verr("error on opening basefile=%s", fin);
            }
            else if (nzs == 1) { // There is only a single file
                vmess("1 file detected");
            }
            else { // Stop after the final file has been detected
                vmess("%li files detected",nzs);
                file_det = 0;
                break;
            }
        }
        fclose(fp_in);
        nzs++;
		if (nzmax!=0 && nzs == nzmax) { // Stop if the amount of files exceed the indicated maximum
			vmess("%li files detected",nzs);
            file_det = 0;
            break;
		}
    }

    /*----------------------------------------------------------------------------*
    *   Read in the first two files and determine the header values
    *   of the output
    *----------------------------------------------------------------------------*/
	sprintf(fins,"%li",numb);
    sprintf(fin2,"%s%s%s",fbegin,fins,fend);
	nshots = 0;
    getFileInfo3D(fin2, &nt, &nx, &ny, &nshots, &dt, &dx, &dy, &t0, &x0, &y0, &sclsxgx, &ntraces);

	sprintf(fins,"%li",numb+dnumb);
    sprintf(fin2,"%s%s%s",fbegin,fins,fend);
    nshots = 0;
    getFileInfo3D(fin2, &nt2, &nx2, &ny2, &nshots2, &dt2, &dx2, &dy2, &t02, &x02, &y02, &sclsxgx2, &ntraces2);

	dxrcv=dx*1000;
    dyrcv=dy*1000;
	dzrcv=t02-t0;

	if (nshots==0) nshots=1;
	nxs = ntraces;

	/*----------------------------------------------------------------------------*
    *   Read in a single file to determine if the header values match
    *   and allocate the data
    *----------------------------------------------------------------------------*/
	hdr_bin      = (segy *)calloc(nxs,sizeof(segy));
    indata    	= (float *)calloc(nxs*nt,sizeof(float));

	readSnapData3D(fin2, &indata[0], &hdr_bin[0], nshots, nxs, ny, nt, 0, nxs, 0, ny, 0, nt);
	nshots 	= hdr_bin[nxs-1].fldr;
	nxs		= hdr_bin[nxs-1].tracf;

	nshot_out = nshots/2;
	free(indata);

	hdr_out     = (segy *)calloc(nshot_out*nxs,sizeof(segy));	
	outdata		= (float *)calloc(nshot_out*nxs*nt,sizeof(float));

    /*----------------------------------------------------------------------------*
    *   Write out the file info in case of verbose
    *----------------------------------------------------------------------------*/
    if (verbose) vmess("Number of virtual receivers: %li, nx=%li, ny=%li, nz=%li",nx*ny*nt,nx,ny,nt);

    /*----------------------------------------------------------------------------*
    *   Parallel loop for reading in and combining the various shots
    *----------------------------------------------------------------------------*/
#pragma omp parallel default(shared) \
  private(indata,iz,hdr_in,fins,fin2,fp_in,is,ix,iy,it,ishift)
{
	indata     = (float *)calloc(ntraces*nt,sizeof(float));
	hdr_in      = (segy *)calloc(ntraces,sizeof(segy));

#pragma omp for
	for (iz = 0; iz < nzs; iz++) {
		if (verbose) vmess("Depth:%li out of %li",iz+1,nzs);
		sprintf(fins,"%li",iz*dnumb+numb);
       	sprintf(fin2,"%s%s%s",fbegin,fins,fend);
       	fp_in = fopen(fin2, "r");
		if (fp_in == NULL) {
			verr("Error opening file");
		}
		fclose(fp_in);
		readSnapData3D(fin2, &indata[0], &hdr_in[0], nshots, nxs, ny, nt, 0, nxs, 0, ny, 0, nt);
		if (iz==0) fz=hdr_in[0].f1; fx=hdr_in[0].f2;
		if (iz==1) dzrcv=hdr_in[0].f1-fz;
		
		ishift = nshift*iz;
		if (verbose) vmess("Shifting %li timesteps for a total of %.3f seconds",ishift,shift+(dtshift*((float)iz)));
		for (is = ishift; is < nshot_out; is++) {
			for (ix = 0; ix < nxs; ix++) {
				for (it = 0; it < nt; it++) {
					outdata[is*nxs*nt+ix*nt+it] += indata[(is-ishift+(nshots/2))*nxs*nt+ix*nt+it];
				}
			}
		}
	}
	free(indata);free(hdr_in);
}

    /*----------------------------------------------------------------------------*
    *   Write out the data to file
    *----------------------------------------------------------------------------*/
	fp_out = fopen(fout, "w+");

	for (is = 0; is < nshot_out; is++) {
		for (ix = 0; ix < nxs; ix++) {
           	hdr_out[ix].fldr	= is+1;
           	hdr_out[ix].tracl	= is*nxs+ix+1;
           	hdr_out[ix].tracf	= ix+1;
			hdr_out[ix].scalco  = -1000;
   			hdr_out[ix].scalel	= -1000;
			hdr_out[ix].sdepth	= hdr_bin[0].sdepth;
			hdr_out[ix].trid	= 1;
			hdr_out[ix].ns		= nt;
			hdr_out[ix].trwf	= nxs;
			hdr_out[ix].ntr		= hdr_out[ix].fldr*hdr_out[ix].trwf;
			hdr_out[ix].f1		= fz;
			hdr_out[ix].f2		= fx;
			hdr_out[ix].dt      = dt_time*(1E6);
			hdr_out[ix].d1      = dzrcv;
           	hdr_out[ix].d2      = dxrcv;
			hdr_out[ix].sx      = (int)roundf(fx + (ix*hdr_out[ix].d2));
			hdr_out[ix].gx      = (int)roundf(fx + (ix*hdr_out[ix].d2));
           	hdr_out[ix].offset	= (hdr_out[ix].gx - hdr_out[ix].sx)/1000.0;
		}
		ret = writeData3D(fp_out, &outdata[is*nxs*nt], hdr_out, nt, nxs);
		if (ret < 0 ) verr("error on writing output file.");
	}
	
	fclose(fp_out);
	return 0;
}

