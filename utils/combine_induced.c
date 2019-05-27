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
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);
int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);
double wallclock_time(void);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);
int readSnapData(char *filename, float *data, segy *hdr, int ngath, int nx, int ntfft, int sx, int ex, int sz, int ez);

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
	FILE *fp_in, *fp_out;
	char *fin, *fout, *ptr, fbegin[100], fend[100], fins[100], fin2[100], numb1[100];
	float *indata, *outdata, *rtrace, fz, fx, shift, dtshift, dt_time;
	float dt, dx, t0, x0, xmin, xmax, sclsxgx, dt2, dx2, t02, x02, xmin2, xmax2, sclsxgx2, dxrcv, dzrcv;
	int nshots, nt, nx, ntraces, nshots2, nt2, nx2, ntraces2, ix, it, is, iz, pos, ifile, file_det, nxs, nzs;
	int xcount, numb, dnumb, ret, nzmax, verbose, nshot_out, ishift, nshift;
	segy *hdr_in, *hdr_bin, *hdr_out;

	initargs(argc, argv);
	requestdoc(1);

	if (!getparstring("file_in", &fin)) fin = NULL;
    if (!getparstring("file_out", &fout)) fout = "out.su";
	if (!getparint("numb", &numb)) numb=0;
	if (!getparint("dnumb", &dnumb)) dnumb=0;
	if (!getparint("nzmax", &nzmax)) nzmax=0;
	if (!getparint("verbose", &verbose)) verbose=0;
	if (!getparint("nshift", &nshift)) nshift=0;
	if (!getparfloat("shift", &shift)) shift=0.0;
	if (!getparfloat("dtshift", &dtshift)) dtshift=0.0;
	if (!getparfloat("dt_time", &dt_time)) dt_time=0.004;
	if (fin == NULL) verr("Incorrect downgoing input");

	if (dnumb < 1) dnumb = 1;

	sprintf(numb1,"%d",numb);

	ptr  = strstr(fin,numb1);
    pos = ptr - fin + 1;

    sprintf(fbegin,"%*.*s", pos-1, pos-1, fin);
   	sprintf(fend,"%s", fin+pos);

	file_det = 1;
	nzs=0;

	while (file_det) {
        sprintf(fins,"%d",nzs*dnumb+numb);
        sprintf(fin,"%s%s%s",fbegin,fins,fend);
        fp_in = fopen(fin, "r");
        if (fp_in == NULL) {
            if (nzs == 0) {
                verr("error on opening basefile=%s", fin);
            }
            else if (nzs == 1) {
                vmess("1 file detected");
            }
            else {
                vmess("%d files detected",nzs);
                file_det = 0;
                break;
            }
        }
        fclose(fp_in);
        nzs++;
		if (nzmax!=0 && nzs == nzmax) {
			vmess("%d files detected",nzs);
            file_det = 0;
            break;
		}
    }

	sprintf(fins,"%d",numb);
    sprintf(fin2,"%s%s%s",fbegin,fins,fend);
	nshots = 0;
    getFileInfo(fin2, &nt, &nx, &nshots, &dt, &dx, &t0, &x0, &xmin, &xmax, &sclsxgx, &ntraces);

	sprintf(fins,"%d",numb+dnumb);
    sprintf(fin2,"%s%s%s",fbegin,fins,fend);
    nshots = 0;
    getFileInfo(fin2, &nt2, &nx2, &nshots2, &dt2, &dx2, &t02, &x02, &xmin2, &xmax2, &sclsxgx2, &ntraces2);

	dxrcv=dx*1000;
	dzrcv=t02-t0;

	if (nshots==0) nshots=1;
	nxs = ntraces;

	// ngath zijn het aantal schoten
	hdr_bin      = (segy *)calloc(nxs,sizeof(segy));
    indata    	= (float *)calloc(nxs*nt,sizeof(float));

	readSnapData(fin2, &indata[0], &hdr_bin[0], nshots, nxs, nt, 0, nxs, 0, nt);
	nshots 	= hdr_bin[nxs-1].fldr;
	nxs		= hdr_bin[nxs-1].tracf;

	nshot_out = nshots/2;
	free(indata);

	hdr_out     = (segy *)calloc(nshot_out*nxs,sizeof(segy));	
	outdata		= (float *)calloc(nshot_out*nxs*nt,sizeof(float));

	vmess("nshot_out:%d nxs=%d nt:%d",nshot_out,nxs,nt);

#pragma omp parallel default(shared) \
  private(indata,iz,hdr_in,fins,fin2,fp_in,is,ix,it,ishift)
{
	indata     = (float *)calloc(ntraces*nt,sizeof(float));
	hdr_in      = (segy *)calloc(ntraces,sizeof(segy));

#pragma omp for
	for (iz = 0; iz < nzs; iz++) {
		if (verbose) vmess("Depth:%d out of %d",iz+1,nzs);
		sprintf(fins,"%d",iz*dnumb+numb);
       	sprintf(fin2,"%s%s%s",fbegin,fins,fend);
       	fp_in = fopen(fin2, "r");
		if (fp_in == NULL) {
			verr("Error opening file");
		}
		fclose(fp_in);
		readSnapData(fin2, &indata[0], &hdr_in[0], nshots, nxs, nt, 0, nxs, 0, nt);
		if (iz==0) fz=hdr_in[0].f1; fx=hdr_in[0].f2;
		if (iz==1) dzrcv=hdr_in[0].f1-fz;
		//ishift = (int)((shift+(dtshift*((float)iz)))/dt_time);
		ishift = nshift*iz;
		if (verbose) vmess("Shifting %d timesteps for a total of %.3f seconds",ishift,shift+(dtshift*((float)iz)));
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
		ret = writeData(fp_out, &outdata[is*nxs*nt], hdr_out, nt, nxs);
		if (ret < 0 ) verr("error on writing output file.");
	}
	
	fclose(fp_out);
	return 0;
}

