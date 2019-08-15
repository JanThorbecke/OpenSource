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

void name_ext(char *filename, char *extension);

long readSnapData3D(char *filename, float *data, segy *hdrs, long nsnaps, long nx, long ny, long nz,
    long sx, long ex, long sy, long ey, long sz, long ez);

char *sdoc[] = {
" ",
" snap2shot - Reshape snapshot data to receiver data",
" ",
" authors  : Joeri Brackenhoff  (J.A.Brackenhoff@tudelft.nl)",
"		   : Jan Thorbecke      (janth@xs4all.nl)",
" ",
" Required parameters: ",
"",
"   file_snap= ............... First filename that needs to be reshaped to receiver data",
" ",
" Optional parameters: ",
" ",
"   file_rcv= ................ File containing the snapshot data that will determine the mute window",
"   numb=0 ................... Starting number in the first file",
"   dnumb=1 .................. Increment per file name",
"   nxmax=0 .................. Maximum number of files that are to be reshaped",
"   numb=0 ................... Starting number in the first file",
"   nzstart=0 ................ Starting depth number to be written out",
"   nzend=end ................ Final depth number to be written out",
NULL};

int main (int argc, char **argv)
{
	FILE *fp_snap, *fp_rcv;
	char    *file_snap, *file_rcv, file_tmp[150], ins[100], fbegin[100], fend[100], fins[100], fin2[100], numb1[100], *ptr;
	float   *rcvdata, *snapdata;
    float   dxs, dys, dzs, fxs, fys, fzs, sclr;
    long    nxs, nys, nzs, nts, ntrs, ret, file_det;
	long	it, ix, iy, iz, ixr, nxr, dnumb, numb, pos, nxmax;
	long 	nzstart, nzend, dnz;
	int 	*sx, *sy, pf;
	segy    *hdr_snap, *hdr_rcv;

	initargs(argc, argv);
	requestdoc(1);

    /*----------------------------------------------------------------------------*
    *   Get the parameters passed to the function 
    *----------------------------------------------------------------------------*/
	if (!getparstring("file_rcv", &file_rcv)) file_rcv = "rcv.su";
	if (!getparstring("file_snap", &file_snap)) file_snap = NULL;
	if (file_snap == NULL) verr("No input data for the snap file given");
	if (!getparlong("numb", &numb)) numb=0;
	if (!getparlong("dnumb", &dnumb)) dnumb=1;
	if (!getparlong("nxmax", &nxmax)) nxmax=0;
	if (!getparlong("nzstart", &nzstart)) nzstart=0;
	if (!getparlong("nzend", &nzend)) nzend=0;

    /*----------------------------------------------------------------------------*
    *   Split the filename so the number can be changed
    *----------------------------------------------------------------------------*/
	if (dnumb < 1) dnumb = 1;
	sprintf(numb1,"%li",numb);
	ptr  = strstr(file_snap,numb1);
    pos = ptr - file_snap + 1;
	pf = pos-1;
    sprintf(fbegin,"%*.*s", pf, pf, file_snap);
   	sprintf(fend,"%s", file_snap+pos);

    /*----------------------------------------------------------------------------*
    *   Determine the amount of files to be read
    *----------------------------------------------------------------------------*/
	file_det = 1;
	nxr=0;
	while (file_det) {
        sprintf(fins,"%li",nxr*dnumb+numb);
        sprintf(fin2,"%s%s%s",fbegin,fins,fend);
        fp_snap = fopen(fin2, "r");
        if (fp_snap == NULL) {
            if (nxr == 0) {
                verr("error on opening basefile=%s", fin2);
            }
            else if (nxr == 1) {
                vmess("1 file detected");
				file_det = 0;
                break;
            }
            else {
                vmess("%d files detected",nxr);
                file_det = 0;
                break;
            }
        }
        fclose(fp_snap);
        nxr++;
		if (nxr!=0 && nxr == nxmax) {
			vmess("%li files detected",nxr);
            file_det = 0;
            break;
		}
    }

    /*----------------------------------------------------------------------------*
    *   Get the info from the files
    *----------------------------------------------------------------------------*/
	sprintf(fins,"%li",numb);
    sprintf(fin2,"%s%s%s",fbegin,fins,fend);

	getFileInfo3D(fin2, &nzs, &nxs, &nys, &nts, &dzs, &dxs, &dys, &fzs, &fxs, &fys, &sclr, &ntrs);

	if (nzend > nzs || nzend==0) nzend=nzs;
	if (nzstart < 0) nzstart=0;
	if (nzstart > nzend) verr("The value of nzstart (%li) is greater than nzend (%li)",nzstart,nzend);

	dnz = nzend-nzstart;

    /*----------------------------------------------------------------------------*
    *   Allocate the data
    *----------------------------------------------------------------------------*/
	rcvdata		= (float *)malloc(nxr*dnz*nxs*nys*nts*sizeof(float));
	snapdata    = (float *)malloc(dnz*nxs*nys*nts*sizeof(float));
	hdr_snap    = (segy *)calloc(nxs*nys*nts,sizeof(segy));
	sx 			= (int *)malloc(nxr*sizeof(int));
	sy 			= (int *)malloc(nxr*sizeof(int));

    /*----------------------------------------------------------------------------*
    *   Reshape the data
    *----------------------------------------------------------------------------*/
	for (ixr=0; ixr<nxr; ixr++) {
		vmess("Reshaping %li out of %li files",ixr+1,nxr);
		sprintf(fins,"%li",ixr*dnumb+numb);
    	sprintf(fin2,"%s%s%s",fbegin,fins,fend);
		readSnapData3D(fin2, snapdata, hdr_snap, nts, nxs, nys, nzs, 0, nxs, 0, nys, nzstart, nzend);
		sx[ixr] = hdr_snap[0].sx;
		sy[ixr] = hdr_snap[0].sy;
		if (ixr==0) dzs = hdr_snap[0].d1;
		for (iz=0; iz<dnz; iz++) {
			for (iy=0; iy<nys; iy++) {
				for (ix=0; ix<nxs; ix++) {
					for (it=0; it<nts; it++) {
						rcvdata[iz*nys*nxs*nxr*nts+iy*nxs*nxr*nts+ix*nxr*nts+ixr*nts+it] = snapdata[it*nys*nxs*dnz+iy*nxs*dnz+ix*dnz+iz];
					}
				}
			}
		}
	}
	free(snapdata);

    /*----------------------------------------------------------------------------*
    *   Write out the data to new files
    *----------------------------------------------------------------------------*/
	hdr_rcv = (segy *)calloc(nxs*nys*nxr,sizeof(segy));
	for (iz=nzstart; iz<nzend; iz++) {
		vmess("Writing depth %li out of %li",iz-nzstart+1,dnz);
		strcpy(file_tmp, file_rcv);
		sprintf(ins,"_z%li", iz);
		name_ext(file_tmp, ins);
		fp_rcv = fopen(file_tmp, "w+");
		if (fp_rcv==NULL) verr("error on creating output file %s", file_rcv);
		// Set headers
		for (iy=0; iy<nys; iy++){
			for (ix=0; ix<nxs; ix++){
				for (ixr=0; ixr<nxr; ixr++) {
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].ns      = nts;
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].fldr    = iy*nxs+ix+1;
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].tracf   = ixr+1;
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].tracl   = iy*nxs*nxr+ix*nxr+ixr+1;
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].trid    = 1;
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].scalco  = -1000;
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].scalel  = -1000;
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].sx      = (int)((fxs+dxs*ix)*(1e3));
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].sy      = (int)((fys+dys*iy)*(1e3));
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].gx      = sx[ixr];
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].gy      = sy[ixr];
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].sdepth  = (int)((fzs+dzs*iz)*(1e3));
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].f1      = 0.0;
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].f2      = sx[0]/(1e3);
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].d1      = (float)(hdr_snap[0].dt/1e6);
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].d2      = dxs;
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].dt      = hdr_snap[0].dt;
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].trwf    = nxs*nys*nxr;
					hdr_rcv[iy*nxs*nxr+ix*nxr+ixr].ntr     = nxs*nys*nxr;
				}
			}
		}
		// Write out homogeneous Green's function
		ret = writeData3D(fp_rcv, (float *)&rcvdata[(iz-nzstart)*nys*nxs*nts*nxr], hdr_rcv, nts, nxr*nxs*nys);
		if (ret < 0 ) verr("error on writing output file.");
		fclose(fp_rcv);
	}

	free(rcvdata); free(hdr_rcv); free(hdr_snap);

	return 0;
}
