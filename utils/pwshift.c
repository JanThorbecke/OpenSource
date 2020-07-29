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

void timeShift(float *data, long nsam, long nrec, float dt, float *time, float *amp, float *delay, float fmin, float fmax);
void pad_data(float *data, long nsam, long nrec, long nsamout, float *datout);
void scl_data(float *data, long nsam, long nrec, float scl, float *datout, long nsamout);

char *sdoc[] = {
" ",
" combine - Combine results into a single result ",
" ",
" authors  : Joeri Brackenhoff  : (J.A.Brackenhoff@tudelft.nl)",
"          : Jan Thorbecke      : (janth@xs4all.nl)",
" ",
" Required parameters: ",
"",
"   file_gmin= ................. File containing the G- data",
"   file_ray= .................. File containing the ray time data",
"   file_amp= .................. File containing the ray amplitude data",
" ",
" Optional parameters: ",
" ",
"   file_out=out.su .......... Filename of the output",
"   fmin=0.0 ................. Minimum of the frequency range",
"   fmax=70.0 ................ Maximum of the frequency range",
"   src_velox=1500.0 ......... Velocity of the x-angle of the plane wave",
"   src_veloy=1500.0 ......... Velocity of the y-angle of the plane wave",
"   src_anglex=0.0 ........... x-angle of the plane wave",
"   src_angley=0.0 ........... y-angle of the plane wave",
"   numb=0 ................... Starting number of the level for the image",
"   dnumb=1 .................. Increment number of the level for the image",
"   verbose=1 ................ Give detailed information of process",
NULL};

int main (int argc, char **argv)
{
	FILE	*fp_gmin, *fp_ray, *fp_amp, *fp_out;
	char	*file_gmin, *file_ray, *file_amp, *file_out, fbr[100], fer[100], fba[100], fea[100], fins[100], fin2[100], numb1[100], *ptr;
	float	*gmin, *conv, *time, *amp, *image, fmin, fmax;
	float	dt, dy, dx, dz, t0, y0, x0, scl, *shift, *block;
	float	dt_ray, dy_ray, dx_ray, t0_ray, y0_ray, x0_ray, scl_ray, px, py, src_velox, src_veloy, src_anglex, src_angley, grad2rad;
	long	nshots, nt, ny, nx, ntr, delay;
	long	nray, nt_ray, ny_ray, nx_ray, ntr_ray;
    long    verbose, ix, iy, it, iz, is, *gx, *gy, *gz, numb, dnumb, pos, nzs, file_det;
    size_t  ret;
	segy	*hdr_gmin, *hdr_time, *hdr_amp, *hdr_out;

	initargs(argc, argv);
	requestdoc(1);

    /*----------------------------------------------------------------------------*
    *   Get the parameters passed to the function 
    *----------------------------------------------------------------------------*/
	if (!getparstring("file_gmin", &file_gmin)) file_gmin = NULL;
	if (!getparstring("file_ray", &file_ray)) file_ray = NULL;
	if (!getparstring("file_amp", &file_amp)) file_amp = NULL;
    if (!getparstring("file_out", &file_out)) file_out = "out.su";
	if (!getparlong("verbose", &verbose)) verbose=1;
	if (!getparfloat("fmin", &fmin)) fmin=0.0;
	if (!getparfloat("fmax", &fmax)) fmax=70.0;
	if (!getparfloat("src_velox", &src_velox)) src_velox=1500.0;
	if (!getparfloat("src_veloy", &src_veloy)) src_veloy=1500.0;
	if (!getparfloat("src_anglex", &src_anglex)) src_anglex=0.0;
	if (!getparfloat("src_angley", &src_angley)) src_angley=0.0;
	if (!getparlong("numb", &numb)) numb=0;
	if (!getparlong("dnumb", &dnumb)) dnumb=1;
	if (file_gmin == NULL) verr("Incorrect gmin input");
	if (file_ray == NULL) verr("Incorrect ray time input");
	if (file_amp == NULL) verr("Incorrect ray amplitude input");

	/*----------------------------------------------------------------------------*
    *   Determine the position of the number in the string
    *   and split the file into beginning, middle and end
    *----------------------------------------------------------------------------*/
	if (dnumb < 1) dnumb = 1;
	sprintf(numb1,"z%li",numb);

	ptr  = strstr(file_ray,numb1);
    pos = ptr - file_ray + 1;
    sprintf(fbr,"%*.*s", pos-1, pos-1, file_ray);
   	sprintf(fer,"%s", file_ray+pos+1);

	ptr  = strstr(file_amp,numb1);
    pos = ptr - file_amp + 1;
    sprintf(fba,"%*.*s", pos-1, pos-1, file_ray);
   	sprintf(fea,"%s", file_ray+pos+1);

    /*----------------------------------------------------------------------------*
    *   Determine the amount of files that are present
    *----------------------------------------------------------------------------*/
	file_det = 1;
	nzs=0;
	while (file_det) { // Check for a file with the filename
        sprintf(fins,"z%li",nzs*dnumb+numb);
        sprintf(file_ray,"%s%s%s",fbr,fins,fer);
        fp_ray = fopen(file_ray, "r");
        if (fp_ray == NULL) { // If the filename does not exist
            if (nzs == 0) { // The filename is wrong to begin with
                verr("error on opening basefile=%s", file_ray);
            }
            else if (nzs == 1) { // There is only a single file
                vmess("1 file detected");
                file_det = 0;
                break;
            }
            else { // Stop after the final file has been detected
                vmess("%li files detected",nzs);
                file_det = 0;
                break;
            }
        }
        fclose(fp_ray);
        nzs++;
    }

    /*----------------------------------------------------------------------------*
    *   Read in the first two files and determine the header values
    *   of the output
    *----------------------------------------------------------------------------*/
    getFileInfo3D(file_gmin, &nt, &nx, &ny, &nshots, &dt, &dx, &dy, &t0, &x0, &y0, &scl, &ntr);

	if (verbose) {
        vmess("************************ Gmin info ************************");
		vmess("Number of depth levels : %li",nshots);
		vmess("Number of samples     x: %li,  y: %li,  t: %li",nx,ny,nt);
		vmess("Starting distance for x: %.3f, y: %.3f, t: %.3f",x0,y0,t0);
		vmess("Sampling distance for x: %.3f, y: %.3f, t: %.3f",dx,dy,dt);
        vmess("***********************************************************");
	}
   
    nray = 0;
    sprintf(fins,"z%li",0);
    sprintf(file_ray,"%s%s%s",fbr,fins,fer);
    getFileInfo3D(file_ray, &nx_ray, &nt_ray, &ny_ray, &nray, &dx_ray, &dt_ray, &dy_ray, &x0_ray, &t0_ray, &y0_ray, &scl_ray, &ntr_ray);

	if (verbose) {
        vmess("************************ ray info *************************");
		vmess("Number of depth levels : %li",nzs);
		vmess("Number of focal points : %li",nray);
		vmess("Number of samples     x: %li,  y: %li,  t: %li",nx_ray,ny_ray,nt_ray);
		vmess("Starting distance for x: %.3f, y: %.3f, t: %.3f",x0_ray,y0_ray,t0_ray);
		vmess("Sampling distance for x: %.3f, y: %.3f, t: %.3f",dx_ray,dy_ray,dt_ray);
        vmess("***********************************************************");
	}

	if (nshots!=nzs) verr("The depth levels of the rays (%li) do not match those of the plane waves (%li)",nzs,nshots);

	/*----------------------------------------------------------------------------*
    *   Read in a single file to determine if the header values match
    *   and allocate the data
    *----------------------------------------------------------------------------*/
	hdr_gmin    = (segy *)calloc(nshots*nx*ny,sizeof(segy));	
	gmin        = (float *)calloc(nshots*nx*ny*nt,sizeof(float));
	hdr_out     = (segy *)calloc(nray*nshots,sizeof(segy));	

	image       = (float *)calloc(nray*nshots,sizeof(float));
	gx          = (long *)calloc(nray*nshots,sizeof(long));
	gy       	= (long *)calloc(nray*nshots,sizeof(long));
	gz       	= (long *)calloc(nray*nshots,sizeof(long));

	block       = (float *)calloc(nray*nshots*nt,sizeof(float));

	readSnapData3D(file_gmin, gmin, hdr_gmin, nshots, nx, ny, nt, 0, nx, 0, ny, 0, nt);
	if (verbose) vmess("Read in Gmin data");

	/*----------------------------------------------------------------------------*
    *   Add the delay in case the plane wave is at an angle
    *----------------------------------------------------------------------------*/

    shift = (float *)calloc(nx*ny,sizeof(float));
	grad2rad = 17.453292e-3;
	px = sin(src_anglex*grad2rad)/src_velox;
	py = sin(src_angley*grad2rad)/src_veloy;
	if (verbose) vmess("px value is %f and py value is %f",px,py);

	// if (py < 0.0) {
	// 	for (iy=0; iy<ny; iy++) {
	// 		if (px < 0.0) {
	// 			for (ix=0; ix<nx; ix++) {
	// 				shift[iy*nx+ix] = fabsf((nx-1-ix)*dx*px) + fabsf((ny-1-iy)*dy*py);
	// 			}
	// 		}
	// 		else {
	// 			for (ix=0; ix<nx; ix++) {
	// 				shift[iy*nx+ix] = ix*dx*px + fabsf((ny-1-iy)*dy*py);
	// 			}
	// 		}
	// 	}
	// }
	// else {
	// 	for (iy=0; iy<ny; iy++) {
	// 		if (px < 0.0) {
	// 			for (ix=0; ix<nx; ix++) {
	// 				shift[iy*nx+ix] = fabsf((nx-1-ix)*dx*px) + iy*dy*py;
	// 			}
	// 		}
	// 		else {
	// 			for (ix=0; ix<nx; ix++) {
	// 				shift[iy*nx+ix] = ix*dx*px + iy*dy*py;
	// 			}
	// 		}
	// 	}
	// }

	/*----------------------------------------------------------------------------*
    *   Apply the imaging condition
    *----------------------------------------------------------------------------*/
#pragma omp parallel for schedule(static,1) default(shared) \
  private(is,it,ix,fins,fin2,time,amp,hdr_time,hdr_amp,conv)
    for (iy = 0; iy < nshots; iy++) {

		hdr_time    = (segy *)calloc(ny*nray,sizeof(segy));	
		time        = (float *)calloc(nx*ny*nray,sizeof(float));
		hdr_amp     = (segy *)calloc(ny*nray,sizeof(segy));	
		amp         = (float *)calloc(nx*ny*nray,sizeof(float));
		conv        = (float *)calloc(nx*ny*nt,sizeof(float));

		vmess("Depth level %li out of %li",iy+1,nshots);
		sprintf(fins,"z%li",iy*dnumb+numb);
        sprintf(fin2,"%s%s%s",fbr,fins,fer);
		readSnapData3D(fin2, time, hdr_time, nray, 1, ny, nx, 0, 1, 0, ny, 0, nx);
        sprintf(fin2,"%s%s%s",fba,fins,fea);
		readSnapData3D(file_amp, amp,  hdr_amp,  nray, 1, ny, nx, 0, 1, 0, ny, 0, nx);

		for (is = 0; is < nray; is++) {
			gx[is*nshots+iy] = hdr_time[is*ny].sx;
			gy[is*nshots+iy] = hdr_time[is*ny].sy;
			gz[is*nshots+iy] = hdr_time[is*ny].sdepth;

			x0_ray = ((float)gx[is*nshots+iy])/1000.0;
			y0_ray = ((float)gy[is*nshots+iy])/1000.0;

			if (py < 0.0) {
				if (px < 0.0) {
					delay = NINT(fabsf((x0-x0_ray)*px)/dt) + NINT(fabsf((y0-y0_ray)*py)/dt);
				}
				else {
					delay = NINT((x0_ray-x0)*px/dt) + NINT(fabsf((y0-y0_ray)*py)/dt);
				}
			}
			else {
				if (px < 0.0) {
					delay = NINT(fabsf((x0-x0_ray)*px)/dt) + NINT((y0_ray-y0)*py/dt);
				}
				else {
					delay = NINT((x0_ray-x0)*px/dt) + NINT((y0_ray-y0)*py/dt);
				}
			}

			if (delay > nt-1) delay = nt-1;

			for (it = 0; it < nx*ny*nt; it++) {
				conv[it] = gmin[iy*nx*ny*nt+it];
			}
			timeShift(&conv[0],nt,nx*ny,dt,&time[is*ny*nx],&amp[is*ny*nx],shift,fmin,fmax);
			for (ix = 0; ix < ny*nx; ix++) {
				image[is*nshots+iy] += conv[ix*nt+delay+2]*dx*dy*dt;
				for (it=0; it<nt; it++) {
					block[is*nshots*nt+iy*nt+it] += conv[ix*nt+it];
				}
			}
		}

		free(hdr_time); free(time); free(hdr_amp); free(amp); free(conv);
	}
	free(gmin); free(hdr_gmin);

	/*----------------------------------------------------------------------------*
    *	Write out the data
    *----------------------------------------------------------------------------*/
	fp_out = fopen(file_out, "w+");

	if (nshots>1) dz = ((float)(gz[1] - gz[0]))/1000.0;
	else dz = 1.0;

	for (it = 0; it < nray; it++) {
		hdr_out[it].fldr		= 1;
		hdr_out[it].tracl		= it+1;
		hdr_out[it].tracf		= it+1;
		hdr_out[it].scalco		= -1000;
		hdr_out[it].scalel		= -1000;
		hdr_out[it].trid		= 1;
		hdr_out[it].ns			= nshots;
		hdr_out[it].trwf		= nray;
		hdr_out[it].ntr			= nray;
		hdr_out[it].f1			= (((float)gz[0])/1000.0);
		hdr_out[it].f2			= (((float)gx[0])/1000.0);
		hdr_out[it].dt			= ((int)(dt*1E6));
		hdr_out[it].d1			= roundf(dz*1000.0)/1000.0;
		hdr_out[it].d2			= roundf(dx*1000.0)/1000.0;
		hdr_out[it].gx			= gx[it*nshots];
		hdr_out[it].gy			= gy[it*nshots];
	}

	ret = writeData3D(fp_out, &image[0], hdr_out, nshots, nray);
	if (ret < 0 ) verr("error on writing output file.");

	fclose(fp_out);

	fp_out = fopen("block.su", "w+");

	if (nshots>1) dz = ((float)(gz[1] - gz[0]))/1000.0;
	else dz = 1.0;

	for (is = 0; is < nshots; is++) {
		for (it = 0; it < nray; it++) {
			hdr_out[it*nshots+is].fldr		= it+1;
			hdr_out[it*nshots+is].tracl		= it+1;
			hdr_out[it*nshots+is].tracf		= it+1;
			hdr_out[it*nshots+is].scalco	= -1000;
			hdr_out[it*nshots+is].scalel	= -1000;
			hdr_out[it*nshots+is].trid		= 1;
			hdr_out[it*nshots+is].ns		= nt;
			hdr_out[it*nshots+is].trwf		= nray;
			hdr_out[it*nshots+is].ntr		= nray;
			hdr_out[it*nshots+is].f1		= (((float)gz[0])/1000.0);
			hdr_out[it*nshots+is].f2		= (((float)gx[0])/1000.0);
			hdr_out[it*nshots+is].dt		= ((int)(dt*1E6));
			hdr_out[it*nshots+is].d1		= roundf(dz*1000.0)/1000.0;
			hdr_out[it*nshots+is].d2		= roundf(dx*1000.0)/1000.0;
			hdr_out[it*nshots+is].gx		= gx[it*nshots+is];
			hdr_out[it*nshots+is].gy		= gy[it*nshots+is];
			hdr_out[it*nshots+is].sdepth	= gz[it*nshots+is];
		}
	}

	ret = writeData3D(fp_out, &block[0], hdr_out, nt, nray*nshots);
	if (ret < 0 ) verr("error on writing output file.");

	fclose(fp_out);

	free(image); free(hdr_out); free(block);

	vmess("Wrote data");

	return 0;
}

void timeShift(float *data, long nsam, long nrec, float dt, float *time, float *amp, float *delay, float fmin, float fmax)
{
	long 	optn, iom, iomin, iomax, nfreq, ix, sign;
	float	omin, omax, deltom, om, tom, df, *rdata, scl;
	complex *cdata, *cdatascl;

	optn = optncr(nsam);
	nfreq = optn/2+1;
	df    = 1.0/(optn*dt);

	cdata = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");

	rdata = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");

	/* pad zeroes until Fourier length is reached */
	pad_data(data,nsam,nrec,optn,rdata);

	/* Forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata[0], &cdata[0], optn, nrec, optn, nfreq, sign);

	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;
	iomin  = (long)MIN((omin/deltom), (nfreq));
	iomax  = MIN((long)(omax/deltom), (nfreq));

	cdatascl = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");

	for (ix = 0; ix < nrec; ix++) {
		for (iom = 0; iom < iomin; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		for (iom = iomax; iom < nfreq; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		for (iom = iomin ; iom < iomax ; iom++) {
			om = deltom*iom;
			tom = om*-1.0*(time[ix]-delay[ix]);
			cdatascl[ix*nfreq+iom].r = (cdata[ix*nfreq+iom].r*cos(-tom) - cdata[ix*nfreq+iom].i*sin(-tom))/(amp[ix]*amp[ix]);
			cdatascl[ix*nfreq+iom].i = (cdata[ix*nfreq+iom].i*cos(-tom) + cdata[ix*nfreq+iom].r*sin(-tom))/(amp[ix]*amp[ix]);
		}
	}
	free(cdata);

	/* Inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&cdatascl[0], &rdata[0], optn, nrec, nfreq, optn, sign);
	scl_data(rdata,optn,nrec,scl,data,nsam);

	free(cdatascl);
	free(rdata);

	return;
}

void pad_data(float *data, long nsam, long nrec, long nsamout, float *datout)
{
	long it,ix;
	for (ix=0;ix<nrec;ix++) {
		for (it=0;it<nsam;it++)
			datout[ix*nsamout+it]=data[ix*nsam+it];
		for (it=nsam;it<nsamout;it++)
			datout[ix*nsamout+it]=0.0;
	}
}

void scl_data(float *data, long nsam, long nrec, float scl, float *datout, long nsamout)
{
	long it,ix;
	for (ix = 0; ix < nrec; ix++) {
		for (it = 0 ; it < nsamout ; it++)
			datout[ix*nsamout+it] = scl*data[ix*nsam+it];
	}
}