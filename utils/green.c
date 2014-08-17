#include <genfft.h>
#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)


#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

void getrecpos(float *xi, float *zi, int nx, float *xrcv, float *zrcv, int verbose);
int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);
int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);
int readData(FILE *fp, float *data, segy *hdrs, int n1);

void xwgreen(float *data, int nt, int nx, float dt, float fmin, float fmax, float *xi, float xsrc, 
			 float dx, float *zi, float zsrc, float c, float cs, float rho, float *wavelet, 
             float dipx, float maxdip, int far, int p_vz, int dip, int verbose);

/*********************** self documentation **********************/
char *sdoc[] = {
" 								",
" green - calculation of 2D Greens function in homogenoeus medium based one exact expressions",
" 								",
" green c= zsrc1= [optional parameters]",
" 							        ",
" Required parameters:",
" ",
"   c= ....................... P-wave velocity",
"   cs=0.7*c ................. S-wave velocity",
"   zsrc1= ................... depth of source",
" 							        ",
" Optional parameters:",
" ",
"   file_out= ................ output file (default SU-pipe)",
" RECEIVER POSITIONS ",
"   xrcv=-1500,1500 .......... x-position's of receivers (array)",
"   zrcv=0,0 ................. z-position's of receivers (array)",
"   dxrcv=15 ................. step in receiver x-direction",
"   var=0 .................... variance for irregular sampling (dxrcv +- var)",
"   seed=0 ................... seed for random generator",
"   lint=1 ................... linear interpolate between the rcv points",
"   rrcv= .................... radius for receivers on a circle ",
"   oxrcv=0.0 ................ x-center position of circle",
"   ozrcv=0.0 ................ z-center position of circle",
"   dphi=2 ................... angle between receivers on circle ",
" SOURCE POSITIONS ",
"   xsrc1=0.0 ................ x-position of first source",
"   xsrc2=xsrc1 .............. x-position of last source",
"   dxsrc=0.0 ................ step in source x-direction",
"   zsrc2=zsrc1 .............. depth position (z) of last source",
"   dzsrc=0.0 ................ step in source z-direction",
" SAMPLING AND SOURCE DEFINITION ",
"   file_src=spike ........... source wavelet (overrules dt)",
"   nt=256 ................... number of samples",
"   dt=0.004 ................. stepsize in time-direction ",
"   fmin=0 ................... minimum frequency",
"   fmax=70 .................. maximum frequency",
"   dipx=0 ................... local dip of the dipole",
"   dip=1 .................... 1; dipole 0; monopole source",
"   rho=1000 ................. density",
" FIELD DEFINITION ",
"   far=0 .................... farfield approximation 0=off)",
"   p_vz=0  .................. P or Vz field (0 = P field, 1 = Vz field)",
"   Fz=0  .................... Force source in z with Vz receivers",
"   Fx=0  .................... Force source in x with Vz receivers",
"   maxdip=90 ................ maximum angle (degrees) to be computed ",
"   sum=0 .................... sum all sources",
"   verbose=0 ................ silent option; >0 display info",
"",
"  The P or Vz field of a dipole source at depth z below the receivers",
"  in a homogeneous 2-D medium is calculated.",
"   ",
" author  : Jan Thorbecke : 23-03-1995 (janth@xs4all.nl)",
"                         : revision 2010",
" ",
NULL};
/**************** end self doc ***********************************/

int main(int argc, char **argv)
{
	FILE	*fp_in, *fp_out;
	int     n1, n2, i, ret, j, nkeys, nrx, nrz, dip, key_idx;
	int     far, p_vz, nt, nx, Ns, is, sum, lint, verbose;
	int     size, ntraces, ngath, Fz, Fx;
	float   scl, xmin, xmax;
	float   dx, dt, d1, d2, fmin, fmax, f1, f2, c, cs, rho;
	float 	*data, *wavelet, *tmpdata, dipx, xsrc1, xsrc2;
	float 	*xrcv, *zrcv, *xi, *zi, x0, maxdip;
    float   rrcv, dphi, oxrcv, ozrcv;
	float	zsrc1, zsrc2, dxsrc, dzsrc, xsrc, zsrc, dxrcv;
	char    *file_src, *file_out;
	size_t  nwrite;
	segy	*hdrs;

/* ========================= Reading parameters ====================== */

	initargs(argc, argv);
	requestdoc(1);

	if(!getparint("verbose", &verbose)) verbose = 0;
	if(!getparstring("file_out", &file_out)){
		if (verbose) vwarn("parameter file_out not found, assume pipe");
		file_out = NULL;
	}
	if(!getparstring("file_src", &file_src)) file_src = NULL;
	if(!getparfloat("c", &c)) verr("velocity must be specified.");
    if(!getparfloat("cs", &cs)) cs=0.7*c;
	if(!getparfloat("zsrc1", &zsrc1)) verr("zsrc1(depth) must be specified.");
	if(!getparint("lint", &lint)) lint=1;
	if(!getparfloat("maxdip", &maxdip)) maxdip=90.0;

	nrx  = countparval("xrcv");
	nrz  = countparval("zrcv");
	if(!getparfloat("dxrcv",&dxrcv)) dxrcv = 15;

    if (getparfloat("rrcv", &rrcv)) {
		if (!getparfloat("dphi",&dphi)) dphi=2.0;
		if (!getparfloat("oxrcv",&oxrcv)) oxrcv=0.0;
		if (!getparfloat("ozrcv",&ozrcv)) ozrcv=0.0;
		nx = NINT(360.0/dphi);
        xi = (float *)malloc(nx*sizeof(float));
		zi = (float *)malloc(nx*sizeof(float));
        getrecpos(xi, zi, nx, xrcv, zrcv, verbose);
    }
	else if (nrx == 0 && nrz == 0) {
		nx = NINT((3000)/dxrcv) + 1;
		xi = (float *)malloc(nx*sizeof(float));
		zi = (float *)malloc(nx*sizeof(float));
		x0 = -1500;
		for (i = 0; i < nx; i++) {
			xi[i] = x0 + i*dxrcv;
			zi[i] = 0;
		}
	}
	else if (nrx != 0 && nrz == 0) {
		if (nrx != 2) verr("xrcv should have only two values");
		xrcv = (float *)malloc(nrx*sizeof(float));
		getparfloat("xrcv",xrcv);
		nx = NINT((xrcv[1] - xrcv[0])/dxrcv) + 1;
		xi = (float *)malloc(nx*sizeof(float));
		zi = (float *)malloc(nx*sizeof(float));
		x0 = xrcv[0];
		for (i = 0; i < nx; i++) {
			xi[i] = x0 + i*dxrcv;
			zi[i] = 0;
		}
	}
	else if (nrx == nrz) {
		xrcv = (float *)malloc(nrx*sizeof(float));
		zrcv = (float *)malloc(nrz*sizeof(float));
		getparfloat("xrcv",xrcv);
		getparfloat("zrcv",zrcv);
		if (lint) nx = NINT((xrcv[nrx-1] - xrcv[0])/dxrcv) + 1;
		else nx = nrx;
		xi = (float *)malloc(nx*sizeof(float));
		zi = (float *)malloc(nx*sizeof(float));
		x0 = xrcv[0];
		getrecpos(xi, zi, nx, xrcv, zrcv, verbose);
	}
	else verr("Number of xrcv and zrcv values are not equal");

	if (verbose) vmess("number of receivers = %d", nx);
	if (verbose == 13) {
		for (i = 0; i < nx; i++) {
			vmess("i = %d x = %f z = %f", i, xi[i], zi[i]);
		}
	}

	if(!getparfloat("xsrc1", &xsrc1)) xsrc1=0;
	if(!getparfloat("xsrc2", &xsrc2)) xsrc2=xsrc1;
	if(!getparfloat("dxsrc", &dxsrc)) dxsrc=0.0;
	if(!getparfloat("zsrc2", &zsrc2)) zsrc2=zsrc1;
	if(!getparfloat("dzsrc", &dzsrc)) dzsrc=0;
	if(!getparint("nt", &nt)) nt = 256;
	if(!getparfloat("fmin", &fmin)) fmin = 0.0;
	if(!getparfloat("fmax", &fmax)) fmax = 70.0;
	if(!getparfloat("dipx", &dipx)) dipx = 0.0;
	if(!getparfloat("rho", &rho)) rho = 1000.0;
	if(!getparint("far", &far)) far = 0;
	if(!getparint("p_vz", &p_vz)) p_vz = 0;
    if(!getparint("Fz", &Fz)) Fz = 0;
    if(!getparint("Fx", &Fx)) Fx = 0;
	if(!getparint("dip", &dip)) dip = 1;
	if(!getparint("sum", &sum)) sum = 0;
    if(Fz) p_vz=2;
    if(Fx) p_vz=3;

/* ========================= Opening wavelet file ====================== */

	if (file_src == NULL){
		if(!getparfloat("dt", &dt)) dt = 0.004;
		wavelet = (float *)calloc(nt,sizeof(float));
		wavelet[0] = 1.0;
	}
	else {
		if (verbose) vmess("Reading wavelet from file %s.", file_src);
		ngath = 1;
		getFileInfo(file_src, &n1, &n2, &ngath, &d1, &d2, &f1, &f2, &xmin, &xmax, &scl, &ntraces);
		
		fp_in = fopen(file_src, "r");
		if (fp_in == NULL) verr("error on opening input file_src=%s", file_src);
		
		tmpdata = (float *)calloc(n1*n2,sizeof(float));
		hdrs = (segy *) calloc(n2,sizeof(segy));
		
		n2 = readData(fp_in, tmpdata, hdrs, n1);
		fclose(fp_in);
		if (verbose) {
			disp_fileinfo(file_src, n1, n2, f1,  f2,  d1,  d2, hdrs);
		}
		dt = d1;
		wavelet = (float *)calloc(nt,sizeof(float));

		if (n1 <= nt) {
			for (i = 0; i < n1; i++) wavelet[i] = tmpdata[i];
			for (i = n1; i < nt; i++) wavelet[i] = 0.0;
		}
		else {
			vwarn("file_src has more samples than output");
			for (i = 0; i < nt; i++) wavelet[i] = tmpdata[i];
		}
		if( tmpdata ) free(tmpdata);
		if( hdrs ) free( (void *) hdrs);
	}

/* ============ INITIALIZE AND CHECK PARAMETERS =============== */

	if (dxsrc == 0 && dzsrc == 0) 
		Ns = 1;
	else if (dxsrc == 0 && dzsrc != 0)
		Ns = NINT((zsrc2 - zsrc1)/dzsrc) + 1;
	else if (dzsrc == 0 && dxsrc != 0)
		Ns = NINT((xsrc2 - xsrc1)/dxsrc) + 1;
	else if (dzsrc != 0 && dxsrc != 0)
		Ns = MAX(NINT((xsrc2 - xsrc1)/dxsrc), NINT((zsrc2 - zsrc1)/dzsrc))+1;

	if (verbose) vmess("Number of shot records to generate = %d", Ns);
	if (Ns > 1) {
		dxsrc = (xsrc2-xsrc1)/(Ns-1);
		dzsrc = (zsrc2-zsrc1)/(Ns-1);
		if (verbose) {
			vmess("dxsrc = %f", dxsrc);
			vmess("dzsrc = %f", dzsrc);
		}
	}

	size = nt * nx;
	dx   = dxrcv;
	tmpdata = (float *)calloc(size,sizeof(float));
	data = (float *)calloc(size,sizeof(float));
	hdrs = (segy *) calloc(nx,sizeof(segy));
	for(j = 0; j < nx; j++) {
		hdrs[j].f1= 0.0;
		hdrs[j].f2= x0;
		hdrs[j].d1= dt;
		hdrs[j].d2= dx;
		hdrs[j].ns= nt;
		hdrs[j].dt= (int)1000000*dt;
		hdrs[j].trwf= nx;
		hdrs[j].tracl= j+1;
		hdrs[j].tracf= j+1;
		hdrs[j].gx = (x0 + j*dx)*1000;
		hdrs[j].scalco = -1000;
		hdrs[j].trid = TREAL;
	}
	if (file_out==NULL) fp_out=stdout;
	else fp_out = fopen(file_out,"w");
	if (fp_out == NULL) verr("error in creating output file");

	for (is = 0; is < Ns; is++) {
		xsrc = xsrc1 + is*dxsrc;
		zsrc = zsrc1 + is*dzsrc;
		if (verbose) vmess("xsrc = %f zsrc = %f", xsrc, zsrc);

		xwgreen(data,nt,nx,dt,fmin,fmax,xi,xsrc,dx,zi,zsrc,c,cs,rho,wavelet,
			dipx, maxdip, far, p_vz, dip, verbose);

		if (sum == 1) {
			for (i = 0; i < nx; i++) {
				for (j = 0; j < nt; j++) tmpdata[i*nt+j] += data[i*nt+j];
			}
		}
		else {
			for (i = 0; i < nx; i++) {
				for (j = 0; j < nt; j++) tmpdata[i*nt+j] = data[i*nt+j];
				hdrs[i].sx = NINT(xsrc*1000);
				hdrs[i].scalco = -1000;
				hdrs[i].offset = xi[i]-xsrc;
				hdrs[i].gx = NINT(xi[i]*1000);
				hdrs[i].fldr = is+1;
				hdrs[i].trwf = nx;
				nwrite = fwrite( &hdrs[i], 1, TRCBYTES, fp_out);
				assert(nwrite == TRCBYTES);
				nwrite = fwrite( &tmpdata[i*nt], sizeof(float), nt, fp_out);
				assert(nwrite == nt);
			}
		}
	}

	if( xi ) free(xi);
	if( zi ) free(zi);
	if( wavelet ) free( wavelet );

	if (sum == 1) {
		for (i = 0; i < nx; i++) {
			hdrs[i].sx = NINT(xsrc1*1000);
			hdrs[i].scalco = -1000;
			hdrs[i].offset = x0-xsrc1;
			hdrs[i].gx = NINT(x0 + i*dx)*1000;
			hdrs[i].fldr = 1;
			hdrs[i].trwf = nx;
			nwrite = fwrite( &hdrs[i], 1, TRCBYTES, fp_out);
			assert(nwrite == TRCBYTES);
			nwrite = fwrite( &tmpdata[i*nt], sizeof(float), nt, fp_out);
			assert(nwrite == nt);
		}
	}

    fclose(fp_out);

	if( data ) free(data);
	if( tmpdata ) free(tmpdata);
	if( hdrs ) free( hdrs);

	exit ( 0 );
}

/***************************************************************************
*  
*   Calculation of pulse response in homogeneous medium
*
*
***************************************************************************/

void xwgreen(float *data, int nt, int nx, float dt, float fmin, float fmax, float *xi, float xsrc, float dx, float *zi, float zsrc, float c, float cs, float rho, float *wavelet, float dipx, float maxdip, int far, int p_vz, int dip, int verbose)
{
	int    	iomin, iomax, iom, ix, nfreq, i, sign, optn;
	float  	df, deltom, om, k, r, x, invr, phi, phi2, cosphi;
	float	*rwave, *rdata, cos2, scl, z, kp, ks, sclr;
	complex	*cwave, *cdata, tmp, tmp2, sum;
    complex H02p, H12p, H02s, H12s, Gp, Gs;

	optn	= optncr(nt);
	nfreq	= 1+(optn/2);
	df		= 1.0/(dt*optn);
	deltom	= 2.*M_PI*df;
	iomin	= (int)MIN((fmin*dt*optn), (nfreq-1));
	iomin	= MAX(iomin, 1);
	iomax	= MIN((int)(fmax*dt*optn), (nfreq-1));

	rdata = (float *)calloc(optn*nx,sizeof(float));
	cdata = (complex *)calloc(nfreq*nx,sizeof(complex));
	rwave = (float *)calloc(optn,sizeof(float));
	cwave = (complex *)calloc(nfreq,sizeof(complex));

	for (i = 0; i < nt; i++) rwave[i] = wavelet[i]*dt;
	for (i = nt; i < optn; i++) rwave[i] = 0.0;
	
	sign = -1;
	rc1fft(rwave, cwave, optn, sign);

	for (ix = 0; ix < nx; ix++) {
		for (iom = 0; iom < iomin; iom++) {
			cdata[ix*nfreq+iom].r = 0.0;
			cdata[ix*nfreq+iom].i = 0.0;
		}
	}
	for (ix = 0; ix < nx; ix++) {
		for (iom = iomax; iom < nfreq; iom++) {
			cdata[ix*nfreq+iom].r = 0.0;
			cdata[ix*nfreq+iom].i = 0.0;
		}
	}

	if (p_vz == 0) {
		if (far == 0 && dip == 1) {
			if (verbose) vmess("near and far P field of dipole");
			for (ix = 0; ix < nx; ix++) {
				x      = xi[ix] - xsrc;
				z      = fabs(zi[ix] - zsrc);
				r      = sqrt(x*x + z*z);
				if (r != 0) phi = acos(z/r);
				else phi = M_PI/2;
				phi2   = SGN(x)*phi - (dipx*M_PI/180.0);
				cosphi = 0.25*cos(phi2)*rho;
				if (fabs(phi) < maxdip*M_PI/180.0) {
					for (iom = iomin; iom <= iomax; iom++) {
						om = iom*deltom;
						k = om/c;
						tmp.r = -k*cosphi*y1(k*r);
						tmp.i = -k*cosphi*j1(k*r);
						cdata[ix*nfreq+iom].r = tmp.r*cwave[iom].r - 
										   tmp.i*cwave[iom].i;
						cdata[ix*nfreq+iom].i = tmp.r*cwave[iom].i +
										   tmp.i*cwave[iom].r;
					}
				}
				else {
					for (iom = iomin; iom <= iomax; iom++) {
						cdata[ix*nfreq+iom].r = 0.0;
						cdata[ix*nfreq+iom].i = 0.0;
					}
				}
			}
		}
		else if (far == 1 && dip == 1){
			if (verbose) vmess("far P field of dipole");
			for (ix = 0; ix < nx; ix++) {
				x = xi[ix] - xsrc;
				z = fabs(zi[ix] - zsrc);
				r = sqrt(x*x + z*z);
				if (r != 0) phi = acos(z/r);
				else phi = M_PI/2;
				phi2   = SGN(x)*phi - (dipx*M_PI/180.0);
				cosphi = 0.5*cos(phi2)*rho/sqrt(r);
				if (fabs(phi) < maxdip*M_PI/180.0) {
					for (iom = iomin; iom <= iomax; iom++) {
						om = iom*deltom;
						k = om/c;
						tmp.r = sqrt(k/(2.0*M_PI))*cosphi*cos(k*r-M_PI/4.0);
						tmp.i = -sqrt(k/(2.0*M_PI))*cosphi*sin(k*r-M_PI/4.0);

						cdata[ix*nfreq+iom].r = tmp.r*cwave[iom].r -
										   tmp.i*cwave[iom].i;
						cdata[ix*nfreq+iom].i = tmp.r*cwave[iom].i +
										   tmp.i*cwave[iom].r;
					}
				}
				else {
					for (iom = iomin; iom <= iomax; iom++) {
						cdata[ix*nfreq+iom].r = 0.0;
						cdata[ix*nfreq+iom].i = 0.0;
					}
				}
			}
		}
		else if (far == 0 && dip == 0){
			if (verbose) vmess("near and far P field of monopole");
			for (ix = 0; ix < nx; ix++) {
				x = xi[ix] - xsrc;
				z = fabs(zi[ix] - zsrc);
				r = sqrt(x*x + z*z);
				if (r != 0) phi = acos(z/r);
				else phi = M_PI/2;
				scl = 0.25*rho;
				if (fabs(phi) < maxdip*M_PI/180.0) {
					for (iom = iomin; iom <= iomax; iom++) {
						om = iom*deltom;
						k  = om/c;
						tmp.r = -scl*y0(k*r);
						tmp.i = -scl*j0(k*r);

						cdata[ix*nfreq+iom].r = tmp.r*cwave[iom].r - 
										   tmp.i*cwave[iom].i;
						cdata[ix*nfreq+iom].i = tmp.r*cwave[iom].i +
										   tmp.i*cwave[iom].r;
					}
				}
				else {
					for (iom = iomin; iom <= iomax; iom++) {
						cdata[ix*nfreq+iom].r = 0.0;
						cdata[ix*nfreq+iom].i = 0.0;
					}
				}
			}
		}
		else if (far == 1 && dip == 0){
			if (verbose) vmess("far P field of monopole");
			for (ix = 0; ix < nx; ix++) {
				x = xi[ix] - xsrc;
				z = fabs(zi[ix] - zsrc);
				r = sqrt(x*x + z*z);
				if (r != 0) phi = acos(z/r);
				else phi = M_PI*0.5;
				scl = 0.5*rho/sqrt(r);
				if (fabs(phi) <= M_PI*(maxdip/180.0)) {
					for (iom = iomin; iom <= iomax; iom++) {
						om = iom*deltom;
						k = om/c;
						tmp.r = -sqrt(1.0/(2.0*M_PI*k))*scl*sin(k*r-M_PI/4.0);
						tmp.i = -sqrt(1.0/(2.0*M_PI*k))*scl*cos(k*r-M_PI/4.0);

						cdata[ix*nfreq+iom].r = tmp.r*cwave[iom].r -
										   tmp.i*cwave[iom].i;
						cdata[ix*nfreq+iom].i = tmp.r*cwave[iom].i +
										   tmp.i*cwave[iom].r;
					}
				}
				else {
					for (iom = iomin; iom <= iomax; iom++) {
						cdata[ix*nfreq+iom].r = 0.0;
						cdata[ix*nfreq+iom].i = 0.0;
					}
				}
			}
		}
	}
	else if (p_vz == 1) {
		if (dip == 1) {
	    	if (far == 0) {
	    		if (verbose) vmess("near and far Vz field of dipole");
	    		for (ix = 0; ix < nx; ix++) {
	    			x = xi[ix] - xsrc;
	    			z = fabs(zi[ix] - zsrc);
	    			r = sqrt(x*x + z*z);
	    			invr   = -0.25/(c);
	    			if (r != 0) phi = acos(z/r);
	    			else phi = M_PI/2;
	    			phi2   = SGN(x)*phi - (dipx*M_PI/180.0);
	    			cosphi = cos(phi2);
	    			cos2 = cosphi*cosphi;
	    			if (fabs(phi) < maxdip*M_PI/180.0) {
	    				for (iom = iomin; iom <= iomax; iom++) {
	    					om = iom*deltom;
	    					k = om/c;
	    					tmp.r = k*cos2*invr*j0(k*r);
	    					tmp.i = -k*cos2*invr*y0(k*r);
	    					tmp2.r = k*(1-2*cos2)*invr*j1(k*r)/(k*r);
	    					tmp2.i = -k*(1-2*cos2)*invr*y1(k*r)/(k*r);
	    					sum.r = tmp.r + tmp2.r;
	    					sum.i = tmp.i + tmp2.i;

	    					cdata[ix*nfreq+iom].r = sum.r*cwave[iom].r -
	    									   sum.i*cwave[iom].i;
	    					cdata[ix*nfreq+iom].i = sum.r*cwave[iom].i +
	    									   sum.i*cwave[iom].r;
	    				}
	    			}
	    			else {
	    				for (iom = iomin; iom <= iomax; iom++) {
	    					cdata[ix*nfreq+iom].r = 0.0;
	    					cdata[ix*nfreq+iom].i = 0.0;
	    				}
	    			}
	    		}
	    	}
	    	else {
	    		if (verbose) vmess("far Vz field of dipole");
	    		for (ix = 0; ix < nx; ix++) {
	    			x = xi[ix] - xsrc;
	    			z = fabs(zi[ix] - zsrc);
	    			r = sqrt(x*x + z*z);
	    			invr   = -0.5/(c*sqrt(r));
	    			if (r != 0) phi = acos(z/r);
	    			else phi = M_PI/2;
	    			phi2   = SGN(x)*phi - (dipx*M_PI/180.0);
	    			cosphi = cos(phi2);
	    			cos2 = cosphi*cosphi;
	    			if (fabs(phi) < maxdip*M_PI/180.0) {
	    				for (iom = iomin; iom <= iomax; iom++) {
	    					om = iom*deltom;
	    					k = om/c;
	    					tmp.r = sqrt(k/(2.0*M_PI))*invr*cos2*cos(k*r-M_PI/4.0);
	    					tmp.i = -sqrt(k/(2.0*M_PI))*invr*cos2*sin(k*r-M_PI/4.0);

	    					cdata[ix*nfreq+iom].r = tmp.r*cwave[iom].r -
	    									   tmp.i*cwave[iom].i;
	    					cdata[ix*nfreq+iom].i = tmp.r*cwave[iom].i +
	    									   tmp.i*cwave[iom].r;
	    				}
	    			}
	    			else {
	    				for (iom = iomin; iom <= iomax; iom++) {
	    					cdata[ix*nfreq+iom].r = 0.0;
	    					cdata[ix*nfreq+iom].i = 0.0;
	    				}
	    			}
	    		}
	    	}
		}
		else {
	    	if (verbose) vmess("near and far Vz field of monopole");

			for (ix = 0; ix < nx; ix++) {
				x      = xi[ix] - xsrc;
				z      = fabs(zi[ix] - zsrc);
				r      = sqrt(x*x + z*z);
				if (r != 0) phi = acos(z/r);
				else phi = M_PI/2;
				phi2   = SGN(x)*phi - (dipx*M_PI/180.0);
				cosphi = cos(phi2);
				if (fabs(phi) < maxdip*M_PI/180.0) {
					for (iom = iomin; iom <= iomax; iom++) {
						om = iom*deltom;
						k = om/c;
						tmp.i = -cosphi*y1(k*r)/(4.0*c);
						tmp.r = cosphi*j1(k*r)/(4.0*c);
						cdata[ix*nfreq+iom].r = tmp.r*cwave[iom].r - 
										   tmp.i*cwave[iom].i;
						cdata[ix*nfreq+iom].i = tmp.r*cwave[iom].i +
										   tmp.i*cwave[iom].r;
					}
				}
				else {
					for (iom = iomin; iom <= iomax; iom++) {
						cdata[ix*nfreq+iom].r = 0.0;
						cdata[ix*nfreq+iom].i = 0.0;
					}
				}
			}
		}
	}
    else if (p_vz == 2) { /* Fz source with Vz receivers Fz=1 == p_vz=2 */
        for (ix = 0; ix < nx; ix++) {
            x      = xi[ix] - xsrc;
            z      = fabs(zi[ix] - zsrc);
            r      = sqrt(x*x + z*z);

            if (r != 0) phi = acos(z/r);
            else phi = M_PI/2;
            phi2   = SGN(x)*phi - (dipx*M_PI/180.0);
            cosphi = cos(phi2);
            sclr = (z*z-x*x)/(r);
            if (fabs(phi) < maxdip*M_PI/180.0) {
                for (iom = iomin; iom <= iomax; iom++) {
                    om = iom*deltom;
                    kp = om/c;
                    ks = om/cs;
                    H02p.r = j0(kp*r);
                    H02p.i = -y0(kp*r);
                    H12p.r = j1(kp*r);
                    H12p.i = -y1(kp*r);
                    
                    H02s.r = j0(ks*r);
                    H02s.i = -y0(ks*r);
                    H12s.r = j1(ks*r);
                    H12s.i = -y1(ks*r);

                    Gp.r = kp/(4*om*rho*r*r)*(-z*z*kp*H02p.r + sclr*H12p.r);
                    Gp.i = kp/(4*om*rho*r*r)*(-z*z*kp*H02p.i + sclr*H12p.i);

                    Gs.r = ks/(4*om*rho*r*r)*(-z*z*ks*H02s.r + sclr*H12s.r);
                    Gs.i = ks/(4*om*rho*r*r)*(-z*z*ks*H02s.i + sclr*H12s.i);

                    tmp.i = (-1.0/om)*((om/(4*rho*cs*cs))*(H02s.r) - Gp.r + Gs.r);
                    tmp.r = ( 1.0/om)*((om/(4*rho*cs*cs))*(H02s.i) - Gp.i + Gs.i);

                    cdata[ix*nfreq+iom].r = tmp.r*cwave[iom].r - 
                    tmp.i*cwave[iom].i;
                    cdata[ix*nfreq+iom].i = tmp.r*cwave[iom].i +
                    tmp.i*cwave[iom].r;
                }
            }
            else {
                for (iom = iomin; iom <= iomax; iom++) {
                    cdata[ix*nfreq+iom].r = 0.0;
                    cdata[ix*nfreq+iom].i = 0.0;
                }
            }
        }

    }
    else if (p_vz == 3) { /* Fx source with Vz receivers Fx=1 == p_vz=3 */
        for (ix = 0; ix < nx; ix++) {
            x      = xi[ix] - xsrc;
            z      = fabs(zi[ix] - zsrc);
            r      = sqrt(x*x + z*z);

            if (r != 0) phi = acos(z/r);
            else phi = M_PI/2;
            phi2   = SGN(x)*phi - (dipx*M_PI/180.0);
            cosphi = cos(phi2);
            scl = (z*x)/(4.0*r*r*rho);
            if (fabs(phi) < maxdip*M_PI/180.0) {
                for (iom = iomin; iom <= iomax; iom++) {
                    om = iom*deltom;
                    kp = om/c;
                    ks = om/cs;
                    H02p.r = kp*kp*j0(kp*r);
                    H02p.i = -kp*kp*y0(kp*r);
                    H12p.r = 2.0*kp*j1(kp*r)/r;
                    H12p.i = -2.0*kp*y1(kp*r)/r;
                    
                    H02s.r = ks*ks*j0(ks*r);
                    H02s.i = -ks*ks*y0(ks*r);
                    H12s.r = 2.0*ks*j1(ks*r)/r;
                    H12s.i = -2.0*ks*y1(ks*r)/r;
                    
                    tmp.i = (scl/(om*om))*((H02p.r-H12p.r) - (H02s.r-H12s.r));
                    tmp.r = -(scl/(om*om))*((H02p.i-H12p.i) - (H02s.i-H12s.i));
                    
                    cdata[ix*nfreq+iom].r = tmp.r*cwave[iom].r - 
                    tmp.i*cwave[iom].i;
                    cdata[ix*nfreq+iom].i = tmp.r*cwave[iom].i +
                    tmp.i*cwave[iom].r;
                }
            }
            else {
                for (iom = iomin; iom <= iomax; iom++) {
                    cdata[ix*nfreq+iom].r = 0.0;
                    cdata[ix*nfreq+iom].i = 0.0;
                }
            }
        }

    }


	scl  = df;
	sign = 1;
	crmfft(&cdata[0], &rdata[0], optn, nx, nfreq, optn, sign);
	for (ix = 0; ix < nx; ix++) {
		for (i = 0; i < nt; i++) {
			data[ix*nt+i] = scl*rdata[ix*optn+i];
		}
	}

	free(cdata);
	free(cwave);
	free(rdata);
	free(rwave);

	return;
}
