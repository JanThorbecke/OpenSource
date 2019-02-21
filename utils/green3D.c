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

int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);
int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);
int readData(FILE *fp, float *data, segy *hdrs, int n1);

void xwgreen3D(float *data, int nt, int nx, int ny, float dt, float fmin, float fmax, float *xi, float xsrc,
			float dx, float *yi, float ysrc, float dy, float *zi, float zsrc, float c, float cs, float rho,
			float *wavelet, float dipx, float maxdip, int far, int p_vz, int dip, int verbose);


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
"   yrcv=-1500,1500 .......... y-position's of receivers (array)",
"   zrcv=0,0 ................. z-position's of receivers (array)",
"   dxrcv=15 ................. step in receiver x-direction",
"   dyrcv=15 ................. step in receiver y-direction",
"   var=0 .................... variance for irregular sampling (dxrcv +- var)",
"   seed=0 ................... seed for random generator",
"   lint=1 ................... linear interpolate between the rcv points",
"   rrcv= .................... radius for receivers on a circle ",
"   oxrcv=0.0 ................ x-center position of circle",
"   oyrcv=0.0 ................ y-center position of circle",
"   ozrcv=0.0 ................ z-center position of circle",
"   dphi=2 ................... angle between receivers on circle ",
" SOURCE POSITIONS ",
"   xsrc1=0.0 ................ x-position of first source",
"   xsrc2=xsrc1 .............. x-position of last source",
"   dxsrc=0.0 ................ step in source x-direction",
"   ysrc1=0.0 ................ y-position of first source",
"   ysrc2=ysrc1 .............. y-position of last source",
"   dysrc=0.0 ................ step in source y-direction",
"   zsrc2=zsrc1 .............. depth position (z) of last source",
"   dzsrc=0.0 ................ step in source z-direction",
" SAMPLING AND SOURCE DEFINITION ",
"   file_src=spike ........... source wavelet (overrules dt)",
"   nt=256 ................... number of samples",
"   dt=0.004 ................. stepsize in time-direction ",
"   fmin=0 ................... minimum frequency",
"   fmax=70 .................. maximum frequency",
"   dipx=0 ................... local dip of the dipole in x-direction",
"   dipy=0 ................... local dip of the dipole in y-direction",
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
	int     n1, n2, n3, i, j, l, nrx, nry, nrz, dip;
	int     far, p_vz, nt, nx, ny, Nsx, Nsy, is, isy, sum, lint, verbose;
	int     size, ntraces, ngath, Fz, Fx;
	float   scl, xmin, xmax, ymin, ymax;
	float   dx, dy, dt, d1, d2, d3, fmin, fmax, f1, f2, f3, c, cs, rho;
	float 	*data, *wavelet, *tmpdata, dipx, dipy, xsrc1, xsrc2, ysrc1, ysrc2;
	float 	*xrcv, *yrcv, *zrcv, *xi, *yi, *zi, x0, y0, maxdip;
    float   rrcv, dphi, oxrcv, ozrcv;
	float	zsrc1, zsrc2, dxsrc, dysrc, dzsrc, xsrc, ysrc, zsrc, dxrcv, dyrcv;
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
    nry  = countparval("yrcv");
	// nrz  = countparval("zrcv");
	nrz = 0;
	if(!getparfloat("dxrcv",&dxrcv)) dxrcv = 15;
    if(!getparfloat("dyrcv",&dyrcv)) dyrcv = 15;

	if (nrx != 0 && nry != 0 && nrz == 0) {
		if (nrx != 2) verr("xrcv should have only two values");
        if (nry != 2) verr("yrcv should have only two values");
		xrcv = (float *)malloc(nrx*sizeof(float));
        yrcv = (float *)malloc(nry*sizeof(float));
		getparfloat("xrcv",xrcv);
        getparfloat("yrcv",yrcv);
		nx = NINT((xrcv[1] - xrcv[0])/dxrcv) + 1;
        ny = NINT((yrcv[1] - yrcv[0])/dyrcv) + 1;
		xi = (float *)malloc(nx*ny*sizeof(float));
        yi = (float *)malloc(nx*ny*sizeof(float));
		zi = (float *)malloc(nx*ny*sizeof(float));
		x0 = xrcv[0];
        y0 = yrcv[0];
		for (i = 0; i < ny; i++) {
            for (j = 0; j < nx; j++) {
                xi[i*nx+j] = x0 + j*dxrcv;
                yi[i*nx+j] = y0 + i*dyrcv;
			    zi[i*nx+j] = 0;
            }
		}
	}
	else if (nrx == 0 && nry == 0 && nrz == 0) {
		nx = NINT((3000)/dxrcv) + 1;
		ny = NINT((3000)/dyrcv) + 1;
		xi = (float *)malloc(nx*ny*sizeof(float));
		yi = (float *)malloc(nx*ny*sizeof(float));
		zi = (float *)malloc(nx*ny*sizeof(float));
		x0 = -1500;
		y0 = -1500;
		for (i = 0; i < ny; i++) {
            for (j = 0; j < nx; j++) {
                xi[i*nx+j] = x0 + j*dxrcv;
                yi[i*nx+j] = y0 + i*dyrcv;
			    zi[i*nx+j] = 0;
            }
		}
	}
	else verr("Number of xrcv and yrcv values are not equal");

	if (verbose) vmess("number of receivers nx = %d, ny = %d total = %d", nx, ny, nx*ny);
	if (verbose == 13) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				vmess("xi = %d yi = %d x = %f y=%f z = %f", i, j, xi[j*nx+i], yi[j*nx+i], zi[j*nx+i]);
			}
		}
	}

	if(!getparfloat("xsrc1", &xsrc1)) xsrc1=0;
	if(!getparfloat("xsrc2", &xsrc2)) xsrc2=xsrc1;
	if(!getparfloat("dxsrc", &dxsrc)) dxsrc=0.0;
    if(!getparfloat("ysrc1", &ysrc1)) ysrc1=0;
	if(!getparfloat("ysrc2", &ysrc2)) ysrc2=ysrc1;
	if(!getparfloat("dysrc", &dysrc)) dysrc=0.0;
	if(!getparfloat("zsrc2", &zsrc2)) zsrc2=zsrc1;
	if(!getparfloat("dzsrc", &dzsrc)) dzsrc=0;
	if(!getparint("nt", &nt)) nt = 256;
	if(!getparfloat("fmin", &fmin)) fmin = 0.0;
	if(!getparfloat("fmax", &fmax)) fmax = 70.0;
	if(!getparfloat("dipx", &dipx)) dipx = 0.0;
    if(!getparfloat("dipy", &dipy)) dipy = 0.0;
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

	if (xsrc2==xsrc1) Nsx = 1;
	else Nsx = NINT((xsrc2 - xsrc1)/dxsrc) + 1;
	if (ysrc2==ysrc1) Nsy = 1;
	else Nsy = NINT((ysrc2 - ysrc1)/dysrc) + 1;

	if (verbose) vmess("Number of shot records to generate x = %d y = %d", Nsx, Nsy);
	if (Nsx > 1 && Nsy > 1) {
		dxsrc = (xsrc2-xsrc1)/(Nsx-1);
		dysrc = (ysrc2-ysrc1)/(Nsy-1);
		dzsrc = (zsrc2-zsrc1)/(Nsx-1);
		if (verbose) {
			vmess("dxsrc = %f", dxsrc);
			vmess("dysrc = %f", dysrc);
			vmess("dzsrc = %f", dzsrc);
		}
	}

	size = nt * nx *ny;
	dx   = dxrcv;
	dy   = dyrcv;
	tmpdata = (float *)calloc(size,sizeof(float));
	data = (float *)calloc(size,sizeof(float));
	hdrs = (segy *) calloc(nx*ny,sizeof(segy));
	for (i = 0; i < ny; i++) {
		for(j = 0; j < nx; j++) {
			hdrs[i*nx+j].f1= 0.0;
			hdrs[i*nx+j].f2= x0;
			hdrs[i*nx+j].d1= dt;
			hdrs[i*nx+j].d2= dx;
			hdrs[i*nx+j].ns= nt;
			hdrs[i*nx+j].dt= (int)1000000*dt;
			hdrs[i*nx+j].trwf= nx*ny;
			hdrs[i*nx+j].tracl= i*nx+j+1;
			hdrs[i*nx+j].tracf= i*nx+j+1;
			hdrs[i*nx+j].gx = (x0 + j*dx)*1000;
			hdrs[i*nx+j].gy = (y0 + i*dy)*1000;
			hdrs[i*nx+j].scalco = -1000;
			hdrs[i*nx+j].trid = TREAL;
		}
	}
	if (file_out==NULL) fp_out=stdout;
	else fp_out = fopen(file_out,"w");
	if (fp_out == NULL) verr("error in creating output file");

	for (isy = 0; isy < Nsy; isy++) {
		for (is = 0; is < Nsx; is++) {
			xsrc = xsrc1 + is*dxsrc;
			ysrc = ysrc1 + isy*dysrc;
			zsrc = zsrc1 + is*dzsrc;
			if (verbose) vmess("xsrc = %f ysrc=%f zsrc = %f", xsrc, ysrc, zsrc);

			xwgreen3D(data,nt,nx,ny,dt,fmin,fmax,xi,xsrc,dx,yi,ysrc,dy,zi,zsrc,c,cs,rho,wavelet,
				dipx, maxdip, far, p_vz, dip, verbose);

			for (l = 0; l < ny; l++) {
				for (i = 0; i < nx; i++) {
					for (j = 0; j < nt; j++) tmpdata[l*nx*nt+i*nt+j] = data[l*nx*nt+i*nt+j];
					hdrs[l*nx+i].sx = NINT(xsrc*1000);
					hdrs[l*nx+i].sy = NINT(ysrc*1000);
					hdrs[l*nx+i].scalco = -1000;
					hdrs[l*nx+i].offset = xi[l*nx+i]-xsrc;
					hdrs[l*nx+i].gx = NINT(xi[l*nx+i]*1000);
					hdrs[l*nx+i].gy = NINT(yi[l*nx+i]*1000);
					hdrs[l*nx+i].fldr = isy*Nsx+is+1;
					hdrs[l*nx+i].trwf = nx*ny;
					nwrite = fwrite( &hdrs[l*nx+i], 1, TRCBYTES, fp_out);
					assert(nwrite == TRCBYTES);
					nwrite = fwrite( &tmpdata[l*nx*nt+i*nt], sizeof(float), nt, fp_out);
					assert(nwrite == nt);
				}
			}
		}
	}

	if( xi ) free(xi);
	if( yi ) free(yi);
	if( zi ) free(zi);
	if( wavelet ) free( wavelet );

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

void xwgreen3D(float *data, int nt, int nx, int ny, float dt, float fmin, float fmax, float *xi, float xsrc, float dx, float *yi, float ysrc, float dy, float *zi, float zsrc, float c, float cs, float rho, float *wavelet, float dipx, float maxdip, int far, int p_vz, int dip, int verbose)
{
	int    	iomin, iomax, iom, ix, iy, nfreq, i, sign, optn;
	float  	df, deltom, om, k, r, x, y, invr, phi, phi2, cosphi;
	float	*rwave, *rdata, cos2, scl, z, kp, ks, sclr;
	complex	*cwave, *cdata, tmp, tmp2, ekr, sum;
    complex H02p, H12p, H02s, H12s, Gp, Gs;

	optn	= optncr(nt);
	nfreq	= 1+(optn/2);
	df		= 1.0/(dt*optn);
	deltom	= 2.*M_PI*df;
	iomin	= (int)MIN((fmin*dt*optn), (nfreq-1));
	iomin	= MAX(iomin, 1);
	iomax	= MIN((int)(fmax*dt*optn), (nfreq-1));

	rdata = (float *)calloc(optn*nx*ny,sizeof(float));
	cdata = (complex *)calloc(nfreq*nx*ny,sizeof(complex));
	rwave = (float *)calloc(optn,sizeof(float));
	cwave = (complex *)calloc(nfreq,sizeof(complex));

	for (i = 0; i < nt; i++) rwave[i] = wavelet[i]*dt;
	for (i = nt; i < optn; i++) rwave[i] = 0.0;
	
	sign = -1;
	rc1fft(rwave, cwave, optn, sign);

	for (iy = 0; iy < ny; iy++) {
		for (ix = 0; ix < nx; ix++) {
			for (iom = 0; iom < iomin; iom++) {
				cdata[iy*nx*nfreq+ix*nfreq+iom].r = 0.0;
				cdata[iy*nx*nfreq+ix*nfreq+iom].i = 0.0;
			}
		}
	}
	for (iy = 0; iy < ny; iy++) {
		for (ix = 0; ix < nx; ix++) {
			for (iom = iomax; iom < nfreq; iom++) {
				cdata[iy*nx*nfreq+ix*nfreq+iom].r = 0.0;
				cdata[iy*nx*nfreq+ix*nfreq+iom].i = 0.0;
			}
		}
	}

	if (p_vz == 0) {
		if (dip == 1) {
			if (verbose) vmess("P field of dipole");
			for (iy = 0; iy < ny; iy++) {
				for (ix = 0; ix < nx; ix++) {
					x      = xi[iy*nx+ix] - xsrc;
					y      = yi[iy*nx+ix] - ysrc;
					z      = fabs(zi[iy*nx+ix] - zsrc);
					r      = sqrt(x*x + y*y + z*z);
					if (r != 0) phi = acos(z/r);
					else phi = M_PI/2;
					phi2   = SGN(x)*phi - (dipx*M_PI/180.0);
					cosphi = 0.25*cos(phi2)*rho;
					if (fabs(phi) < maxdip*M_PI/180.0) {
						/* exp(-jkr) = cos(kr) - j*sin(kr) */
						for (iom = iomin; iom <= iomax; iom++) {
							om = iom*deltom;
							k = om/c;
							ekr.r = cos(k*r)/(r*r);
							ekr.i = sin(-k*r)/(r*r);
							tmp.r = cosphi*(ekr.r - k*r*ekr.i);
							tmp.i = cosphi*(ekr.i + k*r*ekr.r);
							cdata[iy*nx*nfreq+ix*nfreq+iom].r = tmp.r*cwave[iom].r - 
											tmp.i*cwave[iom].i;
							cdata[iy*nx*nfreq+ix*nfreq+iom].i = tmp.r*cwave[iom].i +
											tmp.i*cwave[iom].r;
						}
					}
					else {
						for (iom = iomin; iom <= iomax; iom++) {
							cdata[iy*nx*nfreq+ix*nfreq+iom].r = 0.0;
							cdata[iy*nx*nfreq+ix*nfreq+iom].i = 0.0;
						}
					}
				}
			}
		}
/* There is no far-field definition 'needed' in 3D 
		else if (far == 1 && dip == 1){
			if (verbose) vmess("far P field of dipole");
			for (iy = 0; iy < ny; iy++) {
				for (ix = 0; ix < nx; ix++) {
					x = xi[ix] - xsrc;
					y = yi[iy*nx+ix] - ysrc;
					z = fabs(zi[iy*nx+ix] - zsrc);
					r = sqrt(x*x + y*y + z*z);
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

							cdata[iy*nx*nfreq+ix*nfreq+iom].r = tmp.r*cwave[iom].r -
											tmp.i*cwave[iom].i;
							cdata[iy*nx*nfreq+ix*nfreq+iom].i = tmp.r*cwave[iom].i +
											tmp.i*cwave[iom].r;
						}
					}
					else {
						for (iom = iomin; iom <= iomax; iom++) {
							cdata[iy*nx*nfreq+ix*nfreq+iom].r = 0.0;
							cdata[iy*nx*nfreq+ix*nfreq+iom].i = 0.0;
						}
					}
				}
			}
		}
*/
		else if (dip == 0){
			if (verbose) vmess("P field of monopole");
			for (iy = 0; iy < ny; iy++) {
				for (ix = 0; ix < nx; ix++) {
					x = xi[iy*nx+ix] - xsrc;
					y = yi[iy*nx+ix] - ysrc;
					z = fabs(zi[iy*nx+ix] - zsrc);
					r = sqrt(x*x + y*y + z*z);
					if (r != 0) phi = acos(z/r);
					else phi = M_PI/2;
					scl = 0.25*rho;
					if (fabs(phi) < maxdip*M_PI/180.0) {
						for (iom = iomin; iom <= iomax; iom++) {
							om = iom*deltom;
							k  = om/c;
							tmp.r = scl*cos(k*r)/(r);
							tmp.i = scl*sin(-k*r)/(r);

							cdata[iy*nx*nfreq+ix*nfreq+iom].r = tmp.r*cwave[iom].r - 
											tmp.i*cwave[iom].i;
							cdata[iy*nx*nfreq+ix*nfreq+iom].i = tmp.r*cwave[iom].i +
											tmp.i*cwave[iom].r;
						}
					}
					else {
						for (iom = iomin; iom <= iomax; iom++) {
							cdata[iy*nx*nfreq+ix*nfreq+iom].r = 0.0;
							cdata[iy*nx*nfreq+ix*nfreq+iom].i = 0.0;
						}
					}
				}
			}
		}
/* There is no far-field definition 'needed' in 3D 
		else if (far == 1 && dip == 0){
			if (verbose) vmess("far P field of monopole");
			for (iy = 0; iy < ny; iy++) {
				for (ix = 0; ix < nx; ix++) {
					x = xi[iy*nx+ix] - xsrc;
					y = yi[iy*nx+ix] - ysrc;
					z = fabs(zi[iy*nx+ix] - zsrc);
					r = sqrt(x*x + y*y + z*z);
					if (r != 0) phi = acos(z/r);
					else phi = M_PI*0.5;
					scl = 0.5*rho/sqrt(r);
					if (fabs(phi) <= M_PI*(maxdip/180.0)) {
						for (iom = iomin; iom <= iomax; iom++) {
							om = iom*deltom;
							k = om/c;
							tmp.r = -sqrt(1.0/(2.0*M_PI*k))*scl*sin(k*r-M_PI/4.0);
							tmp.i = -sqrt(1.0/(2.0*M_PI*k))*scl*cos(k*r-M_PI/4.0);

							cdata[iy*nx*nfreq+ix*nfreq+iom].r = tmp.r*cwave[iom].r -
											tmp.i*cwave[iom].i;
							cdata[iy*nx*nfreq+ix*nfreq+iom].i = tmp.r*cwave[iom].i +
											tmp.i*cwave[iom].r;
						}
					}
					else {
						for (iom = iomin; iom <= iomax; iom++) {
							cdata[iy*nx*nfreq+ix*nfreq+iom].r = 0.0;
							cdata[iy*nx*nfreq+ix*nfreq+iom].i = 0.0;
						}
					}
				}
			}
		}
*/
	}
	else if (p_vz == 1) {
		if (dip == 1) {
	    		if (verbose) vmess("Vz field of dipole, this is not yet implemented");
				for (iy = 0; iy < ny; iy++) {
					for (ix = 0; ix < nx; ix++) {
						x = xi[iy*nx+ix] - xsrc;
						y = yi[iy*nx+ix] - ysrc;
						z = fabs(zi[iy*nx+ix] - zsrc);
						r = sqrt(x*x + y*y + z*z);
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

								cdata[iy*nx*nfreq+ix*nfreq+iom].r = sum.r*cwave[iom].r -
												sum.i*cwave[iom].i;
								cdata[iy*nx*nfreq+ix*nfreq+iom].i = sum.r*cwave[iom].i +
												sum.i*cwave[iom].r;
							}
						}
						else {
							for (iom = iomin; iom <= iomax; iom++) {
								cdata[iy*nx*nfreq+ix*nfreq+iom].r = 0.0;
								cdata[iy*nx*nfreq+ix*nfreq+iom].i = 0.0;
							}
						}
					}
				}
		}
		else {
	    	if (verbose) vmess("Vz field of monopole");

			for (iy = 0; iy < ny; iy++) {
				for (ix = 0; ix < nx; ix++) {
					x      = xi[iy*nx+ix] - xsrc;
					y      = yi[iy*nx+ix] - ysrc;
					z      = fabs(zi[iy*nx+ix] - zsrc);
					r      = sqrt(x*x + y*y + z*z);
					if (r != 0) phi = acos(z/r);
					else phi = M_PI/2;
					phi2   = SGN(x)*phi - (dipx*M_PI/180.0);
					cosphi = 0.25*cos(phi2)/c;
					if (fabs(phi) < maxdip*M_PI/180.0) {
						/* exp(-jkr) = cos(kr) - j*sin(kr) */
						for (iom = iomin; iom <= iomax; iom++) {
							om = iom*deltom;
							k = om/c;
							ekr.r = cos(k*r)/(r*r);
							ekr.i = sin(-k*r)/(r*r);
							tmp.r = cosphi*(ekr.r - k*r*ekr.i);
							tmp.i = cosphi*(ekr.i + k*r*ekr.r);
							cdata[iy*nx*nfreq+ix*nfreq+iom].r = tmp.r*cwave[iom].r - 
											tmp.i*cwave[iom].i;
							cdata[iy*nx*nfreq+ix*nfreq+iom].i = tmp.r*cwave[iom].i +
											tmp.i*cwave[iom].r;
						}
					}
					else {
						for (iom = iomin; iom <= iomax; iom++) {
							cdata[iy*nx*nfreq+ix*nfreq+iom].r = 0.0;
							cdata[iy*nx*nfreq+ix*nfreq+iom].i = 0.0;
						}
					}
				}
			}
		}
	}
    else if (p_vz == 2) { /* Fz source with Vz receivers Fz=1 == p_vz=2 */
        for (iy = 0; iy < ny; iy++) {
			for (ix = 0; ix < nx; ix++) {
				x = xi[iy*nx+ix] - xsrc;
				y = yi[iy*nx+ix] - ysrc;
				z = fabs(zi[iy*nx+ix] - zsrc);
				r = sqrt(x*x + y*y + z*z);

				if (r != 0) phi = acos(z/r);
				else phi = M_PI/2;
				phi2   = SGN(x)*phi - (dipx*M_PI/180.0);
				cosphi = cos(phi2);
				sclr = (z*z-x*x-y*y)/(r);
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

						cdata[iy*nx*nfreq+ix*nfreq+iom].r = tmp.r*cwave[iom].r - 
						tmp.i*cwave[iom].i;
						cdata[iy*nx*nfreq+ix*nfreq+iom].i = tmp.r*cwave[iom].i +
						tmp.i*cwave[iom].r;
					}
				}
				else {
					for (iom = iomin; iom <= iomax; iom++) {
						cdata[iy*nx*nfreq+ix*nfreq+iom].r = 0.0;
						cdata[iy*nx*nfreq+ix*nfreq+iom].i = 0.0;
					}
				}
			}
        }

    }
    else if (p_vz == 3) { /* Fx source with Vz receivers Fx=1 == p_vz=3 */
        for (iy = 0; iy < ny; iy++) {
			for (ix = 0; ix < nx; ix++) {
				x = xi[iy*nx+ix] - xsrc;
				y = yi[iy*nx+ix] - ysrc;
				z = fabs(zi[iy*nx+ix] - zsrc);
				r = sqrt(x*x + y*y + z*z);

				if (r != 0) phi = acos(z/r);
				else phi = M_PI/2;
				phi2   = SGN(x)*phi - (dipx*M_PI/180.0);
				cosphi = cos(phi2);
				scl = (z*x*y)/(4.0*r*r*rho);
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
						
						cdata[iy*nx*nfreq+ix*nfreq+iom].r = tmp.r*cwave[iom].r - 
						tmp.i*cwave[iom].i;
						cdata[iy*nx*nfreq+ix*nfreq+iom].i = tmp.r*cwave[iom].i +
						tmp.i*cwave[iom].r;
					}
				}
				else {
					for (iom = iomin; iom <= iomax; iom++) {
						cdata[iy*nx*nfreq+ix*nfreq+iom].r = 0.0;
						cdata[iy*nx*nfreq+ix*nfreq+iom].i = 0.0;
					}
				}
			}
        }

    }


	scl  = df;
	sign = 1;
	crmfft(&cdata[0], &rdata[0], optn, nx*ny, nfreq, optn, sign);
	for (iy = 0; iy < ny; iy++) {
		for (ix = 0; ix < nx; ix++) {
			for (i = 0; i < nt; i++) {
				data[iy*nx*nt+ix*nt+i] = scl*rdata[iy*nx*optn+ix*optn+i];
			}
		}
	}

	free(cdata);
	free(cwave);
	free(rdata);
	free(rwave);

	return;
}
