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
int readData(FILE *fp, float *data, segy *hdrs, int n1);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);
int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);
double wallclock_time(void);

void synthesis(float *shotdata, float *syndata, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float xsrc, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int di, float fmin, float fmax, int mode, int reci, int k, int off, int nmo, int nshots, int verbose);

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" SYN2D - synthesize shots records to CFP gathers in frequency domain",
" ",
" syn2d file_syn= file_shot= nshots= [optional parameters]",
" ",
" Required parameters: ",
" ",
"   file_syn= ................ synthesis operator(s)",
"   file_shot= ............... shot records (can be input pipe) ",
"   nshots= .................. number of shot records",
" ",
" Optional parameters: ",
" ",
" INPUT DEFINITION ",
"   nxmax=512 ................ maximum number of traces in input files",
"   ntmax=1024 ............... maximum number of samples/trace in input files",
" SYNTHESIS ",
"   ixa=0 .................... number of traces after focus point",
"   ixb=ixa .................. number of traces before focus point",
"   tap=0 .................... lateral taper synthesis(1), shot(2) or both(3)",
"   ntap=0 ................... number of taper points at boundaries",
"   nmo=0 .................... 1; corrects CFP-gather with operator times",
"   reci=0 ................... 1; add focusing in emission 2; emission only",
"   mode=-1 .................. -1; inverse 1; forward",
"   off=0 .................... trace offset to start synthesis",
"   alpha=0 .................. Laplace factor (good default = -2)",
"   fmin=0 ................... minimum frequency",
"   fmax=70 .................. maximum frequency",
" OUTPUT DEFINITION ",
"   file_cfp= ................ output file with areal shot records",
"   file_fresnel= ............ output file with shot records convolved with operator",
"   verbose=0 ................ silent option; >0 displays info",
" ",
" The nmo option determines the display of the calculated CFP-gather.",
"     nmo = 0: displays CFP-gather at operator times,",
"     nmo = 1: displays CFP-gather corrected with operator times,",
" ",
"  Note that if ixa=0 and ixb=0 all shots are used.",
" ",
" author  : Jan Thorbecke : 4-9-1995 (jan@delphi.tn.tudelft.nl)",
    " product : Originates from DELPHI software",
    "                         : revision 2012",
" ",
NULL};
/**************** end self doc ***********************************/

int main (int argc, char **argv)
{
    FILE	*fp_syn, *fp_shot, *fp_out, *fp_fres;
	int		i, j, k, l, ret, nshots, Nsyn, nt, nx, nts, nxs, more, ngath;
	int		size, n1, n2, ntap, tap, di, ixrcv, ixsrc, off, optn, ntraces;
	int		nxmax, ntmax, reci, mode, nmo, ixa, ixb, n2out, verbose;
	float	fmin, fmax, *taper, fxf, dxf, fxs2, xsrc, *xrcv;
    double  t0, t1, t2, t3, tsyn, tread;
	float	*shotdata, d1, d2, f1, f2, fts, fxs, ft, fx, *etap, *xsyn, dxsrc;
	float   *syndata, *tmpdata, dt, dx, dts, dxs, scl, alpha, mem;
	float   max, *zsyn, scel, xmin, xmax;
	char	*file_syn, *file_shot, *file_cfp, *file_fresnel;
	segy	*hdrs, *hdrs_in, *hdrs_out;

	initargs(argc, argv);
	requestdoc(1);

	if (!getparstring("file_shot", &file_shot)) file_shot = NULL;
	if (!getparstring("file_syn", &file_syn)) file_syn = NULL;
	if (!getparstring("file_fresnel", &file_fresnel)) file_fresnel = NULL;
	if (!getparint("verbose", &verbose)) verbose = 0;
	if (file_syn == NULL && file_shot == NULL) 
		verr("file_syn and file_shot cannot be both input pipe");
	if (!getparstring("file_cfp", &file_cfp)) {
		if (verbose) vwarn("parameter file_cfp not found, assume pipe");
		file_cfp = NULL;
	}
	if (!getparint("nshots", &nshots)) verr("nshots should be given");
	if (!getparfloat("fmin", &fmin)) fmin = 0.0;
	if (!getparfloat("fmax", &fmax)) fmax = 70.0;
	if (!getparint("ixa", &ixa)) ixa = 0;
	if (!getparint("ixb", &ixb)) ixb = ixa;
	if (!getparint("reci", &reci)) reci = 0;
	if (!getparint("off", &off)) off = 0;
	if (!getparint("mode", &mode)) mode = -1;
	if (!getparfloat("alpha", &alpha)) alpha = 0.0;
	if (!getparint("tap", &tap)) tap = 0;
	if (!getparint("ntap", &ntap)) ntap = 0;
	if (!getparint("nmo", &nmo)) nmo = 0;

	if(mode >= 0) mode = 1;
	if(mode < 0) mode = -1;
	if (nmo && mode > 0) verr("nmo = 1 cannot be used with mode = 1");
	if (reci && ntap) vwarn("tapering influences the reciprocal result");

/*================ Reading all synthesis operator(s) ================*/

	ngath = 1;
	ret = getFileInfo(file_syn, &n1, &n2, &ngath, &d1, &d2, &f1, &f2, &xmin, &xmax, &scl, &ntraces);
    if (ret == 0) {
        if (!getparint("ntmax", &ntmax)) ntmax = n1;
        if (!getparint("nxmax", &nxmax)) nxmax = n2;
        if (verbose>=2 && (ntmax!=n1 || nxmax!=n2))
            vmess("dimensions overruled: %d x %d",ntmax,nxmax);
    }
    else {
        if (verbose>=2) vmess("dimensions used: %d x %d",ntmax,nxmax);
    }
    
	size    = ntmax * nxmax;
	tmpdata = (float *)malloc(size*sizeof(float));
	hdrs = (segy *) calloc(nxmax,sizeof(segy));

    
    fp_syn = fopen(file_syn, "r");
	if (fp_syn == NULL) verr("error on opening input file_syn=%s", file_syn);
    
	n2 = readData(fp_syn, tmpdata, hdrs, n1);
	if (n2 == 0) {
		fclose(fp_syn);
		if (verbose) vmess("end of file_syn data reached");
	}
	if (verbose) {
		disp_fileinfo(file_syn, n1, n2, f1, f2, d1, d2, hdrs);
	}
	nxs = n2; 
    nts = n1;

                             
    if (hdrs[0].scalco < 0) scl = 1.0/fabs(hdrs[0].scalco);
	else if (hdrs[0].scalco == 0) scl = 1.0;
	else scl = hdrs[0].scalco;

	if (hdrs[0].scalel < 0) scel = 1.0/fabs(hdrs[0].scalel);
	else if (hdrs[0].scalel == 0) scel = 1.0;
	else scel = hdrs[0].scalel;


	taper = (float *)malloc(n2*sizeof(float));
	if (tap == 1 || tap == 3) {
		for (j = 0; j < ntap; j++)
			taper[j] = (cos(PI*(j-ntap)/ntap)+1)/2.0;
		for (j = ntap; j < n2-ntap; j++)
			taper[j] = 1.0;
		for (j = n2-ntap; j < n2; j++)
			taper[j] =(cos(PI*(j-(n2-ntap))/ntap)+1)/2.0;
	}
	else {
		for (j = 0; j < n2; j++) taper[j] = 1.0;
	}

    etap = (float *)malloc(n1*sizeof(float));
	for (j = 0; j < n1; j++) etap[j] = exp(alpha*j*d1);

	syndata = (float *)malloc(n1*n2*sizeof(float));
	xsyn    = (float *)malloc(sizeof(float));
	zsyn    = (float *)malloc(sizeof(float));
	Nsyn    = 0;
    
	while (n2 > 0) {
		syndata = (float *)realloc(syndata, (Nsyn+1)*n1*n2*sizeof(float));
		if (syndata == NULL) verr("Out of memory for syndata array!");
		xsyn = (float *)realloc(xsyn, (Nsyn+1)*sizeof(float));
		if (xsyn == NULL) verr("Out of memory for xsyn array!");
		zsyn = (float *)realloc(zsyn, (Nsyn+1)*sizeof(float));
		if (zsyn == NULL) verr("Out of memory for zsyn array!");

		xsyn[Nsyn] = hdrs[0].sx*scl;
		zsyn[Nsyn] = hdrs[0].selev*scel;
		for (i = 0; i < n2; i++) {
			for (j = 0; j < n1; j++) 
				syndata[Nsyn*n1*n2+i*n1+j] = tmpdata[i*n1+j]*taper[i]*etap[j];
		}
		Nsyn += 1;
        n2 = readData(fp_syn, tmpdata, hdrs, n1);
		if (verbose>=2 && n2!=0) disp_fileinfo(file_syn, n1, n2, f1, f2, d1, d2, hdrs);
	}
    fclose(fp_syn);

    
	dxs = d2; dts = d1;
	fxs = f2; fts = f1;
	if (hdrs[0].gx != 0 || hdrs[1].gx != 0 ) fxs = hdrs[0].gx*scl;
	fxs2 = fxs + (float)(nxs-1)*dxs;
	dxf = (hdrs[nxs-1].gx - hdrs[0].gx)*scl/(float)(nxs-1);
	if (NINT(dxs*1e3) != NINT(fabs(dxf)*1e3)) {
		vmess("dx in hdr.d1 (%.3f) and hdr.gx (%.3f) not equal",d2, dxf);
		if (dxf != 0) dxs = fabs(dxf);
		vmess("dx in operator => %f", dxs);
	}

	free(tmpdata);
	free(taper);
	free(etap);
	free(hdrs);


/*================ Reading first shot record ================*/

    ngath = 1;
	ret = getFileInfo(file_shot, &nt, &nx, &ngath, &d1, &dx, &ft, &fx, &xmin, &xmax, &scl, &ntraces);
    if (ret == 0) {
        if (!getparint("ntmax", &ntmax)) ntmax = nt;
        if (!getparint("nxmax", &nxmax)) nxmax = nx;
        if (verbose>=2 && (ntmax!=nt || nxmax!=nx))
            vmess("dimensions overruled: %d x %d",ntmax,nxmax);
    }
    else {
        if (verbose>=2) vmess("dimensions used: %d x %d",ntmax,nxmax);
    }
	if (!getparfloat("dt", &dt)) dt = d1;
	fprintf(stderr,"dt=%e\n", dt);
    
    size     = ntmax * nxmax;
	xrcv     = (float *)malloc(nxmax*sizeof(float));
	tmpdata  = (float *)malloc(size*sizeof(float));
	hdrs_in  = (segy *) calloc(nxmax,sizeof(segy));
	if (tmpdata == NULL || hdrs_in == NULL)
		verr("memory allocation error for input data");

    fp_shot = fopen(file_shot, "r");
	if (fp_shot == NULL) verr("error on opening input file_shot=%s", file_shot);
    
	nx = readData(fp_shot, tmpdata, hdrs_in, nt);
	if (nx == 0) {
		fclose(fp_shot);
		if (verbose) verr("end of file_shot data reached");
	}
	if (verbose) {
		disp_fileinfo(file_shot, nt, nx, ft, fx, dt, dx, hdrs_in);
	}

	taper = (float *)malloc(nx*sizeof(float));
	if (tap == 2 || tap == 3) {
		for (j = 0; j < ntap; j++)
			taper[j] = (cos(PI*(j-ntap)/ntap)+1)/2.0;
		for (j = ntap; j < nx-ntap; j++)
			taper[j] = 1.0;
		for (j = nx-ntap; j < nx; j++)
			taper[j] =(cos(PI*(j-(nx-ntap))/ntap)+1)/2.0;
	}
	else {
		for (j = 0; j < nx; j++) taper[j] = 1.0;
	}

	optn = optncr(MAX(nt, nts)); 
    etap = (float *)malloc(optn*sizeof(float));
	for (j = 0; j < optn; j++) etap[j] = exp(alpha*j*dt);

	shotdata = (float *)malloc(optn*nxmax*sizeof(float));
	for (i = 0; i < nx; i++) {
		for (j = 0; j < nt; j++) 
			shotdata[i*optn+j] = tmpdata[i*nt+j]*taper[i]*etap[j];
		for (j = nt; j < optn; j++) shotdata[i*optn+j] = 0.0;
	}

	if (hdrs_in[0].scalco < 0) scl = 1.0/fabs(hdrs_in[0].scalco);
	else if (hdrs_in[0].scalco == 0) scl = 1.0;
	else scl = hdrs_in[0].scalco;

	fxf = (float)hdrs_in[0].sx*scl;
	if (nx > 1) dxf = (hdrs_in[nx-1].gx - hdrs_in[0].gx)*scl/(float)(nx-1);
	else dxf = d2;
	if (NINT(dx*1e3) != NINT(fabs(dxf)*1e3)) {
		vmess("dx in hdr.d1 (%.3f) and hdr.gx (%.3f) not equal",dx, dxf);
		if (dxf != 0) dx = fabs(dxf);
		else verr("gx hdrs not set");
		vmess("dx used => %f", dx);
	}
    

/*=========== Reading second shot record (for source sampling) ===========*/

	hdrs    = (segy *) calloc(nxmax,sizeof(segy));
	if (hdrs == NULL)
		verr("memory allocation error for input data");

    n2 = readData(fp_shot, tmpdata, hdrs, nt);
	if (n2 == 0) {
		fclose(fp_shot);
        vwarn("only one shot record is available");
		dxsrc = dx;
		more = 0;
	}
    else {
        if (verbose>=3) {
            disp_fileinfo(file_shot, nt, n2, ft, fx, dt, dx, hdrs);
        }  
        dxsrc = (float)hdrs[0].sx*scl - fxf;
		if (dxsrc == 0) {
			vwarn("sx hdrs are not filled in!!");
			dxsrc = dx;
		}
		more = 1;
        nx = n2;
    }


/*================ Check the size of the files ================*/

	if (NINT(dxsrc/dx)*dx != NINT(dxsrc)) {
		vwarn("source (%.2f) and receiver step (%.2f) don't match",dxsrc,dx);
		if (reci == 2) vwarn("step used from operator (%.2f) ",dxs);
	}
	di = NINT(dxf/dxs);
	if ((NINT(di*dxs) != NINT(dxf)) && verbose) 
		vwarn("dx in receiver (%.2f) and operator (%.2f) don't match",dx,dxs);
	if (nt != nts) 
		vmess("Time samples in shot (%d) and synthesis operator (%d) are not equal",nt, nts);
	if (verbose) {
		vmess("Number of synthesis operators  = %d", Nsyn);
		vmess("number of shots                = %d", nshots);
		vmess("first model position           = %.2f", fxs);
		vmess("last model position            = %.2f", fxs2);
		vmess("first source position fxf      = %.2f", fxf);
		vmess("source distance dxsrc          = %.2f", dxsrc);
		vmess("last source position           = %.2f", fxf+(nshots-1)*dxsrc);
		vmess("receiver distance     dxf      = %.2f", dxf);
		vmess("direction of increasing traces = %d", di);
		vmess("number of time samples(fft)  s = %d %d", nt, optn);
	}

/*================ initializations ================*/

	if (ixa || ixb) n2out = ixa + ixb + 1;
	else if (reci) n2out = nxs;
	else n2out = nshots;
	mem = Nsyn*n2out*nts*sizeof(float)/1048576.0;
	if (verbose) {
		vmess("number of output traces        = %d", n2out);
		vmess("number of output samples       = %d", nts);
		vmess("Size of output data            = %.1f Mb", mem);
	}

	t0   = wallclock_time();
	tsyn = tread = 0.0;
	ret  = 0;
	k    = 0;
	size = ntmax*nxmax;

	if (file_fresnel) {
        fp_fres = fopen(file_fresnel, "w+");
	}

/*================ loop over all shot records ================*/

	while (nx > 0) {
		t1    = wallclock_time();
		xsrc  = (float)hdrs_in[0].sx*scl;
		for (i = 0; i < nx; i++) xrcv[i] = (float)hdrs_in[i].gx*scl;
		if (verbose>=2) {
			vmess("source position:     %.2f", xsrc);
			vmess("receiver positions:  %.2f <--> %.2f", xrcv[0], xrcv[nx-1]);
		}

		if ((NINT(xsrc-fxs2) > 0) || (NINT(xrcv[nx-1]-fxs2) > 0) ||
			(NINT(xrcv[nx-1]-fxs) < 0) || (NINT(xsrc-fxs) < 0) || 
			(NINT(xrcv[0]-fxs) < 0) || (NINT(xrcv[0]-fxs2) > 0) ) {
			vwarn("source/receiver positions are outside synthesis model");
			vwarn("CFP calculation is stopped at gather %d", k);
			vmess("xsrc = %.2f xrcv_1 = %.2f xrvc_N = %.2f", xsrc, xrcv[0], xrcv[nx-1]);
			fclose(fp_shot);
			break;
		}

		synthesis(shotdata, syndata, nx, nt, nxs, nts, dt, xsyn, Nsyn, 
			xrcv, xsrc, fxs, dxs, dxsrc, dx, ixa, ixb, di, fmin, fmax, mode, 
			reci, k, off, nmo, nshots, verbose);

		t3 = wallclock_time();
		tsyn +=  t3 - t1;

		if (file_fresnel) {
			for (i = 0; i < nx; i++) {
				hdrs_in[i].ns   = optn;
			}
            ret = writeData(fp_fres, (float *)&shotdata[0], hdrs_in, optn, nx);
            if (ret < 0 ) verr("error on writing output file.");
		}

		if (k == 0 && more) {
			nx = n2;
			for (i = 0; i < nx; i++) {
				hdrs_in[i] = hdrs[i];
				for (j = 0; j < nt; j++) 
					shotdata[i*optn+j] = tmpdata[i*nt+j];
				for (j = nt; j < optn; j++) 
					shotdata[i*optn+j] = 0.0;
			}
			free(hdrs);

			t2 = wallclock_time();
			if (verbose>=3) {
				vmess("CPU-time one gather      = %.3f s.", t2-t1);
				vmess("with CPU-time synthesis  = %.3f", tsyn);
				vmess("and CPU-time read data   = %.3f", t2-t3);
			}
		}
		else {
            nx = readData(fp_shot, tmpdata, hdrs_in, nt);

            for (i = 0; i < nx; i++) {
                for (j = 0; j < nt; j++) shotdata[i*optn+j] = tmpdata[i*nt+j];
                for (j = nt; j < optn; j++) shotdata[i*optn+j] = 0.0;
            }
		}
		k++;
		if (verbose) vmess("*** Shot gather %d processed ***", k);

		if (nx <= 0 ) {
			ret = fclose(fp_shot);
			if (ret < 0) vwarn("err %d on closing input file",ret);
			if (verbose) vmess("end of data reached");
			break;
		}
		for (i = 0; i < nx; i++) {
			for (j = 0; j < nt; j++) shotdata[i*optn+j] *= (taper[i]*etap[j]);
			for (j = nt; j < optn; j++) shotdata[i*optn+j] = 0.0;
		}

        if (verbose>=3) disp_fileinfo(file_shot, nt, nx, ft, fx, dt, dx, hdrs_in);
		tread += wallclock_time() - t3;
	}

	if (file_fresnel) {
		ret = fclose(fp_fres);
		if (ret < 0) verr("err %d on closing fresnel output file",ret);
	}

	t2 = wallclock_time();
	if (verbose) {
		vmess("Total CPU-time synthesis = %.3f", t2-t0);
		vmess("with CPU-time synthesis  = %.3f", tsyn);
		vmess("and CPU-time read data   = %.3f", tread);
	}

	if (k != nshots) {
		vwarn("number of shots given not equal to shots read in (%d)", k);
		if (ixa == 0 && ixb == 0 && reci == 0) n2out = k;
	}

/*================ write output files and free memory ================*/

	n1 = nts; n2 = n2out;
	f1 = ft; f2 = fxs;
	d1 = dt;
	if (reci == 0) d2 = dxsrc;
	else if (reci == 1) d2 = dxs;
	else if (reci == 2) d2 = dx;

	hdrs_out = (segy *) calloc(n2,sizeof(segy));
	if (hdrs_out == NULL) verr("allocation for hdrs_out");
	size  = nxs*nts;

	free(etap);
    etap = (float *)malloc(n1*sizeof(float));
	for (j = 0; j < n1; j++) etap[j] = exp(-alpha*j*dt);

    fp_out = fopen(file_cfp, "w+");
    if (fp_out==NULL) verr("error on ceating output file");

	for (l = 0; l < Nsyn; l++) {
		if (ixa || ixb) f2 = xsyn[l]-ixb*d2;
		else {
			if (reci) f2 = fxs;
			else f2 = fxf;
		}

		for (i = 0; i < n2; i++) {
            hdrs_out[i].ns     = n1;
            hdrs_out[i].trid   = hdrs_in[i].trid;
            hdrs_out[i].f1     = f1;
            hdrs_out[i].f2     = f2;
            hdrs_out[i].d1     = d1;
            hdrs_out[i].d2     = d2;
			hdrs_out[i].fldr   = l+1;
			hdrs_out[i].trwf   = n2out;
			hdrs_out[i].scalco = -1000;
			hdrs_out[i].sx = NINT(xsyn[l]*1000);
			hdrs_out[i].gx = NINT(1000*(f2+i*d2));
			hdrs_out[i].offset = (long)NINT((f2+i*d2) - xsyn[l]);
			hdrs_out[i].scalel = -1000;
			hdrs_out[i].selev  = NINT(zsyn[l]*1000);
			hdrs_out[i].sdepth = NINT(zsyn[l]*1000);
			for (j = 0; j < n1; j++) syndata[l*size+i*n1+j] *= etap[j];
		}

        ret = writeData(fp_out, (float *)&syndata[l*size], hdrs_out, n1, n2);
        if (ret < 0 ) verr("error on writing output file.");
	}
	ret = fclose(fp_out);
	if (ret < 0) verr("err %d on closing output file",ret);

	if (verbose) {
		t1 = wallclock_time();
		vmess("and CPU-time write data  = %.3f", t1-t2);
	}

	free(hdrs_in);
	free(hdrs_out);
	free(tmpdata);

	exit(0);
}

void synthesis(float *shotdata, float *syndata, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float xsrc, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int di, float fmin, float fmax, int mode, int reci, int k, int off, int nmo, int nshots, int verbose)
{
	static int first=0, iomin, iomax, nfreq, optn, Niom, size, iox, inx;
	static complex *syncdata;
	static float scl;
	int 	i, j, l, m, ixsrc, ixsyn, ix, ixrcv, dosrc;
	float	*rdata, *p, **dum, x0, x1;
	static double t0, t1, tfft, t;
	complex *sum, *cdata, *shotcdata, tmp, ts, to;
	int      npe;


/*============= Transform synthesis operator(s), only one time =============*/

	if (first == 0) {
		t     = 0.0;
		tfft  = 0.0;
		size  = nxs*nts;
		optn  = optncr(MAX(nt, nts));
		scl   = 1.0/optn;
		nfreq = (optn+2)/2;
		iomin = (int)MIN((fmin*optn*dt), nfreq-1);
		iomin = MAX(iomin, 1);
		iomax = MIN((int)(fmax*optn*dt), nfreq-1);
		Niom  = iomax - iomin + 1;
		if (verbose>=3) 
			vmess("iomin = %d iomax =%d Niom = %d",iomin,iomax,Niom);
		mode *= -1;

		syncdata = (complex *)malloc(nxs*Niom*Nsyn*sizeof(complex));
		cdata    = (complex *)malloc(nfreq*nxs*sizeof(complex));
		rdata    = (float *)malloc(nxs*optn*sizeof(float));

		for (l = 0; l < Nsyn; l++) {
			for(i = 0; i < nxs; i++) {
				for (j = 0; j < nts; j++) 
					rdata[i*optn+j] = syndata[l*size+i*nts+j];
				for (j = nts; j < optn; j++)
					 rdata[i*optn+j] = 0.0;
			}

			rcmfft(rdata, cdata, optn, nxs, optn, nfreq, -1);
			for(i = 0; i < nxs; i++) {
				for (j = iomin, m = 0; j < iomax; j++, m++) {
					syncdata[l*Niom*nxs+m*nxs+i].r = cdata[i*nfreq+j].r;
					syncdata[l*Niom*nxs+m*nxs+i].i = mode*cdata[i*nfreq+j].i;
				}
			}
		}
		free(cdata);
		free(rdata);

		p = (float *) &syndata[0];
		for (i = 0; i < size*Nsyn; i++) *p++ = 0.0;

		first = 1;
		if (verbose >= 2) {
			vmess("Operators are transformed to x-w");
			vmess("nt = %d nts = %d, optn = %d", nt, nts, optn);
		}
	}

/*============= FFT of shot to frequency domain =============*/

	t1 = wallclock_time();

	shotcdata = (complex *)calloc(nx*nfreq,sizeof(complex));
#pragma omp parallel default(none) \
 shared(shotdata, nx, optn) \
 shared(shotcdata, nfreq) \
 private(j, i, cdata)
	{ /* start of parallel region */
	cdata     = (complex *)malloc(nfreq*sizeof(complex));
#pragma omp for
	for(i = 0; i < nx; i++) {
#pragma omp critical 
{
		rc1fft(&shotdata[i*optn], cdata, optn, -1);
}
		for (j = 0; j < nfreq; j++) {
			shotcdata[j*nx+i] = cdata[j];
		}
	}
	free(cdata);
	} /* end of parallel region */

	tfft += wallclock_time() - t1;
//	xt2wx(shotdata, shotcdata, optn, nx, optn, nx);

/*================ SYNTHESIS ================*/

	ixsrc = NINT((xsrc - fxs)/dxs);

	if (abs(xrcv[0]-xsrc) > 0.5*nx*dx) { iox = 0; inx = nx-off; }
	else { iox = off; inx = nx; }

	t0 = wallclock_time();

#pragma omp parallel default(none) \
 shared(syndata, dx, npe, Niom) \
 shared(shotcdata, Nsyn, reci, xrcv, xsrc, xsyn, fxs, nxs, dxs) \
 shared(nx, ixa, ixb, dxsrc, nmo, iox, inx, k, nfreq, iomin, iomax) \
 shared(syncdata, size, nts, optn, scl, ixsrc) \
 private(l, ixsyn, x0, x1, ix, dosrc, j, m, i, ixrcv, sum, rdata, tmp, ts, to)
	{ /* start of parallel region */
	sum   = (complex *)malloc(nfreq*sizeof(complex));
	rdata = (float *)calloc(optn,sizeof(float));
#pragma omp for 
	for (l = 0; l < Nsyn; l++) {
		ixsyn = NINT((xsyn[l] - fxs)/dxs);

		if (ixa || ixb) { 
			if (reci == 0) {
				x0 = xsyn[l]-ixb*dxsrc; 
				x1 = xsyn[l]+ixa*dxsrc; 
				if ((xsrc < x0) || (xsrc > x1)) continue;
				ix = NINT((xsrc-x0)/dxsrc);
				dosrc = 1;
			}
			else if (reci == 1) {
				x0 = xsyn[l]-ixb*dxs; 
				x1 = xsyn[l]+ixa*dxs; 
				if (((xsrc < x0) || (xsrc > x1)) && 
					(xrcv[0] < x0) && (xrcv[nx-1] < x0)) continue;
				if (((xsrc < x0) || (xsrc > x1)) && 
					(xrcv[0] > x1) && (xrcv[nx-1] > x1)) continue;
				if ((xsrc < x0) || (xsrc > x1)) dosrc = 0;
				else dosrc = 1;
				ix = NINT((xsrc-x0)/dxs);
			}
			else if (reci == 2) {
				if (NINT(dxsrc/dx)*dx != NINT(dxsrc)) dx = dxs;
				x0 = xsyn[l]-ixb*dx; 
				x1 = xsyn[l]+ixa*dx; 
				if ((xrcv[0] < x0) && (xrcv[nx-1] < x0)) continue;
				if ((xrcv[0] > x1) && (xrcv[nx-1] > x1)) continue;
			}
		}
		else { 
			ix = k; 
			x0 = fxs; 
			x1 = fxs+dxs*nxs;
			dosrc = 1;
		}
		if (reci == 1 && dosrc) ix = NINT((xsrc-x0)/dxs);

		if (reci < 2 && dosrc) {
			for (j = 0; j < nfreq; j++) sum[j].r = sum[j].i = 0.0;
			for (j = iomin, m = 0; j < iomax; j++, m++) {
				for (i = iox; i < inx; i++) {
					ixrcv = NINT((xrcv[i]-fxs)/dxs);
					tmp = syncdata[l*Niom*nxs+m*nxs+ixrcv];
					sum[j].r += shotcdata[j*nx+i].r*tmp.r +
								shotcdata[j*nx+i].i*tmp.i;
					sum[j].i += shotcdata[j*nx+i].i*tmp.r -
								shotcdata[j*nx+i].r*tmp.i;
				}
				if (nmo) {
					ts = syncdata[l*Niom*nxs+m*nxs+ixsrc];
					to = syncdata[l*Niom*nxs+m*nxs+ixsyn];
					tmp.r = sum[j].r*to.r - sum[j].i*to.i;
					tmp.i = sum[j].i*to.r + sum[j].r*to.i;
					sum[j].r = tmp.r*ts.r + tmp.i*ts.i;
					sum[j].i = tmp.i*ts.r - tmp.r*ts.i;
				}
			}
#pragma omp critical
{
			cr1fft(sum, rdata, optn, 1);
}
/* dx = receiver distance */
			for (j = 0; j < nts; j++) 
				syndata[l*size+ix*nts+j] += rdata[j]*scl*dx;
		}

		if (reci == 1 || reci == 2) {
			for (j = 0; j < nfreq; j++) sum[j].r = sum[j].i = 0.0;
			for (i = iox; i < inx; i++) {
				if ((xrcv[i] < x0) || (xrcv[i] > x1)) continue;
				if (reci == 1) ix = NINT((xrcv[i]-x0)/dxs);
				else ix = NINT((xrcv[i]-x0)/dx);

				for (j = iomin, m = 0; j < iomax; j++, m++) {
					tmp = syncdata[l*Niom*nxs+m*nxs+ixsrc];
					sum[j].r = shotcdata[j*nx+i].r*tmp.r +
							   shotcdata[j*nx+i].i*tmp.i;
					sum[j].i = shotcdata[j*nx+i].i*tmp.r -
							   shotcdata[j*nx+i].r*tmp.i;
				}
#pragma omp critical
{
				cr1fft(sum, rdata, optn, 1);
}
/* dxsrc = source distance */
				for (j = 0; j < nts; j++) 
					syndata[l*size+ix*nts+j] += rdata[j]*scl*dxsrc;
			}
		}
	} /* end of Nsyn loop */

	free(sum);
	free(rdata);

    } /* end of parallel region */
	t += wallclock_time() - t0;
	if (k == nshots-1) {
#pragma omp single 
{ 
#ifdef __OPENMP
    npe   = omp_get_num_threads();
#endif
}
		vmess("OMP: parallel region = %f seconds (%d threads)", t, npe);
		vmess("FFT: parallel region = %f seconds (%d threads)", tfft, npe);
	}


	free(shotcdata);

	return;
}
