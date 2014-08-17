#include <DELPHI_IOc.h>
static char rcsid[] = "$Id: syn2d.c,v 1.9 1998/10/06 12:05:54 seistool Exp $";

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
"   key=fldr ................. input data sorting key",
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
"   key=fldr ................. input data sorting key",
"   verbose=0 ................ silent option; >0 displays info",
" ",
" The nmo option determines the display of the calculated CFP-gather.",
"     nmo = 0: displays CFP-gather at operator times,",
"     nmo = 1: displays CFP-gather corrected with operator times,",
" ",
"  Note that if ixa=0 and ixb=0 all shots are used.",
" ",
" author  : Jan Thorbecke : 4-9-1995 (jan@delphi.tn.tudelft.nl)",
" product : CFP",
" ",
NULL};
/**************** end self doc ***********************************/

void main (int argc, char **argv)
{
	intn	seqnr[MAX_KEYS];
	int	type;
	int		i, j, k, l, ret, nshots, Nsyn, nt, nx, nts, nxs, more;
	int		size, n1, n2, nkeys, ntap, tap, di, ixrcv, ixsrc, off, optn;
	int		nxmax, ntmax, reci, mode, nmo, ixa, ixb, n2out, verbose;
	float	fmin, fmax, *taper, t0, t1, t2, t3, fxf, dxf, fxs2, xsrc, *xrcv;
	float	*shotdata, d1, d2, f1, f2, fts, fxs, ft, fx, *etap, *xsyn, dxsrc;
	float   *syndata, *tmpdata, dt, dx, dts, dxs, scl, tsyn, tread, alpha, mem;
	float max, *zsyn, scel;
	char	*file_syn, *file_shot, *file_cfp, *keys[MAX_KEYS], *key, *file_fresnel;
	segyhdr	*hdrs, *hdrs_in, *hdrs_out;

	initargs(argc, argv);
	requestdoc(1);

	if (!getparstring("file_shot", &file_shot)) file_shot = NULL;
	if (!getparstring("file_syn", &file_syn)) file_syn = NULL;
	if (!getparstring("file_fresnel", &file_fresnel)) file_fresnel = NULL;
	if (!getparint("verbose", &verbose)) verbose = 0;
	if (file_syn == NULL && file_shot == NULL) 
		saerr("file_syn and file_shot cannot be both input pipe");
	if (!getparstring("file_cfp", &file_cfp)) {
		if (verbose) sawarn("parameter file_cfp not found, assume pipe");
		file_cfp = NULL;
	}
	if (!getparint("nshots", &nshots)) saerr("nshots should be given");
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
	if (nmo && mode > 0) saerr("nmo = 1 cannot be used with mode = 1");
	if (reci && ntap) sawarn("tapering influences the reciprocal result");

/*================ Reading synthesis operator(s) ================*/

	ret = open_file(file_syn, GUESS_TYPE, DELPHI_READ);
	if (ret < 0) saerr("error in opening file %s", file_syn);
	ret = get_dims(file_syn, &n1, &n2, &type);
	if (ret >= 0) {
		if (!getparint("ntmax", &ntmax)) ntmax = n1;
		if (!getparint("nxmax", &nxmax)) nxmax = n2;
		if (verbose>=2 && (ntmax!=n1 || nxmax!=n2))
		    samess("dimensions overruled: %d x %d",ntmax,nxmax);
	}
	else {
		if (!getparint("ntmax", &ntmax)) ntmax = 1024;
		if (!getparint("nxmax", &nxmax)) nxmax = 512;
		if (verbose>=2) samess("dimensions used: %d x %d",ntmax,nxmax);
	}
	size    = ntmax * nxmax;
	tmpdata = alloc1float(size);
	hdrs    = (segyhdr *) malloc(nxmax*sizeof(segyhdr));

	ret = read_data(file_syn,tmpdata,size,&n1,&n2,&f1,&f2,&d1,&d2,&type,hdrs);
	if (ret < 0) saerr("error in reading data from file %s", file_syn);
	if (hdrs[0].scalco < 0) scl = 1.0/fabs(hdrs[0].scalco);
	else if (hdrs[0].scalco == 0) scl = 1.0;
	else scl = hdrs[0].scalco;

	if (hdrs[0].scalel < 0) scel = 1.0/fabs(hdrs[0].scalel);
	else if (hdrs[0].scalel == 0) scel = 1.0;
	else scel = hdrs[0].scalel;

	taper = alloc1float(n2);
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

	etap = alloc1float(n1);
	for (j = 0; j < n1; j++) etap[j] = exp(alpha*j*d1);

	syndata = (float *)malloc(n1*n2*sizeof(float));
	xsyn    = (float *)malloc(sizeof(float));
	zsyn    = (float *)malloc(sizeof(float));
	Nsyn    = 0;
	while (ret >= 0) {
		if (verbose>=2) disp_info(file_syn,n1,n2,f1,f2,d1,d2,type);
		syndata = (float *)realloc(syndata, (Nsyn+1)*n1*n2*sizeof(float));
		if (syndata == NULL) saerr("Out of memory for syndata array!");
		xsyn = (float *)realloc(xsyn, (Nsyn+1)*sizeof(float));
		if (xsyn == NULL) saerr("Out of memory for xsyn array!");
		zsyn = (float *)realloc(zsyn, (Nsyn+1)*sizeof(float));
		if (zsyn == NULL) saerr("Out of memory for zsyn array!");

		xsyn[Nsyn] = hdrs[0].sx*scl;
		zsyn[Nsyn] = hdrs[0].selev*scel;
		for (i = 0; i < n2; i++) {
			for (j = 0; j < n1; j++) 
				syndata[Nsyn*n1*n2+i*n1+j] = tmpdata[i*n1+j]*taper[i]*etap[j];
		}
		Nsyn += 1;
		ret = read_data(file_syn, tmpdata, size, &n1, &n2, &f1, &f2, 
			&d1, &d2, &type, hdrs);
	}
	ret = close_file(file_syn);
	if (ret < 0) sawarn("err %d on closing input file",ret);

	nxs = n2; nts = n1;
	dxs = d2; dts = d1;
	fxs = f2; fts = f1;
	if (hdrs[0].gx != 0 || hdrs[1].gx != 0 ) fxs = hdrs[0].gx*scl;
	fxs2 = fxs + (float)(nxs-1)*dxs;
	dxf = (hdrs[nxs-1].gx - hdrs[0].gx)*scl/(float)(nxs-1);
	if (NINT(dxs*1e3) != NINT(fabs(dxf)*1e3)) {
		samess("dx in hdr.d1 (%.3f) and hdr.gx (%.3f) not equal",d2, dxf);
		if (dxf != 0) dxs = fabs(dxf);
		samess("dx in operator => %f", dxs);
	}

	free1float(tmpdata);
	free1float(taper);
	free1float(etap);
	free(hdrs);

/*================ Reading first shot record ================*/

	ret = open_file(file_shot, GUESS_TYPE, DELPHI_READ);
	if (ret < 0 ) saerr("error in opening file %s", file_shot);
	ret = get_dims(file_shot, &n1, &n2, &type);
	if (ret >= 0) {
		if (!getparint("ntmax", &ntmax)) ntmax = n1;
		if (!getparint("nxmax", &nxmax)) nxmax = n2;
		if (verbose>=2 && (ntmax!=n1 || nxmax!=n2))
		    samess("dimensions overruled: %d x %d",ntmax,nxmax);
	}
	else {
		if (!getparint("ntmax", &ntmax)) ntmax = 1024;
		if (!getparint("nxmax", &nxmax)) nxmax = 512;
		if (verbose>=2) samess("dimensions used: %d x %d",ntmax,nxmax);
	}
	size     = ntmax * nxmax;
	xrcv     = alloc1float(nxmax);
	tmpdata  = alloc1float(size);
	hdrs_in  = (segyhdr *) malloc(nxmax*sizeof(segyhdr));
	if (tmpdata == NULL || hdrs_in == NULL)
		saerr("memory allocation error for input data");

	if (!getparstring("key", &key)) {
		ret = get_sort(file_shot);
		if (ret < 0) key = "fldr";
		else key = getkey(ret);
	}
	if (verbose) samess("input sorting key is %s",key);
	set_sukey(key);

	read_data(file_shot,tmpdata,size,&nt,&nx,&ft,&fx,&dt,&dx,&type,hdrs_in);
	if (verbose>=2) disp_info(file_shot,nt,nx,ft,fx,dt,dx,type);

	taper = alloc1float(nx);
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
	etap = alloc1float(optn);
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
		samess("dx in hdr.d1 (%.3f) and hdr.gx (%.3f) not equal",dx, dxf);
		if (dxf != 0) dx = fabs(dxf);
		else saerr("gx hdrs not set");
		samess("dx used => %f", dx);
	}

/*=========== Reading second shot record (for source sampling) ===========*/

	hdrs    = (segyhdr *) malloc(nxmax*sizeof(segyhdr));
	if (hdrs == NULL)
		saerr("memory allocation error for input data");

	ret = read_data(file_shot,tmpdata,size,&n1,&n2,&f1,&f2,&d1,&d2,&type,hdrs);
	if (ret < 0) {
		sawarn("only one shot record is available");
		dxsrc = dx;
		more = 0;
	}
	else {
		if (verbose>=3) disp_info(file_shot,n1,n2,f1,f2,d1,d2,type);
		dxsrc = (float)hdrs[0].sx*scl - fxf;
		if (dxsrc == 0) {
			sawarn("sx hdrs are not filled in!!");
			dxsrc = dx;
		}
		more = 1;
	}

/*================ Check the size of the files ================*/

	if (NINT(dxsrc/dx)*dx != NINT(dxsrc)) {
		sawarn("source (%.2f) and receiver step (%.2f) don't match",dxsrc,dx);
		if (reci == 2) sawarn("step used from operator (%.2f) ",dxs);
	}
	di = NINT(dxf/dxs);
	if ((NINT(di*dxs) != NINT(dxf)) && verbose) 
		sawarn("dx in receiver (%.2f) and operator (%.2f) don't match",dx,dxs);
	if (nt != nts) 
		samess("Time samples in shot and synthesis operator are not equal");
	if (verbose) {
		samess("Number of synthesis operators  = %d", Nsyn);
		samess("number of shots                = %d", nshots);
		samess("first model position           = %.2f", fxs);
		samess("last model position            = %.2f", fxs2);
		samess("first source position fxf      = %.2f", fxf);
		samess("source distance dxsrc          = %.2f", dxsrc);
		samess("last source position           = %.2f", fxf+(nshots-1)*dxsrc);
		samess("receiver distance     dxf      = %.2f", dxf);
		samess("direction of increasing traces = %d", di);
	}

/*================ initializations ================*/

	if (ixa || ixb) n2out = ixa + ixb + 1;
	else if (reci) n2out = nxs;
	else n2out = nshots;
	mem = Nsyn*n2out*nts*sizeof(float)/1048576.0;
	if (verbose) {
		samess("number of output traces        = %d", n2out);
		samess("number of output samples       = %d", nts);
		samess("Size of output data            = %.1f Mb", mem);
	}

	t0   = cputime();
	tsyn = tread = 0.0;
	ret  = 0;
	k    = 0;
	size = ntmax*nxmax;
	keys[0]  = (char *) malloc(MAX_KEY_LENGTH);
	nkeys    = 1;
	keys[0]  = SA_ASG;

	if (file_fresnel) {
		ret = open_file(file_fresnel, GUESS_TYPE, DELPHI_CREATE);
	}

/*================ loop over all shot records ================*/

	while (ret >= 0) {
		t1    = cputime();
		xsrc  = (float)hdrs_in[0].sx*scl;
		for (i = 0; i < nx; i++) xrcv[i] = (float)hdrs_in[i].gx*scl;
		if (verbose>=2) {
			samess("source position:     %.2f", xsrc);
			samess("receiver positions:  %.2f <--> %.2f", xrcv[0], xrcv[nx-1]);
		}

		if ((NINT(xsrc-fxs2) > 0) || (NINT(xrcv[nx-1]-fxs2) > 0) ||
			(NINT(xrcv[nx-1]-fxs) < 0) || (NINT(xsrc-fxs) < 0) || 
			(NINT(xrcv[0]-fxs) < 0) || (NINT(xrcv[0]-fxs2) > 0) ) {
			sawarn("source/receiver positions are outside synthesis model");
			sawarn("CFP calculation is stopped at gather %d", k);
			samess("xsrc = %.2f xrcv_1 = %.2f xrvc_N = %.2f", 
			xsrc, xrcv[0], xrcv[nx-1]);
			ret = close_file(file_shot);
			if (ret < 0) sawarn("err %d on closing input file",ret);
			break;
		}

		synthesis(shotdata, syndata, nx, nt, nxs, nts, dt, xsyn, Nsyn, 
			xrcv, xsrc, fxs, dxs, dxsrc, dx, ixa, ixb, di, fmin, fmax, mode, 
			reci, k, off, nmo, nshots, verbose);

		t3 = cputime();
		tsyn +=  t3 - t1;

		if (file_fresnel) {
			for (i = 0; i < nx; i++) {
				hdrs_in[i].ns   = optn;
			}
			ret = write_data(file_fresnel,(float *)&shotdata[0], optn, nx, 
				ft, fx, dt, dx, type, hdrs_out);
			if (ret < 0 ) saerr("error on writing output file.");
		}

		if (k == 0 && more) {
			nt = n1; nx = n2;
			for (i = 0; i < nx; i++) {
				hdrs_in[i] = hdrs[i];
				for (j = 0; j < nt; j++) 
					shotdata[i*optn+j] = tmpdata[i*n1+j];
				for (j = nt; j < optn; j++) 
					shotdata[i*optn+j] = 0.0;
			}
			free(hdrs);

			t2 = cputime();
			if (verbose>=3) {
				samess("CPU-time one gather      = %.3f s.", t2-t1);
				samess("with CPU-time synthesis  = %.3f", tsyn);
				samess("and CPU-time read data   = %.3f", t2-t3);
			}
		}
		else {
			ret = read_data(file_shot, tmpdata, size, 
				&nt, &nx, &ft, &fx, &dt, &d2, &type, hdrs_in);

			for (i = 0; i < nx; i++) {
				for (j = 0; j < nt; j++) shotdata[i*optn+j] = tmpdata[i*nt+j];
				for (j = nt; j < optn; j++) shotdata[i*optn+j] = 0.0;
			}
		}
		k++;
		if (verbose) samess("*** Shot gather %d processed ***", k);

		if (ret < 0 ) {
			ret = close_file(file_shot);
			if (ret < 0) sawarn("err %d on closing input file",ret);
			if (verbose) samess("end of data reached");
			break;
		}
		for (i = 0; i < nx; i++) {
			for (j = 0; j < nt; j++) shotdata[i*optn+j] *= (taper[i]*etap[j]);
			for (j = nt; j < optn; j++) shotdata[i*optn+j] = 0.0;
		}

		if (verbose>=3) disp_info(file_shot,nt,nx,f1,f2,d1,d2,type);
		tread += cputime() - t3;
	}

	if (file_fresnel) {
		ret = close_file(file_fresnel);
		if (ret < 0) saerr("err %d on closing fresnel output file",ret);
	}

	t2 = cputime();
	if (verbose) {
		samess("Total CPU-time synthesis = %.3f", t2-t0);
		samess("with CPU-time synthesis  = %.3f", tsyn);
		samess("and CPU-time read data   = %.3f", tread);
	}

	if (k != nshots) {
		sawarn("number of shots given not equal to shots read in (%d)", k);
		if (ixa == 0 && ixb == 0 && reci == 0) n2out = k;
	}

/*================ write output files and free memory ================*/

	n1 = nts; n2 = n2out;
	f1 = ft; f2 = fxs;
	d1 = dt;
	if (reci == 0) d2 = dxsrc;
	else if (reci == 1) d2 = dxs;
	else if (reci == 2) d2 = dx;

	hdrs_out = (segyhdr *) malloc(n2*sizeof(segyhdr));
	if (hdrs_out == NULL) saerr("allocation for hdrs_out");
	gen_hdrs(hdrs_out,n1,n2,f1,f2,d1,d2,TRID_TIME);
	size  = nxs*nts;

	free1float(etap);
	etap = alloc1float(n1);
	for (j = 0; j < n1; j++) etap[j] = exp(-alpha*j*dt);

	ret = open_file(file_cfp, GUESS_TYPE, DELPHI_CREATE);
	if (ret < 0 ) saerr("error on creating output file");
	for (l = 0; l < Nsyn; l++) {
		seqnr[0] = l+1;
		if (ixa || ixb) f2 = xsyn[l]-ixb*d2;
		else {
			if (reci) f2 = fxs;
			else f2 = fxf;
		}

		for (i = 0; i < n2; i++) {
			hdrs_out[i].fldr   = l+1;
			hdrs_out[i].trwf   = n2out;
			hdrs_out[i].timbas = getindex("fldr")+1;;
			hdrs_out[i].scalco = -1000;
			hdrs_out[i].sx = NINT(xsyn[l]*1000);
			hdrs_out[i].gx = NINT(1000*(f2+i*d2));
			hdrs_out[i].offset = (long)NINT((f2+i*d2) - xsyn[l]);
			hdrs_out[i].scalel = -1000;
			hdrs_out[i].selev  = NINT(zsyn[l]*1000);
			hdrs_out[i].sdepth = NINT(zsyn[l]*1000);
			for (j = 0; j < n1; j++) syndata[l*size+i*n1+j] *= etap[j];
		}
		ret = set_keys(keys, seqnr, nkeys);
		if (ret < 0 ) sawarn("error on writing keys.");
		ret = set_axis(SA_AXIS_TIME, SA_AXIS_X);
		if (ret < 0 ) saerr("error on writing axis.");
		ret = write_data(file_cfp,(float *)&syndata[l*size], n1, n2, 
			f1, f2, d1, d2, type, hdrs_out);
		if (ret < 0 ) saerr("error on writing output file.");
	}
	ret = close_file(file_cfp);
	if (ret < 0) saerr("err %d on closing output file",ret);

	free(hdrs_in);
	free(hdrs_out);
	free1float(tmpdata);

	exit(0);
}

void synthesis(float *shotdata, float *syndata, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float xsrc, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int di, float fmin, float fmax, int mode, int reci, int k, int off, int nmo, int nshots, int verbose)
{
	static int first=0, iomin, iomax, nfreq, optn, Niom, size, iox, inx;
	static complex *syncdata;
	static float scl;
	int 	i, j, l, m, ixsrc, ixsyn, ix, ixrcv, dosrc;
	float	*rdata, *p, **dum, x0, x1;
	static float t0, t;
	complex *sum, *cdata, *shotcdata, tmp, ts, to;
#ifdef SGI
	int     mype, np = mp_numthreads(), npe;
#endif


/*============= Transform synthesis operator(s), only one time =============*/

	if (first == 0) {
		size  = nxs*nts;
		optn  = optncr(MAX(nt, nts));
		scl   = 1.0/optn;
		nfreq = (optn+2)/2;
		iomin = (int)MIN((fmin*optn*dt), nfreq-1);
		iomin = MAX(iomin, 1);
		iomax = MIN((int)(fmax*optn*dt), nfreq-1);
		Niom  = iomax - iomin + 1;
		if (verbose>=3) 
			samess("iomin = %d iomax =%d Niom = %d",iomin,iomax,Niom);
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
				for (j = iomin, m = 0; j <= iomax; j++, m++) {
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
			samess("Operators are transformed to x-w");
			samess("nt = %d nts = %d, optn = %d", nt, nts, optn);
		}
	}

/*============= FFT of shot to frequency domain =============*/

	shotcdata = (complex *)malloc(nx*nfreq*sizeof(complex));
	xt2wx(shotdata, shotcdata, optn, nx, optn, nx);

/*================ SYNTHESIS ================*/

	ixsrc = NINT((xsrc - fxs)/dxs);

	if (abs(xrcv[0]-xsrc) > 0.5*nx*dx) { iox = 0; inx = nx-off; }
	else { iox = off; inx = nx; }

#ifdef SGI
	t0 = cputime();
	npe = MIN(np,Nsyn);
	mpc_set_numthreads(npe);
	mp_setup();
#pragma parallel
#pragma shared(syndata, dx)
#pragma byvalue(shotcdata, Nsyn, reci, xrcv, xsrc, xsyn, fxs, nxs, dxs)
#pragma byvalue(nx, ixa, ixb, dxsrc, nmo, iox, inx, k, nfreq, iomin, iomax)
#pragma byvalue(syncdata, size, nts, optn, scl, ixsrc)
#pragma local(l, ixsyn, x0, x1, ix, dosrc, j, m, i, ixrcv, sum, rdata, tmp, ts, to, mype)
	{ /* start of parallel region */
#endif
	sum   = (complex *)malloc(nfreq*sizeof(complex));
	rdata = (float *)calloc(optn,sizeof(float));
#ifdef SGI
#pragma pfor iterate(l=0;Nsyn;1)
#endif
	for (l = 0; l < Nsyn; l++) {
#ifdef SGI
		mype = mp_my_threadnum();
#endif
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
			for (j = iomin, m = 0; j <= iomax; j++, m++) {
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
#ifdef SGI
#pragma critical
    {
#endif
			cr1fft(sum, rdata, optn, 1);
			for (j = 0; j < nts; j++) 
				syndata[l*size+ix*nts+j] += rdata[j]*scl;
#ifdef SGI
    }
#endif
		}

		if (reci == 1 || reci == 2) {
			for (j = 0; j < nfreq; j++) sum[j].r = sum[j].i = 0.0;
			for (i = iox; i < inx; i++) {
				if ((xrcv[i] < x0) || (xrcv[i] > x1)) continue;
				if (reci == 1) ix = NINT((xrcv[i]-x0)/dxs);
				else ix = NINT((xrcv[i]-x0)/dx);

				for (j = iomin, m = 0; j <= iomax; j++, m++) {
					tmp = syncdata[l*Niom*nxs+m*nxs+ixsrc];
					sum[j].r = shotcdata[j*nx+i].r*tmp.r +
							   shotcdata[j*nx+i].i*tmp.i;
					sum[j].i = shotcdata[j*nx+i].i*tmp.r -
							   shotcdata[j*nx+i].r*tmp.i;
				}
#ifdef SGI
#pragma critical
    {
#endif
				cr1fft(sum, rdata, optn, 1);
				for (j = 0; j < nts; j++) 
					syndata[l*size+ix*nts+j] += rdata[j]*scl;
#ifdef SGI
    }
#endif
			}
		}
	}

	free(sum);
	free(rdata);

#ifdef SGI
    } /* end of parallel region */
	t += cputime() - t0;
	if (k == nshots-1) fprintf(stderr,"SGI: parallel region = %f seconds (%d CPUS's)\n", t, npe);
#endif


	free(shotcdata);

	return;
}
