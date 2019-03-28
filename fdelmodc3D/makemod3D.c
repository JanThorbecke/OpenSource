#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif
#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

float **alloc2float (size_t n1, size_t n2);
void free2float (float **p);

void grid3D(float **gridcp, float **gridcs, float **gridro, long *zp, 
	  float **cp, float **cs, float **ro, long minx, long maxx, long miny, long maxy,
      long optgrad, float gradlen, float gradcp, float gradcs, float gradro, float dx, 
	  float dy, float dz, long nz);

void gridabove(float **gridcp, float **gridcs, float **gridro, long *zp, float **cp, float **cs, float **ro, long minx, long maxx, long optgrad, float gradlen, float gradcp, float gradcs, float gradro, float dx, float dz, long nz);

void plotexample();

void name_ext(char *filename, char *extension);

void interpolation(float *x, float *z, long nxp, long nx, long poly, long *pminx, long *pmaxx, float dx, float **cp, float **cs, float **ro, long nvel, float *interface);

void linearint(long *zp, long minx, long maxx, float dz, float *interface);

void sinusint(long *zp, long minx, long maxx, float dz, float *interface, float dx, float ampl, float wavel);

void roughint(long *zp, long minx, long maxx, float dz, float *interface, float ampl, float beta, float seed);

void fractint(long *zp, long minx, long maxx, float dx, float dz, float *interface, float Nsin, float ampl, float D, float k0, float b, float seed);

void elipse(float *x, float *z, long nxp, float dx, float dz, float **gridcp, float **gridcs, float **gridro, float **cp, float **cs, float **ro, float *interface, long *zp, long nz, long nx, float r1, float r2, float gradcp, float gradcs, float gradro);

void diffraction(float *x, float *z, long nxp, float dx, float dz, float **gridcp, float **gridcs, float **gridro, float **cp, float **cs, float **ro, float *interface, long *zp, long nx, long diffrwidth, long type);

void randdf(float *x, float *z, long nxp, float dx, float dz, float **gridcp, float **gridcs, float **gridro, float **cp, float **cs, float **ro, float *interface, long *zp, long nx, float sizex, float sizez, long ndiff, long diffrwidth, long type);

/*********************** self documentation **********************/
char *sdoc[] = {
  " 								",
  " makemod - gridded subsurface model builder",
  " 								",
  " makemod file_base= cp0= sizex= sizez= dx= dz= [optional parameters] ",
  " 							        ",
  " Required parameters:",
  " ",
  "   file_base= ............... base name of the output file(s).",
  "   cp0= ..................... Cp for background medium.",
  "   sizex= ................... x-size of the model in meters.",
  "   sizez= ................... z-size of the model in meters.",
  "   dx= ...................... grid distance in x in meters.",
  "   dz= ...................... grid distance in z in meters.",
  " ",
  " Optional parameters:",
  " ",
  "   orig=0,0 (x,z) ........... origin at the top left corner of the model.",
  " MEDIUM ",
  "   cs0=0 .................... Cs for background medium (0 is none).",
  "   ro0=0 .................... Rho for background medium (0 is none).",
  "   gradt=1 .................. type of boundary gradient (1=linear, 2=cos)",
  "   cp=none .................. P-wave velocities below the interface",
  "   cs=none .................. S-wave velocities below the interface",
  "   ro=none .................. Density below the interface",
  "   above=0 .................. define model below interface",
  "                              =1: define model above interface",
  " INTERFACE ",
  "   intt=none ................ Type of interface",
  "   var=none ................. variables to describe the interface",
  "   grad=0.0 ................. gradient(m) of the boundary",
  "   gradunit=0 ............... gradient unit (m/s per gradunit)",
  "   gradcp=0.0 ............... gradient(m/s per grad-unit) in the layer",
  "   gradcs=0.0 ............... gradient(m/s per grad-unit) in the layer",
  "   gradro=0.0 ............... gradient(kg/m3 per grad-unit) in the layer",
  "   poly=0 ................... polynominal interpolation through (x,z) points",
  "   x=none ................... x-positions for the interface",
  "   z=none ................... z-positions for the interface",
  "   dtype=0 .................. diffractor type for diffr and randdf",
  " OUTPUT ",
  "   writeint=0 ............... interfaces as function of x (ext: _int)",
  "   supersmooth=0 ............ samples at dx/4 and apply 5x5 smoothing operator(1)",
  "   sigma=0.8 ................ smoothing value, higher values gives more smoothing",
  "   reflectivity=0 ........... compute reflectivity (ext: _rfl)",
  "   rayfile=0 ................ interfaces as function of x in ASCII file.mod",
  "   skip=5 ................... number of x position to skip in file_ray",
  "   example=0 ................ makes an example parameter file",
  "   verbose=0 ................ silent option; >0 display info",
  " ",
  "   Options for intt:",
  "         - def       = default interface through the points(Xi, Zi)",
  "         - sin       = sinus shaped interface",
  "         - rough     = rough interface with beta(smoothness)",
  "         - fract     = cosinus fractal shaped interface",
  "         - random    = define random velocities in layer",
  "         - elipse    = define elipse shaped body",
  "         - diffr     = point diffractions",
  "         - randdf    = define random diffractors ",
  "   Options for var in case of intt =:",
  "         - sin(2)    = wavelength,amplitude",
  "         - rough(3)  = amplitude,beta,seed",
  "         - fract(6)  = Nsinus,amplitude,dim,k0,freqscale,seed",
  "         - random(1) = min-max variation around cp",
  "         - elipse(2) = r1, r2: vertical and horizontal radius",
  "         - diffr(1)  = width of each point, type(optional)",
  "         - randdf(2) = number of points, width of each point",
  "   Options for poly in default interface:",
  "         - 0         = linear",
  "         - 1         = polynomal",
  "         - 2         = cubic spline",
  "   Options for dtype value in var=width,dtype for diffr:",
  "         - -1        = random (0, 1, or 2) diffractor type",
  "         - 0         = cubic diffractor",
  "         - 1         = diamond diffractor",
  "         - 2         = circular diffractor",
  "   Option for gradunit, gradient unit per layer:",
  "         - 0         = gradient unit per layer is m/s per dz (default)",
  "         - 1         = gradient unit per layer is m/s per m",
  "  (1) Reference: Zeng and West: Geophysics 1996 ",
  " ",
  " makemod builds a gridded subsurface file which can be used in a migration",
  " or finite difference program. The gridded model is stored in files with",
  " extensions _cp, _cs, _ro. The extensions _int and .mod are used for the ",
  " interface files. The output format of the file(s) depends on the .(dot)",
  " extension of file_base.",
  " ",
  " author  : Jan Thorbecke : 18-01-1994 (janth@xs4all.nl)",
  " product : Originates from DELPHI software",
  "                         : revision 2010",
  " ",
  NULL};
/**************** end self doc ***********************************/

int main(int argc, char **argv)
{
  FILE *fpint, *fpcp, *fpcs, *fpro;
  long   example, verbose, writeint, nb;
  long	above, diffrwidth, dtype, reflectivity, supersmooth;
  long   Ngp, Ngs, Ngr, Np, Ns, Nr, Ng, Ni, Nv, Nvi, No, Noi;
  long   jint, jcount, j, l, ix, iy, iz, nx, ny, nz, nxp, nyp, nzp, *zp, nmaxx, nminx, nmaxy, nminy, optgrad, poly, gradt; 
  long	ncp, nro, ncs, nvel, skip, rayfile, store_int;
  long lseed;
  size_t  nwrite;
  float *data_out, orig[3], cp0, cs0, ro0, Z1, Z2, sigma;
  float	*x, *y, *z, *var, *interface, **inter;
  float back[3], sizex, sizey, sizez, dx, dy, dz;
  float **cp ,**cs, **ro, aver, gradlen, gradcp, gradcs, gradro;
  
  /* Gradient unit flag    */
  /* ------------------    */
  /* - 0 Unit : m/s per dz */
  /* - 1 Unit : m/s per m  */
  long gradunit;
  /* Number of Z-reference points (one per layer) */
  long Nzr=0;
  
  float   **gridcp, **gridcs, **gridro;
  segy    *hdrs;
  char    *file, intt[10], *file_base, filename[150];
  
  initargs(argc, argv);
  requestdoc(1);
  
  if (!getparint("example", &example)) example=0;
  else { plotexample(); exit(0); }
  
  if(getparstring("file",&file_base)) 
    vwarn("parameters file is changed to file_base");
  else {
    if(!getparstring("file_base",&file_base))
      verr("file_base not specified.");
  }
  if(getparfloat("back", back)) {
    vwarn("parameters back is not used anymore");
    vwarn("it has changed into cp0 (ro0,cs0 are optional)");
    nb = countparval("back");
    if (nb == 1) {
      vwarn("The new call should be cp0=%.1f",back[0]);
      cp0 = back[0];
    }
    if (nb == 2) {
      vwarn("The new call should be cp0=%.1f",back[0]);
      vwarn("                       ro0=%.1f",back[1]);
      cp0 = back[0];
      ro0 = back[1];
    }
    if (nb == 3) {
      vwarn("The new call should be cp0=%.1f",back[0]);
      vwarn("                       cs0=%.1f",back[1]);
      vwarn("                       ro0=%.1f",back[2]);
      cp0 = back[0];
      cs0 = back[1];
      ro0 = back[2];
    }
    vmess("Don't worry everything still works fine");
  }
  else {
    if(!getparfloat("cp0", &cp0)) verr("cp0 not specified.");
    if(!getparfloat("cs0", &cs0)) cs0 = -1;
    if(!getparfloat("ro0", &ro0)) ro0 = -1;
  }
  if(!getparfloat("sizex", &sizex)) verr("x-model size not specified.");
  if(!getparfloat("sizey", &sizey)) verr("y-model size not specified.");
  if(!getparfloat("sizez", &sizez)) verr("z-model size not specified.");
  if(!getparfloat("dx", &dx)) verr("grid distance dx not specified.");
  if(!getparfloat("dy", &dy)) verr("grid distance dy not specified.");
  if(!getparfloat("dz", &dz)) verr("grid distance dz not specified.");
  if(!getparfloat("orig", orig)) orig[0] = orig[1] = orig[2] = 0.0;
  if(!getparlong("gradt", &gradt)) gradt = 1;
  if(!getparlong("gradunit", &gradunit)) gradunit = 0;
  if(!getparlong("writeint", &writeint)) writeint = 0;
  if(!getparlong("reflectivity", &reflectivity)) reflectivity = 0;
  if(!getparlong("supersmooth", &supersmooth)) supersmooth = 0;
  if(!getparfloat("sigma", &sigma)) sigma = 0.8;
  if(!getparlong("rayfile", &rayfile)) rayfile = 0;
  if(!getparlong("skip", &skip)) skip = 5;
  if(!getparlong("above", &above)) above=0;
  if(!getparlong("verbose", &verbose)) verbose=0;
  if(!getparlong("dtype", &dtype)) dtype = 0;
  
  if (supersmooth==1) {
    dx = dx/4;
    dy = dy/4;
    dz = dz/4;
  }
  if ((writeint == 1) || (rayfile == 1)) store_int = 1;
  else store_int = 0;
  
  /*=================== check parameters =================*/
  
  Np = countparname("cp");
  Ns = countparname("cs");
  Nr = countparname("ro");
  Ng = countparname("grad");
  No = countparname("poly");
  Ni = countparname("intt");
  Nv = countparname("var");
  Ngp = countparname("gradcp");
  Ngs = countparname("gradcs");
  Ngr = countparname("gradro");
  
  Nvi = 0;
  for (jint = 1; jint <= Ni; jint++) {
    getnparstring(jint,"intt", &file);
    strcpy(intt, file);
    if (strstr(intt,"sin") != NULL) Nvi++;
    if (strstr(intt,"rough") != NULL) Nvi++;
    if (strstr(intt,"fract") != NULL) Nvi++;
    if (strstr(intt,"elipse") != NULL) Nvi++;
	if (strstr(intt,"random") != NULL) Nvi++;
//	if (strstr(intt,"randdf") != NULL) Nvi++;
    if (strstr(intt,"diffr") != NULL || strstr(intt,"randdf") != NULL) {
	  Nvi++;
//      if (Ng != 0) Ng++; 
//      if (No != 0) No++;
    }
  }
//  fprintf(stderr,"Nvi=%li ng=%li No=%li np=%li,", Nvi,Ng,No,Np);
  
  if (Np != Nr && ro0 > 0) verr("number of cp and ro not equal.");
  if (Np != Ni) verr("number of cp and interfaces not equal.");
  if (Nvi != Nv) verr("number of interface variables(var) not correct.");
  if (Ns == 0 && Nr == 0) if (verbose>=2) vmess("Velocity model.");
  if (Ns == 0) { if (verbose>=2) vmess("Acoustic model."); }
  else {
    if (verbose>=2) vmess("Elastic model.");
    if (Np != Ns) verr("number of cp and cs not equal");
  }
  
  if (Ng == 0) {
    if (verbose>=2) vmess("No boundary gradients are defined.");
  }
  else if (Ng != Np) {
    verr("number of boundary gradients and interfaces are not equal.");
  }
  if (Ngp == 0) {
    if (verbose>=2) vmess("No interface gradients for cp defined.");
  }
  else if (Ngp != Np) {
    verr("gradcp gradients and interfaces are not equal.");
  }
  if (Ngs == 0) {
    if (verbose>=2) vmess("No interface gradients for cs defined.");
  }
  else if (Ngs != Np) {
    verr("gradcs gradients and interfaces are not equal.");
  }
  if (Ngr == 0) {
    if (verbose>=2) vmess("No interface gradients for rho defined.");
  }
  else if (Ngr != Np) {
    verr("gradro gradients and interfaces are not equal.");
  }
  if (No == 0) {
    if (verbose>=2) vmess("All interfaces are linear.");
  }
//  else if (No != Np) {
//    verr("number of poly variables and interfaces are not equal.");
//  }
  
  if (Np > 0) {
    if (countparname("x") != Np)
      verr("a x array must be specified for each interface.");
    if (countparname("y") != Np)
      verr("a y array must be specified for each interface.");
    if (countparname("z") != Np)
      verr("a z array must be specified for each interface.");
  } 
  else Np = 1;

  if (Nzr != Np && Nzr !=0) {
    verr("number of zref gradients not equal to number of interfaces");
  }
  
  /*=================== initialization of arrays =================*/
  
  nz = NINT(sizez/dz)+1;
  ny = NINT(sizey/dy)+1;
  nx = NINT(sizex/dx)+1;
  
  zp = (long *)malloc(nx*ny*sizeof(long));
  interface = (float *)malloc(nx*ny*sizeof(float));
  var = (float *)malloc(8*sizeof(float));
  gridcp = alloc2float(nz, nx*ny);
  if(gridcp == NULL) verr("memory allocation error gridcp");
  if (Ns || (NINT(cs0*1e3) >= 0)) {
    gridcs = alloc2float(nz, nx*ny);
    if(gridcs == NULL) verr("memory allocation error gridcs");
  }
  else gridcs = NULL;
  if (Nr || (NINT(ro0*1e3) >= 0)) {
    gridro = alloc2float(nz, nx*ny);
    if(gridro == NULL) verr("memory allocation error gridro");
  }
  else gridro = NULL;
  
  cp = alloc2float(nx*ny,3);
  cs = alloc2float(nx*ny,3);
  ro = alloc2float(nx*ny,3);
  if (store_int == 1) inter = alloc2float(nx*ny, 2*Np);
  
  if (verbose) {
    vmess("Origin top left (x,y,z) . = %.1f, %.1f, %.1f", orig[0], orig[1], orig[2]);
    vmess("Base name ............. = %s", file_base);
    vmess("Number of interfaces .. = %li", Np);
    vmess("length in x ........... = %f (=%li)", sizex, nx);
    vmess("length in y ........... = %f (=%li)", sizey, ny);
    vmess("length in z ........... = %f (=%li)", sizez, nz);
    vmess("delta x ............... = %f", dx);
    vmess("delta y ............... = %f", dy);
    vmess("delta z ............... = %f", dz);
    vmess("cp0 ................... = %f", cp0);
    if (Ns) vmess("cs0 ................... = %f", cs0);
    if (Nr) vmess("ro0 ................... = %f", ro0);
    vmess("write interfaces ...... = %li", writeint);
    vmess("store interfaces ...... = %li", store_int);
    if (above) vmess("define model above interface");
    else vmess("define model below interface");
  }

  /*========== initializing for homogeneous background =============*/
  
  nminx = 0;
  nmaxx = nx;
  nminy = 0;
  nmaxy = ny;
  for (l = nminy; l < nmaxy; l++) {
    for (j = nminx; j < nmaxx; j++) {
        cp[0][l*nmaxx+j] = cp0;
        cs[0][l*nmaxx+j] = cs0;
        ro[0][l*nmaxx+j] = ro0;
        zp[l*nmaxx+j] = 0;
        cp[1][l*nmaxx+j] = cp0;
        cs[1][l*nmaxx+j] = cs0;
        ro[1][l*nmaxx+j] = ro0;
    }
  }
  
  gradlen = 0.0;
  gradcp = gradcs = gradro = 0.0;
  optgrad = 3;
  if (above == 0) {
    Nvi = 1; Noi = 1;
  } else {
    Nvi = Ngp; Noi = Ngp;
  }
  grid3D(gridcp, gridcs, gridro, zp, cp, cs, ro, nminx, nmaxx, 
       nminy, nmaxy, optgrad, gradlen, gradcp, gradcs, gradro, dx, dy, dz, nz);
  
  nxp = nyp = nzp = 2;
  x = (float *)malloc(nxp*sizeof(float));
  y = (float *)malloc(nyp*sizeof(float));
  z = (float *)malloc(nzp*sizeof(float));
  
  if (Ni == 0) {
    if (verbose) vmess("No interfaces are defined, Homogeneous model.");
    Np = 0;
  }
  
  /*========== filling gridded model with interfaces =============*/
  
  for (jcount = 1; jcount <= Np; jcount++) {
    
    /* for above interface definition reverse	*/
    /* order of interfaces to scan 			*/
    if (above == 0) jint=jcount;
    else jint=Np+1-jcount;
    
    if (verbose) vmess("***** Interface number %li *****", jint);
    
    getnparstring(jint,"intt", &file);
    strcpy(intt, file);
    
    nxp = countnparval(jint,"x");
    nyp = countnparval(jint,"y");
    nzp = countnparval(jint,"z");
    if ((nxp != nzp) || (nxp != nyp)) {
      vmess("nxp = %li nyp = %li nzp =%li for interface %li",nxp, nyp, nzp, jint);
      verr("Number of x, y and z values not equal for interface %li",jint);
    }
    ncp = countnparval(jint,"cp");
    nro = countnparval(jint,"ro");
    ncs = countnparval(jint,"cs");
    
    if (ncp == 1) {
      if (verbose>=2) vmess("No lateral gradient in P velocity");
    }
    else if (ncp == 2) {
      if (verbose) vmess("lateral P-gradient from begin to end");
    }
    else if (ncp != nxp) {
      vmess("ncp = %li nxp =%li for interface %li",ncp, nxp, jint);
      verr("Number of cp's and x not equal for interface %li",jint);
    }
    else if (ncp != nyp) {
      vmess("ncp = %li nyp =%li for interface %li",ncp, nyp, jint);
      verr("Number of cp's and y not equal for interface %li",jint);
    }
    if (nro <= 1) {
      if (verbose>=2) vmess("No lateral gradient in density");
    }
    else if (nro == 2) {
      if (verbose) vmess("lateral Rho-gradient from begin to end");
    }
    else if (nro != nxp) {
      vmess("nro = %li nxp =%li for interface %li",nro, nxp, jint);
      verr("Number of ro's and x not equal for interface %li",jint);
    }
    else if (nro != nyp) {
      vmess("nro = %li nyp =%li for interface %li",nro, nyp, jint);
      verr("Number of ro's and y not equal for interface %li",jint);
    }
    if (ncs <= 1) {
      if (verbose>=2) vmess("No lateral gradient in S velocity");
    }
    else if (ncs == 2) {
      if (verbose) vmess("lateral S-gradient from begin to end");
    }
    else if (ncs != nxp) {
      vmess("ncs = %li nxp =%li for interface %li",ncs, nxp, jint);
      verr("Number of cs's and x not equal for interface %li",jint);
    }
    else if (ncs != nyp) {
      vmess("ncs = %li nyp =%li for interface %li",ncs, nyp, jint);
      verr("Number of cs's and y not equal for interface %li",jint);
    }
    
    nvel = MAX(ncp, MAX(nro, ncs));
    
    free(x);
    free(y);
    free(z);
    x = (float *)malloc(nxp*sizeof(float));
    y = (float *)malloc(nyp*sizeof(float));
    z = (float *)malloc(nzp*sizeof(float));
	memset(interface, 0, nx*ny*sizeof(float));
    
    getnparfloat(jint,"x",x);
    getnparfloat(jint,"y",y);
    getnparfloat(jint,"z",z);
    getnparfloat(jint,"cp",cp[2]);
    if (Nr == 0) ro[2][0] = 0.0;
    else getnparfloat(jint,"ro", ro[2]);
    if (Ns == 0) cs[2][0] = 0.0;
    else getnparfloat(jint,"cs", cs[2]);
    if (Ng == 0) gradlen = 0.0;
    else getnparfloat(Noi,"grad", &gradlen);
    if (No == 0) poly = 0;
    else getnparint(Noi,"poly", &poly);
    if (Ngp == 0) gradcp = 0.0;
    else getnparfloat(Noi,"gradcp", &gradcp);
    if (Ngs == 0) gradcs = 0.0;
    else getnparfloat(Noi,"gradcs", &gradcs);
    if (Ngr == 0) gradro = 0.0;
    else getnparfloat(Noi,"gradro", &gradro);
    /* if gradunit is in meters, recalculate gradcp,gradcs and gradro */
    if (gradunit > 0) {
      gradcs = gradcs * dz;
      gradcp = gradcp * dz;
      gradro = gradro * dz;
    }
    
    if (nvel != 1) {
      if (ncp == 1) {
	for (j = 1; j < nvel; j++) cp[2][j] = cp[2][0];
      }
      if (ncs == 1) {
	for (j = 1; j < nvel; j++) cs[2][j] = cs[2][0];
      }
      if (nro == 1) {
	for (j = 1; j < nvel; j++) ro[2][j] = ro[2][0];
      }
    }
    
    if (verbose) {
      vmess("Interface type .......... = %s", intt);
      vmess("Boundary gradient ....... = %f", gradlen);
      vmess("Interface gradient cp ... = %f", gradcp);
      if (Ns) vmess("Interface gradient cs ... = %f", gradcs);
      if (Nr) vmess("Interface gradient ro ... = %f", gradro);
      if (verbose>=2) {
	vmess("Polynomal ............... = %li", poly);
	vmess("Number of (x,z) points... = %li", nxp);
	vmess("P-wave velocities ....... = %li", ncp);
	if (Ns) vmess("S-wave velocities ....... = %li", ncs);
	if (Nr) vmess("Densities ............... = %li", nro);
      }
      for (j = 0; j < nxp; j++) {
	vmess("x = %6.1f \t z = %6.1f", x[j], z[j]);
	if (nvel != 1) {
	  vmess("    cp = %f", cp[2][j]);
	  if (Ns) vmess("    cs = %f", cs[2][j]);
	  if (Nr) vmess("   rho = %f", ro[2][j]);
	}
      }
      if (nvel == 1) {
	vmess("    cp = %f", cp[2][0]);
	if (Ns) vmess("    cs = %f", cs[2][0]);
	if (Nr) vmess("    rho = %f", ro[2][0]);
      }
    }
    
    for (j = 0; j < nxp; j++) {
      x[j] -= orig[0];
      y[j] -= orig[1];
      z[j] -= orig[2];
    }
    for (j = 0; j < nxp; j++) {
      if(x[j] > sizex) verr("x coordinate bigger than model");
      if(y[j] > sizey) verr("y coordinate bigger than model");
      if(z[j] > sizez) verr("z coordinate bigger than model");
    }
    if (gradlen > 0) optgrad = gradt;
    else optgrad = 3;
    
	if (strstr(intt,"random") != NULL) {
		Nv = countnparval(Nvi,"var");
		if (Nv != 1) verr("Random interface must have 1 variables.");
		getnparfloat(Nvi,"var", var);
		lseed = (long)var[0];
		srand48(lseed);
		gradcp=gradcs=gradro=var[0];
		optgrad = 4;
        if (above == 0) Noi++; else Noi--;
        if (above == 0) Nvi++; else Nvi--;
	}
	  
    if ((strstr(intt,"diffr") == NULL) && (strstr(intt,"randdf") == NULL)) {
      interpolation(x, z, nxp, nx, poly, &nminx, &nmaxx, dx, 
		    cp, cs, ro, nvel, interface);
    }

    if ( (strstr(intt,"def") != NULL) || (strstr(intt,"random") != NULL) ) {
      linearint(zp, nminx, nmaxx, dz, interface);
      if (above == 0) Noi++; else Noi--;
    }

    if (strstr(intt,"sin") != NULL) {
      Nv = countnparval(Nvi,"var");
      if (Nv != 2) verr("Sinus interface must have 2 variables.");
      getnparfloat(Nvi,"var", var);
      sinusint(zp, nminx, nmaxx, dz, interface, dx, var[0], var[1]);
      if (above == 0) Noi++; else Noi--;
      if (above == 0) Nvi++; else Nvi--;
    }
    else if (strstr(intt,"rough") != NULL) {
      Nv = countnparval(Nvi,"var");
      if (Nv != 3) verr("Rough interface must have 3 variables.");
      getnparfloat(Nvi,"var", var);
      roughint(zp, nminx, nmaxx, dz, interface, var[0],var[1],var[2]);
      if (above == 0) Noi++; else Noi--;
      if (above == 0) Nvi++; else Nvi--;
    }
    else if (strstr(intt,"fract") != NULL) {
      Nv = countnparval(Nvi, "var");
      if (Nv != 6) verr("Fractal interface must have 6 variables.");
      getnparfloat(Nvi,"var", var);
      fractint(zp, nminx, nmaxx, dx, dz, interface, var[0], var[1], 
	       var[2], var[3], var[4], var[5]);
      if (above == 0) Noi++; else Noi--;
      if (above == 0) Nvi++; else Nvi--;
    }
    if ((strstr(intt,"elipse") != NULL) || (strstr(intt,"diffr") != NULL) || (strstr(intt,"randdf") != NULL)) {
   	    if (strstr(intt,"randdf") != NULL) {
        	Nv = countnparval(Nvi, "var");
			if (Nv != 2) verr("randdf interface must have 2 variables: number of points, width.");
			getnparfloat(Nvi,"var", var);
        	if(!getparlong("dtype", &dtype)) dtype = -1;
        
        	randdf(x, z, nxp, dx, dz, gridcp, gridcs, gridro, cp, cs, ro,
           	   interface, zp, nx, sizex, sizez, var[0], (long)var[1], dtype);

			if (above == 0) Noi++; else Noi--;
			if (above == 0) Nvi++; else Nvi--;
		}
		if (strstr(intt,"elipse") != NULL) {
			Nv = countnparval(Nvi, "var");
			if (Nv != 2) verr("Elipse interface must have 2 variables.");
			getnparfloat(Nvi,"var", var);
			elipse(x, z, nxp, dx, dz, gridcp, gridcs, gridro, 
				cp, cs, ro, interface, zp, nz, nx, var[0], var[1], gradcp, gradcs, gradro);
			if (above == 0) Noi++; else Noi--;
			if (above == 0) Nvi++; else Nvi--;
		}
		else if ((strstr(intt,"diffr") != NULL)) {
			Nv = countnparval(Nvi, "var");
        	if (Nv == 2 || Nv == 1) {
           		getnparfloat(Nvi,"var", var);
           		diffrwidth=(long)var[0];
           		if (Nv==1) {
           			if(!getparint("dtype", &dtype)) dtype = 0;
			}
           		else dtype=(long)var[1];
        	}
        	else {
           		verr("diffr interface must have 1 or 2 variables: width,type.");
        	}
        	
			diffraction(x, z, nxp, dx, dz, gridcp, gridcs, gridro,
				cp, cs, ro, interface, zp, nx, diffrwidth, dtype);
			if (above == 0) Noi++; else Noi--;
			if (above == 0) Nvi++; else Nvi--;
		}
	}
    else {
		if (above == 0) {
			grid(gridcp, gridcs, gridro, zp, cp, cs, ro, nminx, nmaxx, 
				optgrad, gradlen, gradcp, gradcs, gradro, dx, dz, nz);
		} 
		else {
			gridabove(gridcp, gridcs, gridro, zp, cp, cs, ro, nminx, nmaxx, 
				optgrad, gradlen, gradcp, gradcs, gradro, dx, dz, nz);
		}
	}
    
    if (store_int == 1) {
      for(j = 0; j < nminx; j++) inter[jint-1][j] = 0.0;
      for(j = nminx; j < nmaxx; j++) inter[jint-1][j] = interface[j];
      for(j = nmaxx; j < nx; j++) inter[jint-1][j] = 0.0;
      
      for(j = 0; j < nminx; j++) inter[jint-1+Np][j] = 0.0;
      for(j = nminx; j < nmaxx; j++) inter[jint-1+Np][j] = zp[j]*dz;
      for(j = nmaxx; j < nx; j++) inter[jint-1+Np][j] = 0.0;
    }
  } /* end of loop over interfaces */
  
/* apply supersampled data smoothing operator according to Zeng and West (1996) */
  
  if (supersmooth==1) {
	float Wsmooth[5][5], C, iC, xx, xz, *dataS, smooth;
	long ixi, izi, nxout, nzout;

	C=0;
	sigma = -1.0*log(sigma)/(dx*(powf(0.25*2.0,2.0)));
    for(ixi = -2; ixi < 3; ixi++) {
      for(izi = -2; izi < 3; izi++) {
        xx = 0.25*ixi;
        xz = 0.25*izi;
        Wsmooth[ixi+2][izi+2] = exp(-sigma*dx*(xx*xx+xz*xz) );
		//fprintf(stderr,"Wsmooth[%li][%li] = %f\n", ixi, izi, Wsmooth[ixi+2][izi+2]);
		C += Wsmooth[ixi+2][izi+2];
	  }
    }
    iC = 1.0/C;
	//fprintf(stderr,"sigma=%f %f C=%f\n",sigma, exp(-sigma*dx*powf(0.25*2.0,2.0)),C);

	nxout = (nx-1)/4+1;
	nzout = (nz-1)/4+1;
	//fprintf(stderr,"nxout=%li nzout=%li nx=%li nz=%li \n",nxout,nzout, nx, nz);
    dataS = malloc(nxout*nzout*sizeof(float));
    for(ix = 0; ix < nxout-1; ix++) {
      for(iz = 0; iz < nzout-1; iz++) {
        smooth = 0.0;
        for(ixi = -2; ixi < 3; ixi++) {
          for(izi = -2; izi < 3; izi++) {
            smooth += Wsmooth[ixi+2][izi+2]/gridcp[4*ix+2+ixi][4*iz+2+izi];
	      }
        }
        dataS[ix*nzout+iz] = 1.0/(iC*smooth);
	  }
      iz = nzout-1;
      dataS[ix*nzout+iz] = dataS[ix*nzout+iz-1];
    }
    ix = nxout-1;
    for(iz = 0; iz < nzout-1; iz++) {
      dataS[ix*nzout+iz] = dataS[(ix-1)*nzout+iz];
	}
    ix = nxout-1; iz = nzout-1;
    dataS[ix*nzout+iz] = dataS[(ix-1)*nzout+iz-1];

    for(ix = 0; ix < nxout; ix++) {
      for(iz = 0; iz < nzout; iz++) {
        gridcp[ix][iz] = dataS[ix*nzout+iz];
      }
    }

    /* smoothin densities */
    if (Nr > 0 || getparfloat("ro0", &ro0)) {
      for(ix = 0; ix < nxout-1; ix++) {
        for(iz = 0; iz < nzout-1; iz++) {
          smooth = 0.0;
          for(ixi = -2; ixi < 3; ixi++) {
            for(izi = -2; izi < 3; izi++) {
              smooth += Wsmooth[ixi+2][izi+2]*gridro[4*ix+2+ixi][4*iz+2+izi];
	        }
          }
          dataS[ix*nzout+iz] = iC*smooth;
	    }
        iz = nzout-1;
        dataS[ix*nzout+iz] = dataS[ix*nzout+iz-1];
      }
      ix = nxout-1;
      for(iz = 0; iz < nzout-1; iz++) {
        dataS[ix*nzout+iz] = dataS[(ix-1)*nzout+iz];
	  }
      ix = nxout-1; iz = nzout-1;
      dataS[ix*nzout+iz] = dataS[(ix-1)*nzout+iz-1];

      for(ix = 0; ix < nxout; ix++) {
        for(iz = 0; iz < nzout; iz++) {
          gridro[ix][iz] = dataS[ix*nzout+iz];
        }
      }
    }
    nx = nxout;
	nz = nzout;
    dx = dx*4;
    dz = dz*4;
	free(dataS);
  }

  if (verbose) vmess("Writing data to disk.");
  
  hdrs = (segy *) calloc(nx,sizeof(segy));
  for(j = 0; j < nx; j++) {
	hdrs[j].f1= orig[1];
	hdrs[j].f2= orig[0];
	hdrs[j].d1= dz;
	hdrs[j].d2= dx;
    hdrs[j].ns= nz;
	hdrs[j].trwf= nx;
	hdrs[j].tracl= j;
	hdrs[j].tracf= j;
    hdrs[j].gx = (orig[0] + j*dx)*1000;
    hdrs[j].scalco = -1000;
    hdrs[j].timbas = 25;
    hdrs[j].trid = TRID_DEPTH;
  }

  /* due to bug in SU, int-file has to be opened before other files are closed */
	if (writeint == 1) {
		strcpy(filename, file_base);
		name_ext(filename, "_int");
		fpint = fopen(filename,"w");
		assert(fpint != NULL);
	}

 
  /* write P-velocities in file */
  strcpy(filename, file_base);
  name_ext(filename, "_cp");
  fpcp = fopen(filename,"w");
  assert(fpcp != NULL);
  
  data_out = (float *)malloc(nx*nz*sizeof(float));
  for(ix = 0; ix < nx; ix++) {
    for(iz = 0; iz < nz; iz++) {
      data_out[ix*nz+iz] = gridcp[ix][iz];
	}
    nwrite = fwrite( &hdrs[ix], 1, TRCBYTES, fpcp);
    assert(nwrite == TRCBYTES);
    nwrite = fwrite( &data_out[ix*nz], sizeof(float), nz, fpcp);
    assert(nwrite == nz);
  }
  fclose(fpcp);

  /* write reflectivity of the interfaces */
  if (reflectivity == 1) {
    strcpy(filename, file_base);
    name_ext(filename, "_rfl");
    fpro = fopen(filename,"w");
    assert(fpro != NULL);
    
    for(ix = 0; ix < nx-1; ix++) {
      for(iz = 0; iz < nz-1; iz++) {
         Z1=gridro[ix][iz]*gridcp[ix][iz];
         Z2=gridro[ix][iz+1]*gridcp[ix][iz+1];
	     data_out[ix*nz+iz] = (Z2-Z1)/(Z2+Z1);
      }
	  data_out[ix*nz+nz-1] = 0.0;
      nwrite = fwrite( &hdrs[ix], 1, TRCBYTES, fpro);
      assert(nwrite == TRCBYTES);
      nwrite = fwrite( &data_out[ix*nz], sizeof(float), nz, fpro);
      assert(nwrite == nz);
    }
    ix = nx-1;
    for(iz = 0; iz < nz; iz++) {
	   data_out[ix*nz+iz] = 0.0;
    }
    nwrite = fwrite( &hdrs[ix], 1, TRCBYTES, fpro);
    assert(nwrite == TRCBYTES);
    nwrite = fwrite( &data_out[ix*nz], sizeof(float), nz, fpro);
    assert(nwrite == nz);
    fclose(fpro);
  } /* end of writing reflectivity file */

  free2float(gridcp);
  
  /* write S-velocities in file */
  if (Ns > 0 || getparfloat("cs0", &cs0)) {
    strcpy(filename, file_base);
    name_ext(filename, "_cs");
    fpcs = fopen(filename,"w");
    assert(fpcs != NULL);
    
    for(ix = 0; ix < nx; ix++) {
      for(iz = 0; iz < nz; iz++) {
	     data_out[ix*nz+iz] = gridcs[ix][iz];
      }
      nwrite = fwrite( &hdrs[ix], 1, TRCBYTES, fpcs);
      assert(nwrite == TRCBYTES);
      nwrite = fwrite( &data_out[ix*nz], sizeof(float), nz, fpcs);
      assert(nwrite == nz);
    }
    fclose(fpcs);
    free2float(gridcs);
    
  } /* end of writing S-velocity file */
  
  /* write densities in file */
  if (Nr > 0 || getparfloat("ro0", &ro0)) {
    strcpy(filename, file_base);
    name_ext(filename, "_ro");
    fpro = fopen(filename,"w");
    assert(fpro != NULL);
    
    for(ix = 0; ix < nx; ix++) {
      for(iz = 0; iz < nz; iz++) {
	     data_out[ix*nz+iz] = gridro[ix][iz];
      }
      nwrite = fwrite( &hdrs[ix], 1, TRCBYTES, fpro);
      assert(nwrite == TRCBYTES);
      nwrite = fwrite( &data_out[ix*nz], sizeof(float), nz, fpro);
      assert(nwrite == nz);
    }
    fclose(fpro);
    free2float(gridro);
  } /* end of writing density file */
  
  
  /* write depths of the interfaces */
  if (writeint == 1) {
    free(hdrs);
    hdrs = (segy *) calloc(Np,sizeof(segy));
    for(j = 0; j < Np; j++) {
      hdrs[j].fldr = 1;
      hdrs[j].timbas = 25;
	  hdrs[j].f1= orig[0];
	  hdrs[j].f2= 0.0;
	  hdrs[j].d1= dx;
	  hdrs[j].d2= dz;
	  hdrs[j].ns= nx;
	  hdrs[j].trwf= Np;
	  hdrs[j].tracl= j;
	  hdrs[j].tracf= j;
	  hdrs[j].trid= TRID_DEPTH;
    }
    
    /* note that due to bug in SU, interface file has already been opened */
    strcpy(filename, file_base);
    name_ext(filename, "_int");

    free(data_out);
    data_out = (float *)malloc(nx*Np*sizeof(float));
    for(jint = 0; jint < Np; jint++) {
      for(j = 0; j < nx; j++) {
	    data_out[jint*nx+j] = inter[jint][j]+orig[1];
      }
      nwrite = fwrite( &hdrs[jint], 1, TRCBYTES, fpint);
      assert(nwrite == TRCBYTES);
      nwrite = fwrite( &data_out[jint*nx], sizeof(float), nx, fpint);
      assert(nwrite == nx);
    }
    
	for(j = 0; j < Np; j++) hdrs[j].fldr = 2;
    for(jint = 0; jint < Np; jint++) {
      for(j = 0; j < nx; j++) {
	     data_out[jint*nx+j] = inter[jint+Np][j]+orig[1];
      }
      nwrite = fwrite( &hdrs[jint], 1, TRCBYTES, fpint);
      assert(nwrite == TRCBYTES);
      nwrite = fwrite( &data_out[jint*nx], sizeof(float), nx, fpint);
      assert(nwrite == nx);
    }
    fclose(fpint);
    
  } /* end of writing interface file */
  
  if (rayfile == 1) {
    FILE    *fpout, *fpcurve;
    strcpy(filename, file_base);
    strcpy(strrchr(filename, '.'), ".mod");
    
    fpout = fopen(filename, "w+");
    fprintf(fpout,"RAYTRACE MODEL FILE\n");
    fprintf(fpout,"# ASCII file for ray-tracer\n\n");
    fprintf(fpout,"# Top interface\n\n");
    fprintf(fpout,"x=0,%.1f\n", sizex);
    fprintf(fpout,"z=0.,0.\n");
    
    /*		for(i = 1; i <= Np; i++) {
		fprintf(fpout,"\n# %li th interface\n\nx=",i);
		nxp = countnparval(i,"x");
		nzp = countnparval(i,"z");
		free(x);
		free(z);
		x = (float *)malloc(nxp*sizeof(float));
		z = (float *)malloc(nzp*sizeof(float));
		getnparfloat(i,"x",x);
		getnparfloat(i,"z",z);
		for(j = 0; j < (nxp-1); j ++) fprintf(fpout,"%.1f,", x[j]);
		fprintf(fpout,"%.1f\nz=", x[nxp-1]);
		for(j = 0; j < (nxp-1); j ++) fprintf(fpout,"%.1f,", z[j]);
		fprintf(fpout,"%.1f\n", z[nxp-1]);
		}
    */
	for(jint = 0; jint < Np; jint++) {
		fprintf(fpout,"\n# %li th interface\n\nx=0",jint+1);
		sprintf(filename,"curve%li",jint+1);
		fpcurve = fopen(filename, "w+");
		if (verbose) vmess("writing %li th interface to ascii file %s for usage in su(ps/x)image ",jint+1, filename);
		for(j = 0; j < nx; j += skip)  {
			fprintf(fpout,",%.1f", (float)j*dx);
			fprintf(fpcurve,"%.1f %.1f\n", inter[jint][j], orig[0]+(float)j*dx);
		}
		fclose(fpcurve);
      
		fprintf(fpout,"\nz=%.1f", inter[jint][0]);
		for(j = skip; j < nx; j += skip) 
			fprintf(fpout,",%.1f", inter[jint][j]);
		fprintf(fpout,"\n");
	}
	fprintf(fpout,"\n# %li th interface\n\nx=0",jint+1);
    for(j = skip; j < nx; j += skip) 
      fprintf(fpout,",%.1f", (float)j*dx);
    
    fprintf(fpout,"\nz=%.1f", sizez);
    for(j = skip; j < nx; j += skip) 
      fprintf(fpout,",%.1f", sizez);
    fprintf(fpout,"\n");
    
    /**/
    fprintf(fpout,"\n\n");
    fprintf(fpout,"cp=%.1f,", back[0]);
    for(jint = 1; jint <= Np; jint++) {
      aver = 0.0;
      ncp = countnparval(jint,"cp");
      getnparfloat(jint,"cp",cp[2]);
      for(j = 0; j < ncp; j++) 
	aver += cp[2][j];
      aver = aver/((float)ncp);
      if (jint == Np ) fprintf(fpout,"%.1f", aver);
      else fprintf(fpout,"%.1f,", aver);
    }
    
    fclose(fpout);
    free2float(inter);
  }
  
  free(hdrs);
  
  return 0;
}
 
