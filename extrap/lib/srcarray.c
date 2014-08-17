#include "optim.h"
#include "par.h"

float setzsrc(int nb, int *boundary, float *inter, int nx, int ni, float zsrc1, float dzsrc, float h, float oz, int nz, float xsrc, float ox, int id, int verbose);

void srcarray(complex *source, int nx, int nz, int nb, int *Ns, int ik, 
	int *boundary, float *inter, int ni, float zsrc1, float dzsrc, 
	float dz, float oz, float xsrc1, float dxsrc, float dx, float ox, 
	int *izmax, int *izmin, int add, int wnx, float *latwav, int *izsrc, int *ixsrc,
	int verbose)
{
	int 	izs, ixs, ix, iz, is, ib, i, hnx, io, zero=0;
	float	xsrc, zsrc;
	static int todo=0;

	if (ISODD(wnx)) hnx = (wnx+1)/2;
	else hnx = wnx/2;

	*izmax = 0;
	*izmin = nz;
	for (iz = 0; iz < nz; iz++) {
		for (ix = 0; ix < nx; ix++) {
			source[iz*nx+ix].r = 0.0;
			source[iz*nx+ix].i = 0.0;
		}
	}

	if (add) {
		for (is = 0; is < *Ns; is++) {
			xsrc = xsrc1 + is*dxsrc - ox;
			ixs = NINT(xsrc/dx);
			if (nb) {
				for (ib = 0; ib < nb; ib++) {
					io = ixs-hnx+1;
					if (io < 0) {
						if (verbose) vwarn("lateral wavelet outside model");
						ix = 0;
					}
					else ix = io;
					if (io+wnx > nx && verbose) 
						vwarn("lateral wavelet outside model");
					for (i = ix; i < MIN(io+wnx, nx); i++) {
						xsrc += (i-ixs)*dx;
						zsrc = setzsrc(nb,boundary,inter,nx,ni,zsrc1,dzsrc,dx,
							oz,nz,xsrc,ox,ib,verbose);
						izs = NINT(zsrc/dz);
						*izmax = MAX(izs, *izmax);
						*izmin = MIN(izs, *izmin);
						source[izs*nx+i].r += latwav[i-ix];
					}
					if (verbose>=2) 
						vmess("xsrc = %f, zsrc = %f", xsrc+ox, zsrc+oz);
				}
			}
			else {
				zsrc = zsrc1 + is*dzsrc - oz;
				izs = NINT(zsrc/dz);
				*izmax = MAX(izs, *izmax);
				*izmin = MIN(izs, *izmin);
				io = ixs-hnx+1;
				if (io < 0) {
					if (verbose) vwarn("lateral wavelet outside model");
					ix = 0;
				}
				else ix = io;
				if (io+wnx > nx && verbose) 
					vwarn("lateral wavelet outside model");
				for (i = ix; i < MIN(io+wnx, nx); i++) {
					source[izs*nx+i].r += latwav[i-ix];
				}
				if (verbose) 
					vmess("Adding xsrc = %f(%d), zsrc = %f(%d)", xsrc+ox ,ixs, zsrc+oz, izs);
			}
		}
		*Ns = 1;
	}
	else {
		if (ik < *Ns) {
			xsrc = xsrc1 + ik*dxsrc - ox;
			ixs  = NINT(xsrc/dx);
			if (nb) {
				if (todo < nb) {
					zsrc = setzsrc(nb,boundary,inter,nx,ni,zsrc1,dzsrc,dx,
						oz,nz,xsrc,ox,todo,zero);
					izs = NINT(zsrc/dz);
					todo++;
				}
				else {
					todo=0;
					zsrc = setzsrc(nb,boundary,inter,nx,ni,zsrc1,dzsrc,dx,
						oz,nz,xsrc,ox,todo,zero);
					izs = NINT(zsrc/dz);
					todo++;
				}
			}
			else {
				zsrc = zsrc1 + ik*dzsrc - oz;
				izs = NINT(zsrc/dz);
				*izmax = izs;
				*izmin = izs;
			}
			if (verbose) vmess("xsrc = %f (%d), zsrc = %f (%d)", 
				xsrc+ox, ixs, zsrc+oz, izs);
			ixs = NINT(xsrc/dx);

			if ( (NINT((ixs*dx)*1000)-NINT(xsrc*1000)) || 
				 (NINT((izs*dz)*1000)-NINT(zsrc*1000)) ) {
				if (verbose) vmess("source position not on a gridpoint");

/* To Do implement a source distribution to simulate the correct position */

				/*izs = floor(zsrc/dz);*/
				io = ixs-hnx+1;
				if (io < 0) ix = 0;
				else ix = io;
				for (i = ix; i < MIN(io+wnx, nx); i++) {
					if (nb) {
						xsrc += (i-ixs)*dx;
						zsrc = setzsrc(nb,boundary,inter,nx,ni,zsrc1,dzsrc,dx,
							oz,nz,xsrc,ox,(todo-1),verbose);
						izs = NINT(zsrc/dz);
						*izmax = MAX(izs, *izmax);
						*izmin = MIN(izs, *izmin);
					}
					source[izs*nx+i].r = latwav[i-ix];
				}

			}
			else {
				io = ixs-hnx+1;
				if (io < 0) ix = 0;
				else ix = io;
				for (i = ix; i < MIN(io+wnx, nx); i++) {
					if (nb) {
						xsrc += (i-ixs)*dx;
						zsrc = setzsrc(nb,boundary,inter,nx,ni,zsrc1,dzsrc,dx,
							oz,nz,xsrc,ox,(todo-1),verbose);
						izs = NINT(zsrc/dz);
						*izmax = MAX(izs, *izmax);
						*izmin = MIN(izs, *izmin);
					}
					source[izs*nx+i].r = latwav[i-ix];
				}
			}

		}	
	}
	*ixsrc = ixs;
	*izsrc = izs;

	return;
}


