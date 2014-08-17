#include <math.h>
#include <stdlib.h>
#include "Area.h"

void spline3(float x1, float x2, float z1, float z2, float dzdx1, float dzdx2, float *a, float *b, float *c, float *d);

void Woper3d(float k, float dx, float dy, float dz, int nkx, int nky, complex *oper);

void Woper3d_ph(float k, float dx, float dy, float dz, int nkx2, int nky2, float alpha, complex *oper);

void W3d(float kp, float dx, float dy, float dz, int nkx, int nky, float alpha, complex *oper, float *wfacto)
{


	if (kp <= 0.85*M_PI/dx) 
		Woper3d(kp, dx, dy, dz, nkx, nky, oper);
	else 
		Woper3d_ph(kp, dx, dy, dz, nkx, nky, alpha, oper);

	if (NINT(dz) == 5) {
		if (kp < 0.20*M_PI/dx) *wfacto = 9e-5;
		else if (kp < 0.30*M_PI/dx) *wfacto = 3e-5;
		else if (kp < 0.40*M_PI/dx) *wfacto = 1e-4;
		else if (kp < 0.80*M_PI/dx) *wfacto = 7e-5;
		else if (kp < 1.00*M_PI/dx) *wfacto = 1e-4;
		else if (kp < 1.15*M_PI/dx) *wfacto = 3e-5;
		else if (kp < 1.40*M_PI/dx) *wfacto = 2e-5;
		else *wfacto = 5e-6;
	}
	else if (NINT(dz) == 6) {
		if (kp < 0.30*M_PI/dx) *wfacto = 6e-5;
		else if (kp < 0.80*M_PI/dx) *wfacto = 8e-5;
		else if (kp < 1.00*M_PI/dx) *wfacto = 1e-4;
		else if (kp < 1.15*M_PI/dx) *wfacto = 6e-5;
		else if (kp < 1.45*M_PI/dx) *wfacto = 2e-5;
		else *wfacto = 5e-6;
	}
/*	else if (NINT(dz) == 7) {
		if (kp < 0.30*M_PI/dx) *wfacto = 4e-5;
		else if (kp < 0.50*M_PI/dx) *wfacto = 7e-5;
		else if (kp < 0.88*M_PI/dx) *wfacto = 6e-5;
		else if (kp < 1.00*M_PI/dx) *wfacto = 2e-4;
		else if (kp < 1.15*M_PI/dx) *wfacto = 7e-5;
		else if (kp < 1.30*M_PI/dx) *wfacto = 3e-5;
		else if (kp < 1.40*M_PI/dx) *wfacto = 1e-5;
		else *wfacto = 7e-6;
	}*/
	else if (NINT(dz) == 10) {
		if (kp < 0.25*M_PI/dx) *wfacto = 9e-5;
		else if (kp < 0.85*M_PI/dx) *wfacto = 4e-5;
		else if (kp < 1.10*M_PI/dx) *wfacto = 3e-4;
		else if (kp < 1.20*M_PI/dx) *wfacto = 1e-4;
		else if (kp < 1.30*M_PI/dx) *wfacto = 5e-5;
		else if (kp < 1.50*M_PI/dx) *wfacto = 2e-5;
		else *wfacto = 1e-5;
	}
/* optimw oper=1 fmin=5 fmax=100 dx=20 dy=20 nt=630 oplx=19 cp=1500 
** df=0.22 dz=7.5 verbose=1 alpha=50 */
	else if (NINT(dx) == 20) {
		if (kp < 0.70*M_PI/dx) *wfacto = 1e-4;
		else if (kp < 0.85*M_PI/dx) *wfacto = 1e-4;
		else if (kp < 1.00*M_PI/dx) *wfacto = 1e-4;
		else if (kp < 1.15*M_PI/dx) *wfacto = 6e-5;
		else if (kp < 1.50*M_PI/dx) *wfacto = 1e-5;
		else *wfacto = 1e-5;
	}
	else {
		if (kp < 0.35*M_PI/dx) *wfacto = 1e-4;
		else if (kp < 0.90*M_PI/dx) *wfacto = 1e-4;
		else if (kp < 1.10*M_PI/dx) *wfacto = 5e-4;
		else if (kp < 1.30*M_PI/dx) *wfacto = 7e-5;
		else *wfacto = 5e-5;
	}


	return;
}


void Woper3d(float k, float dx, float dy, float dz, int nkx, int nky, complex *oper)
{
	int 	ikx, iky;
	float 	kx, kx2, kz2, ky, ky2, k2;
	float 	dkx, dky, kz;

	k2 	= k*k;
	dkx = M_PI/((nkx-1)*dx);
	dky = M_PI/((nky-1)*dy);

	for (iky = 0; iky < nky; iky++) {
		ky   = iky*dky;
		ky2  = ky*ky;
		for (ikx = 0; ikx < nkx; ikx++) {
			kx   = ikx*dkx;
			kx2  = kx*kx;
			kz2 = k2 - (kx2 + ky2);
			if (kz2 > 0) {
				kz = dz*sqrt(kz2);
				oper[iky*nkx+ikx].r = cos(kz);
				oper[iky*nkx+ikx].i = sin(kz);
			}
			else {
				kz = dz*sqrt(-kz2);
				oper[iky*nkx+ikx].r = exp(-kz);
				oper[iky*nkx+ikx].i = 0.0;
			} 
		}
	}

	return;
}

void Woper3d_ph(float k, float dx, float dy, float dz, int nkx2, int nky2, float alpha, complex *oper)
{
	int 	ikx, iky;
	float 	k2, kr, kr2, x1, dzdx1, dzdx2, x2, z1, z2;
	float 	kx, kx2, ky, ky2, a, b, c, d;
	float 	dkx, dky, kmax, *phase, *window;

	k2 	= k*k;
	dkx = M_PI/((nkx2-1)*dx);
	dky = M_PI/((nky2-1)*dy);

	kmax = k*sin(M_PI*alpha/180.0);
	if (kmax > 0.85*M_PI/dx) kmax = 0.85*M_PI/dx;

/* Calculation of cubic spline */

	x1 = dkx*NINT(kmax/dkx);
	z1 = sqrt(k2-x1*x1)*dz;
	dzdx1 = -1.0*x1*dz / sqrt(k2-x1*x1);
	x2 = MIN(M_PI/dx, M_PI/dy);
	z2 = 0.0;
	dzdx2 = 0.0;
	spline3(x1, x2, z1, z2, dzdx1, dzdx2, &a, &b, &c, &d);

/* Calculation of the phase and the window function*/

	phase  = (float *)malloc(nkx2*nky2*sizeof(float));
	window = (float *)malloc(nkx2*nky2*sizeof(float));

	for (iky = 0; iky < nky2; iky++) {
		ky   = iky*dky;
		ky2  = ky*ky;
		for (ikx = 0; ikx < nkx2; ikx++) {
			kx   = ikx*dkx;
			kx2  = kx*kx;
			kr2  = kx2 + ky2;
			kr   = sqrt(kr2);
			if (kr <= kmax) {
				phase[iky*nkx2+ikx] = sqrt(k2-kr2)*dz;
				window[iky*nkx2+ikx] = 1.0;
			}
			else if ((kr > kmax) && (kr <= x2))  {
				phase[iky*nkx2+ikx] = a*kr2*kr + b*kr2 + c*kr + d;
				window[iky*nkx2+ikx] = (cos(M_PI*(kr-kmax)/(x2-kmax))+1.0)/2.0;
			} 
			else {
				phase[iky*nkx2+ikx]  = 0.0;
				window[iky*nkx2+ikx] = 0.0;
			}
		}
	}

	for (iky = 0; iky < nky2; iky++) {
		for (ikx = 0; ikx < nkx2; ikx++) {
			oper[iky*nkx2+ikx].r = cos(phase[iky*nkx2+ikx])*window[iky*nkx2+ikx];
			oper[iky*nkx2+ikx].i = sin(phase[iky*nkx2+ikx])*window[iky*nkx2+ikx];
		}
	}

	free(phase);
	free(window);
	return;
}

