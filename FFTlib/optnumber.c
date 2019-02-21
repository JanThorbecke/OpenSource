#include "genfft.h"
#ifdef SGI
#include "sgintab.h"
int optnfft(int n);
long loptnfft(long n);
#endif

#if defined(CRAY_MPP_64)
static int list[756], first=0;
#endif

int npfar (int nmin);
int npfa (int nmin);
/**
* NAME:        optncc
* 
* DESCRIPTION: function to determine the nearest valid number of points
* 	     in a Complex to Complex fourier transform
* 
* USAGE:	     int optncc(int n)
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands
* ----------------------------------------------------------------------*/

int optncc(int n)
{
#ifdef SGI
	return optnfft(n);
#else
	int n2, n3;

	n2 = pow(2.0, 1.0*(int)(log((float)n)/log(2.0)+0.9999));
	if (n2 != n) {
		n3 = npfa(n);
		if((n3-n) < (n2-n)) return npfa(n);
		else return n2;
	}
	else return n;
#endif

}

/*----------------------------------------------------------------------
NAME:        optncr

DESCRIPTION: function to determine the nearest valid number of points
	     in a Real to Complex and Complex to Real fourier transform

USAGE:	     int optncr(int n)
----------------------------------------------------------------------*/

int optncr(int n)
{

#ifdef SGI
	return optnfft(n);
#else
	int n2, n3;

	n2 = pow(2.0, 1.0*(int)(log((float)n)/log(2.0)+0.9999));
	if (n2 != n) {
		n3 = npfar(n);
		if((n3-n) < (n2-n)) return npfar(n);
		else return n2;
	}
	else return n;
#endif
}

	
#ifdef SGI
int optnfft(int n)
{
	int i,j, nmax;

	if (n > NTAB+3) {
		i=13;
		nmax=NTAB+3;
		while (nmax < n) { nmax=(int)pow(2.0,(double)++i); }
		return nmax;
	}

	nmax = NTAB;
	for (i=0; i<NTAB-1 && ntab[i].n<n; ++i);
	for (j=i+1; j<NTAB-1 && ntab[j].n<=n+nmax; ++j)
	if (ntab[j].c<ntab[i].c) i = j;
	return ntab[i].n;
}
#endif

long loptncr(long n)
{

#ifdef SGI
	return loptnfft(n);
#else
	long n2, n3;

	n2 = pow(2.0, 1.0*(long)(log((float)n)/log(2.0)+0.9999));
	if (n2 != n) {
		n3 = npfar(n);
		if((n3-n) < (n2-n)) return npfar(n);
		else return n2;
	}
	else return n;
#endif
}

	
#ifdef SGI
long loptnfft(long n)
{
	long i,j, nmax;

	if (n > NTAB+3) {
		i=13;
		nmax=NTAB+3;
		while (nmax < n) { nmax=(long)pow(2.0,(double)++i); }
		return nmax;
	}

	nmax = NTAB;
	for (i=0; i<NTAB-1 && ntab[i].n<n; ++i);
	for (j=i+1; j<NTAB-1 && ntab[j].n<=n+nmax; ++j)
	if (ntab[j].c<ntab[i].c) i = j;
	return ntab[i].n;
}
#endif

#if defined(CRAY_MPP_64)
int factorized(int n)
{
 	int m, i, j, k;

	if (!first) {
		m = 0;
		for (i=0; i<=5; i++) {
			for (j=0; j<=8; j++) {
				for (k=0; k<=12; k++) {
					list[m++] = pow(2.0,k)*pow(3.0,j)*pow(5.0,i);
				}
			}
		}
		first = 1;
	}

	for (i=0; i<756; i++) {
		if (n == list[i]) return 1;
	}
	return 0;
}
#endif


/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define noptncc	FNAME(OPTNCCF)
#else
#define noptncc	FNAME(optnccf)
#endif

int noptncc(int *n)
{
	int optn;

	optn = optncc(*n);
	return optn;
}


#ifdef DF_CAPFNAMES
#define noptncr	FNAME(OPTNCRF)
#else
#define noptncr	FNAME(optncrf)
#endif

int noptncr(int *n)
{
	int optn;

	optn = optncr(*n);
	return optn;
}

