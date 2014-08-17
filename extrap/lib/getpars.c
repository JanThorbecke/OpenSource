/* This file is property of the Colorado School of Mines.
 
 Copyright (C) 2007, Colorado School of Mines,
 All rights reserved.
 
 
 Redistribution and use in source and binary forms, with or 
 without modification, are permitted provided that the following 
 conditions are met:
 
 *  Redistributions of source code must retain the above copyright 
 notice, this list of conditions and the following disclaimer.
 *  Redistributions in binary form must reproduce the above 
 copyright notice, this list of conditions and the following 
 disclaimer in the documentation and/or other materials provided 
 with the distribution.
 *  Neither the name of the Colorado School of Mines nor the names of
 its contributors may be used to endorse or promote products 
 derived from this software without specific prior written permission.
 
 Warranty Disclaimer:
 THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS 
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
 COLORADO SCHOOL OF MINES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
 STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
 IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 POSSIBILITY OF SUCH DAMAGE.
 
 
 Export Restriction Disclaimer:
 We believe that CWP/SU: Seismic Un*x is a low technology product that does
 not appear on the Department of Commerce CCL list of restricted exports.
 Accordingly, we believe that our product meets the qualifications of
 an ECCN (export control classification number) of EAR99 and we believe
 it fits the qualifications of NRR (no restrictions required), and
 is thus not subject to export restrictions of any variety.
 
 Approved Reference Format:
 In publications, please refer to SU as per the following example:
 Cohen, J. K. and Stockwell, Jr. J. W., (200_), CWP/SU: Seismic Un*x 
 Release No. __: an open source software  package for seismic 
 research and processing, 
 Center for Wave Phenomena, Colorado School of Mines.
 
 Articles about SU in peer-reviewed journals:
 Saeki, T., (1999), A guide to Seismic Un*x (SU)(2)---examples of data processing (part 1), data input and preparation of headers, Butsuri-Tansa (Geophysical Exploration), vol. 52, no. 5, 465-477.
 Stockwell, Jr. J. W. (1999), The CWP/SU: Seismic Un*x Package, Computers and Geosciences, May 1999.
 Stockwell, Jr. J. W. (1997), Free Software in Education: A case study of CWP/SU: Seismic Un*x, The Leading Edge, July 1997.
 Templeton, M. E., Gough, C.A., (1998), Web Seismic Un*x: Making seismic reflection processing more accessible, Computers and Geosciences.
 
 Acknowledgements:
 SU stands for CWP/SU:Seismic Un*x, a processing line developed at Colorado 
 School of Mines, partially based on Stanford Exploration Project (SEP) 
 software.
 */

/*********************** self documentation **********************/
/*****************************************************************************
GETPARS - Functions to GET PARameterS from the command line. Numeric
	parameters may be single values or arrays of int, uint,
	short, ushort, long, ulong, float, or double.  Single character
	strings (type string or char *) may also be gotten. 
	Arrays of strings, delimited by, but not containing
        commas are permitted.

The functions are:

initargs 	Makes command line args available to subroutines (re-entrant).
		Every par program starts with this call!
getparint		get integers
getparuint		get unsigned integers
getparshort		get short integers
getparushort		get unsigned short integers
getparlong		get long integers 
getparulong		get unsigned long integers
getparfloat		get float
getpardouble		get double
getparstring		get a single string
getparstringarray	get string array (fields delimited by commas) 
getpar			get parameter by type
getnparint		get n'th occurrence of integer
getnparuint		get n'th occurrence of unsigned int
getnparshort		get n'th occurrence of short integer
getnparushort		get n'th occurrence of unsigned short int
getnparlong		get n'th occurrence of long integer
getnparulong		get n'th occurrence of unsigned long int
getnparfloat		get n'th occurrence of float integer
getnpardouble		get n'th occurrence of double integer
getnparstring		get n'th occurrence of string integer
getnparstringarray	get n'th occurrence of string integer array
getnpar			get n'th occurrence by type
countparname		return the number of times a parameter names is used
countparval		return the number of values in the last occurrence
				of a parameter
countnparval		return the number of values in the n'th occurrence
				of a parameter
getPar			Promax compatible version of getpar

******************************************************************************
Function Prototypes:
void initargs (int argc, char **argv);
int getparint (char *name, int *p);
int getparuint (char *name, unsigned int *p);
int getparshort (char *name, short *p);
int getparushort (char *name, unsigned short *p);
int getparlong (char *name, long *p);
int getparulong (char *name, unsigned long *p);
int getparfloat (char *name, float *p);
int getpardouble (char *name, double *p);
int getparstring (char *name, char **p);
int getparstringarray (char *name, char **p);
int getnparint (int n, char *name, int *p);
int getnparuint (int n, char *name, unsigned int *p);
int getnparshort (int n, char *name, short *p);
int getnparushort (int n, char *name, unsigned short *p);
int getnparlong (int n, char *name, long *p);
int getnparulong (int n, char *name, unsigned long *p);
int getnparfloat (int n, char *name, float *p);
int getnpardouble (int n, char *name, double *p);
int getnparstring (int n, char *name, char **p);
int getnparstringarray (int n, char *name, char **p);
int getnpar (int n, char *name, char *type, void *ptr);
int countparname (char *name);
int countparval (char *name);
int countnparval (int n, char *name);
void getPar(char *name, char *type, void *ptr);

******************************************************************************
Notes:
Here are some usage examples:

	... if integer n not specified, then default to zero. 
	if (!getparint("n", &n)) n = 0;
	
	... if array of floats vx is specified, then
	if (nx=countparval("vx")) {
		... allocate space for array
		vx = (float *)malloc(nx*sizeof(float));
		... and get the floats
		getparfloat("vx",vx);
	}
	
The command line for the above examples might look like:
	progname n=35 vx=3.21,4,9.5
	Every par program starts with this call!

More examples are provided in the DTEST code at the end of this file.

The functions: eatoh, eatou, eatol, eatov, eatoi, eatop used
below are versions of atoi that check for overflow.  The source
file for these functions is atopkge.c.

******************************************************************************	
Authors:
Rob Clayton & Jon Claerbout, Stanford University, 1979-1985
Shuki Ronen & Jack Cohen, Colorado School of Mines, 1985-1990
Dave Hale, Colorado School of Mines, 05/29/90
Credit to John E. Anderson for re-entrant initargs 03/03/94
*****************************************************************************/	
/**************** end self doc ********************************/

#include "par.h"

#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

/* parameter table */
typedef struct {
	char *name;		/* external name of parameter	*/
	char *asciival;		/* ascii value of parameter	*/
} pointer_table;

/* global variables declared and used internally */
static pointer_table *argtbl;	/* parameter table		*/
static int nargs;		/* number of args that parse	*/
static int tabled = FALSE;	/* true when parameters tabled 	*/
static int targc;		/* total number of args		*/
static char **targv;		/* pointer to arg strings	*/
static char *argstr;		/* storage for command line	*/

/* functions declared and used internally */
static int getparindex (int n, char *name);
static void getparinit(void);
static void tabulate (int argc, char **argv);
static char *getpfname (void);
static int white2null (char *str, int len);
static int ccount (char c, char *s);
static void strchop(char *s, char *t);

/* make command line args available to subroutines -- re-entrant version */
void initargs(int argc, char **argv)
{
	xargc = argc; xargv = argv;
	if(tabled==TRUE){
		free(argstr);
		free(targv);
		free(argtbl);
	}
	tabled =  FALSE;
	return;
}

/* functions to get values for the last occurrence of a parameter name */
int getparint (char *name, int *ptr)
{
	return getnpar(0,name,"i",ptr);
}
int getparuint (char *name, unsigned int *ptr)
{
	return getnpar(0,name,"p",ptr);
}
int getparshort (char *name, short *ptr)
{
	return getnpar(0,name,"h",ptr);
}
int getparushort (char *name, unsigned short *ptr)
{
	return getnpar(0,name,"u",ptr);
}
int getparlong (char *name, long *ptr)
{
	return getnpar(0,name,"l",ptr);
}
int getparulong (char *name, unsigned long *ptr)
{
	return getnpar(0,name,"v",ptr);
}
int getparfloat (char *name, float *ptr)
{
	return getnpar(0,name,"f",ptr);
}
int getpardouble (char *name, double *ptr)
{
	return getnpar(0,name,"d",ptr);
}
int getparstring (char *name, char **ptr)
{
	return getnpar(0,name,"s",ptr);
}
int getparstringarray (char *name, char **ptr)
{
	return getnpar(0,name,"a",ptr);
}
int getpar (char *name, char *type, void *ptr)
{
	return getnpar(0,name,type,ptr);
}

/* functions to get values for the n'th occurrence of a parameter name */
int getnparint (int n, char *name, int *ptr)
{
	return getnpar(n,name,"i",ptr);
}
int getnparuint (int n, char *name, unsigned int *ptr)
{
	return getnpar(n,name,"p",ptr);
}
int getnparshort (int n, char *name, short *ptr)
{
	return getnpar(n,name,"h",ptr);
}
int getnparushort (int n, char *name, unsigned short *ptr)
{
	return getnpar(n,name,"u",ptr);
}
int getnparlong (int n, char *name, long *ptr)
{
	return getnpar(n,name,"l",ptr);
}
int getnparulong (int n, char *name, unsigned long *ptr)
{
	return getnpar(n,name,"v",ptr);
}
int getnparfloat (int n, char *name, float *ptr)
{
	return getnpar(n,name,"f",ptr);
}
int getnpardouble (int n, char *name, double *ptr)
{
	return getnpar(n,name,"d",ptr);
}
int getnparstring (int n, char *name, char **ptr)
{
	return getnpar(n,name,"s",ptr);
}
int getnparstringarray (int n, char *name, char **ptr)
{
	return getnpar(n,name,"a",ptr);
}
int getnpar (int n, char *name, char *type, void *ptr)
{
	int i;			/* index of name in symbol table	*/
	int nval;		/* number of parameter values found	*/
	char *aval;		/* ascii field of symbol		*/

	if (xargc == 1) return 0;
	if (!tabled) getparinit();/* Tabulate command line and parfile */
	i = getparindex(n,name);/* Get parameter index */
	if (i < 0) return 0;	/* Not there */
	
	/* 
	 * handle string type as a special case, since a string 
	 * may contain commas. 
	 */
	if (type[0]=='s') {
		*((char**)ptr) = argtbl[i].asciival;
		return 1;
	} 

	/* convert vector of ascii values to numeric values */
	for (nval=0,aval=argtbl[i].asciival; *aval; nval++) {
		switch (type[0]) {
			case 'i':
				*(int*)ptr = eatoi(aval);
				ptr = (int*)ptr+1;
				break;
			case 'p':
				*(unsigned int*)ptr = eatop(aval);
				ptr = (unsigned int*)ptr+1;
				break;
			case 'h':
				*(short*)ptr = eatoh(aval);
				ptr = (short*)ptr+1;
				break;
			case 'u':
				*(unsigned short*)ptr = eatou(aval);
				ptr = (unsigned short*)ptr+1;
				break;
			case 'l':
				*(long*)ptr = eatol(aval);
				ptr = (long*)ptr+1;
				break;
			case 'v':
				*(unsigned long*)ptr = eatov(aval);
				ptr = (unsigned long*)ptr+1;
				break;
			case 'f':
				*(float*)ptr = eatof(aval);
				ptr = (float*)ptr+1;
				break;
			case 'd':
				*(double*)ptr = eatod(aval);
				ptr = (double*)ptr+1;
				break;
			case 'a':
				{ char *tmpstr="";
				   tmpstr = (char *)calloc(strlen(aval),1);

				   strchop(aval,tmpstr);
				   *(char**)ptr = tmpstr;
				   ptr=(char **)ptr + 1;
				}
				   break;
			default:
				err("%s: invalid parameter type = %s",
					__FILE__,type);
		}
		while (*aval++ != ',') {
			if (!*aval) break;
		}
	}
	return nval;
}
/* Promax compatible version of getnpar */
void getPar(char *name, char *type, void *ptr)
{
	(void) getnpar(0,name,type,ptr);
	return;
}

/* return number of occurrences of parameter name */
int countparname (char *name)
{
	int i,nname;

	if (xargc == 1) return 0;
	if (!tabled) getparinit();
	for (i=0,nname=0; i<nargs; ++i)
		if (!strcmp(name,argtbl[i].name)) ++nname;
	return nname;
}

/* return number of values in n'th occurrence of parameter name */
int countnparval (int n, char *name)
{
	int i;

	if (xargc == 1) return 0;
	if (!tabled) getparinit();
	i = getparindex(n,name);
	if (i>=0) 
		return ccount(',',argtbl[i].asciival) + 1;
	else
		return 0;
}

/* return number of values in last occurrence of parameter name */
int countparval (char *name)
{
	return countnparval(0,name);
}



/*
 * Return the index of the n'th occurrence of a parameter name, 
 * except if n==0, return the index of the last occurrence.
 * Return -1 if the specified occurrence does not exist.
 */
static int getparindex (int n, char *name)
{
	int i;
	if (n==0) {
		for (i=nargs-1; i>=0; --i)
			if (!strcmp(name,argtbl[i].name)) break;
		return i;
	} else {
		for (i=0; i<nargs; ++i)
			if (!strcmp(name,argtbl[i].name))
				if (--n==0) break;
		if (i<nargs)
			return i;
		else
			return -1;
	}
}

/* Initialize getpar */
static void getparinit (void)
{
	static char *pfname;	/* name of parameter file		*/
	FILE *pffd=NULL;	/* file id of parameter file		*/
	int pflen;		/* length of parameter file in bytes	*/ 
	static int pfargc;	/* arg count from parameter file	*/
	int parfile;		/* parfile existence flag		*/
	int argstrlen;
	char *pargstr;		/* storage for parameter file args	*/
	int nread;		/* bytes fread				*/
	int i, j;		/* counters				*/


	tabled = TRUE;		/* remember table is built		*/

	/* Check if xargc was initiated */
	if(!xargc)
		err("%s: xargc=%d -- not initiated in main", __FILE__, xargc);

	/* Space needed for command lines */
	for (i = 1, argstrlen = 0; i < xargc; i++) {
		argstrlen += strlen(xargv[i]) + 1;
	}

	/* Get parfile name if there is one */
	/* parfile = (pfname = getpfname()) ? TRUE : FALSE; */
	if ((pfname = getpfname())) {
		parfile = TRUE;
	} else {
		parfile = FALSE;
	}

	if (parfile) {
	 	pffd = fopen(pfname, "r");

		/* Get the length */
		fseek(pffd, 0, SEEK_END);
		pflen = ftell(pffd);
		rewind(pffd);
		argstrlen += pflen;
	} else {
		pflen = 0;
	}

	/* Allocate space for command line and parameter file
		plus nulls at the ends to help with parsing. */
	/* argstr = (char *) calloc((size_t) (1+argstrlen+1), 1); */
	/*argstr = (char *) ealloc1(1+argstrlen+1, 1);*/
	argstr = (char *) calloc((size_t) (1+argstrlen+1), 1);

	if (parfile) {
		/* Read the parfile */
		nread = fread(argstr + 1, 1, pflen, pffd);
  		if (nread != pflen) {
  	 	    err("%s: fread only %d bytes out of %d from %s",
  					__FILE__,  nread, pflen, pfname);
		}
		fclose(pffd);

		/* Zap whites in parfile to help in parsing */
		pfargc = white2null(argstr, pflen);

	} else {
		pfargc = 0;
	}

	/* Total arg count */
	targc = pfargc + xargc - 1;

	/* Allocate space for total arg pointers */
	targv = (char **) calloc(targc, sizeof(char*));

	if (parfile) {
		/* Parse the parfile.  Skip over multiple NULLs */
		for (j = 1, i = 0; j < pflen; j++) {
			if (argstr[j] && !argstr[j-1]) {
			       targv[i++] = argstr + j;
			}
		}
	} else {
		i = 0;
	}

	/* Copy command line arguments */
	for (j = 1, pargstr = argstr + pflen + 2; j < xargc; j++) {
		strcpy(pargstr,xargv[j]);
		targv[i++] = pargstr;
		pargstr += strlen(xargv[j]) + 1;
	}

	/* Allocate space for the pointer table */
	argtbl = (pointer_table*) calloc(targc, sizeof(pointer_table));

	/* Tabulate targv */
	tabulate(targc, targv);
	
	return;
}

#define PFNAME "par="
/* Get name of parameter file */
static char *getpfname (void)
{
	int i;
	int pfnamelen;

	pfnamelen = strlen(PFNAME);
	for (i = xargc-1 ; i > 0 ; i--) {
		if(!strncmp(PFNAME, xargv[i], pfnamelen)
		    && strlen(xargv[i]) != pfnamelen) {
			return xargv[i] + pfnamelen;
		}	
	}
	return NULL;
}

#define iswhite(c)	((c) == ' ' || (c) == '\t' || (c) == '\n')

/* 
 * Replace the whites by (possibly multiple) nulls.  If we see a non-white
 * and the previous char is a null, this signals the start of a string
 * and we bump the count.  This routine returns a count of the strings.
 */
static int white2null (char *str, int len)
{
	int i;
	int count;
	int inquote = FALSE;

	str[0] = '\0'; /* This line added by Dave Hale, 1/30/96. */
	for (i = 1, count = 0; i < len; i++) {
		if (str[i]=='"') inquote=(inquote==TRUE)?FALSE:TRUE;
		if (!inquote) {
			if (iswhite(str[i])) { /* Is this a new word ? */
				str[i] = '\0';
			} else if (!str[i-1]) { /* multiple whites */
				count++;
			}
		}
	}
	for (i = 1, inquote=FALSE; i < len; i++) {
		if (str[i]=='"') inquote=(inquote==TRUE)?FALSE:TRUE;
		if (inquote) {
			if (str[i+1]!='"') {
				str[i] = str[i+1];
			} else {
				str[i] = '\0';
				str[i+1] = '\0';
				inquote = FALSE;
			}
		}
	}
	str[len] = '\0';
	return count;
}

/* Install symbol table */
static void tabulate (int argc, char **argv)
{
	int i;
	char *eqptr;

	for (i = 0, nargs = 0 ; i < argc; i++) {
		eqptr = (char *)strchr(argv[i], '=');
		if (eqptr) {
			argtbl[nargs].name = argv[i];
			argtbl[nargs].asciival = eqptr + 1;
			*eqptr = (char)0;

			/* Debugging dump */
/* 			fprintf(stderr, */
/* 			"argtbl[%d]: name=%s asciival=%s\n", */
/* 			nargs,argtbl[nargs].name,argtbl[nargs].asciival); */

			nargs++;
		}
	}
	return;
}

/* Count characters in a string */
static int ccount (char c, char *s)
{
	int i, count;
	for (i = 0, count = 0; s[i] != 0; i++)
		if(s[i] == c) count++;
	return count;
}

static void strchop(char *s, char *t)
/***********************************************************************
strchop - chop off the tail end of a string "s" after a "," returning
          the front part of "s" as "t".
************************************************************************
Notes:
Based on strcpy in Kernighan and Ritchie's C [ANSI C] book, p. 106.
************************************************************************
Author: CWP: John Stockwell and Jack K. Cohen, July 1995
***********************************************************************/
{

	while ( (*s != ',') && (*s != '\0') ) {
		 *t++ = *s++;
	}
	*t='\0';
}


#ifdef TEST
#define N 100
main(int argc, char **argv)
{
	char *s;
	short h, vh[N];
	unsigned short u, vu[N];
	long l, vl[N];
	unsigned long v, vv[N];
	int i, vi[N], ipar, npar, nval;
	unsigned int p, vp[N];
	float f, vf[N];
	double d, vd[N];

	initargs(argc, argv);

	/* int parameters */
	npar = countparname("i");
	printf("\nnumber of i pars = %d\n",npar);
	for (ipar=1; ipar<=npar; ++ipar) {
		getnparint(ipar,"i",&i);
		printf("occurrence %d of i=%d\n",ipar,i);
	}
	if (getparint("i", &i))	
		printf("last occurrence of i=%d\n",i);
	npar = countparname("vi");
	printf("number of vi pars = %d\n",npar);
	for (ipar=1; ipar<=npar; ++ipar) {
		nval = countnparval(ipar,"vi");
		printf("occurrence %d has %d values\n",ipar,nval);
		nval = getnparint(ipar,"vi",vi);
		printf("vi=");
		for (i=0; i<nval; i++)
			printf("%d%c",vi[i],i==nval-1?'\n':',');
	}
	if (npar>0) {
		nval = countparval("vi");
		printf("last occurrence has %d values\n",nval);
		getparint("vi",vi);
		printf("vi=");
		for (i=0; i<nval; i++)
			printf("%d%c",vi[i],i==nval-1?'\n':',');
	}

	/* float parameters */
	npar = countparname("f");
	printf("\nnumber of f pars = %d\n",npar);
	for (ipar=1; ipar<=npar; ++ipar) {
		getnparfloat(ipar,"f",&f);
		printf("occurrence %d of f=%g\n",ipar,f);
	}
	if (getparfloat("f", &f))	
		printf("last occurrence of f=%g\n",f);
	npar = countparname("vf");
	printf("number of vf pars = %d\n",npar);
	for (ipar=1; ipar<=npar; ++ipar) {
		nval = countnparval(ipar,"vf");
		printf("occurrence %d has %d values\n",ipar,nval);
		nval = getnparfloat(ipar,"vf",vf);
		printf("vf=");
		for (i=0; i<nval; i++)
			printf("%g%c",vf[i],i==nval-1?'\n':',');
	}
	if (npar>0) {
		nval = countparval("vf");
		printf("last occurrence has %d values\n",nval);
		getparfloat("vf",vf);
		printf("vf=");
		for (i=0; i<nval; i++)
			printf("%g%c",vf[i],i==nval-1?'\n':',');
	}

	/* string parameters */
	npar = countparname("s");
	printf("\nnumber of s pars = %d\n",npar);
	for (ipar=1; ipar<=npar; ++ipar) {
		getnparstring(ipar,"s",&s);
		printf("occurrence %d of s=%s\n",ipar,s);
	}
	if (getparstring("s", &s))	
		printf("last occurrence of s=%s\n",s);
	
	return EXIT_SUCCESS;
}
#endif

