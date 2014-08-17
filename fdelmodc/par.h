/* This file is property of the Colorado School of Mines.
 
 Copyright © 2007, Colorado School of Mines,
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
 .                       */
/* par.h - include file for getpar, selfdoc, and error handling functions */

#ifndef PAR_H
#define PAR_H
void verr(char *fmt, ...);
void vmess(char *fmt, ...);
void vwarn(char *fmt, ...);
void vsyserr(char *fmt, ...);

/* TYPEDEFS */
typedef union { /* storage for arbitrary type */
    char s[8];
    short h;
    unsigned short u;
    long l;
    unsigned long v;
    int i;
    unsigned int p;
    float f;
    double d;
    unsigned int U:16;
    unsigned int P:32;
} Value;

/* INCLUDES */

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>	/* non-ANSI */
#include <unistd.h>	/* non-ANSI */
#include <sys/types.h>	/* non-ANSI */
#include<string.h>


/* GLOBAL DECLARATIONS */
extern int xargc; extern char **xargv;


/* TYPEDEFS */
typedef char *cwp_String;

typedef enum {BADFILETYPE = -1,
	TTY, DISK, DIRECTORY, TAPE, PIPE, FIFO, SOCKET, SYMLINK} FileType;
	
/* DEFINES */

/* getpar macros */
#define MUSTGETPARINT(x,y) if(!getparint(x,y)) err("must specify %s=",x)
#define MUSTGETPARFLOAT(x,y) if(!getparfloat(x,y)) err("must specify %s=",x)
#define MUSTGETPARSTRING(x,y) if(!getparstring(x,y)) err("must specify %s=",x)

#define STDIN (0)
#define STDOUT (1)
#define STDERR (2)

/* FUNCTION PROTOTYPES */

#ifdef __cplusplus  /* if C++, specify external C linkage */
extern "C" {
#endif

/* getpar parameter parsing */
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

/* For ProMAX */
void getPar(char *name, char *type, void *ptr);

/* errors and warnings */
void err (char *fmt, ...);
void syserr (char *fmt, ...);
void warn (char *fmt, ...);

/* self documentation */
void pagedoc (void);
void requestdoc (int i);

/* system calls with error trapping */
int ecreat(char *path, int perms);
int efork(void);
int eopen(char *path, int flags, int perms);
int eclose(int fd);
int eunlink(char *path);
long elseek(int fd, long offset, int origin);
int epipe(int fd[2]);
int eread(int fd, char *buf, int nbytes);
int ewrite(int fd, char *buf, int nbytes);

/* system subroutine calls with error trapping */
FILE *efopen(const char *file, const char *mode);
FILE *efreopen(const char *file, const char *mode, FILE *stream1);
FILE *efdopen(int fd, const char *mode);
FILE *epopen(char *command, char *type);
int efclose(FILE *stream);
int epclose(FILE *stream);
int efflush(FILE *stream);
int eremove(const char *file);
int erename(const char *oldfile, const char* newfile);
int efseek(FILE *stream, long offset, int origin);
void erewind(FILE *stream);
long eftell(FILE *stream);
FILE *etmpfile(void);
char *etmpnam(char *namebuffer);
void *emalloc(size_t size);
void *erealloc(void *memptr, size_t size);
void *ecalloc(size_t count, size_t size);
size_t efread(void *bufptr, size_t size, size_t count, FILE *stream);
size_t efwrite(void *bufptr, size_t size, size_t count, FILE *stream);

#ifndef SUN_A
int efgetpos(FILE *stream, fpos_t *position);
int efsetpos(FILE *stream, const fpos_t *position);
#endif

/* string to numeric conversion with error checking */
short eatoh(char *s);
unsigned short eatou(char *s);
int eatoi(char *s);
unsigned int eatop(char *s);
long eatol(char *s);
unsigned long eatov(char *s);
float eatof(char *s);
double eatod(char *s);

/* file type checking */
FileType filestat(int fd);
char *printstat(int fd);

#ifdef __cplusplus  /* if C++ (external C linkage is being specified) */
}
#endif

#endif /* PAR_H */
