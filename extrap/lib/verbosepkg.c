#include <stdio.h>
#include <stdarg.h>
#include "par.h"
#include <string.h>
#ifdef _CRAYMPP
#include <intrinsics.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

void verr(char *fmt, ...)
{
	va_list args;

	if (EOF == fflush(stdout)) {
		fprintf(stderr, "\nverr: fflush failed on stdout");
	}
	fprintf(stderr, "    Error in %s: ", xargv[0]);
#ifdef _OPENMP
        fprintf(stderr, "PE %d: ", omp_get_thread_num());
#endif
	va_start(args,fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");

	exit(EXIT_FAILURE);
}

void vwarn(char *fmt, ...)
{
	va_list args;

	if (EOF == fflush(stdout)) {
		fprintf(stderr, "\nvwarn: fflush failed on stdout");
	}
	fprintf(stderr, "    Warning in %s: ", xargv[0]);
#ifdef _OPENMP
        fprintf(stderr, "PE %d: ", omp_get_thread_num());
#endif
	va_start(args,fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	return;
}

void vmess(char *fmt, ...)
{
	va_list args;

	if (EOF == fflush(stdout)) {
		fprintf(stderr, "\nvmess: fflush failed on stdout");
	}
	fprintf(stderr, "    %s: ", xargv[0]);
#ifdef _OPENMP
        fprintf(stderr, "PE %d: ", omp_get_thread_num());
#endif
	va_start(args,fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	return;
}
