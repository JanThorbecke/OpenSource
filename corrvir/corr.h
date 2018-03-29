#include<stdlib.h>
#include<stdio.h>
#include<math.h>

typedef struct _tracePos { /* Type */
	int gx;
	int gy;
	int gelev;
} tracePos;

typedef struct _traceCoord { /* Type */
	int x;
	int y;
	int peg;
	size_t fpos;
} traceCoord;

typedef struct _crgPos { /* Type */
	int gx;
	int gy;
	int rpeg;
	int gelev;
	int nsrc;
	traceCoord *src;
} crgPos;

