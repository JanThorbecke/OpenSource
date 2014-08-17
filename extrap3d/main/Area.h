
typedef struct area {
	float xmin;
	float ymin;
	float zmin;
	int ixmin;
	int ixmax;
	int iymin;
	int iymax;
	int nx;
	int ny;
	int nz;
	float dx;
	float dy;
	float dz;
	int sxy;
} Area;

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#define COMPLEX
#endif/* complex */

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define ISODD(n) ((n) & 01)



