#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include"par.h"

#define SQR2 1.414213562
#define SQR3 1.732050808
#define SQR6 2.449489743
// #define t0(x,y,z)   time0[nxz*(y) + nz*(x) + (z)]
// #define s0(x,y,z)   slow0[nxz*(y) + nz*(x) + (z)]
#define t0(x,y,z)   time0[nxy*(z) + nx*(y) + (x)]
#define s0(x,y,z)   slow0[nxy*(z) + nx*(y) + (x)]

/* definitions from verbose.c */
extern void verr(char *fmt, ...);
extern void vwarn(char *fmt, ...);
extern void vmess(char *fmt, ...);

struct sorted
	{ float time; long i1, i2;};

// int compar(struct sorted *a,struct sorted *b);
int compar(const void * a, const void * b);

float fdhne(float t1,float t2,float t3,float t4,float t5,float ss0,float s1,float s2,float s3);
float fdh3d(float t1,float t2,float t3,float t4,float t5,float t6,float t7,float ss0,float s1,float s2,float s3,float s4,float s5,float s6,float s7);
float fdh2d(float t1,float t2,float t3,float ss0,float s1,float s2,float s3);
float fdhnf(float t1,float t2,float t3,float t4,float t5,float ss0,float s1);

void vidale3d(float *slow0, float *time0, long nz, long nx, long ny, float h, long xs, long ys, long zs, long NCUBE)
{
	long
		srctype=1,	/* if 1, source is a point;
						2, source is on the walls of the data volume;
						3, source on wall, time field known; */
		srcwall,	/* if 1, source on x=0 wall, if 2, on x=nx-1 wall
						if 3, source on y=0 wall, if 4, on y=ny-1 wall
						if 5, source on z=0 wall, if 6, on z=nz-1 wall */
		iplus=1,	/* rate of expansion of "cube" in the */
		iminus=1,	/*    plus/minus x/y/z direction */
		jplus=1,
		jminus=1,
		kplus=1,
		kminus=1,
		igrow,		/* counter for "cube" growth */
		X1, X2, lasti, index, ii, i, j, k, radius, 
		nxy, nyz, nxz, nxyz, nwall,
		/* counters for the position of the sides of current cube */
		x1, x2, y1, y2, z1, z2,
		/* flags set to 1 until a side has been reached */
		dx1=1, dx2=1, dy1=1, dy2=1, dz1=1, dz2=1, rad0=1,
		maxrad,		/* maximum radius to compute */
		reverse=1,	/* will automatically do up to this number of
						reverse propagation steps to fix waves that travel 
						back into expanding cube */
		headpref=6,	/* if headpref starts > 0, will determine 
						model wall closest to source and will prefer to start
						reverse calculations on opposite wall */
		/* counters for detecting head waves on sides of current cube */
		head,headw[7], verbose;
	float
		*wall,
		guess, try,
		/* maximum offset (real units) to compute */
		maxoff = -1.,
		/* used to detect head waves:  if headwave operator decreases 
		   the previously-computed traveltime by at least 
		   headtest*<~time_across_cell> then the headwave counter is 
		   triggered */
		fhead,headtest=1.e-3;

	/* ARRAY TO ORDER SIDE FOR SOLUTION IN THE RIGHT ORDER */
	struct sorted *sort;

	if(!getparlong("verbose",&verbose)) verbose=0;
	if(!getparfloat("maxoff",&maxoff)) maxoff = -1.;
	if(!getparlong("iminus",&iminus)) iminus=1;
	if(!getparlong("iplus",&iplus)) iplus=1;
	if(!getparlong("jminus",&jminus)) jminus=1;
	if(!getparlong("jplus",&jplus)) jplus=1;
	if(!getparlong("kminus",&kminus)) kminus=1;
	if(!getparlong("kplus",&kplus)) kplus=1;
	if(!getparlong("reverse",&reverse)) reverse=0;
	if(!getparlong("headpref",&headpref)) headpref=6;
	if(!getparlong("NCUBE",&NCUBE)) NCUBE=2;

	/* SET MAXIMUM RADIUS TO COMPUTE */
	if (maxoff > 0.) {
		maxrad = maxoff/h + 1;
		vwarn("WARNING: Computing only to max radius = %li",maxrad);
	}
	else maxrad = 99999999;

	nxy = nx * ny;
	nyz = ny * nz;
	nxz = nx * nz;
	nxyz = nx * ny * nz;

	/* MAKE ARRAY SORT LARGE ENOUGH FOR ANY SIDE */
	if(nx <= ny && nx <= nz)  {
		sort = (struct sorted *) malloc(sizeof(struct sorted)*ny*nz);
		nwall = nyz;
	}
	else if(ny <= nx && ny <= nz)  {
		sort = (struct sorted *) malloc(sizeof(struct sorted)*nx*nz);
		nwall = nxz;
	}
	else  {
		sort = (struct sorted *) malloc(sizeof(struct sorted)*nx*ny);
		nwall = nxy;
	}
	wall = (float *) malloc(4*nwall);
	if(sort == NULL || wall == NULL) 
		verr("error in allocation of arrays sort and wall");

	if(!getparlong("srctype",&srctype)) srctype=1;
	if(srctype==1) {
		/* SETS LOCATION OF THE SIDES OF THE CUBE	*/
		radius = NCUBE;
		if(xs > NCUBE) x1 = xs - (NCUBE + 1);
		else{ x1 = -1; dx1 = 0;}
		if(xs < nx-(NCUBE + 1)) x2 = xs + (NCUBE + 1);
		else{ x2 = nx; dx2 = 0;}
		if(ys > NCUBE) y1 = ys - (NCUBE + 1);
		else{ y1 = -1; dy1 = 0;}
		if(ys < ny-(NCUBE + 1)) y2 = ys + (NCUBE + 1);
		else{ y2 = ny; dy2 = 0;}
		if(zs > NCUBE) z1 = zs - (NCUBE + 1);
		else{ z1 = -1; dz1 = 0;}
		if(zs < nz-(NCUBE + 1)) z2 = zs + (NCUBE + 1);
		else{ z2 = nz; dz2 = 0;}
	}
	else {
		if (!getparlong("srcwall",&srcwall)) verr("srcwall not given");
		/* SET LOCATIONS OF SIDES OF THE CUBE SO THAT CUBE IS A FACE  */
		radius = 1;
		if (srcwall == 1)	x2=1;
		else	{  x2=nx;	dx2=0;  }
		if (srcwall == 2)	x1=nx-2;
		else	{  x1= -1;	dx1=0;  }
		if (srcwall == 3)	y2=1;
		else	{  y2=ny;	dy2=0;  }
		if (srcwall == 4)	y1=ny-2;
		else	{  y1= -1;	dy1=0;  }
		if (srcwall == 5)	z2=1;
		else	{  z2=nz;	dz2=0;  }
		if (srcwall == 6)	z1=nz-2;
		else	{  z1= -1;	dz1=0;  }
	}

	if (headpref>0) {	/* HOLE - PREFERRED REVERSE DIRECTION */
		head = nx*ny*nz;
		if (nx>5 && x2<=head)   {headpref=2;  head=x2;}
		if (nx>5 && (nx-1-x1)<=head)   {headpref=1;  head=nx-1-x1;}
		if (ny>5 && y2<=head)   {headpref=4;  head=y2;}
		if (ny>5 && (ny-1-y1)<=head)   {headpref=3;  head=ny-1-y1;}
		if (nz>5 && z2<=head)   {headpref=6;  head=z2;}
		if (nz>5 && (nz-1-z1)<=head)   {headpref=5;  head=nz-1-z1;}
	}

	/* BIGGER LOOP - HOLE - ALLOWS AUTOMATIC REVERSE PROPAGATION IF 
		HEAD WAVES ARE ENCOUNTERED ON FACES OF EXPANDING CUBE, 
		ALLOWING WAVES TO TRAVEL BACK INTO THE CUBE */

	while ( reverse > -1 )  {

		headw[1]=0; headw[2]=0; headw[3]=0; headw[4]=0;
		headw[5]=0; headw[6]=0;

	/* BIG LOOP */
	while(rad0 && (dx1 || dx2 || dy1 || dy2 || dz1 || dz2))  {
		/* CALCULATE ON PRIMARY (time0) GRID */

		/* TOP SIDE */
      for (igrow=1;igrow<=kminus;igrow++) {  
	if(dz1){
		ii = 0;
		for(j=y1+1; j<=y2-1; j++){
			for(i=x1+1; i<=x2-1; i++){
				sort[ii].time = t0(i,j,z1+1);
				sort[ii].i1 = i;
				sort[ii].i2 = j;
				ii++;
			}
		}
		qsort((char *)sort,ii,sizeof(struct sorted),compar);
		for(i=0;i<ii;i++){
			X1 = sort[i].i1;
			X2 = sort[i].i2;
			index = z1*nxy + X2*nx + X1;
			lasti = (z1+1)*nxy + X2*nx + X1;
			fhead = 0.;
			guess = time0[index];
                        if(time0[index+1] < 1.e9 && time0[index+nx+1] < 1.e9
			   && time0[index+nx] < 1.e9 && X2<ny-1  && X1<nx-1 ) {
			  try = fdh3d(              t0(X1,X2,z1+1),
				      t0(X1+1,X2,z1+1),t0(X1+1,X2+1,z1+1),t0(X1,X2+1,z1+1),
				      t0(X1+1,X2,z1  ),t0(X1+1,X2+1,z1  ),t0(X1,X2+1,z1  ),
				      s0(X1,X2,z1), s0(X1,X2,z1+1),
				      s0(X1+1,X2,z1+1),s0(X1+1,X2+1,z1+1),s0(X1,X2+1,z1+1),
				      s0(X1+1,X2,z1  ),s0(X1+1,X2+1,z1  ),s0(X1,X2+1,z1  ));
			  if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index+nx-1] < 1.e9
			   && time0[index+nx] < 1.e9 && X2<ny-1  && X1>0 ) {
			  try = fdh3d(              t0(X1,X2,z1+1),
				      t0(X1-1,X2,z1+1),t0(X1-1,X2+1,z1+1),t0(X1,X2+1,z1+1),
				      t0(X1-1,X2,z1  ),t0(X1-1,X2+1,z1  ),t0(X1,X2+1,z1  ),
				      s0(X1,X2,z1), s0(X1,X2,z1+1),
				      s0(X1-1,X2,z1+1),s0(X1-1,X2+1,z1+1),s0(X1,X2+1,z1+1),
				      s0(X1-1,X2,z1  ),s0(X1-1,X2+1,z1  ),s0(X1,X2+1,z1  ));
			  if (try<guess) guess = try;
			}
			if(time0[index+1] < 1.e9 && time0[index-nx+1] < 1.e9
			   && time0[index-nx] < 1.e9 && X2>0  && X1<nx-1 ) {
			  try = fdh3d(              t0(X1,X2,z1+1),
				      t0(X1+1,X2,z1+1),t0(X1+1,X2-1,z1+1),t0(X1,X2-1,z1+1),
				      t0(X1+1,X2,z1  ),t0(X1+1,X2-1,z1  ),t0(X1,X2-1,z1  ),
				      s0(X1,X2,z1), s0(X1,X2,z1+1),
				      s0(X1+1,X2,z1+1),s0(X1+1,X2-1,z1+1),s0(X1,X2-1,z1+1),
				      s0(X1+1,X2,z1  ),s0(X1+1,X2-1,z1  ),s0(X1,X2-1,z1  ));
			  if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index-nx-1] < 1.e9
			   && time0[index-nx] < 1.e9 && X2>0  && X1>0 ) {
			  try = fdh3d(              t0(X1,X2,z1+1),
				      t0(X1-1,X2,z1+1),t0(X1-1,X2-1,z1+1),t0(X1,X2-1,z1+1),
				      t0(X1-1,X2,z1  ),t0(X1-1,X2-1,z1  ),t0(X1,X2-1,z1  ),
				      s0(X1,X2,z1), s0(X1,X2,z1+1),
				      s0(X1-1,X2,z1+1),s0(X1-1,X2-1,z1+1),s0(X1,X2-1,z1+1),
				      s0(X1-1,X2,z1  ),s0(X1-1,X2-1,z1  ),s0(X1,X2-1,z1  ));
			  if (try<guess) guess = try;
			}
			if(guess > 1.0e9){ 
			  if(time0[index+1] < 1.e9 && X1<nx-1 && X2>y1+1 && X2<y2-1 )  {
			      try = fdhne(t0(X1,X2,z1+1),t0(X1+1,X2,z1+1),t0(X1+1,X2,z1),
					  t0(X1+1,X2-1,z1+1),t0(X1+1,X2+1,z1+1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1+1,X2,z1+1),s0(X1+1,X2,z1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-1] < 1.e9 && X1>0 && X2>y1+1 && X2<y2-1 )  {
			      try = fdhne(t0(X1,X2,z1+1),t0(X1-1,X2,z1+1),t0(X1-1,X2,z1),
					  t0(X1-1,X2-1,z1+1),t0(X1-1,X2+1,z1+1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1-1,X2,z1+1),s0(X1-1,X2,z1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nx] < 1.e9 && X2<ny-1 && X1>x1+1 && X1<x2-1 )  {
			      try = fdhne(t0(X1,X2,z1+1),t0(X1,X2+1,z1+1),t0(X1,X2+1,z1),
					  t0(X1-1,X2+1,z1+1),t0(X1+1,X2+1,z1+1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1,X2+1,z1+1),s0(X1,X2+1,z1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X2>0 && X1>x1+1 && X1<x2-1 )  {
			      try = fdhne(t0(X1,X2,z1+1),t0(X1,X2-1,z1+1),t0(X1,X2-1,z1),
					  t0(X1-1,X2-1,z1+1),t0(X1+1,X2-1,z1+1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1,X2-1,z1+1),s0(X1,X2-1,z1) );
			    if (try<guess)  guess = try;
			  }
		        } 
			  if(time0[index+1] < 1.e9 && X1<nx-1 )  {
			    try = fdh2d(t0(X1,X2,z1+1),t0(X1+1,X2,z1+1),t0(X1+1,X2,z1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1+1,X2,z1+1),s0(X1+1,X2,z1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-1] < 1.e9 && X1>0 )  {
			    try = fdh2d(t0(X1,X2,z1+1),t0(X1-1,X2,z1+1),t0(X1-1,X2,z1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1-1,X2,z1+1),s0(X1-1,X2,z1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nx] < 1.e9 && X2<ny-1 )  {
			    try = fdh2d(t0(X1,X2,z1+1),t0(X1,X2+1,z1+1),t0(X1,X2+1,z1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1,X2+1,z1+1),s0(X1,X2+1,z1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X2>0 )  {
			    try = fdh2d(t0(X1,X2,z1+1),t0(X1,X2-1,z1+1),t0(X1,X2-1,z1),
					  s0(X1,X2,z1),
					  s0(X1,X2,z1+1),s0(X1,X2-1,z1+1),s0(X1,X2-1,z1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+1] < 1.e9 && time0[index+nx+1] < 1.e9
			     && time0[index+nx] < 1.e9 && X2<ny-1  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,X2,z1),t0(X1+1,X2+1,z1),t0(X1,X2+1,z1),
					s0(X1,X2,z1),
					s0(X1+1,X2,z1),s0(X1+1,X2+1,z1),s0(X1,X2+1,z1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index+1] < 1.e9 && time0[index-nx+1] < 1.e9
			     && time0[index-nx] < 1.e9 && X2>0  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,X2,z1),t0(X1+1,X2-1,z1),t0(X1,X2-1,z1),
					s0(X1,X2,z1),
					s0(X1+1,X2,z1),s0(X1+1,X2-1,z1),s0(X1,X2-1,z1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-1] < 1.e9 && time0[index+nx-1] < 1.e9
			     && time0[index+nx] < 1.e9 && X2<ny-1  && X1>0 ) {
			    try = fdh2d(t0(X1-1,X2,z1),t0(X1-1,X2+1,z1),t0(X1,X2+1,z1),
					s0(X1,X2,z1),
					s0(X1-1,X2,z1),s0(X1-1,X2+1,z1),s0(X1,X2+1,z1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-1] < 1.e9 && time0[index-nx-1] < 1.e9
			     && time0[index-nx] < 1.e9 && X2>0  && X1>0 ) {
			    try = fdh2d(t0(X1-1,X2,z1),t0(X1-1,X2-1,z1),t0(X1,X2-1,z1),
					s0(X1,X2,z1),
					s0(X1-1,X2,z1),s0(X1-1,X2-1,z1),s0(X1,X2-1,z1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if(guess > 1.0e9){ 
			  	if ( X1>x1+1 && X1<x2-1 && X2>y1+1 && X2<y2-1 ) {
			    	try = fdhnf(t0(X1,X2,z1+1),
						t0(X1+1,X2,z1+1),t0(X1,X2+1,z1+1),
					  	t0(X1-1,X2,z1+1),t0(X1,X2-1,z1+1),
					  	s0(X1,X2,z1),
					  	s0(X1,X2,z1+1) );
					if (try<guess)  guess = try;
			  	}
			}
            try = t0(X1,X2,z1+1) + .5*(s0(X1,X2,z1)+s0(X1,X2,z1+1));
			if (try<guess)  guess = try;
            if ( time0[index+1]<1.e9 && X1<nx-1 )  {
				try = t0(X1+1,X2,z1) + .5*(s0(X1,X2,z1)+s0(X1+1,X2,z1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			}
			if ( time0[index-1]<1.e9 && X1>0 )  {
			    try = t0(X1-1,X2,z1) + .5*(s0(X1,X2,z1)+s0(X1-1,X2,z1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			}
			if ( time0[index+nx]<1.e9 && X2<ny-1 )  {
				try = t0(X1,X2+1,z1) + .5*(s0(X1,X2,z1)+s0(X1,X2+1,z1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			}
			if ( time0[index-nx]<1.e9 && X2>0 )  {
			    try = t0(X1,X2-1,z1) + .5*(s0(X1,X2,z1)+s0(X1,X2-1,z1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			}
			if (guess<time0[index])  {
				time0[index] = guess;
				if (fhead>headtest)  headw[5]++;
			}
		}
		if(z1 == 1) dz1 = 0;
		z1--;
	}
      }
		/* BOTTOM SIDE */
      for (igrow=1;igrow<=kplus;igrow++) {  
	if(dz2){
		ii = 0;
		for(j=y1+1; j<=y2-1; j++){
			for(i=x1+1; i<=x2-1; i++){
				sort[ii].time = t0(i,j,z2-1);
				sort[ii].i1 = i;
				sort[ii].i2 = j;
				ii++;
			}
		}
		qsort((char *)sort,ii,sizeof(struct sorted),compar);
		for(i=0;i<ii;i++){
			X1 = sort[i].i1;
			X2 = sort[i].i2;
			index = z2*nxy + X2*nx + X1;
			lasti = (z2-1)*nxy + X2*nx + X1;
			fhead = 0.;
			guess = time0[index];
                        if(time0[index+1] < 1.e9 && time0[index+nx+1] < 1.e9
			   && time0[index+nx] < 1.e9 && X2<ny-1  && X1<nx-1 ) {
			  try = fdh3d(              t0(X1,X2,z2-1),
				      t0(X1+1,X2,z2-1),t0(X1+1,X2+1,z2-1),t0(X1,X2+1,z2-1),
				      t0(X1+1,X2,z2  ),t0(X1+1,X2+1,z2  ),t0(X1,X2+1,z2  ),
				      s0(X1,X2,z2), s0(X1,X2,z2-1),
				      s0(X1+1,X2,z2-1),s0(X1+1,X2+1,z2-1),s0(X1,X2+1,z2-1),
				      s0(X1+1,X2,z2  ),s0(X1+1,X2+1,z2  ),s0(X1,X2+1,z2  ));
			  if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index+nx-1] < 1.e9
			   && time0[index+nx] < 1.e9 && X2<ny-1  && X1>0 ) {
			  try = fdh3d(              t0(X1,X2,z2-1),
				      t0(X1-1,X2,z2-1),t0(X1-1,X2+1,z2-1),t0(X1,X2+1,z2-1),
				      t0(X1-1,X2,z2  ),t0(X1-1,X2+1,z2  ),t0(X1,X2+1,z2  ),
				      s0(X1,X2,z2), s0(X1,X2,z2-1),
				      s0(X1-1,X2,z2-1),s0(X1-1,X2+1,z2-1),s0(X1,X2+1,z2-1),
				      s0(X1-1,X2,z2  ),s0(X1-1,X2+1,z2  ),s0(X1,X2+1,z2  ));
			  if (try<guess) guess = try;
			}
			if(time0[index+1] < 1.e9 && time0[index-nx+1] < 1.e9
			   && time0[index-nx] < 1.e9 && X2>0  && X1<nx-1 ) {
			  try = fdh3d(              t0(X1,X2,z2-1),
				      t0(X1+1,X2,z2-1),t0(X1+1,X2-1,z2-1),t0(X1,X2-1,z2-1),
				      t0(X1+1,X2,z2  ),t0(X1+1,X2-1,z2  ),t0(X1,X2-1,z2  ),
				      s0(X1,X2,z2), s0(X1,X2,z2-1),
				      s0(X1+1,X2,z2-1),s0(X1+1,X2-1,z2-1),s0(X1,X2-1,z2-1),
				      s0(X1+1,X2,z2  ),s0(X1+1,X2-1,z2  ),s0(X1,X2-1,z2  ));
			  if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index-nx-1] < 1.e9
			   && time0[index-nx] < 1.e9 && X2>0  && X1>0 ) {
			  try = fdh3d(              t0(X1,X2,z2-1),
				      t0(X1-1,X2,z2-1),t0(X1-1,X2-1,z2-1),t0(X1,X2-1,z2-1),
				      t0(X1-1,X2,z2  ),t0(X1-1,X2-1,z2  ),t0(X1,X2-1,z2  ),
				      s0(X1,X2,z2), s0(X1,X2,z2-1),
				      s0(X1-1,X2,z2-1),s0(X1-1,X2-1,z2-1),s0(X1,X2-1,z2-1),
				      s0(X1-1,X2,z2  ),s0(X1-1,X2-1,z2  ),s0(X1,X2-1,z2  ));
			  if (try<guess) guess = try;
			}
                        if(guess > 1.0e9){ 
			  if(time0[index+1] < 1.e9 && X1<nx-1 && X2>y1+1 && X2<y2-1 )  {
			      try = fdhne(t0(X1,X2,z2-1),t0(X1+1,X2,z2-1),t0(X1+1,X2,z2),
					  t0(X1+1,X2-1,z2-1),t0(X1+1,X2+1,z2-1),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1+1,X2,z2-1),s0(X1+1,X2,z2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-1] < 1.e9 && X1>0 && X2>y1+1 && X2<y2-1 )  {
			      try = fdhne(t0(X1,X2,z2-1),t0(X1-1,X2,z2-1),t0(X1-1,X2,z2),
					  t0(X1-1,X2-1,z2-1),t0(X1-1,X2+1,z2-1),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1-1,X2,z2-1),s0(X1-1,X2,z2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nx] < 1.e9 && X2<ny-1 && X1>x1+1 && X1<x2-1 )  {
			      try = fdhne(t0(X1,X2,z2-1),t0(X1,X2+1,z2-1),t0(X1,X2+1,z2),
					  t0(X1-1,X2+1,z2-1),t0(X1+1,X2+1,z2-1),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1,X2+1,z2-1),s0(X1,X2+1,z2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X2>0 && X1>x1+1 && X1<x2-1 )  {
			      try = fdhne(t0(X1,X2,z2-1),t0(X1,X2-1,z2-1),t0(X1,X2-1,z2),
					  t0(X1-1,X2-1,z2-1),t0(X1+1,X2-1,z2-1),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1,X2-1,z2-1),s0(X1,X2-1,z2) );
			    if (try<guess)  guess = try;
			  }
		        }
			  if(time0[index+1] < 1.e9 && X1<nx-1 )  {
			    try = fdh2d(t0(X1,X2,z2-1),t0(X1+1,X2,z2-1),t0(X1+1,X2,z2),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1+1,X2,z2-1),s0(X1+1,X2,z2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-1] < 1.e9 && X1>0 )  {
			    try = fdh2d(t0(X1,X2,z2-1),t0(X1-1,X2,z2-1),t0(X1-1,X2,z2),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1-1,X2,z2-1),s0(X1-1,X2,z2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nx] < 1.e9 && X2<ny-1 )  {
			    try = fdh2d(t0(X1,X2,z2-1),t0(X1,X2+1,z2-1),t0(X1,X2+1,z2),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1,X2+1,z2-1),s0(X1,X2+1,z2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X2>0 )  {
			    try = fdh2d(t0(X1,X2,z2-1),t0(X1,X2-1,z2-1),t0(X1,X2-1,z2),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1),s0(X1,X2-1,z2-1),s0(X1,X2-1,z2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+1] < 1.e9 && time0[index+nx+1] < 1.e9
			     && time0[index+nx] < 1.e9 && X2<ny-1  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,X2,z2),t0(X1+1,X2+1,z2),t0(X1,X2+1,z2),
					s0(X1,X2,z2),
					s0(X1+1,X2,z2),s0(X1+1,X2+1,z2),s0(X1,X2+1,z2) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index+1] < 1.e9 && time0[index-nx+1] < 1.e9
			     && time0[index-nx] < 1.e9 && X2>0  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,X2,z2),t0(X1+1,X2-1,z2),t0(X1,X2-1,z2),
					s0(X1,X2,z2),
					s0(X1+1,X2,z2),s0(X1+1,X2-1,z2),s0(X1,X2-1,z2) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-1] < 1.e9 && time0[index+nx-1] < 1.e9
			     && time0[index+nx] < 1.e9 && X2<ny-1  && X1>0 ) {
			    try = fdh2d(t0(X1-1,X2,z2),t0(X1-1,X2+1,z2),t0(X1,X2+1,z2),
					s0(X1,X2,z2),
					s0(X1-1,X2,z2),s0(X1-1,X2+1,z2),s0(X1,X2+1,z2) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-1] < 1.e9 && time0[index-nx-1] < 1.e9
			     && time0[index-nx] < 1.e9 && X2>0  && X1>0 ) {
			    try = fdh2d(t0(X1-1,X2,z2),t0(X1-1,X2-1,z2),t0(X1,X2-1,z2),
					s0(X1,X2,z2),
					s0(X1-1,X2,z2),s0(X1-1,X2-1,z2),s0(X1,X2-1,z2) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if(guess > 1.0e9){ 
			  if ( X1>x1+1 && X1<x2-1 && X2>y1+1 && X2<y2-1 ) {
			    try = fdhnf(t0(X1,X2,z2-1),
					  t0(X1+1,X2,z2-1),t0(X1,X2+1,z2-1),
					  t0(X1-1,X2,z2-1),t0(X1,X2-1,z2-1),
					  s0(X1,X2,z2),
					  s0(X1,X2,z2-1) );
			    if (try<guess)  guess = try;
			  }
			} 
			  try = t0(X1,X2,z2-1) + .5*(s0(X1,X2,z2)+s0(X1,X2,z2-1));
			  if (try<guess)  guess = try;
                          if ( time0[index+1]<1.e9 && X1<nx-1 )  {
			    try = t0(X1+1,X2,z2) + .5*(s0(X1,X2,z2)+s0(X1+1,X2,z2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-1]<1.e9 && X1>0 )  {
			    try = t0(X1-1,X2,z2) + .5*(s0(X1,X2,z2)+s0(X1-1,X2,z2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index+nx]<1.e9 && X2<ny-1 )  {
			    try = t0(X1,X2+1,z2) + .5*(s0(X1,X2,z2)+s0(X1,X2+1,z2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-nx]<1.e9 && X2>0 )  {
			    try = t0(X1,X2-1,z2) + .5*(s0(X1,X2,z2)+s0(X1,X2-1,z2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if (guess<time0[index]) {
				time0[index] = guess;
				if (fhead>headtest)  headw[6]++;
			}
		}
		if(z2 == nz-2) dz2 = 0;
		z2++;
	}
      }
		/* FRONT SIDE */
      for (igrow=1;igrow<=jminus;igrow++) {  
	if(dy1){
		ii = 0;
		for(k=z1+1; k<=z2-1; k++){
			for(i=x1+1; i<=x2-1; i++){
				sort[ii].time = t0(i,y1+1,k);
				sort[ii].i1 = i;
				sort[ii].i2 = k;
				ii++;
			}
		}
		qsort((char *)sort,ii,sizeof(struct sorted),compar);
		for(i=0;i<ii;i++){
			X1 = sort[i].i1;
			X2 = sort[i].i2;
			index = X2*nxy + y1*nx + X1;
			lasti = X2*nxy + (y1+1)*nx + X1;
			fhead = 0.;
			guess = time0[index];
			if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) {
			  	try = fdh3d(              t0(X1,y1+1,X2),
				      t0(X1+1,y1+1,X2),t0(X1+1,y1+1,X2+1),t0(X1,y1+1,X2+1),
				      t0(X1+1,y1  ,X2),t0(X1+1,y1  ,X2+1),t0(X1,y1  ,X2+1),
				      s0(X1,y1,X2), s0(X1,y1+1,X2),
				      s0(X1+1,y1+1,X2),s0(X1+1,y1+1,X2+1),s0(X1,y1+1,X2+1),
				      s0(X1+1,y1  ,X2),s0(X1+1,y1  ,X2+1),s0(X1,y1  ,X2+1));
			  	if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			  	try = fdh3d(              t0(X1,y1+1,X2),
				      t0(X1-1,y1+1,X2),t0(X1-1,y1+1,X2+1),t0(X1,y1+1,X2+1),
				      t0(X1-1,y1  ,X2),t0(X1-1,y1  ,X2+1),t0(X1,y1  ,X2+1),
				      s0(X1,y1,X2), s0(X1,y1+1,X2),
				      s0(X1-1,y1+1,X2),s0(X1-1,y1+1,X2+1),s0(X1,y1+1,X2+1),
				      s0(X1-1,y1  ,X2),s0(X1-1,y1  ,X2+1),s0(X1,y1  ,X2+1));
			  	if (try<guess) guess = try;
			}
			if(time0[index+1] < 1.e9 && time0[index-nxy+1] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1<nx-1 ) {
			  	try = fdh3d(              t0(X1,y1+1,X2),
				      t0(X1+1,y1+1,X2),t0(X1+1,y1+1,X2-1),t0(X1,y1+1,X2-1),
				      t0(X1+1,y1  ,X2),t0(X1+1,y1  ,X2-1),t0(X1,y1  ,X2-1),
				      s0(X1,y1,X2), s0(X1,y1+1,X2),
				      s0(X1+1,y1+1,X2),s0(X1+1,y1+1,X2-1),s0(X1,y1+1,X2-1),
				      s0(X1+1,y1  ,X2),s0(X1+1,y1  ,X2-1),s0(X1,y1  ,X2-1));
			  	if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index-nxy-1] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			  	try = fdh3d(              t0(X1,y1+1,X2),
				      t0(X1-1,y1+1,X2),t0(X1-1,y1+1,X2-1),t0(X1,y1+1,X2-1),
				      t0(X1-1,y1  ,X2),t0(X1-1,y1  ,X2-1),t0(X1,y1  ,X2-1),
				      s0(X1,y1,X2), s0(X1,y1+1,X2),
				      s0(X1-1,y1+1,X2),s0(X1-1,y1+1,X2-1),s0(X1,y1+1,X2-1),
				      s0(X1-1,y1  ,X2),s0(X1-1,y1  ,X2-1),s0(X1,y1  ,X2-1));
			  	if (try<guess) guess = try;
			}
			if(guess > 1.0e9){ 
			  	if(time0[index+1] < 1.e9 && X1<nx-1 && X2>z1+1 && X2<z2-1 )  {
			      	try = fdhne(t0(X1,y1+1,X2),t0(X1+1,y1+1,X2),t0(X1+1,y1,X2),
					  t0(X1+1,y1+1,X2-1),t0(X1+1,y1+1,X2+1),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1+1,y1+1,X2),s0(X1+1,y1,X2) );
			    	if (try<guess)  guess = try;
				}
				if(time0[index-1] < 1.e9 && X1>0 && X2>z1+1 && X2<z2-1 )  {
			    	try = fdhne(t0(X1,y1+1,X2),t0(X1-1,y1+1,X2),t0(X1-1,y1,X2),
					  t0(X1-1,y1+1,X2-1),t0(X1-1,y1+1,X2+1),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1-1,y1+1,X2),s0(X1-1,y1,X2) );
			    	if (try<guess)  guess = try;
				}
				if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>x1+1 && X1<x2-1 )  {
			    	try = fdhne(t0(X1,y1+1,X2),t0(X1,y1+1,X2+1),t0(X1,y1,X2+1),
					  t0(X1-1,y1+1,X2+1),t0(X1+1,y1+1,X2+1),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1,y1+1,X2+1),s0(X1,y1,X2+1) );
			    	if (try<guess)  guess = try;
				}
				if(time0[index-nxy] < 1.e9 && X2>0 && X1>x1+1 && X1<x2-1 )  {
			    	try = fdhne(t0(X1,y1+1,X2),t0(X1,y1+1,X2-1),t0(X1,y1,X2-1),
					  t0(X1-1,y1+1,X2-1),t0(X1+1,y1+1,X2-1),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1,y1+1,X2-1),s0(X1,y1,X2-1) );
			    	if (try<guess)  guess = try;
				}
		    } 
			if(time0[index+1] < 1.e9 && X1<nx-1 )  {
			    try = fdh2d(t0(X1,y1+1,X2),t0(X1+1,y1+1,X2),t0(X1+1,y1,X2),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1+1,y1+1,X2),s0(X1+1,y1,X2) );
			    if (try<guess)  guess = try;
			}
			if(time0[index-1] < 1.e9 && X1>0 )  {
			    try = fdh2d(t0(X1,y1+1,X2),t0(X1-1,y1+1,X2),t0(X1-1,y1,X2),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1-1,y1+1,X2),s0(X1-1,y1,X2) );
			    if (try<guess)  guess = try;
			}
			if(time0[index+nxy] < 1.e9 && X2<nz-1 )  {
			    try = fdh2d(t0(X1,y1+1,X2),t0(X1,y1+1,X2+1),t0(X1,y1,X2+1),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1,y1+1,X2+1),s0(X1,y1,X2+1) );
			    if (try<guess)  guess = try;
			}
			if(time0[index-nxy] < 1.e9 && X2>0 )  {
			    try = fdh2d(t0(X1,y1+1,X2),t0(X1,y1+1,X2-1),t0(X1,y1,X2-1),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2),s0(X1,y1+1,X2-1),s0(X1,y1,X2-1) );
			    if (try<guess)  guess = try;
			}
			if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,y1,X2),t0(X1+1,y1,X2+1),t0(X1,y1,X2+1),
					s0(X1,y1,X2),
					s0(X1+1,y1,X2),s0(X1+1,y1,X2+1),s0(X1,y1,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			}
			if(time0[index+1] < 1.e9 && time0[index-nxy+1] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,y1,X2),t0(X1+1,y1,X2-1),t0(X1,y1,X2-1),
					s0(X1,y1,X2),
					s0(X1+1,y1,X2),s0(X1+1,y1,X2-1),s0(X1,y1,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			}
			if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			    try = fdh2d(t0(X1-1,y1,X2),t0(X1-1,y1,X2+1),t0(X1,y1,X2+1),
					s0(X1,y1,X2),
					s0(X1-1,y1,X2),s0(X1-1,y1,X2+1),s0(X1,y1,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			}
			if(time0[index-1] < 1.e9 && time0[index-nxy-1] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			    try = fdh2d(t0(X1-1,y1,X2),t0(X1-1,y1,X2-1),t0(X1,y1,X2-1),
					s0(X1,y1,X2),
					s0(X1-1,y1,X2),s0(X1-1,y1,X2-1),s0(X1,y1,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			}
			if(guess > 1.0e9){ 
			  	if ( X1>x1+1 && X1<x2-1 && X2>z1+1 && X2<z2-1 ) {
			    	try = fdhnf(t0(X1,y1+1,X2),
					  t0(X1+1,y1+1,X2),t0(X1,y1+1,X2+1),
					  t0(X1-1,y1+1,X2),t0(X1,y1+1,X2-1),
					  s0(X1,y1,X2),
					  s0(X1,y1+1,X2) );
			    	if (try<guess)  guess = try;
			  	}
			} 
			try = t0(X1,y1+1,X2) + .5*(s0(X1,y1,X2)+s0(X1,y1+1,X2));
			if (try<guess)  guess = try;
            if ( time0[index+1]<1.e9 && X1<nx-1 )  {
			    try = t0(X1+1,y1,X2) + .5*(s0(X1,y1,X2)+s0(X1+1,y1,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			}
			if ( time0[index-1]<1.e9 && X1>0 )  {
			    try = t0(X1-1,y1,X2) + .5*(s0(X1,y1,X2)+s0(X1-1,y1,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			}
			if ( time0[index+nxy]<1.e9 && X2<nz-1 )  {
			    try = t0(X1,y1,X2+1) + .5*(s0(X1,y1,X2)+s0(X1,y1,X2+1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			}
			if ( time0[index-nxy]<1.e9 && X2>0 )  {
			    try = t0(X1,y1,X2-1) + .5*(s0(X1,y1,X2)+s0(X1,y1,X2-1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			}
			if (guess<time0[index]) {
				time0[index] = guess;
				if (fhead>headtest)  headw[3]++;
			}
		}
		if(y1 == 1) dy1 = 0;
		y1--;
	}
      }
		/* BACK SIDE */
      for (igrow=1;igrow<=jplus;igrow++) {  
	if(dy2){
		ii = 0;
		for(k=z1+1; k<=z2-1; k++){
			for(i=x1+1; i<=x2-1; i++){
				sort[ii].time = t0(i,y2-1,k);
				sort[ii].i1 = i;
				sort[ii].i2 = k;
				ii++;
			}
		}
		qsort((char *)sort,ii,sizeof(struct sorted),compar);
		for(i=0;i<ii;i++){
			X1 = sort[i].i1;
			X2 = sort[i].i2;
			index = X2*nxy + y2*nx + X1;
			lasti = X2*nxy + (y2-1)*nx + X1;
			fhead = 0.;
			guess = time0[index];
			if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) {
			  try = fdh3d(              t0(X1,y2-1,X2),
				      t0(X1+1,y2-1,X2),t0(X1+1,y2-1,X2+1),t0(X1,y2-1,X2+1),
				      t0(X1+1,y2  ,X2),t0(X1+1,y2  ,X2+1),t0(X1,y2  ,X2+1),
				      s0(X1,y2,X2), s0(X1,y2-1,X2),
				      s0(X1+1,y2-1,X2),s0(X1+1,y2-1,X2+1),s0(X1,y2-1,X2+1),
				      s0(X1+1,y2  ,X2),s0(X1+1,y2  ,X2+1),s0(X1,y2  ,X2+1));
			  if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			  try = fdh3d(              t0(X1,y2-1,X2),
				      t0(X1-1,y2-1,X2),t0(X1-1,y2-1,X2+1),t0(X1,y2-1,X2+1),
				      t0(X1-1,y2  ,X2),t0(X1-1,y2  ,X2+1),t0(X1,y2  ,X2+1),
				      s0(X1,y2,X2), s0(X1,y2-1,X2),
				      s0(X1-1,y2-1,X2),s0(X1-1,y2-1,X2+1),s0(X1,y2-1,X2+1),
				      s0(X1-1,y2  ,X2),s0(X1-1,y2  ,X2+1),s0(X1,y2  ,X2+1));
			  if (try<guess) guess = try;
			}
			if(time0[index+1] < 1.e9 && time0[index-nxy+1] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1<nx-1 ) {
			  try = fdh3d(              t0(X1,y2-1,X2),
				      t0(X1+1,y2-1,X2),t0(X1+1,y2-1,X2-1),t0(X1,y2-1,X2-1),
				      t0(X1+1,y2  ,X2),t0(X1+1,y2  ,X2-1),t0(X1,y2  ,X2-1),
				      s0(X1,y2,X2), s0(X1,y2-1,X2),
				      s0(X1+1,y2-1,X2),s0(X1+1,y2-1,X2-1),s0(X1,y2-1,X2-1),
				      s0(X1+1,y2  ,X2),s0(X1+1,y2  ,X2-1),s0(X1,y2  ,X2-1));
			  if (try<guess) guess = try;
			}
			if(time0[index-1] < 1.e9 && time0[index-nxy-1] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			  try = fdh3d(              t0(X1,y2-1,X2),
				      t0(X1-1,y2-1,X2),t0(X1-1,y2-1,X2-1),t0(X1,y2-1,X2-1),
				      t0(X1-1,y2  ,X2),t0(X1-1,y2  ,X2-1),t0(X1,y2  ,X2-1),
				      s0(X1,y2,X2), s0(X1,y2-1,X2),
				      s0(X1-1,y2-1,X2),s0(X1-1,y2-1,X2-1),s0(X1,y2-1,X2-1),
				      s0(X1-1,y2  ,X2),s0(X1-1,y2  ,X2-1),s0(X1,y2  ,X2-1));
			  if (try<guess) guess = try;
			}
			if(guess > 1.0e9){ 
			  if(time0[index+1] < 1.e9 && X1<nx-1 && X2>z1+1 && X2<z2-1 )  {
			      try = fdhne(t0(X1,y2-1,X2),t0(X1+1,y2-1,X2),t0(X1+1,y2,X2),
					  t0(X1+1,y2-1,X2-1),t0(X1+1,y2-1,X2+1),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1+1,y2-1,X2),s0(X1+1,y2,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-1] < 1.e9 && X1>0 && X2>z1+1 && X2<z2-1 )  {
			      try = fdhne(t0(X1,y2-1,X2),t0(X1-1,y2-1,X2),t0(X1-1,y2,X2),
					  t0(X1-1,y2-1,X2-1),t0(X1-1,y2-1,X2+1),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1-1,y2-1,X2),s0(X1-1,y2,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>x1+1 && X1<x2-1 )  {
			      try = fdhne(t0(X1,y2-1,X2),t0(X1,y2-1,X2+1),t0(X1,y2,X2+1),
					  t0(X1-1,y2-1,X2+1),t0(X1+1,y2-1,X2+1),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1,y2-1,X2+1),s0(X1,y2,X2+1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nxy] < 1.e9 && X2>0 && X1>x1+1 && X1<x2-1 )  {
			      try = fdhne(t0(X1,y2-1,X2),t0(X1,y2-1,X2-1),t0(X1,y2,X2-1),
					  t0(X1-1,y2-1,X2-1),t0(X1+1,y2-1,X2-1),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1,y2-1,X2-1),s0(X1,y2,X2-1) );
			    if (try<guess)  guess = try;
			  }
		        } 
			  if(time0[index+1] < 1.e9 && X1<nx-1 )  {
			    try = fdh2d(t0(X1,y2-1,X2),t0(X1+1,y2-1,X2),t0(X1+1,y2,X2),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1+1,y2-1,X2),s0(X1+1,y2,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-1] < 1.e9 && X1>0 )  {
			    try = fdh2d(t0(X1,y2-1,X2),t0(X1-1,y2-1,X2),t0(X1-1,y2,X2),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1-1,y2-1,X2),s0(X1-1,y2,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 )  {
			    try = fdh2d(t0(X1,y2-1,X2),t0(X1,y2-1,X2+1),t0(X1,y2,X2+1),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1,y2-1,X2+1),s0(X1,y2,X2+1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nxy] < 1.e9 && X2>0 )  {
			    try = fdh2d(t0(X1,y2-1,X2),t0(X1,y2-1,X2-1),t0(X1,y2,X2-1),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2),s0(X1,y2-1,X2-1),s0(X1,y2,X2-1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+1] < 1.e9 && time0[index+nxy+1] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,y2,X2),t0(X1+1,y2,X2+1),t0(X1,y2,X2+1),
					s0(X1,y2,X2),
					s0(X1+1,y2,X2),s0(X1+1,y2,X2+1),s0(X1,y2,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index+1] < 1.e9 && time0[index-nxy+1] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1<nx-1 ) {
			    try = fdh2d(t0(X1+1,y2,X2),t0(X1+1,y2,X2-1),t0(X1,y2,X2-1),
					s0(X1,y2,X2),
					s0(X1+1,y2,X2),s0(X1+1,y2,X2-1),s0(X1,y2,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-1] < 1.e9 && time0[index+nxy-1] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			    try = fdh2d(t0(X1-1,y2,X2),t0(X1-1,y2,X2+1),t0(X1,y2,X2+1),
					s0(X1,y2,X2),
					s0(X1-1,y2,X2),s0(X1-1,y2,X2+1),s0(X1,y2,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-1] < 1.e9 && time0[index-nxy-1] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			    try = fdh2d(t0(X1-1,y2,X2),t0(X1-1,y2,X2-1),t0(X1,y2,X2-1),
					s0(X1,y2,X2),
					s0(X1-1,y2,X2),s0(X1-1,y2,X2-1),s0(X1,y2,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if(guess > 1.0e9){ 
			  if ( X1>x1+1 && X1<x2-1 && X2>z1+1 && X2<z2-1 ) {
			    try = fdhnf(t0(X1,y2-1,X2),
					  t0(X1+1,y2-1,X2),t0(X1,y2-1,X2+1),
					  t0(X1-1,y2-1,X2),t0(X1,y2-1,X2-1),
					  s0(X1,y2,X2),
					  s0(X1,y2-1,X2) );
			    if (try<guess)  guess = try;
			  }
			} 
			  try = t0(X1,y2-1,X2) + .5*(s0(X1,y2,X2)+s0(X1,y2-1,X2));
			  if (try<guess)  guess = try;
                          if ( time0[index+1]<1.e9 && X1<nx-1 )  {
			    try = t0(X1+1,y2,X2) + .5*(s0(X1,y2,X2)+s0(X1+1,y2,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-1]<1.e9 && X1>0 )  {
			    try = t0(X1-1,y2,X2) + .5*(s0(X1,y2,X2)+s0(X1-1,y2,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index+nxy]<1.e9 && X2<nz-1 )  {
			    try = t0(X1,y2,X2+1) + .5*(s0(X1,y2,X2)+s0(X1,y2,X2+1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-nxy]<1.e9 && X2>0 )  {
			    try = t0(X1,y2,X2-1) + .5*(s0(X1,y2,X2)+s0(X1,y2,X2-1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if (guess<time0[index]) {
				time0[index] = guess;
				if (fhead>headtest)  headw[4]++;
			}
		}
		if(y2 == ny-2) dy2 = 0;
		y2++;
	}
      }
		/* LEFT SIDE */
      for (igrow=1;igrow<=iminus;igrow++) {  
	if(dx1){
		ii = 0;
		for(k=z1+1; k<=z2-1; k++){
			for(j=y1+1; j<=y2-1; j++){
				sort[ii].time = t0(x1+1,j,k);
				sort[ii].i1 = j;
				sort[ii].i2 = k;
				ii++;
			}
		}
		qsort((char *)sort,ii,sizeof(struct sorted),compar);
		for(i=0;i<ii;i++){
			X1 = sort[i].i1;
			X2 = sort[i].i2;
			index = X2*nxy + X1*nx + x1;
			lasti = X2*nxy + X1*nx + (x1+1);
			fhead = 0.;
			guess = time0[index];
			if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) {
			  try = fdh3d(              t0(x1+1,X1,X2),
				      t0(x1+1,X1+1,X2),t0(x1+1,X1+1,X2+1),t0(x1+1,X1,X2+1),
				      t0(x1  ,X1+1,X2),t0(x1  ,X1+1,X2+1),t0(x1  ,X1,X2+1),
				      s0(x1,X1,X2), s0(x1+1,X1,X2),
				      s0(x1+1,X1+1,X2),s0(x1+1,X1+1,X2+1),s0(x1+1,X1,X2+1),
				      s0(x1  ,X1+1,X2),s0(x1  ,X1+1,X2+1),s0(x1  ,X1,X2+1));
			  if (try<guess) guess = try;
			}
			if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			  try = fdh3d(              t0(x1+1,X1,X2),
				      t0(x1+1,X1-1,X2),t0(x1+1,X1-1,X2+1),t0(x1+1,X1,X2+1),
				      t0(x1  ,X1-1,X2),t0(x1  ,X1-1,X2+1),t0(x1  ,X1,X2+1),
				      s0(x1,X1,X2), s0(x1+1,X1,X2),
				      s0(x1+1,X1-1,X2),s0(x1+1,X1-1,X2+1),s0(x1+1,X1,X2+1),
				      s0(x1  ,X1-1,X2),s0(x1  ,X1-1,X2+1),s0(x1  ,X1,X2+1));
			  if (try<guess) guess = try;
			}
			if(time0[index+nx] < 1.e9 && time0[index-nxy+nx] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1<ny-1 ) {
			  try = fdh3d(              t0(x1+1,X1,X2),
				      t0(x1+1,X1+1,X2),t0(x1+1,X1+1,X2-1),t0(x1+1,X1,X2-1),
				      t0(x1  ,X1+1,X2),t0(x1  ,X1+1,X2-1),t0(x1  ,X1,X2-1),
				      s0(x1,X1,X2), s0(x1+1,X1,X2),
				      s0(x1+1,X1+1,X2),s0(x1+1,X1+1,X2-1),s0(x1+1,X1,X2-1),
				      s0(x1  ,X1+1,X2),s0(x1  ,X1+1,X2-1),s0(x1  ,X1,X2-1));
			  if (try<guess) guess = try;
			}
			if(time0[index-nx] < 1.e9 && time0[index-nxy-nx] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			  try = fdh3d(              t0(x1+1,X1,X2),
				      t0(x1+1,X1-1,X2),t0(x1+1,X1-1,X2-1),t0(x1+1,X1,X2-1),
				      t0(x1  ,X1-1,X2),t0(x1  ,X1-1,X2-1),t0(x1  ,X1,X2-1),
				      s0(x1,X1,X2), s0(x1+1,X1,X2),
				      s0(x1+1,X1-1,X2),s0(x1+1,X1-1,X2-1),s0(x1+1,X1,X2-1),
				      s0(x1  ,X1-1,X2),s0(x1  ,X1-1,X2-1),s0(x1  ,X1,X2-1));
			  if (try<guess) guess = try;
			}
			if(guess > 1.0e9){ 
			  if(time0[index+nx] < 1.e9 && X1<ny-1 && X2>z1+1 && X2<z2-1 )  {
			      try = fdhne(t0(x1+1,X1,X2),t0(x1+1,X1+1,X2),t0(x1,X1+1,X2),
					  t0(x1+1,X1+1,X2-1),t0(x1+1,X1+1,X2+1),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1+1,X2),s0(x1,X1+1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X1>0 && X2>z1+1 && X2<z2-1 )  {
			      try = fdhne(t0(x1+1,X1,X2),t0(x1+1,X1-1,X2),t0(x1,X1-1,X2),
					  t0(x1+1,X1-1,X2-1),t0(x1+1,X1-1,X2+1),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1-1,X2),s0(x1,X1-1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>y1+1 && X1<y2-1 )  {
			      try = fdhne(t0(x1+1,X1,X2),t0(x1+1,X1,X2+1),t0(x1,X1,X2+1),
					  t0(x1+1,X1-1,X2+1),t0(x1+1,X1+1,X2+1),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1,X2+1),s0(x1,X1,X2+1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nxy] < 1.e9 && X2>0 && X1>y1+1 && X1<y2-1 )  {
			      try = fdhne(t0(x1+1,X1,X2),t0(x1+1,X1,X2-1),t0(x1,X1,X2-1),
					  t0(x1+1,X1-1,X2-1),t0(x1+1,X1+1,X2-1),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1,X2-1),s0(x1,X1,X2-1) );
			    if (try<guess)  guess = try;
			  }
		        } 
			  if(time0[index+nx] < 1.e9 && X1<ny-1 )  {
			    try = fdh2d(t0(x1+1,X1,X2),t0(x1+1,X1+1,X2),t0(x1,X1+1,X2),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1+1,X2),s0(x1,X1+1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X1>0 )  {
			    try = fdh2d(t0(x1+1,X1,X2),t0(x1+1,X1-1,X2),t0(x1,X1-1,X2),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1-1,X2),s0(x1,X1-1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 )  {
			    try = fdh2d(t0(x1+1,X1,X2),t0(x1+1,X1,X2+1),t0(x1,X1,X2+1),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1,X2+1),s0(x1,X1,X2+1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nxy] < 1.e9 && X2>0 )  {
			    try = fdh2d(t0(x1+1,X1,X2),t0(x1+1,X1,X2-1),t0(x1,X1,X2-1),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2),s0(x1+1,X1,X2-1),s0(x1,X1,X2-1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) {
			    try = fdh2d(t0(x1,X1+1,X2),t0(x1,X1+1,X2+1),t0(x1,X1,X2+1),
					s0(x1,X1,X2),
					s0(x1,X1+1,X2),s0(x1,X1+1,X2+1),s0(x1,X1,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index+nx] < 1.e9 && time0[index-nxy+nx] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1<ny-1 ) {
			    try = fdh2d(t0(x1,X1+1,X2),t0(x1,X1+1,X2-1),t0(x1,X1,X2-1),
					s0(x1,X1,X2),
					s0(x1,X1+1,X2),s0(x1,X1+1,X2-1),s0(x1,X1,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			    try = fdh2d(t0(x1,X1-1,X2),t0(x1,X1-1,X2+1),t0(x1,X1,X2+1),
					s0(x1,X1,X2),
					s0(x1,X1-1,X2),s0(x1,X1-1,X2+1),s0(x1,X1,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-nx] < 1.e9 && time0[index-nxy-nx] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			    try = fdh2d(t0(x1,X1-1,X2),t0(x1,X1-1,X2-1),t0(x1,X1,X2-1),
					s0(x1,X1,X2),
					s0(x1,X1-1,X2),s0(x1,X1-1,X2-1),s0(x1,X1,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if(guess > 1.0e9){ 
			  if ( X1>y1+1 && X1<y2-1 && X2>z1+1 && X2<z2-1 ) {
			    try = fdhnf(t0(x1+1,X1,X2),
					  t0(x1+1,X1+1,X2),t0(x1+1,X1,X2+1),
					  t0(x1+1,X1-1,X2),t0(x1+1,X1,X2-1),
					  s0(x1,X1,X2),
					  s0(x1+1,X1,X2) );
			    if (try<guess)  guess = try;
			  }
			} 
			  try = t0(x1+1,X1,X2) + .5*(s0(x1,X1,X2)+s0(x1+1,X1,X2));
			  if (try<guess)  guess = try;
                          if ( time0[index+nx]<1.e9 && X1<ny-1 )  {
			    try = t0(x1,X1+1,X2) + .5*(s0(x1,X1,X2)+s0(x1,X1+1,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-nx]<1.e9 && X1>0 )  {
			    try = t0(x1,X1-1,X2) + .5*(s0(x1,X1,X2)+s0(x1,X1-1,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index+nxy]<1.e9 && X2<nz-1 )  {
			    try = t0(x1,X1,X2+1) + .5*(s0(x1,X1,X2)+s0(x1,X1,X2+1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-nxy]<1.e9 && X2>0 )  {
			    try = t0(x1,X1,X2-1) + .5*(s0(x1,X1,X2)+s0(x1,X1,X2-1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if (guess<time0[index]) {
				time0[index] = guess;
				if (fhead>headtest)  headw[1]++;
			}
		}
		if(x1 == 1) dx1 = 0;
		x1--;
	}
      }
		/* RIGHT SIDE */
      for (igrow=1;igrow<=iplus;igrow++) {  
	if(dx2){
		ii = 0;
		for(k=z1+1; k<=z2-1; k++){
			for(j=y1+1; j<=y2-1; j++){
				sort[ii].time = t0(x2-1,j,k);
				sort[ii].i1 = j;
				sort[ii].i2 = k;
				ii++;
			}
		}
		qsort((char *)sort,ii,sizeof(struct sorted),compar);
		for(i=0;i<ii;i++){
			X1 = sort[i].i1;
			X2 = sort[i].i2;
			index = X2*nxy + X1*nx + x2;
			lasti = X2*nxy + X1*nx + (x2-1);
			fhead = 0.;
			guess = time0[index];
			if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) {
			  try = fdh3d(              t0(x2-1,X1,X2),
				      t0(x2-1,X1+1,X2),t0(x2-1,X1+1,X2+1),t0(x2-1,X1,X2+1),
				      t0(x2  ,X1+1,X2),t0(x2  ,X1+1,X2+1),t0(x2  ,X1,X2+1),
				      s0(x2,X1,X2), s0(x2-1,X1,X2),
				      s0(x2-1,X1+1,X2),s0(x2-1,X1+1,X2+1),s0(x2-1,X1,X2+1),
				      s0(x2  ,X1+1,X2),s0(x2  ,X1+1,X2+1),s0(x2  ,X1,X2+1));
			  if (try<guess) guess = try;
			}
			if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9
			   && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			  try = fdh3d(              t0(x2-1,X1,X2),
				      t0(x2-1,X1-1,X2),t0(x2-1,X1-1,X2+1),t0(x2-1,X1,X2+1),
				      t0(x2  ,X1-1,X2),t0(x2  ,X1-1,X2+1),t0(x2  ,X1,X2+1),
				      s0(x2,X1,X2), s0(x2-1,X1,X2),
				      s0(x2-1,X1-1,X2),s0(x2-1,X1-1,X2+1),s0(x2-1,X1,X2+1),
				      s0(x2  ,X1-1,X2),s0(x2  ,X1-1,X2+1),s0(x2  ,X1,X2+1));
			  if (try<guess) guess = try;
			}
			if(time0[index+nx] < 1.e9 && time0[index-nxy+nx] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1<ny-1 ) {
			  try = fdh3d(              t0(x2-1,X1,X2),
				      t0(x2-1,X1+1,X2),t0(x2-1,X1+1,X2-1),t0(x2-1,X1,X2-1),
				      t0(x2  ,X1+1,X2),t0(x2  ,X1+1,X2-1),t0(x2  ,X1,X2-1),
				      s0(x2,X1,X2), s0(x2-1,X1,X2),
				      s0(x2-1,X1+1,X2),s0(x2-1,X1+1,X2-1),s0(x2-1,X1,X2-1),
				      s0(x2  ,X1+1,X2),s0(x2  ,X1+1,X2-1),s0(x2  ,X1,X2-1));
			  if (try<guess) guess = try;
			}
			if(time0[index-nx] < 1.e9 && time0[index-nxy-nx] < 1.e9
			   && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			  try = fdh3d(              t0(x2-1,X1,X2),
				      t0(x2-1,X1-1,X2),t0(x2-1,X1-1,X2-1),t0(x2-1,X1,X2-1),
				      t0(x2  ,X1-1,X2),t0(x2  ,X1-1,X2-1),t0(x2  ,X1,X2-1),
				      s0(x2,X1,X2), s0(x2-1,X1,X2),
				      s0(x2-1,X1-1,X2),s0(x2-1,X1-1,X2-1),s0(x2-1,X1,X2-1),
				      s0(x2  ,X1-1,X2),s0(x2  ,X1-1,X2-1),s0(x2  ,X1,X2-1));
			  if (try<guess) guess = try;
			}
			if(guess > 1.0e9){ 
			  if(time0[index+nx] < 1.e9 && X1<ny-1 && X2>z1+1 && X2<z2-1 )  {
			      try = fdhne(t0(x2-1,X1,X2),t0(x2-1,X1+1,X2),t0(x2,X1+1,X2),
					  t0(x2-1,X1+1,X2-1),t0(x2-1,X1+1,X2+1),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1+1,X2),s0(x2,X1+1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X1>0 && X2>z1+1 && X2<z2-1 )  {
			      try = fdhne(t0(x2-1,X1,X2),t0(x2-1,X1-1,X2),t0(x2,X1-1,X2),
					  t0(x2-1,X1-1,X2-1),t0(x2-1,X1-1,X2+1),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1-1,X2),s0(x2,X1-1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 && X1>y1+1 && X1<y2-1 )  {
			      try = fdhne(t0(x2-1,X1,X2),t0(x2-1,X1,X2+1),t0(x2,X1,X2+1),
					  t0(x2-1,X1-1,X2+1),t0(x2-1,X1+1,X2+1),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1,X2+1),s0(x2,X1,X2+1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nxy] < 1.e9 && X2>0 && X1>y1+1 && X1<y2-1 )  {
			      try = fdhne(t0(x2-1,X1,X2),t0(x2-1,X1,X2-1),t0(x2,X1,X2-1),
					  t0(x2-1,X1-1,X2-1),t0(x2-1,X1+1,X2-1),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1,X2-1),s0(x2,X1,X2-1) );
			    if (try<guess)  guess = try;
			  }
		        } 
			  if(time0[index+nx] < 1.e9 && X1<ny-1 )  {
			    try = fdh2d(t0(x2-1,X1,X2),t0(x2-1,X1+1,X2),t0(x2,X1+1,X2),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1+1,X2),s0(x2,X1+1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nx] < 1.e9 && X1>0 )  {
			    try = fdh2d(t0(x2-1,X1,X2),t0(x2-1,X1-1,X2),t0(x2,X1-1,X2),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1-1,X2),s0(x2,X1-1,X2) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nxy] < 1.e9 && X2<nz-1 )  {
			    try = fdh2d(t0(x2-1,X1,X2),t0(x2-1,X1,X2+1),t0(x2,X1,X2+1),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1,X2+1),s0(x2,X1,X2+1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index-nxy] < 1.e9 && X2>0 )  {
			    try = fdh2d(t0(x2-1,X1,X2),t0(x2-1,X1,X2-1),t0(x2,X1,X2-1),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2),s0(x2-1,X1,X2-1),s0(x2,X1,X2-1) );
			    if (try<guess)  guess = try;
			  }
			  if(time0[index+nx] < 1.e9 && time0[index+nxy+nx] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1<ny-1 ) {
			    try = fdh2d(t0(x2,X1+1,X2),t0(x2,X1+1,X2+1),t0(x2,X1,X2+1),
					s0(x2,X1,X2),
					s0(x2,X1+1,X2),s0(x2,X1+1,X2+1),s0(x2,X1,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index+nx] < 1.e9 && time0[index-nxy+nx] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1<ny-1 ) {
			    try = fdh2d(t0(x2,X1+1,X2),t0(x2,X1+1,X2-1),t0(x2,X1,X2-1),
					s0(x2,X1,X2),
					s0(x2,X1+1,X2),s0(x2,X1+1,X2-1),s0(x2,X1,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-nx] < 1.e9 && time0[index+nxy-nx] < 1.e9
			     && time0[index+nxy] < 1.e9 && X2<nz-1  && X1>0 ) {
			    try = fdh2d(t0(x2,X1-1,X2),t0(x2,X1-1,X2+1),t0(x2,X1,X2+1),
					s0(x2,X1,X2),
					s0(x2,X1-1,X2),s0(x2,X1-1,X2+1),s0(x2,X1,X2+1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if(time0[index-nx] < 1.e9 && time0[index-nxy-nx] < 1.e9
			     && time0[index-nxy] < 1.e9 && X2>0  && X1>0 ) {
			    try = fdh2d(t0(x2,X1-1,X2),t0(x2,X1-1,X2-1),t0(x2,X1,X2-1),
					s0(x2,X1,X2),
					s0(x2,X1-1,X2),s0(x2,X1-1,X2-1),s0(x2,X1,X2-1) );
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if(guess > 1.0e9){ 
			  if ( X1>y1+1 && X1<y2-1 && X2>z1+1 && X2<z2-1 ) {
			    try = fdhnf(t0(x2-1,X1,X2),
					  t0(x2-1,X1+1,X2),t0(x2-1,X1,X2+1),
					  t0(x2-1,X1-1,X2),t0(x2-1,X1,X2-1),
					  s0(x2,X1,X2),
					  s0(x2-1,X1,X2) );
			    if (try<guess)  guess = try;
			  }
			} 
			  try = t0(x2-1,X1,X2) + .5*(s0(x2,X1,X2)+s0(x2-1,X1,X2));
			  if (try<guess)  guess = try;
                          if ( time0[index+nx]<1.e9 && X1<ny-1 )  {
			    try = t0(x2,X1+1,X2) + .5*(s0(x2,X1,X2)+s0(x2,X1+1,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-nx]<1.e9 && X1>0 )  {
			    try = t0(x2,X1-1,X2) + .5*(s0(x2,X1,X2)+s0(x2,X1-1,X2));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index+nxy]<1.e9 && X2<nz-1 )  {
			    try = t0(x2,X1,X2+1) + .5*(s0(x2,X1,X2)+s0(x2,X1,X2+1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			  if ( time0[index-nxy]<1.e9 && X2>0 )  {
			    try = t0(x2,X1,X2-1) + .5*(s0(x2,X1,X2)+s0(x2,X1,X2-1));
			    if (try<guess)  {fhead=(guess-try)/slow0[index]; guess=try;}
			  }
			if (guess<time0[index]) {
				time0[index] = guess;
				if (fhead>headtest)  headw[2]++;
			}
		}
		if(x2 == nx-2) dx2 = 0;
		x2++;
	}
      }

		/* UPDATE RADIUS */
		radius++;
		if(radius%10 == 0 && verbose>5) vmess("Completed radius = %li",radius);
        if(radius == maxrad) rad0 = 0;

	}	/* END BIG LOOP */


	/* TEST IF REVERSE PROPAGATION IS NEEDED */

	if (headw[1]==0 && headw[2]==0 && headw[3]==0 && headw[4]==0 
		     && headw[5]==0 && headw[6]==0)
		reverse=0;
	else {
		head=0;
		if (headw[1]>0) {
			if(verbose) vmess("Head waves found on left: %li",headw[1]);
			if (headw[1]>head)  {
				head = headw[1];
				srcwall = 1;
			}
		}
		if (headw[2]>0) {
			if(verbose) vmess("Head waves found on right: %li",headw[2]);
			if (headw[2]>head)  {
				head = headw[2];
				srcwall = 2;
			}
		}
		if (headw[3]>0) {
			if(verbose) vmess("Head waves found on front: %li",headw[3]);
			if (headw[3]>head)  {
				head = headw[3];
				srcwall = 3;
			}
		}
		if (headw[4]>0) {
			if(verbose) vmess("Head waves found on back: %li",headw[4]);
			if (headw[4]>head)  {
				head = headw[4];
				srcwall = 4;
			}
		}
		if (headw[5]>0) {
			if(verbose) vmess("Head waves found on top: %li",headw[5]);
			if (headw[5]>head)  {
				head = headw[5];
				srcwall = 5;
			}
		}
		if (headw[6]>0) {
			if(verbose) vmess("Head waves found on bottom: %li",headw[6]);
			if (headw[6]>head)  {
				head = headw[6];
				srcwall = 6;
			}
		}
		if (headpref>0 && headw[headpref]>0) {
			if(verbose) 
				vmess("Preference to restart on wall opposite source");
			srcwall = headpref;
		}
		/* SET LOCATIONS OF SIDES OF THE CUBE SO THAT CUBE IS A FACE */
		dx1=1; dx2=1; dy1=1; dy2=1; dz1=1; dz2=1; rad0=1;
		radius = 1;
		if (srcwall == 1)	{  x2=1;
			vmess("RESTART at left side of model");  }
		else	{  x2=nx;	dx2=0;  }
		if (srcwall == 2)	{ x1=nx-2;
			vmess("RESTART at right side of model");  }
		else	{  x1= -1;	dx1=0;  }
		if (srcwall == 3)	{ y2=1;
			vmess("RESTART at front side of model");  }
		else	{  y2=ny;	dy2=0;  }
		if (srcwall == 4)	{ y1=ny-2;
			vmess("RESTART at back side of model");  }
		else	{  y1= -1;	dy1=0;  }
		if (srcwall == 5)	{ z2=1;
			vmess("RESTART at top side of model");  }
		else	{  z2=nz;	dz2=0;  }
		if (srcwall == 6)	{ z1=nz-2;
			vmess("RESTART at bottom side of model");  }
		else	{  z1= -1;	dz1=0;  }
		if (reverse == 0)  
			vwarn("RESTART CANCELLED by choice of input parameter `reverse`");
	}
	reverse--;

	}	/* END BIGGER LOOP - HOLE */

	free(sort);
	free(wall);
}

// int compar(struct sorted *a,struct sorted *b)
// {
// 	if(a->time > b->time) return(1);
// 	if(b->time > a->time) return(-1);
// 	else return(0);
// }
int compar(const void * a, const void * b)
{
	struct sorted *A = (struct sorted *)a;
	struct sorted *B = (struct sorted *)b;
	if(A->time > B->time) return(1);
	if(B->time > A->time) return(-1);
	else return(0);
}


/* 3D TRANSMISSION STENCIL
   STENCIL FROM VIDALE; CONDITIONS AND OTHER OPTIONS FROM HOLE
   JAH 11/91 */
float fdh3d(float  t1, float t2, float t3, float t4, float t5, float t6, float t7, float ss0, float s1, float s2, float s3, float s4, float s5, float s6, float s7)
     //float  t1,t2,t3,t4,t5,t6,t7,ss0,s1,s2,s3,s4,s5,s6,s7;
     /* ss0 at newpoint; s1,t1 adjacent on oldface;
	s2,t2 and s4,t4 on oldface adjacent to s1;
	s3,t3 on oldface diametrically opposite newpoint;
	s5,t5 on newface adjacent to newpoint AND to s2;
	s6,t6 on newface diagonal to newpoint (adjacent to s3);
	s7,t7 on newface adjacent to newpoint AND to s4
	*/
{
  float x,slo;
  double sqrt();
  slo = .125*(ss0+s1+s2+s3+s4+s5+s6+s7);
  x = 6.*slo*slo - (t4-t2)*(t4-t2) - (t2-t6)*(t2-t6) - (t6-t4)*(t6-t4)
                 - (t7-t5)*(t7-t5) - (t5-t1)*(t5-t1) - (t1-t7)*(t1-t7);
  if (x>=0.)  {
    x = t3 + sqrt(.5*x);
    if ( (x<t1) || (x<t2) || (x<t4) || (x<t5) || (x<t6) || (x<t7) )  
      x = 1.e11;   /* ACAUSAL; ABORT */
  }
  else  x = 1.e11;   /* SQRT IMAGINARY; ABORT */
  return(x);
}

/* 3D STENCIL FOR NEW EDGE
   STENCIL FROM VIDALE; CONDITIONS AND OTHER OPTIONS FROM HOLE
   JAH 11/91 */
float fdhne(float  t1, float t2, float t3, float t4, float t5, float ss0, float s1, float s2, float s3)
     //float  t1,t2,t3,t4,t5,ss0,s1,s2,s3;
     /* ss0 at newpoint; s1,t1 adjacent on oldface;
	s2,t2 diagonal on oldface; s3,t3 adjacent on newface;
	t4,t5 beside t2 on old face opposite each other */
{
  float x,slo;
  double sqrt();
  slo = .25*(ss0+s1+s2+s3);
  x = 2.*slo*slo - (t3-t1)*(t3-t1) - .5*(t5-t4)*(t5-t4);
  if (x>=0.)  {
    x = t2 + sqrt(x);
    if ( (x<t1) || (x<t3) || (x<t4) || (x<t5) )     /* ACAUSAL; ABORT */
      x = 1.e11;
  }
  else  x = 1.e11;   /* SQRT IMAGINARY; ABORT */
  return(x);
}

/* 2D TRANSMISSION STENCIL (FOR HEAD WAVES ON FACES OF GRID CELLS)
   STENCIL FROM VIDALE (1988 2D PAPER); CONDITIONS AND OTHER OPTIONS FROM HOLE
   JAH 11/91 */
float fdh2d(float t1, float t2, float t3, float ss0, float s1, float s2, float s3)
     //float  t1,t2,t3,ss0,s1,s2,s3;
     /* ss0 at newpoint; s1,t1 & s3,t3 adjacent; s2,t2 diagonal
      */
{
  float x,slo;
  double sqrt();
  slo = .25*(ss0+s1+s2+s3);
  x = 2.*slo*slo - (t3-t1)*(t3-t1);
  if (x>=0.)  {
    x = t2 + sqrt(x);
    if ( (x<t1) || (x<t3) )  x = 1.e11;   /* ACAUSAL; ABORT */
  }
  else  x = 1.e11;   /* SQRT IMAGINARY; ABORT */
  return(x);
}

/* 3D STENCIL FOR NEW FACE
   STENCIL FROM VIDALE; CONDITIONS AND OTHER OPTIONS FROM HOLE
   JAH 11/91 */
float fdhnf(float t1, float t2, float t3, float t4, float t5, float ss0, float s1)
     //float  t1,t2,t3,t4,t5,ss0,s1;
     /* ss0 at newpoint; s1,t1 adjacent on old face;
	t2,t4 beside t1 on old face and opposite each other;
	t3,t5 beside t1 on old face and opposite each other
	*/
{
  float x,slo;
  double sqrt();
  slo = .5*(ss0+s1);
  x = slo*slo - .25*( (t4-t2)*(t4-t2) + (t5-t3)*(t5-t3) );
  if (x>=0.)  {
    x = t1 + sqrt(x);
    if ( (x<t2) || (x<t3) || (x<t4) || (x<t5) )     /* ACAUSAL; ABORT */
      x = 1.e11;
  }
  else  x = 1.e11;   /* SQRT IMAGINARY; ABORT */
  return(x);
}


