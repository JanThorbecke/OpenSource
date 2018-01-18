/*  FILE:  vidale.c
 *  AUTHOR:  Joseph R. Matarese
 *  DATE:  Decenber 18, 1992
 *  Copyright (c) 1992 Joseph R. Matarese and Massachusetts Institute of
 *    Technology
 *
 *  Permission is granted to copy or modify this code under the condition
 *  that the above copyright notice is maintained.  The author and MIT
 *  make ABSOLUTELY NO WARRANTY regarding the fitness of this code for
 *  any purpose.
 *
 */

//#include <raytime2.h>
#include <raytime.h>

extern void  near_source(float *ttime, float *slow, icoord *nm, icoord *isrc, fcoord *scale, int order);

extern void  corner(float *ttime, float *slow,icoord *nm, icoord *iin, icoord *iout, fcoord *scale, struct s_ecount *ecount);

extern void  side(float *ttime, float *slow, icoord *nm, icoord *iin, icoord *iout, fcoord *scale, struct s_ecount *ecount);

void rm_head(float *slow, icoord *ndim, icoord *isrc, int mzrcv, float dx, int *nzm);

void vidale(float *ttime, float *slow1, icoord *isrc, icoord grid, float dx, int order, int mzrcv)
{
  int             iz_lo, iz_hi, ix_lo, ix_hi, iz, ix, node_src;
  int             nzm;
  icoord          iin, iout, dnm, *nm;
  fcoord          dscale, *scale;
  short           ascending, finished;
  struct s_ecount ecount;
  float sx, sz, sign, *trueslow, *slow;

  dscale.x = dx; dscale.y = 0.; dscale.z = dx;
  scale = &dscale;
  dnm.x = grid.x; dnm.y = 1; dnm.z = grid.z;
  nm = &dnm;

/* transpose slowness field for usage in vidale */
  trueslow = (float *)calloc(grid.x*grid.z,sizeof(float));
  for (ix=0; ix<grid.x; ix++) {
    for (iz=0; iz<grid.z; iz++) {
      trueslow[iz*grid.x + ix] = slow1[ix*grid.z + iz];
    }
  }
  slow = trueslow;
  rm_head(slow,&grid,isrc,mzrcv,dx,&nzm);
 
  //for actual position of source not on grid: current no-use
  node_src = isrc->z*grid.x + isrc->x;
  sx = isrc->x*dx-isrc->x*dx;
  sz = isrc->z*dx-isrc->z*dx;
  if (sz < 0) sign = 1;
  else sign = -1;
  ttime[node_src] = sign*sqrt(sx*sx+sz*sz)*slow[node_src];

  if(nm->y != 1)
    verr("only 2D models implemented");

  ecount.corner = ecount.corner_min = ecount.side = 0;

/*---------------------------------------------------------------------------*
 * Do near source region.
 *---------------------------------------------------------------------------*/

  near_source(ttime,slow,nm,isrc,scale,order);

  iz_hi = min2(isrc->z+order,nm->z-1);
  iz_lo = max2(isrc->z-order,0);
  ix_hi = min2(isrc->x+order,nm->x-1);
  ix_lo = max2(isrc->x-order,0);

/*---------------------------------------------------------------------------*
 * Loop over outer ring - verticals, then horizontals.
 *---------------------------------------------------------------------------*/

  finished = FALSE;
  while(!finished) {

    finished = TRUE;
    if(ix_hi < nm->x-1) {
      finished = FALSE;
      iin.x = ix_hi;                                     /* right-hand edge */
      iout.x = ix_hi+1;
      
      ascending = FALSE;
      if(ttime[iz_lo*nm->x+iin.x] <= ttime[(iz_lo+1)*nm->x+iin.x]) {
	ttime[iz_lo*nm->x+iout.x] = ttime[iz_lo*nm->x+iin.x] + 0.5 *
	  scale->x * (slow[iz_lo*nm->x+iout.x] + slow[iz_lo*nm->x+iin.x]);
	ecount.corner_min++;
      }
      if(ttime[iz_lo*nm->x+iin.x] < ttime[(iz_lo+1)*nm->x+iin.x]) {
	ascending = TRUE;
        iin.z = iz_lo;
        iout.z = iz_lo+1;
        corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
      }
      for(iz=iz_lo+1; iz<=iz_hi-1; iz++) {
	iin.z = iz;
	if(ttime[iz*nm->x+iin.x] == ttime[(iz+1)*nm->x+iin.x]) {
	  iout.z = iz;
	  side(ttime,slow,nm,&iin,&iout,scale,&ecount);
	}
	if(ttime[iz*nm->x+iin.x] < ttime[(iz+1)*nm->x+iin.x]) {
	  if(!ascending) {
	    iout.z = iz;
	    side(ttime,slow,nm,&iin,&iout,scale,&ecount);
	  }
	  ascending = TRUE;
	  iout.z = iz+1;
	  corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
	} else ascending = FALSE;
      }
      
      if(ttime[iz_hi*nm->x+iin.x] <= ttime[(iz_hi-1)*nm->x+iin.x]) {
	ttime[iz_hi*nm->x+iout.x] = ttime[iz_hi*nm->x+iin.x] + 0.5 *
	  scale->x * (slow[iz_hi*nm->x+iout.x] + slow[iz_hi*nm->x+iin.x]);
	ecount.corner_min++;
      }
      for(iz=iz_hi; iz>=iz_lo+1; iz--) {
	if(ttime[iz*nm->x+iin.x] < ttime[(iz-1)*nm->x+iin.x]) {
	  iin.z = iz;
	  iout.z = iz-1;
	  corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
	}
      }

    }
    
    if(ix_lo > 0) {
      finished = FALSE;
      iin.x = ix_lo;                                     /* left-hand edge */
      iout.x = ix_lo-1;
      
      ascending = FALSE;
      if(ttime[iz_lo*nm->x+iin.x] <= ttime[(iz_lo+1)*nm->x+iin.x]) {
	ttime[iz_lo*nm->x+iout.x] = ttime[iz_lo*nm->x+iin.x] + 0.5 *
	  scale->x * (slow[iz_lo*nm->x+iout.x] + slow[iz_lo*nm->x+iin.x]);
	ecount.corner_min++;
      }
      if(ttime[iz_lo*nm->x+iin.x] < ttime[(iz_lo+1)*nm->x+iin.x]) {
	ascending = TRUE;
        iin.z = iz_lo;
        iout.z = iz_lo+1;
        corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
      }
      for(iz=iz_lo+1; iz<=iz_hi-1; iz++) {
	iin.z = iz;
	if(ttime[iz*nm->x+iin.x] == ttime[(iz+1)*nm->x+iin.x]) {
	  iout.z = iz;
	  side(ttime,slow,nm,&iin,&iout,scale,&ecount);
	}
	if(ttime[iz*nm->x+iin.x] < ttime[(iz+1)*nm->x+iin.x]) {
	  if(!ascending) {
	    iout.z = iz;
	    side(ttime,slow,nm,&iin,&iout,scale,&ecount);
	  }
	  ascending = TRUE;
	  iout.z = iz+1;
	  corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
	} else ascending = FALSE;
      }
      
      if(ttime[iz_hi*nm->x+iin.x] <= ttime[(iz_hi-1)*nm->x+iin.x]) {
	ttime[iz_hi*nm->x+iout.x] = ttime[iz_hi*nm->x+iin.x] + 0.5 *
	  scale->x * (slow[iz_hi*nm->x+iout.x] + slow[iz_hi*nm->x+iin.x]);
	ecount.corner_min++;
      }
      for(iz=iz_hi; iz>=iz_lo+1; iz--) {
	if(ttime[iz*nm->x+iin.x] < ttime[(iz-1)*nm->x+iin.x]) {
	  iin.z = iz;
	  iout.z = iz-1;
	  corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
	}
      }

    }
    
    if(iz_hi < nm->z-1) {
      finished = FALSE;
      iin.z = iz_hi;                                     /* bottom edge */
      iout.z = iz_hi+1;
      
      ascending = FALSE;
      if(ttime[iin.z*nm->x+ix_lo] <= ttime[iin.z*nm->x+ix_lo+1]) {
	ttime[iout.z*nm->x+ix_lo] = ttime[iin.z*nm->x+ix_lo] + 0.5 *
	  scale->x * (slow[iout.z*nm->x+ix_lo] + slow[iin.z*nm->x+ix_lo]);
	ecount.corner_min++;
      }
      if(ttime[iin.z*nm->x+ix_lo] < ttime[iin.z*nm->x+ix_lo+1]) {
	ascending = TRUE;
	iin.x = ix_lo;
	iout.x = ix_lo+1;
        corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
      }
      for(ix=ix_lo+1; ix<=ix_hi-1; ix++) {
	iin.x = ix;
	if(ttime[iin.z*nm->x+ix] == ttime[iin.z*nm->x+ix+1]) {
	  iout.x = ix;
	  side(ttime,slow,nm,&iin,&iout,scale,&ecount);
	}
	if(ttime[iin.z*nm->x+ix] < ttime[iin.z*nm->x+ix+1]) {
	  if(!ascending) {
	    iout.x = ix;
	    side(ttime,slow,nm,&iin,&iout,scale,&ecount);
	  }
	  ascending = TRUE;
	  iout.x = ix+1;
	  corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
	} else ascending = FALSE;
      }
      
      if(ttime[iin.z*nm->x+ix_hi] <= ttime[iin.z*nm->x+ix_hi-1]) {
	ttime[iout.z*nm->x+ix_hi] = ttime[iin.z*nm->x+ix_hi] + 0.5 *
	  scale->x * (slow[iout.z*nm->x+ix_hi] + slow[iin.z*nm->x+ix_hi]);
	ecount.corner_min++;
      }
      for(ix=ix_hi; ix>=ix_lo+1; ix--) {
	if(ttime[iin.z*nm->x+ix] < ttime[iin.z*nm->x+ix-1]) {
	  iin.x = ix;
	  iout.x = ix-1;
	  corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
	}
      }

    }
    
    if(iz_lo > 0) {
      finished = FALSE;
      iin.z = iz_lo;                                     /* top edge */
      iout.z = iz_lo-1;
      
      ascending = FALSE;
      if(ttime[iin.z*nm->x+ix_lo] <= ttime[iin.z*nm->x+ix_lo+1]) {
	ttime[iout.z*nm->x+ix_lo] = ttime[iin.z*nm->x+ix_lo] + 0.5 *
	  scale->x * (slow[iout.z*nm->x+ix_lo] + slow[iin.z*nm->x+ix_lo]);
	ecount.corner_min++;
      }
      if(ttime[iin.z*nm->x+ix_lo] < ttime[iin.z*nm->x+ix_lo+1]) {
	ascending = TRUE;
	iin.x = ix_lo;
	iout.x = ix_lo+1;
        corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
      }
      for(ix=ix_lo+1; ix<=ix_hi-1; ix++) {
	iin.x = ix;
	if(ttime[iin.z*nm->x+ix] == ttime[iin.z*nm->x+ix+1]) {
	  iout.x = ix;
	  side(ttime,slow,nm,&iin,&iout,scale,&ecount);
	}
	if(ttime[iin.z*nm->x+ix] < ttime[iin.z*nm->x+ix+1]) {
	  if(!ascending) {
	    iout.x = ix;
	    side(ttime,slow,nm,&iin,&iout,scale,&ecount);
	  }
	  ascending = TRUE;
	  iout.x = ix+1;
	  corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
	} else ascending = FALSE;
      }
      
      if(ttime[iin.z*nm->x+ix_hi] <= ttime[iin.z*nm->x+ix_hi-1]) {
	ttime[iout.z*nm->x+ix_hi] = ttime[iin.z*nm->x+ix_hi] + 0.5 *
	  scale->x * (slow[iout.z*nm->x+ix_hi] + slow[iin.z*nm->x+ix_hi]);
	ecount.corner_min++;
      }
      for(ix=ix_hi; ix>=ix_lo+1; ix--) {
	if(ttime[iin.z*nm->x+ix] < ttime[iin.z*nm->x+ix-1]) {
	  iin.x = ix;
	  iout.x = ix-1;
	  corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
	}
      }

    }

    if((ix_hi < nm->x-1) || (iz_hi < nm->z-1)) {     /* bottom right corner */
      iin.x = (ix_hi < nm->x-1) ? ix_hi : nm->x-2;
      iout.x = iin.x+1;
      iin.z = (iz_hi < nm->z-1) ? iz_hi : nm->z-2;
      iout.z = iin.z+1;
      corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
    }

    if((ix_lo > 0) || (iz_hi < nm->z-1)) {            /* bottom left corner */
      iin.x = (ix_lo > 0) ? ix_lo : 1;
      iout.x = iin.x-1;
      iin.z = (iz_hi < nm->z-1) ? iz_hi : nm->z-2;
      iout.z = iin.z+1;
      corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
    }

    if((ix_hi < nm->x-1) || (iz_lo > 0)) {              /* top right corner */
      iin.x = (ix_hi < nm->x-1) ? ix_hi : nm->x-2;
      iout.x = iin.x+1;
      iin.z = (iz_lo > 0) ? iz_lo : 1;
      iout.z = iin.z-1;
      corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
    }

    if((ix_lo > 0) || (iz_lo > 0)) {                     /* top left corner */
      iin.x = (ix_lo > 0) ? ix_lo : 1;
      iout.x = iin.x-1;
      iin.z = (iz_lo > 0) ? iz_lo : 1;
      iout.z = iin.z-1;
      corner(ttime,slow,nm,&iin,&iout,scale,&ecount);
    }

    ix_hi = min2(ix_hi+1,nm->x-1);
    ix_lo = max2(ix_lo-1,0);
    iz_hi = min2(iz_hi+1,nm->z-1);
    iz_lo = max2(iz_lo-1,0);

  } /* while(!finished) */

/*
  jm_message("debug",__FILE__,__LINE__,
	     "errors:\n   %s: %d\n   %s: %d\n   %s: %d",
	     "corner minima",ecount.corner_min,
	     "corner negative sqrts",ecount.corner,
	     "side negative sqrts",ecount.side);
 */
  return;
}



/*  FILE:  near_source.c
 *  AUTHOR:  Joseph R. Matarese
 *  DATE:  Decenber 18, 1992
 *  Copyright (c) 1992 Joseph R. Matarese and Massachusetts Institute of
 *    Technology
 *
 *  Permission is granted to copy or modify this code under the condition
 *  that the above copyright notice is maintained.  The author and MIT
 *  make ABSOLUTELY NO WARRANTY regarding the fitness of this code for
 *  any purpose.
 *
 */

/*---------------------------------------------------------------------------*
 *
 *  Calculate near source traveltimes.
 *
 *
 *
 *---------------------------------------------------------------------------*/


void  near_source(float *ttime, float *slow, icoord *nm, icoord *isrc, fcoord *scale, int order)
{
  int                 ix_lo, ix_hi, iz_lo, iz_hi, ix, iz;
  float               slow_0, ttime_0, dist;
  
  if(nm->y != 1)
    verr("only 2D models implemented");

/*---------------------------------------------------------------------------*
 * Boundaries of source region.
 *---------------------------------------------------------------------------*/

  ix_lo = max2(isrc->x-order,0);
  ix_hi = min2(isrc->x+order,nm->x-1);
  iz_lo = max2(isrc->z-order,0);
  iz_hi = min2(isrc->z+order,nm->z-1);

/*---------------------------------------------------------------------------*
 * Calculate traveltimes for source region.
 * We assume this region to nearly homogeneous.
 *---------------------------------------------------------------------------*/

  slow_0 = slow[isrc->z*nm->x+isrc->x];
  ttime_0 = ttime[isrc->z*nm->x+isrc->x];
  for(iz=iz_lo; iz<=iz_hi; iz++) {
    for(ix=ix_lo; ix<=ix_hi; ix++) {
      dist = hypot((ix-isrc->x)*scale->x,(iz-isrc->z)*scale->z);
      ttime[iz*nm->x+ix] = 0.5*dist*(slow[iz*nm->x+ix] + slow_0) + ttime_0;
    }
  }

  return;
}




/*  FILE:  side.c
 *  AUTHOR:  Joseph R. Matarese
 *  DATE:  Decenber 18, 1992
 *  Copyright (c) 1992 Joseph R. Matarese and Massachusetts Institute of
 *    Technology
 *
 *  Permission is granted to copy or modify this code under the condition
 *  that the above copyright notice is maintained.  The author and MIT
 *  make ABSOLUTELY NO WARRANTY regarding the fitness of this code for
 *  any purpose.
 *
 */

/*---------------------------------------------------------------------------*
 *
 * Apply the side stencil.
 *
 *---------------------------------------------------------------------------*/

enum e_orientation { vertical, horizontal };

void  side(float *ttime, float *slow, icoord *nm, icoord *iin, icoord *iout, fcoord *scale, struct s_ecount *ecount)
{
  enum e_orientation  orientation;
  float				dt, slow_0, operand, ttnew, tttmp;
  static double		c_sqrt2;
  static char		initialized = (char)0;

  if(!initialized) {
    if(nm->y != 1)
    	verr("only 2D models implemented");
    if(scale->x != scale->z)
		verr("only uniform grid implemented");
    c_sqrt2 = sqrt(2.);
    initialized = (char)1;
  }

  orientation = (iin->x == iout->x) ? horizontal : vertical;

  if(orientation == vertical) {
    dt = ttime[(iin->z+1)*nm->x+iin->x] - ttime[(iin->z-1)*nm->x+iin->x];
    slow_0 = 0.25*(slow[(iin->z+1)*nm->x+iin->x] +
		   slow[(iin->z-1)*nm->x+iin->x] +
		   slow[iin->z*nm->x+iin->x] + slow[iin->z*nm->x+iout->x]);

    if((operand = scale->x*scale->x*slow_0*slow_0 - 0.25*dt*dt) < 0.) {
      slow_0 = 0.5*(slow[(iin->z-1)*nm->x+iin->x]+slow[iin->z*nm->x+iout->x]);
      ttnew = ttime[(iin->z-1)*nm->x+iin->x] + c_sqrt2*scale->x*slow_0;

      slow_0 = 0.5*(slow[(iin->z+1)*nm->x+iin->x]+slow[iin->z*nm->x+iout->x]);
      tttmp = ttime[(iin->z+1)*nm->x+iin->x] + c_sqrt2*scale->x*slow_0;
      if(tttmp < ttnew) ttnew = tttmp;

      slow_0 = 0.5*(slow[iin->z*nm->x+iin->x] + slow[iin->z*nm->x+iout->x]);
      tttmp = ttime[iin->z*nm->x+iin->x] + scale->x*slow_0;
      if(tttmp < ttnew) ttnew = tttmp;

      ecount->side++;
    } else ttnew = ttime[iin->z*nm->x+iin->x] + sqrt(operand);
    ttime[iin->z*nm->x+iout->x] = ttnew;

  } else {
    dt = ttime[iin->z*nm->x+iin->x+1] - ttime[iin->z*nm->x+iin->x-1];
    slow_0 = 0.25*(slow[iin->z*nm->x+iin->x+1] + slow[iin->z*nm->x+iin->x-1] +
		   slow[iin->z*nm->x+iin->x] + slow[iout->z*nm->x+iin->x]);

    if((operand = scale->x*scale->x*slow_0*slow_0 - 0.25*dt*dt) < 0.) {
      slow_0 = 0.5*(slow[iin->z*nm->x+iin->x-1] + slow[iout->z*nm->x+iin->x]);
      ttnew = ttime[iin->z*nm->x+iin->x-1] + c_sqrt2*scale->x*slow_0;

      slow_0 = 0.5*(slow[iin->z*nm->x+iin->x+1] + slow[iout->z*nm->x+iin->x]);
      tttmp = ttime[iin->z*nm->x+iin->x+1] + c_sqrt2*scale->x*slow_0;
      if(tttmp < ttnew) ttnew = tttmp;

      slow_0 = 0.5*(slow[iin->z*nm->x+iin->x] + slow[iout->z*nm->x+iin->x]);
      tttmp = ttime[iin->z*nm->x+iin->x] + scale->x*slow_0;
      if(tttmp < ttnew) ttnew = tttmp;

      ecount->side++;
    } else ttnew = ttime[iin->z*nm->x+iin->x] + sqrt(operand);
    ttime[iout->z*nm->x+iin->x] = ttnew;

  }  

  return;
}




/*  FILE:  corner.c
 *  AUTHOR:  Joseph R. Matarese
 *  DATE:  Decenber 18, 1992
 *  Copyright (c) 1992 Joseph R. Matarese and Massachusetts Institute of
 *    Technology
 *
 *  Permission is granted to copy or modify this code under the condition
 *  that the above copyright notice is maintained.  The author and MIT
 *  make ABSOLUTELY NO WARRANTY regarding the fitness of this code for
 *  any purpose.
 *
 */

/*---------------------------------------------------------------------------*
 *
 * Apply the corner stencil.
 *
 *---------------------------------------------------------------------------*/


void  corner(float *ttime, float *slow,icoord *nm, icoord *iin, icoord *iout, fcoord *scale, struct s_ecount *ecount)
{
  float               dt, slow_0, operand, ttnew, tttmp;
  static double		  c_sqrt2;
  static char         initialized = (char)0;

  if(!initialized) {
    if(nm->y != 1)
    	verr("only 2D models implemented");
    if(scale->x != scale->z)
		verr("only uniform grid implemented");
    c_sqrt2 = sqrt(2.);
    initialized = (char)1;
  }

  dt = ttime[iout->z*nm->x+iin->x] - ttime[iin->z*nm->x+iout->x];
  slow_0 = 0.25*(slow[iout->z*nm->x+iin->x] + slow[iin->z*nm->x+iout->x] +
		 slow[iout->z*nm->x+iout->x] + slow[iin->z*nm->x+iin->x]);
  
  if((operand = 2*scale->x*scale->x*slow_0*slow_0 - dt*dt) < 0.) {
    slow_0 = 0.5*(slow[iin->z*nm->x+iout->x] + slow[iout->z*nm->x+iout->x]);
    ttnew = ttime[iin->z*nm->x+iout->x] + scale->x*slow_0;

    slow_0 = 0.5*(slow[iout->z*nm->x+iin->x] + slow[iout->z*nm->x+iout->x]);
    tttmp = ttime[iout->z*nm->x+iin->x] + scale->x*slow_0;
    if(tttmp < ttnew) ttnew = tttmp;

    slow_0 = 0.5*(slow[iout->z*nm->x+iout->x] + slow[iin->z*nm->x+iin->x]);
    tttmp = ttime[iin->z*nm->x+iin->x] + c_sqrt2*scale->x*slow_0;
    if(tttmp < ttnew) ttnew = tttmp;

    ecount->corner++;
  } else ttnew = ttime[iin->z*nm->x+iin->x] + sqrt(operand);

/*---------------------------------------------------------------------------*
 * Nonuniform grid?
 *
  scale2->x = scale->x * scale->x;
  scale2->z = scale->z * scale->z;
  scale2_plus = scale2->z + scale2->x;
  scale2_minus = scale2->z + scale2->x;

  ttnew = ttime[iin->z*nm->x+iin->x] +
    (scale2_minus * dt + 2 * scale->z * scale->x * 
     sqrt(scale2_plus * slow_0 * slow_0 - dt * dt)) / scale2_plus;
 *---------------------------------------------------------------------------*/

  if(ttnew < ttime[iout->z*nm->x+iout->x])
    ttime[iout->z*nm->x+iout->x] = ttnew;

  return;
}


void rm_head(float *slow, icoord *ndim, icoord *isrc, int mzrcv, float dx, int *nzm)
{
	int iz, ix, k, l, i, izn, lprev, nz, nx;
	float brd, val, dz;

	iz = isrc->z;
	ix = isrc->x;
	nz = ndim->z;
	nx = ndim->x;
	dz = dx;

	if (iz == 0) { ndim->z = 2; return;}

		if (iz >= mzrcv) { 
			brd = slow[iz*nx+ix];
			val = slow[(iz-1)*nx+ix];
			while ((brd == val) && iz < nz-1) brd = slow[++iz*nx+ix];

			lprev = iz;
			for (l = iz; l > 0; l--) {
				for (k = 0; k < nx; k++) {
					if (slow[l*nx+k] == brd) {
						slow[l*nx+k] = val;
						lprev = l;
					}
				}
				if (lprev-l) l = 0;
			}

			iz = isrc->z;
		}
		else {
			brd = slow[iz*nx+ix];
			val = slow[(iz+1)*nx+ix];
			while (brd == val && iz > 0) brd = slow[--iz*nx+ix];

			if (iz < nz) {
				lprev = iz;
				for (l = iz; l < nz; l++) {
					for (k = 0; k < nx; k++) {
						if (slow[l*nx+k] == brd) {
							slow[l*nx+k] = val;
							lprev = l;
						}
					}
					if (lprev-l) l = nz;
				}
			}
			iz = mzrcv;
		}

	*nzm = iz+1;

	return;
}

