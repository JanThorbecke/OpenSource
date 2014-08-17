#include "genfft.h"
#include "pfaconst.h" 

/* Copyright (c) Colorado School of Mines, 1995.*/
/* All rights reserved.                       */

// stolen by Joerg Arndt from the cwplib ...
//
// Original author:  Dave Hale, Colorado School of Mines, 04/27/89
//
// edited by Joerg Arndt (july 96):
// C++ only 
// changed float to REAL
// replaced 'magic numbers' (P120 etc), they moved to pfaconst.h 
// 
// see pfafft.doc !


int 
npfa (int nmin)
{
    int i;
    for (i=0; i<NTAB-1 && nctab[i].n<nmin; ++i);
    return nctab[i].n;
}

int 
npfao (int nmin, int nmax)
{
    int i,j;
    for (i=0; i<NTAB-1 && nctab[i].n<nmin; ++i);
    for (j=i+1; j<NTAB-1 && nctab[j].n<=nmax; ++j)
    if (nctab[j].c<nctab[i].c) i = j;
    return nctab[i].n;
}

int 
npfar (int nmin)
{
    return 2*npfa((nmin+1)/2);
}

int 
npfaro (int nmin, int nmax)
{
    return 2*npfao((nmin+1)/2,(nmax+1)/2);
}

void 
pfacc (int isign, int n, complex cz[])
{
    static int kfax[] = { 16,13,11,9,8,7,5,4,3,2 };
    REAL *z=(REAL*)cz;
    int j0,j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15,jt;
    int nleft,jfax,ifac,jfac,jinc,jmax,ndiv,m,mm=0,mu=0,l;
    REAL t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,
    t6r,t6i,t7r,t7i,t8r,t8i,t9r,t9i,t10r,t10i,
    t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,
    t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t20r,t20i,
    t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,
    t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t30r,t30i,
    t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,
    t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t40r,t40i,
    t41r,t41i,t42r,t42i,
    y1r,y1i,y2r,y2i,y3r,y3i,y4r,y4i,y5r,y5i,
    y6r,y6i,y7r,y7i,y8r,y8i,y9r,y9i,y10r,y10i,
    y11r,y11i,y12r,y12i,y13r,y13i,y14r,y14i,y15r,y15i,
    c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;

    /* keep track of n left after dividing by factors */
    nleft = n;

    /* begin loop over possible factors (from biggest to smallest) */
    for (jfax=0; jfax<NFAX; jfax++) {

	/* skip if not a mutually prime factor of n */
        ifac = kfax[jfax];
        ndiv = nleft/ifac;
        if (ndiv*ifac!=nleft) continue;
 
	/* update n left and determine n divided by factor */
        nleft = ndiv;
        m = n/ifac;
 
	/* determine rotation factor mu and stride mm */
        for (jfac=1; jfac<=ifac; jfac++) {
	    mu = jfac;
	    mm = jfac*m;
	    if (mm%ifac==1) break;
	}
 
	/* adjust rotation factor for sign of transform */
        if (isign<0) mu = ifac-mu;
 
	/* compute stride, limit, and pointers */
        jinc = 2*mm;
	jmax = 2*n;
        j0 = 0;
        j1 = j0+jinc;

	/* if factor is 2 */
        if (ifac==2) {
	    for (l=0; l<m; l++) {
		t1r = z[j0]-z[j1];
		t1i = z[j0+1]-z[j1+1];
		z[j0] = z[j0]+z[j1];
		z[j0+1] = z[j0+1]+z[j1+1];
		z[j1] = t1r;
		z[j1+1] = t1i;
		jt = j1+2;
		j1 = j0+2;
		j0 = jt;
	    }
	    continue;
	}
        j2 = j1+jinc;
        if (j2>=jmax) j2 = j2-jmax;

	/* if factor is 3 */
        if (ifac==3) {
	    if (mu==1)
	    c1 = P866;
	    else
	    c1 = -P866;
	    for (l=0; l<m; l++) {
		t1r = z[j1]+z[j2];
		t1i = z[j1+1]+z[j2+1];
		y1r = z[j0]-0.5*t1r;
		y1i = z[j0+1]-0.5*t1i;
		y2r = c1*(z[j1]-z[j2]);
		y2i = c1*(z[j1+1]-z[j2+1]);
		z[j0] = z[j0]+t1r;
		z[j0+1] = z[j0+1]+t1i;
		z[j1] = y1r-y2i;
		z[j1+1] = y1i+y2r;
		z[j2] = y1r+y2i;
		z[j2+1] = y1i-y2r;
		jt = j2+2;
		j2 = j1+2;
		j1 = j0+2;
		j0 = jt;
	    }
	    continue;
	}
	j3 = j2+jinc;
	if (j3>=jmax) j3 = j3-jmax;

	/* if factor is 4 */
	if (ifac==4) {
	    if (mu==1)
	    c1 = 1.0;
	    else
	    c1 = -1.0;
	    for (l=0; l<m; l++) {
		t1r = z[j0]+z[j2];
		t1i = z[j0+1]+z[j2+1];
		t2r = z[j1]+z[j3];
		t2i = z[j1+1]+z[j3+1];
		y1r = z[j0]-z[j2];
		y1i = z[j0+1]-z[j2+1];
		y3r = c1*(z[j1]-z[j3]);
		y3i = c1*(z[j1+1]-z[j3+1]);
		z[j0] = t1r+t2r;
		z[j0+1] = t1i+t2i;
		z[j1] = y1r-y3i;
		z[j1+1] = y1i+y3r;
		z[j2] = t1r-t2r;
		z[j2+1] = t1i-t2i;
		z[j3] = y1r+y3i;
		z[j3+1] = y1i-y3r;
		jt = j3+2;
		j3 = j2+2;
		j2 = j1+2;
		j1 = j0+2;
		j0 = jt;
	    }
	    continue;
	}
	j4 = j3+jinc;
	if (j4>=jmax) j4 = j4-jmax;

	/* if factor is 5 */
	if (ifac==5) {
	    if (mu==1) {
		c1 = P559;
		c2 = P951;
		c3 = P587;
	    } else if (mu==2) {
		c1 = -P559;
		c2 = P587;
		c3 = -P951;
	    } else if (mu==3) {
		c1 = -P559;
		c2 = -P587;
		c3 = P951;
	    } else { 
		c1 = P559;
		c2 = -P951;
		c3 = -P587;
	    }
	    for (l=0; l<m; l++) {
		t1r = z[j1]+z[j4];
		t1i = z[j1+1]+z[j4+1];
		t2r = z[j2]+z[j3];
		t2i = z[j2+1]+z[j3+1];
		t3r = z[j1]-z[j4];
		t3i = z[j1+1]-z[j4+1];
		t4r = z[j2]-z[j3];
		t4i = z[j2+1]-z[j3+1];
		t5r = t1r+t2r;
		t5i = t1i+t2i;
		t6r = c1*(t1r-t2r);
		t6i = c1*(t1i-t2i);
		t7r = z[j0]-0.25*t5r;
		t7i = z[j0+1]-0.25*t5i;
		y1r = t7r+t6r;
		y1i = t7i+t6i;
		y2r = t7r-t6r;
		y2i = t7i-t6i;
		y3r = c3*t3r-c2*t4r;
		y3i = c3*t3i-c2*t4i;
		y4r = c2*t3r+c3*t4r;
		y4i = c2*t3i+c3*t4i;
		z[j0] = z[j0]+t5r;
		z[j0+1] = z[j0+1]+t5i;
		z[j1] = y1r-y4i;
		z[j1+1] = y1i+y4r;
		z[j2] = y2r-y3i;
		z[j2+1] = y2i+y3r;
		z[j3] = y2r+y3i;
		z[j3+1] = y2i-y3r;
		z[j4] = y1r+y4i;
		z[j4+1] = y1i-y4r;
		jt = j4+2;
		j4 = j3+2;
		j3 = j2+2;
		j2 = j1+2;
		j1 = j0+2;
		j0 = jt;
	    }
	    continue;
	}
	j5 = j4+jinc;
	if (j5>=jmax) j5 = j5-jmax;
	j6 = j5+jinc;
	if (j6>=jmax) j6 = j6-jmax;

	/* if factor is 7 */
	if (ifac==7) {
	    if (mu==1) {
		c1 = P623;
		c2 = -P222;
		c3 = -P900;
		c4 = P781;
		c5 = P974;
		c6 = P433;
	    } else if (mu==2) {
		c1 = -P222;
		c2 = -P900;
		c3 = P623;
		c4 = P974;
		c5 = -P433;
		c6 = -P781;
	    } else if (mu==3) {
		c1 = -P900;
		c2 = P623;
		c3 = -P222;
		c4 = P433;
		c5 = -P781;
		c6 = P974;
	    } else if (mu==4) {
		c1 = -P900;
		c2 = P623;
		c3 = -P222;
		c4 = -P433;
		c5 = P781;
		c6 = -P974;
	    } else if (mu==5) {
		c1 = -P222;
		c2 = -P900;
		c3 = P623;
		c4 = -P974;
		c5 = P433;
		c6 = P781;
	    } else {
		c1 = P623;
		c2 = -P222;
		c3 = -P900;
		c4 = -P781;
		c5 = -P974;
		c6 = -P433;
	    }
	    for (l=0; l<m; l++) {
		t1r = z[j1]+z[j6];
		t1i = z[j1+1]+z[j6+1];
		t2r = z[j2]+z[j5];
		t2i = z[j2+1]+z[j5+1];
		t3r = z[j3]+z[j4];
		t3i = z[j3+1]+z[j4+1];
		t4r = z[j1]-z[j6];
		t4i = z[j1+1]-z[j6+1];
		t5r = z[j2]-z[j5];
		t5i = z[j2+1]-z[j5+1];
		t6r = z[j3]-z[j4];
		t6i = z[j3+1]-z[j4+1];
		t7r = z[j0]-0.5*t3r;
		t7i = z[j0+1]-0.5*t3i;
		t8r = t1r-t3r;
		t8i = t1i-t3i;
		t9r = t2r-t3r;
		t9i = t2i-t3i;
		y1r = t7r+c1*t8r+c2*t9r;
		y1i = t7i+c1*t8i+c2*t9i;
		y2r = t7r+c2*t8r+c3*t9r;
		y2i = t7i+c2*t8i+c3*t9i;
		y3r = t7r+c3*t8r+c1*t9r;
		y3i = t7i+c3*t8i+c1*t9i;
		y4r = c6*t4r-c4*t5r+c5*t6r;
		y4i = c6*t4i-c4*t5i+c5*t6i;
		y5r = c5*t4r-c6*t5r-c4*t6r;
		y5i = c5*t4i-c6*t5i-c4*t6i;
		y6r = c4*t4r+c5*t5r+c6*t6r;
		y6i = c4*t4i+c5*t5i+c6*t6i;
		z[j0] = z[j0]+t1r+t2r+t3r;
		z[j0+1] = z[j0+1]+t1i+t2i+t3i;
		z[j1] = y1r-y6i;
		z[j1+1] = y1i+y6r;
		z[j2] = y2r-y5i;
		z[j2+1] = y2i+y5r;
		z[j3] = y3r-y4i;
		z[j3+1] = y3i+y4r;
		z[j4] = y3r+y4i;
		z[j4+1] = y3i-y4r;
		z[j5] = y2r+y5i;
		z[j5+1] = y2i-y5r;
		z[j6] = y1r+y6i;
		z[j6+1] = y1i-y6r;
		jt = j6+2;
		j6 = j5+2;
		j5 = j4+2;
		j4 = j3+2;
		j3 = j2+2;
		j2 = j1+2;
		j1 = j0+2;
		j0 = jt;
	    }
	    continue;
	}
	j7 = j6+jinc;
	if (j7>=jmax) j7 = j7-jmax;

	/* if factor is 8 */
	if (ifac==8) {
	    if (mu==1) {
		c1 = 1.0;
		c2 = P707;
	    } else if (mu==3) {
		c1 = -1.0;
		c2 = -P707;
	    } else if (mu==5) {
		c1 = 1.0;
		c2 = -P707;
	    } else {
		c1 = -1.0;
		c2 = P707;
	    }
	    c3 = c1*c2;
	    for (l=0; l<m; l++) {
		t1r = z[j0]+z[j4];
		t1i = z[j0+1]+z[j4+1];
		t2r = z[j0]-z[j4];
		t2i = z[j0+1]-z[j4+1];
		t3r = z[j1]+z[j5];
		t3i = z[j1+1]+z[j5+1];
		t4r = z[j1]-z[j5];
		t4i = z[j1+1]-z[j5+1];
		t5r = z[j2]+z[j6];
		t5i = z[j2+1]+z[j6+1];
		t6r = c1*(z[j2]-z[j6]);
		t6i = c1*(z[j2+1]-z[j6+1]);
		t7r = z[j3]+z[j7];
		t7i = z[j3+1]+z[j7+1];
		t8r = z[j3]-z[j7];
		t8i = z[j3+1]-z[j7+1];
		t9r = t1r+t5r;
		t9i = t1i+t5i;
		t10r = t3r+t7r;
		t10i = t3i+t7i;
		t11r = c2*(t4r-t8r);
		t11i = c2*(t4i-t8i);
		t12r = c3*(t4r+t8r);
		t12i = c3*(t4i+t8i);
		y1r = t2r+t11r;
		y1i = t2i+t11i;
		y2r = t1r-t5r;
		y2i = t1i-t5i;
		y3r = t2r-t11r;
		y3i = t2i-t11i;
		y5r = t12r-t6r;
		y5i = t12i-t6i;
		y6r = c1*(t3r-t7r);
		y6i = c1*(t3i-t7i);
		y7r = t12r+t6r;
		y7i = t12i+t6i;
		z[j0] = t9r+t10r;
		z[j0+1] = t9i+t10i;
		z[j1] = y1r-y7i;
		z[j1+1] = y1i+y7r;
		z[j2] = y2r-y6i;
		z[j2+1] = y2i+y6r;
		z[j3] = y3r-y5i;
		z[j3+1] = y3i+y5r;
		z[j4] = t9r-t10r;
		z[j4+1] = t9i-t10i;
		z[j5] = y3r+y5i;
		z[j5+1] = y3i-y5r;
		z[j6] = y2r+y6i;
		z[j6+1] = y2i-y6r;
		z[j7] = y1r+y7i;
		z[j7+1] = y1i-y7r;
		jt = j7+2;
		j7 = j6+2;
		j6 = j5+2;
		j5 = j4+2;
		j4 = j3+2;
		j3 = j2+2;
		j2 = j1+2;
		j1 = j0+2;
		j0 = jt;
	    }
	    continue;
	}
	j8 = j7+jinc;
	if (j8>=jmax) j8 = j8-jmax;

	/* if factor is 9 */
	if (ifac==9) {
	    if (mu==1) {
		c1 = P866;
		c2 = P766;
		c3 = P642;
		c4 = P173;
		c5 = P984;
	    } else if (mu==2) {
		c1 = -P866;
		c2 = P173;
		c3 = P984;
		c4 = -P939;
		c5 = P342;
	    } else if (mu==4) {
		c1 = P866;
		c2 = -P939;
		c3 = P342;
		c4 = P766;
		c5 = -P642;
	    } else if (mu==5) {
		c1 = -P866;
		c2 = -P939;
		c3 = -P342;
		c4 = P766;
		c5 = P642;
	    } else if (mu==7) {
		c1 = P866;
		c2 = P173;
		c3 = -P984;
		c4 = -P939;
		c5 = -P342;
	    } else {
		c1 = -P866;
		c2 = P766;
		c3 = -P642;
		c4 = P173;
		c5 = -P984;
	    }
	    c6 = c1*c2;
	    c7 = c1*c3;
	    c8 = c1*c4;
	    c9 = c1*c5;
	    for (l=0; l<m; l++) {
		t1r = z[j3]+z[j6];
		t1i = z[j3+1]+z[j6+1];
		t2r = z[j0]-0.5*t1r;
		t2i = z[j0+1]-0.5*t1i;
		t3r = c1*(z[j3]-z[j6]);
		t3i = c1*(z[j3+1]-z[j6+1]);
		t4r = z[j0]+t1r;
		t4i = z[j0+1]+t1i;
		t5r = z[j4]+z[j7];
		t5i = z[j4+1]+z[j7+1];
		t6r = z[j1]-0.5*t5r;
		t6i = z[j1+1]-0.5*t5i;
		t7r = z[j4]-z[j7];
		t7i = z[j4+1]-z[j7+1];
		t8r = z[j1]+t5r;
		t8i = z[j1+1]+t5i;
		t9r = z[j2]+z[j5];
		t9i = z[j2+1]+z[j5+1];
		t10r = z[j8]-0.5*t9r;
		t10i = z[j8+1]-0.5*t9i;
		t11r = z[j2]-z[j5];
		t11i = z[j2+1]-z[j5+1];
		t12r = z[j8]+t9r;
		t12i = z[j8+1]+t9i;
		t13r = t8r+t12r;
		t13i = t8i+t12i;
		t14r = t6r+t10r;
		t14i = t6i+t10i;
		t15r = t6r-t10r;
		t15i = t6i-t10i;
		t16r = t7r+t11r;
		t16i = t7i+t11i;
		t17r = t7r-t11r;
		t17i = t7i-t11i;
		t18r = c2*t14r-c7*t17r;
		t18i = c2*t14i-c7*t17i;
		t19r = c4*t14r+c9*t17r;
		t19i = c4*t14i+c9*t17i;
		t20r = c3*t15r+c6*t16r;
		t20i = c3*t15i+c6*t16i;
		t21r = c5*t15r-c8*t16r;
		t21i = c5*t15i-c8*t16i;
		t22r = t18r+t19r;
		t22i = t18i+t19i;
		t23r = t20r-t21r;
		t23i = t20i-t21i;
		y1r = t2r+t18r;
		y1i = t2i+t18i;
		y2r = t2r+t19r;
		y2i = t2i+t19i;
		y3r = t4r-0.5*t13r;
		y3i = t4i-0.5*t13i;
		y4r = t2r-t22r;
		y4i = t2i-t22i;
		y5r = t3r-t23r;
		y5i = t3i-t23i;
		y6r = c1*(t8r-t12r);
		y6i = c1*(t8i-t12i);
		y7r = t21r-t3r;
		y7i = t21i-t3i;
		y8r = t3r+t20r;
		y8i = t3i+t20i;
		z[j0] = t4r+t13r;
		z[j0+1] = t4i+t13i;
		z[j1] = y1r-y8i;
		z[j1+1] = y1i+y8r;
		z[j2] = y2r-y7i;
		z[j2+1] = y2i+y7r;
		z[j3] = y3r-y6i;
		z[j3+1] = y3i+y6r;
		z[j4] = y4r-y5i;
		z[j4+1] = y4i+y5r;
		z[j5] = y4r+y5i;
		z[j5+1] = y4i-y5r;
		z[j6] = y3r+y6i;
		z[j6+1] = y3i-y6r;
		z[j7] = y2r+y7i;
		z[j7+1] = y2i-y7r;
		z[j8] = y1r+y8i;
		z[j8+1] = y1i-y8r;
		jt = j8+2;
		j8 = j7+2;
		j7 = j6+2;
		j6 = j5+2;
		j5 = j4+2;
		j4 = j3+2;
		j3 = j2+2;
		j2 = j1+2;
		j1 = j0+2;
		j0 = jt;
	    }
	    continue;
	}
	j9 = j8+jinc;
	if (j9>=jmax) j9 = j9-jmax;
	j10 = j9+jinc;
	if (j10>=jmax) j10 = j10-jmax;

	/* if factor is 11 */
	if (ifac==11) {
	    if (mu==1) {
		c1 = P841;
		c2 = P415;
		c3 = -P142;
		c4 = -P654;
		c5 = -P959;
		c6 = P540;
		c7 = P909;
		c8 = P989;
		c9 = P755;
		c10 = P281;
	    } else if (mu==2) {
		c1 = P415;
		c2 = -P654;
		c3 = -P959;
		c4 = -P142;
		c5 = P841;
		c6 = P909;
		c7 = P755;
		c8 = -P281;
		c9 = -P989;
		c10 = -P540;
	    } else if (mu==3) {
		c1 = -P142;
		c2 = -P959;
		c3 = P415;
		c4 = P841;
		c5 = -P654;
		c6 = P989;
		c7 = -P281;
		c8 = -P909;
		c9 = P540;
		c10 = P755;
	    } else if (mu==4) {
		c1 = -P654;
		c2 = -P142;
		c3 = P841;
		c4 = -P959;
		c5 = P415;
		c6 = P755;
		c7 = -P989;
		c8 = P540;
		c9 = P281;
		c10 = -P909;
	    } else if (mu==5) {
		c1 = -P959;
		c2 = P841;
		c3 = -P654;
		c4 = P415;
		c5 = -P142;
		c6 = P281;
		c7 = -P540;
		c8 = P755;
		c9 = -P909;
		c10 = P989;
	    } else if (mu==6) {
		c1 = -P959;
		c2 = P841;
		c3 = -P654;
		c4 = P415;
		c5 = -P142;
		c6 = -P281;
		c7 = P540;
		c8 = -P755;
		c9 = P909;
		c10 = -P989;
	    } else if (mu==7) {
		c1 = -P654;
		c2 = -P142;
		c3 = P841;
		c4 = -P959;
		c5 = P415;
		c6 = -P755;
		c7 = P989;
		c8 = -P540;
		c9 = -P281;
		c10 = P909;
	    } else if (mu==8) {
		c1 = -P142;
		c2 = -P959;
		c3 = P415;
		c4 = P841;
		c5 = -P654;
		c6 = -P989;
		c7 = P281;
		c8 = P909;
		c9 = -P540;
		c10 = -P755;
	    } else if (mu==9) {
		c1 = P415;
		c2 = -P654;
		c3 = -P959;
		c4 = -P142;
		c5 = P841;
		c6 = -P909;
		c7 = -P755;
		c8 = P281;
		c9 = P989;
		c10 = P540;
	    } else {
		c1 = P841;
		c2 = P415;
		c3 = -P142;
		c4 = -P654;
		c5 = -P959;
		c6 = -P540;
		c7 = -P909;
		c8 = -P989;
		c9 = -P755;
		c10 = -P281;
	    }
	    for (l=0; l<m; l++) {
		t1r = z[j1]+z[j10];
		t1i = z[j1+1]+z[j10+1];
		t2r = z[j2]+z[j9];
		t2i = z[j2+1]+z[j9+1];
		t3r = z[j3]+z[j8];
		t3i = z[j3+1]+z[j8+1];
		t4r = z[j4]+z[j7];
		t4i = z[j4+1]+z[j7+1];
		t5r = z[j5]+z[j6];
		t5i = z[j5+1]+z[j6+1];
		t6r = z[j1]-z[j10];
		t6i = z[j1+1]-z[j10+1];
		t7r = z[j2]-z[j9];
		t7i = z[j2+1]-z[j9+1];
		t8r = z[j3]-z[j8];
		t8i = z[j3+1]-z[j8+1];
		t9r = z[j4]-z[j7];
		t9i = z[j4+1]-z[j7+1];
		t10r = z[j5]-z[j6];
		t10i = z[j5+1]-z[j6+1];
		t11r = z[j0]-0.5*t5r;
		t11i = z[j0+1]-0.5*t5i;
		t12r = t1r-t5r;
		t12i = t1i-t5i;
		t13r = t2r-t5r;
		t13i = t2i-t5i;
		t14r = t3r-t5r;
		t14i = t3i-t5i;
		t15r = t4r-t5r;
		t15i = t4i-t5i;
		y1r = t11r+c1*t12r+c2*t13r+c3*t14r+c4*t15r;
		y1i = t11i+c1*t12i+c2*t13i+c3*t14i+c4*t15i;
		y2r = t11r+c2*t12r+c4*t13r+c5*t14r+c3*t15r;
		y2i = t11i+c2*t12i+c4*t13i+c5*t14i+c3*t15i;
		y3r = t11r+c3*t12r+c5*t13r+c2*t14r+c1*t15r;
		y3i = t11i+c3*t12i+c5*t13i+c2*t14i+c1*t15i;
		y4r = t11r+c4*t12r+c3*t13r+c1*t14r+c5*t15r;
		y4i = t11i+c4*t12i+c3*t13i+c1*t14i+c5*t15i;
		y5r = t11r+c5*t12r+c1*t13r+c4*t14r+c2*t15r;
		y5i = t11i+c5*t12i+c1*t13i+c4*t14i+c2*t15i;
		y6r = c10*t6r-c6*t7r+c9*t8r-c7*t9r+c8*t10r;
		y6i = c10*t6i-c6*t7i+c9*t8i-c7*t9i+c8*t10i;
		y7r = c9*t6r-c8*t7r+c6*t8r+c10*t9r-c7*t10r;
		y7i = c9*t6i-c8*t7i+c6*t8i+c10*t9i-c7*t10i;
		y8r = c8*t6r-c10*t7r-c7*t8r+c6*t9r+c9*t10r;
		y8i = c8*t6i-c10*t7i-c7*t8i+c6*t9i+c9*t10i;
		y9r = c7*t6r+c9*t7r-c10*t8r-c8*t9r-c6*t10r;
		y9i = c7*t6i+c9*t7i-c10*t8i-c8*t9i-c6*t10i;
		y10r = c6*t6r+c7*t7r+c8*t8r+c9*t9r+c10*t10r;
		y10i = c6*t6i+c7*t7i+c8*t8i+c9*t9i+c10*t10i;
		z[j0] = z[j0]+t1r+t2r+t3r+t4r+t5r;
		z[j0+1] = z[j0+1]+t1i+t2i+t3i+t4i+t5i;
		z[j1] = y1r-y10i;
		z[j1+1] = y1i+y10r;
		z[j2] = y2r-y9i;
		z[j2+1] = y2i+y9r;
		z[j3] = y3r-y8i;
		z[j3+1] = y3i+y8r;
		z[j4] = y4r-y7i;
		z[j4+1] = y4i+y7r;
		z[j5] = y5r-y6i;
		z[j5+1] = y5i+y6r;
		z[j6] = y5r+y6i;
		z[j6+1] = y5i-y6r;
		z[j7] = y4r+y7i;
		z[j7+1] = y4i-y7r;
		z[j8] = y3r+y8i;
		z[j8+1] = y3i-y8r;
		z[j9] = y2r+y9i;
		z[j9+1] = y2i-y9r;
		z[j10] = y1r+y10i;
		z[j10+1] = y1i-y10r;
		jt = j10+2;
		j10 = j9+2;
		j9 = j8+2;
		j8 = j7+2;
		j7 = j6+2;
		j6 = j5+2;
		j5 = j4+2;
		j4 = j3+2;
		j3 = j2+2;
		j2 = j1+2;
		j1 = j0+2;
		j0 = jt;
	    }
	    continue;
	}
	j11 = j10+jinc;
	if (j11>=jmax) j11 = j11-jmax;
	j12 = j11+jinc;
	if (j12>=jmax) j12 = j12-jmax;

	/* if factor is 13 */
	if (ifac==13) {
	    if (mu==1) {
		c1 = P885;
		c2 = P568;
		c3 = P120;
		c4 = -P354;
		c5 = -P748;
		c6 = -P970;
		c7 = P464;
		c8 = P822;
		c9 = P992;
		c10 = P935;
		c11 = P663;
		c12 = P239;
	    } else if (mu==2) {
		c1 = P568;
		c2 = -P354;
		c3 = -P970;
		c4 = -P748;
		c5 = P120;
		c6 = P885;
		c7 = P822;
		c8 = P935;
		c9 = P239;
		c10 = -P663;
		c11 = -P992;
		c12 = -P464;
	    } else if (mu==3) {
		c1 = P120;
		c2 = -P970;
		c3 = -P354;
		c4 = P885;
		c5 = P568;
		c6 = -P748;
		c7 = P992;
		c8 = P239;
		c9 = -P935;
		c10 = -P464;
		c11 = P822;
		c12 = P663;
	    } else if (mu==4) {
		c1 = -P354;
		c2 = -P748;
		c3 = P885;
		c4 = P120;
		c5 = -P970;
		c6 = P568;
		c7 = P935;
		c8 = -P663;
		c9 = -P464;
		c10 = P992;
		c11 = -P239;
		c12 = -P822;
	    } else if (mu==5) {
		c1 = -P748;
		c2 = P120;
		c3 = P568;
		c4 = -P970;
		c5 = P885;
		c6 = -P354;
		c7 = P663;
		c8 = -P992;
		c9 = P822;
		c10 = -P239;
		c11 = -P464;
		c12 = P935;
	    } else if (mu==6) {
		c1 = -P970;
		c2 = P885;
		c3 = -P748;
		c4 = P568;
		c5 = -P354;
		c6 = P120;
		c7 = P239;
		c8 = -P464;
		c9 = P663;
		c10 = -P822;
		c11 = P935;
		c12 = -P992;
	    } else if (mu==7) {
		c1 = -P970;
		c2 = P885;
		c3 = -P748;
		c4 = P568;
		c5 = -P354;
		c6 = P120;
		c7 = -P239;
		c8 = P464;
		c9 = -P663;
		c10 = P822;
		c11 = -P935;
		c12 = P992;
	    } else if (mu==8) {
		c1 = -P748;
		c2 = P120;
		c3 = P568;
		c4 = -P970;
		c5 = P885;
		c6 = -P354;
		c7 = -P663;
		c8 = P992;
		c9 = -P822;
		c10 = P239;
		c11 = P464;
		c12 = -P935;
	    } else if (mu==9) {
		c1 = -P354;
		c2 = -P748;
		c3 = P885;
		c4 = P120;
		c5 = -P970;
		c6 = P568;
		c7 = -P935;
		c8 = P663;
		c9 = P464;
		c10 = -P992;
		c11 = P239;
		c12 = P822;
	    } else if (mu==10) {
		c1 = P120;
		c2 = -P970;
		c3 = -P354;
		c4 = P885;
		c5 = P568;
		c6 = -P748;
		c7 = -P992;
		c8 = -P239;
		c9 = P935;
		c10 = P464;
		c11 = -P822;
		c12 = -P663;
	    } else if (mu==11) {
		c1 = P568;
		c2 = -P354;
		c3 = -P970;
		c4 = -P748;
		c5 = P120;
		c6 = P885;
		c7 = -P822;
		c8 = -P935;
		c9 = -P239;
		c10 = P663;
		c11 = P992;
		c12 = P464;
	    } else {
		c1 = P885;
		c2 = P568;
		c3 = P120;
		c4 = -P354;
		c5 = -P748;
		c6 = -P970;
		c7 = -P464;
		c8 = -P822;
		c9 = -P992;
		c10 = -P935;
		c11 = -P663;
		c12 = -P239;
	    }
	    for (l=0; l<m; l++) {
		t1r = z[j1]+z[j12];
		t1i = z[j1+1]+z[j12+1];
		t2r = z[j2]+z[j11];
		t2i = z[j2+1]+z[j11+1];
		t3r = z[j3]+z[j10];
		t3i = z[j3+1]+z[j10+1];
		t4r = z[j4]+z[j9];
		t4i = z[j4+1]+z[j9+1];
		t5r = z[j5]+z[j8];
		t5i = z[j5+1]+z[j8+1];
		t6r = z[j6]+z[j7];
		t6i = z[j6+1]+z[j7+1];
		t7r = z[j1]-z[j12];
		t7i = z[j1+1]-z[j12+1];
		t8r = z[j2]-z[j11];
		t8i = z[j2+1]-z[j11+1];
		t9r = z[j3]-z[j10];
		t9i = z[j3+1]-z[j10+1];
		t10r = z[j4]-z[j9];
		t10i = z[j4+1]-z[j9+1];
		t11r = z[j5]-z[j8];
		t11i = z[j5+1]-z[j8+1];
		t12r = z[j6]-z[j7];
		t12i = z[j6+1]-z[j7+1];
		t13r = z[j0]-0.5*t6r;
		t13i = z[j0+1]-0.5*t6i;
		t14r = t1r-t6r;
		t14i = t1i-t6i;
		t15r = t2r-t6r;
		t15i = t2i-t6i;
		t16r = t3r-t6r;
		t16i = t3i-t6i;
		t17r = t4r-t6r;
		t17i = t4i-t6i;
		t18r = t5r-t6r;
		t18i = t5i-t6i;
		y1r = t13r+c1*t14r+c2*t15r+c3*t16r+c4*t17r+c5*t18r;
		y1i = t13i+c1*t14i+c2*t15i+c3*t16i+c4*t17i+c5*t18i;
		y2r = t13r+c2*t14r+c4*t15r+c6*t16r+c5*t17r+c3*t18r;
		y2i = t13i+c2*t14i+c4*t15i+c6*t16i+c5*t17i+c3*t18i;
		y3r = t13r+c3*t14r+c6*t15r+c4*t16r+c1*t17r+c2*t18r;
		y3i = t13i+c3*t14i+c6*t15i+c4*t16i+c1*t17i+c2*t18i;
		y4r = t13r+c4*t14r+c5*t15r+c1*t16r+c3*t17r+c6*t18r;
		y4i = t13i+c4*t14i+c5*t15i+c1*t16i+c3*t17i+c6*t18i;
		y5r = t13r+c5*t14r+c3*t15r+c2*t16r+c6*t17r+c1*t18r;
		y5i = t13i+c5*t14i+c3*t15i+c2*t16i+c6*t17i+c1*t18i;
		y6r = t13r+c6*t14r+c1*t15r+c5*t16r+c2*t17r+c4*t18r;
		y6i = t13i+c6*t14i+c1*t15i+c5*t16i+c2*t17i+c4*t18i;
		y7r = c12*t7r-c7*t8r+c11*t9r-c8*t10r+c10*t11r-c9*t12r;
		y7i = c12*t7i-c7*t8i+c11*t9i-c8*t10i+c10*t11i-c9*t12i;
		y8r = c11*t7r-c9*t8r+c8*t9r-c12*t10r-c7*t11r+c10*t12r;
		y8i = c11*t7i-c9*t8i+c8*t9i-c12*t10i-c7*t11i+c10*t12i;
		y9r = c10*t7r-c11*t8r-c7*t9r+c9*t10r-c12*t11r-c8*t12r;
		y9i = c10*t7i-c11*t8i-c7*t9i+c9*t10i-c12*t11i-c8*t12i;
		y10r = c9*t7r+c12*t8r-c10*t9r-c7*t10r+c8*t11r+c11*t12r;
		y10i = c9*t7i+c12*t8i-c10*t9i-c7*t10i+c8*t11i+c11*t12i;
		y11r = c8*t7r+c10*t8r+c12*t9r-c11*t10r-c9*t11r-c7*t12r;
		y11i = c8*t7i+c10*t8i+c12*t9i-c11*t10i-c9*t11i-c7*t12i;
		y12r = c7*t7r+c8*t8r+c9*t9r+c10*t10r+c11*t11r+c12*t12r;
		y12i = c7*t7i+c8*t8i+c9*t9i+c10*t10i+c11*t11i+c12*t12i;
		z[j0] = z[j0]+t1r+t2r+t3r+t4r+t5r+t6r;
		z[j0+1] = z[j0+1]+t1i+t2i+t3i+t4i+t5i+t6i;
		z[j1] = y1r-y12i;
		z[j1+1] = y1i+y12r;
		z[j2] = y2r-y11i;
		z[j2+1] = y2i+y11r;
		z[j3] = y3r-y10i;
		z[j3+1] = y3i+y10r;
		z[j4] = y4r-y9i;
		z[j4+1] = y4i+y9r;
		z[j5] = y5r-y8i;
		z[j5+1] = y5i+y8r;
		z[j6] = y6r-y7i;
		z[j6+1] = y6i+y7r;
		z[j7] = y6r+y7i;
		z[j7+1] = y6i-y7r;
		z[j8] = y5r+y8i;
		z[j8+1] = y5i-y8r;
		z[j9] = y4r+y9i;
		z[j9+1] = y4i-y9r;
		z[j10] = y3r+y10i;
		z[j10+1] = y3i-y10r;
		z[j11] = y2r+y11i;
		z[j11+1] = y2i-y11r;
		z[j12] = y1r+y12i;
		z[j12+1] = y1i-y12r;
		jt = j12+2;
		j12 = j11+2;
		j11 = j10+2;
		j10 = j9+2;
		j9 = j8+2;
		j8 = j7+2;
		j7 = j6+2;
		j6 = j5+2;
		j5 = j4+2;
		j4 = j3+2;
		j3 = j2+2;
		j2 = j1+2;
		j1 = j0+2;
		j0 = jt;
	    }
	    continue;
	}
	j13 = j12+jinc;
	if (j13>=jmax) j13 = j13-jmax;
	j14 = j13+jinc;
	if (j14>=jmax) j14 = j14-jmax;
	j15 = j14+jinc;
	if (j15>=jmax) j15 = j15-jmax;

	/* if factor is 16 */
	if (ifac==16) {
	    if (mu==1) {
		c1 = 1.0;
		c2 = P923;
		c3 = P382;
		c4 = P707;
	    } else if (mu==3) {
		c1 = -1.0;
		c2 = P382;
		c3 = P923;
		c4 = -P707;
	    } else if (mu==5) {
		c1 = 1.0;
		c2 = -P382;
		c3 = P923;
		c4 = -P707;
	    } else if (mu==7) {
		c1 = -1.0;
		c2 = -P923;
		c3 = P382;
		c4 = P707;
	    } else if (mu==9) {
		c1 = 1.0;
		c2 = -P923;
		c3 = -P382;
		c4 = P707;
	    } else if (mu==11) {
		c1 = -1.0;
		c2 = -P382;
		c3 = -P923;
		c4 = -P707;
	    } else if (mu==13) {
		c1 = 1.0;
		c2 = P382;
		c3 = -P923;
		c4 = -P707;
	    } else {
		c1 = -1.0;
		c2 = P923;
		c3 = -P382;
		c4 = P707;
	    }
	    c5 = c1*c4;
	    c6 = c1*c3;
	    c7 = c1*c2;
	    for (l=0; l<m; l++) {
		t1r = z[j0]+z[j8];
		t1i = z[j0+1]+z[j8+1];
		t2r = z[j4]+z[j12];
		t2i = z[j4+1]+z[j12+1];
		t3r = z[j0]-z[j8];
		t3i = z[j0+1]-z[j8+1];
		t4r = c1*(z[j4]-z[j12]);
		t4i = c1*(z[j4+1]-z[j12+1]);
		t5r = t1r+t2r;
		t5i = t1i+t2i;
		t6r = t1r-t2r;
		t6i = t1i-t2i;
		t7r = z[j1]+z[j9];
		t7i = z[j1+1]+z[j9+1];
		t8r = z[j5]+z[j13];
		t8i = z[j5+1]+z[j13+1];
		t9r = z[j1]-z[j9];
		t9i = z[j1+1]-z[j9+1];
		t10r = z[j5]-z[j13];
		t10i = z[j5+1]-z[j13+1];
		t11r = t7r+t8r;
		t11i = t7i+t8i;
		t12r = t7r-t8r;
		t12i = t7i-t8i;
		t13r = z[j2]+z[j10];
		t13i = z[j2+1]+z[j10+1];
		t14r = z[j6]+z[j14];
		t14i = z[j6+1]+z[j14+1];
		t15r = z[j2]-z[j10];
		t15i = z[j2+1]-z[j10+1];
		t16r = z[j6]-z[j14];
		t16i = z[j6+1]-z[j14+1];
		t17r = t13r+t14r;
		t17i = t13i+t14i;
		t18r = c4*(t15r-t16r);
		t18i = c4*(t15i-t16i);
		t19r = c5*(t15r+t16r);
		t19i = c5*(t15i+t16i);
		t20r = c1*(t13r-t14r);
		t20i = c1*(t13i-t14i);
		t21r = z[j3]+z[j11];
		t21i = z[j3+1]+z[j11+1];
		t22r = z[j7]+z[j15];
		t22i = z[j7+1]+z[j15+1];
		t23r = z[j3]-z[j11];
		t23i = z[j3+1]-z[j11+1];
		t24r = z[j7]-z[j15];
		t24i = z[j7+1]-z[j15+1];
		t25r = t21r+t22r;
		t25i = t21i+t22i;
		t26r = t21r-t22r;
		t26i = t21i-t22i;
		t27r = t9r+t24r;
		t27i = t9i+t24i;
		t28r = t10r+t23r;
		t28i = t10i+t23i;
		t29r = t9r-t24r;
		t29i = t9i-t24i;
		t30r = t10r-t23r;
		t30i = t10i-t23i;
		t31r = t5r+t17r;
		t31i = t5i+t17i;
		t32r = t11r+t25r;
		t32i = t11i+t25i;
		t33r = t3r+t18r;
		t33i = t3i+t18i;
		t34r = c2*t29r-c6*t30r;
		t34i = c2*t29i-c6*t30i;
		t35r = t3r-t18r;
		t35i = t3i-t18i;
		t36r = c7*t27r-c3*t28r;
		t36i = c7*t27i-c3*t28i;
		t37r = t4r+t19r;
		t37i = t4i+t19i;
		t38r = c3*t27r+c7*t28r;
		t38i = c3*t27i+c7*t28i;
		t39r = t4r-t19r;
		t39i = t4i-t19i;
		t40r = c6*t29r+c2*t30r;
		t40i = c6*t29i+c2*t30i;
		t41r = c4*(t12r-t26r);
		t41i = c4*(t12i-t26i);
		t42r = c5*(t12r+t26r);
		t42i = c5*(t12i+t26i);
		y1r = t33r+t34r;
		y1i = t33i+t34i;
		y2r = t6r+t41r;
		y2i = t6i+t41i;
		y3r = t35r+t40r;
		y3i = t35i+t40i;
		y4r = t5r-t17r;
		y4i = t5i-t17i;
		y5r = t35r-t40r;
		y5i = t35i-t40i;
		y6r = t6r-t41r;
		y6i = t6i-t41i;
		y7r = t33r-t34r;
		y7i = t33i-t34i;
		y9r = t38r-t37r;
		y9i = t38i-t37i;
		y10r = t42r-t20r;
		y10i = t42i-t20i;
		y11r = t36r+t39r;
		y11i = t36i+t39i;
		y12r = c1*(t11r-t25r);
		y12i = c1*(t11i-t25i);
		y13r = t36r-t39r;
		y13i = t36i-t39i;
		y14r = t42r+t20r;
		y14i = t42i+t20i;
		y15r = t38r+t37r;
		y15i = t38i+t37i;
		z[j0] = t31r+t32r;
		z[j0+1] = t31i+t32i;
		z[j1] = y1r-y15i;
		z[j1+1] = y1i+y15r;
		z[j2] = y2r-y14i;
		z[j2+1] = y2i+y14r;
		z[j3] = y3r-y13i;
		z[j3+1] = y3i+y13r;
		z[j4] = y4r-y12i;
		z[j4+1] = y4i+y12r;
		z[j5] = y5r-y11i;
		z[j5+1] = y5i+y11r;
		z[j6] = y6r-y10i;
		z[j6+1] = y6i+y10r;
		z[j7] = y7r-y9i;
		z[j7+1] = y7i+y9r;
		z[j8] = t31r-t32r;
		z[j8+1] = t31i-t32i;
		z[j9] = y7r+y9i;
		z[j9+1] = y7i-y9r;
		z[j10] = y6r+y10i;
		z[j10+1] = y6i-y10r;
		z[j11] = y5r+y11i;
		z[j11+1] = y5i-y11r;
		z[j12] = y4r+y12i;
		z[j12+1] = y4i-y12r;
		z[j13] = y3r+y13i;
		z[j13+1] = y3i-y13r;
		z[j14] = y2r+y14i;
		z[j14+1] = y2i-y14r;
		z[j15] = y1r+y15i;
		z[j15+1] = y1i-y15r;
		jt = j15+2;
		j15 = j14+2;
		j14 = j13+2;
		j13 = j12+2;
		j12 = j11+2;
		j11 = j10+2;
		j10 = j9+2;
		j9 = j8+2;
		j8 = j7+2;
		j7 = j6+2;
		j6 = j5+2;
		j5 = j4+2;
		j4 = j3+2;
		j3 = j2+2;
		j2 = j1+2;
		j1 = j0+2;
		j0 = jt;
	    }
	    continue;
	}
    }
}

void pfacr (int isign, int n, complex cz[], REAL rz[])
{
    int i,ir,ii,jr,ji,no2;
    REAL *z,tempr,tempi,sumr,sumi,difr,difi;
    REAL wr,wi,wpr,wpi,wtemp,theta;

    /* copy input to output and fix dc and nyquist */
    z = (REAL*)cz;
    for (i=2; i<n; i++)
    rz[i] = z[i];
    rz[1] = z[0]-z[n];
    rz[0] = z[0]+z[n];
    z = rz;

    /* initialize cosine-sine recurrence */
    theta = 2.0*M_PI/(REAL)n;
    if (isign>0) theta = -theta;
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0+wpr;
    wi = wpi;

    /* twiddle */
    no2 = n/2;
    for (ir=2,ii=3,jr=n-2,ji=n-1; ir<=no2; ir+=2,ii+=2,jr-=2,ji-=2) {
        sumr = z[ir]+z[jr];
        sumi = z[ii]+z[ji];
        difr = z[ir]-z[jr];
        difi = z[ii]-z[ji];
        tempr = wi*difr-wr*sumi;
        tempi = wi*sumi+wr*difr;
        z[ir] = sumr+tempr;
        z[ii] = difi+tempi;
        z[jr] = sumr-tempr;
        z[ji] = tempi-difi;
        wtemp = wr;
        wr += wr*wpr-wi*wpi;
        wi += wi*wpr+wtemp*wpi;
    }

    /* do complex to complex transform */
    pfacc(isign,n/2,(complex*)z);
}



void 
pfarc (int isign, int n, REAL rz[], complex cz[])
{
    int i,ir,ii,jr,ji,no2;
    REAL *z,tempr,tempi,sumr,sumi,difr,difi;
    REAL wr,wi,wpr,wpi,wtemp,theta;

    /* copy input to output while scaling */
    z = (REAL*)cz;
    for (i=0; i<n; i++)
    z[i] = 0.5*rz[i];

    /* do complex to complex transform */
    pfacc(isign,n/2,cz);

    /* fix dc and nyquist */
    z[n] = 2.0*(z[0]-z[1]);
    z[0] = 2.0*(z[0]+z[1]);
    z[n+1] = 0.0;
    z[1] = 0.0;

    /* initialize cosine-sine recurrence */
    theta = 2.0*M_PI/(REAL)n;
    if (isign<0) theta = -theta;
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0+wpr;
    wi = wpi;

    /* twiddle */
    no2 = n/2;
    for (ir=2,ii=3,jr=n-2,ji=n-1; ir<=no2; ir+=2,ii+=2,jr-=2,ji-=2) {
        sumr = z[ir]+z[jr];
        sumi = z[ii]+z[ji];
        difr = z[ir]-z[jr];
        difi = z[ii]-z[ji];
        tempr = wi*difr+wr*sumi;
        tempi = wi*sumi-wr*difr;
        z[ir] = sumr+tempr;
        z[ii] = difi+tempi;
        z[jr] = sumr-tempr;
        z[ji] = tempi-difi;
        wtemp = wr;
        wr += wr*wpr-wi*wpi;
        wi += wi*wpr+wtemp*wpi;
    }
}

void 
pfamcc (int isign, int n, int nt, int k, int kt, complex cz[])
{
    static int kfax[] = { 16,13,11,9,8,7,5,4,3,2 };
    REAL *z=(REAL*)cz;
    int j0,j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15;
    int nleft,jfax,ifac,jfac,iinc,imax,ndiv,m,mm=0,mu=0,l,istep,jstep,
    jt,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,it;
    REAL t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,
    t6r,t6i,t7r,t7i,t8r,t8i,t9r,t9i,t10r,t10i,
    t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,
    t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t20r,t20i,
    t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,
    t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t30r,t30i,
    t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,
    t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t40r,t40i,
    t41r,t41i,t42r,t42i,
    y1r,y1i,y2r,y2i,y3r,y3i,y4r,y4i,y5r,y5i,
    y6r,y6i,y7r,y7i,y8r,y8i,y9r,y9i,y10r,y10i,
    y11r,y11i,y12r,y12i,y13r,y13i,y14r,y14i,y15r,y15i,
    c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;
    
    /* determine step within and between transforms */
    istep = 2*k;
    jstep = 2*kt;
    
    /* keep track of n left after dividing by factors */
    nleft = n;
    
    /* begin loop over possible factors (from biggest to smallest) */
    for (jfax=0; jfax<NFAX; jfax++) {
	
        /* skip if not a mutually prime factor of n */
        ifac = kfax[jfax];
        ndiv = nleft/ifac;
        if (ndiv*ifac!=nleft) continue;
	
        /* update n left and determine n divided by factor */
        nleft = ndiv;
        m = n/ifac;
	
        /* determine rotation factor mu and stride mm */
        for (jfac=1; jfac<=ifac; jfac++) {
            mu = jfac;
            mm = jfac*m;
            if (mm%ifac==1) break;
        }
	
        /* adjust rotation factor for sign of transform */
        if (isign<0) mu = ifac-mu;
	
        /* compute stride, limit, and pointers */
        iinc = istep*mm;
        imax = istep*n;
        i0 = 0;
        i1 = i0+iinc;
	
        /* if factor is 2 */
        if (ifac==2) {
            for (l=0; l<m; l++) {
                j0 = i0;
                j1 = i1;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j0]-z[j1];
                    t1i = z[j0+1]-z[j1+1];
                    z[j0] = z[j0]+z[j1];
                    z[j0+1] = z[j0+1]+z[j1+1];
                    z[j1] = t1r;
                    z[j1+1] = t1i;
                    j0 += jstep;
                    j1 += jstep;
                }
                it = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i2 = i1+iinc;
        if (i2>=imax) i2 = i2-imax;
	
        /* if factor is 3 */
        if (ifac==3) {
            if (mu==1)
	    c1 = P866;
            else
	    c1 = -P866;
            for (l=0; l<m; l++) {
                j0 = i0;
                j1 = i1;
                j2 = i2;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j1]+z[j2];
                    t1i = z[j1+1]+z[j2+1];
                    y1r = z[j0]-0.5*t1r;
                    y1i = z[j0+1]-0.5*t1i;
                    y2r = c1*(z[j1]-z[j2]);
                    y2i = c1*(z[j1+1]-z[j2+1]);
                    z[j0] = z[j0]+t1r;
                    z[j0+1] = z[j0+1]+t1i;
                    z[j1] = y1r-y2i;
                    z[j1+1] = y1i+y2r;
                    z[j2] = y1r+y2i;
                    z[j2+1] = y1i-y2r;
                    j0 += jstep;
                    j1 += jstep;
                    j2 += jstep;
                }
                it = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i3 = i2+iinc;
        if (i3>=imax) i3 = i3-imax;
	
        /* if factor is 4 */
        if (ifac==4) {
            if (mu==1)
	    c1 = 1.0;
            else
	    c1 = -1.0;
            for (l=0; l<m; l++) {
                j0 = i0;
                j1 = i1;
                j2 = i2;
                j3 = i3;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j0]+z[j2];
                    t1i = z[j0+1]+z[j2+1];
                    t2r = z[j1]+z[j3];
                    t2i = z[j1+1]+z[j3+1];
                    y1r = z[j0]-z[j2];
                    y1i = z[j0+1]-z[j2+1];
                    y3r = c1*(z[j1]-z[j3]);
                    y3i = c1*(z[j1+1]-z[j3+1]);
                    z[j0] = t1r+t2r;
                    z[j0+1] = t1i+t2i;
                    z[j1] = y1r-y3i;
                    z[j1+1] = y1i+y3r;
                    z[j2] = t1r-t2r;
                    z[j2+1] = t1i-t2i;
                    z[j3] = y1r+y3i;
                    z[j3+1] = y1i-y3r;
                    j0 += jstep;
                    j1 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                }
                it = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i4 = i3+iinc;
        if (i4>=imax) i4 = i4-imax;
	
        /* if factor is 5 */
        if (ifac==5) {
            if (mu==1) {
                c1 = P559;
                c2 = P951;
                c3 = P587;
            } else if (mu==2) {
                c1 = -P559;
                c2 = P587;
                c3 = -P951;
            } else if (mu==3) {
                c1 = -P559;
                c2 = -P587;
                c3 = P951;
            } else { 
                c1 = P559;
                c2 = -P951;
                c3 = -P587;
            }
            for (l=0; l<m; l++) {
                j0 = i0;
                j1 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j1]+z[j4];
                    t1i = z[j1+1]+z[j4+1];
                    t2r = z[j2]+z[j3];
                    t2i = z[j2+1]+z[j3+1];
                    t3r = z[j1]-z[j4];
                    t3i = z[j1+1]-z[j4+1];
                    t4r = z[j2]-z[j3];
                    t4i = z[j2+1]-z[j3+1];
                    t5r = t1r+t2r;
                    t5i = t1i+t2i;
                    t6r = c1*(t1r-t2r);
                    t6i = c1*(t1i-t2i);
                    t7r = z[j0]-0.25*t5r;
                    t7i = z[j0+1]-0.25*t5i;
                    y1r = t7r+t6r;
                    y1i = t7i+t6i;
                    y2r = t7r-t6r;
                    y2i = t7i-t6i;
                    y3r = c3*t3r-c2*t4r;
                    y3i = c3*t3i-c2*t4i;
                    y4r = c2*t3r+c3*t4r;
                    y4i = c2*t3i+c3*t4i;
                    z[j0] = z[j0]+t5r;
                    z[j0+1] = z[j0+1]+t5i;
                    z[j1] = y1r-y4i;
                    z[j1+1] = y1i+y4r;
                    z[j2] = y2r-y3i;
                    z[j2+1] = y2i+y3r;
                    z[j3] = y2r+y3i;
                    z[j3+1] = y2i-y3r;
                    z[j4] = y1r+y4i;
                    z[j4+1] = y1i-y4r;
                    j0 += jstep;
                    j1 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                }
                it = i4+istep;
                i4 = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i5 = i4+iinc;
        if (i5>=imax) i5 = i5-imax;
        i6 = i5+iinc;
        if (i6>=imax) i6 = i6-imax;
	
        /* if factor is 7 */
        if (ifac==7) {
            if (mu==1) {
                c1 = P623;
                c2 = -P222;
                c3 = -P900;
                c4 = P781;
                c5 = P974;
                c6 = P433;
            } else if (mu==2) {
                c1 = -P222;
                c2 = -P900;
                c3 = P623;
                c4 = P974;
                c5 = -P433;
                c6 = -P781;
            } else if (mu==3) {
                c1 = -P900;
                c2 = P623;
                c3 = -P222;
                c4 = P433;
                c5 = -P781;
                c6 = P974;
            } else if (mu==4) {
                c1 = -P900;
                c2 = P623;
                c3 = -P222;
                c4 = -P433;
                c5 = P781;
                c6 = -P974;
            } else if (mu==5) {
                c1 = -P222;
                c2 = -P900;
                c3 = P623;
                c4 = -P974;
                c5 = P433;
                c6 = P781;
            } else {
                c1 = P623;
                c2 = -P222;
                c3 = -P900;
                c4 = -P781;
                c5 = -P974;
                c6 = -P433;
            }
            for (l=0; l<m; l++) {
                j0 = i0;
                j1 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j1]+z[j6];
                    t1i = z[j1+1]+z[j6+1];
                    t2r = z[j2]+z[j5];
                    t2i = z[j2+1]+z[j5+1];
                    t3r = z[j3]+z[j4];
                    t3i = z[j3+1]+z[j4+1];
                    t4r = z[j1]-z[j6];
                    t4i = z[j1+1]-z[j6+1];
                    t5r = z[j2]-z[j5];
                    t5i = z[j2+1]-z[j5+1];
                    t6r = z[j3]-z[j4];
                    t6i = z[j3+1]-z[j4+1];
                    t7r = z[j0]-0.5*t3r;
                    t7i = z[j0+1]-0.5*t3i;
                    t8r = t1r-t3r;
                    t8i = t1i-t3i;
                    t9r = t2r-t3r;
                    t9i = t2i-t3i;
                    y1r = t7r+c1*t8r+c2*t9r;
                    y1i = t7i+c1*t8i+c2*t9i;
                    y2r = t7r+c2*t8r+c3*t9r;
                    y2i = t7i+c2*t8i+c3*t9i;
                    y3r = t7r+c3*t8r+c1*t9r;
                    y3i = t7i+c3*t8i+c1*t9i;
                    y4r = c6*t4r-c4*t5r+c5*t6r;
                    y4i = c6*t4i-c4*t5i+c5*t6i;
                    y5r = c5*t4r-c6*t5r-c4*t6r;
                    y5i = c5*t4i-c6*t5i-c4*t6i;
                    y6r = c4*t4r+c5*t5r+c6*t6r;
                    y6i = c4*t4i+c5*t5i+c6*t6i;
                    z[j0] = z[j0]+t1r+t2r+t3r;
                    z[j0+1] = z[j0+1]+t1i+t2i+t3i;
                    z[j1] = y1r-y6i;
                    z[j1+1] = y1i+y6r;
                    z[j2] = y2r-y5i;
                    z[j2+1] = y2i+y5r;
                    z[j3] = y3r-y4i;
                    z[j3+1] = y3i+y4r;
                    z[j4] = y3r+y4i;
                    z[j4+1] = y3i-y4r;
                    z[j5] = y2r+y5i;
                    z[j5+1] = y2i-y5r;
                    z[j6] = y1r+y6i;
                    z[j6+1] = y1i-y6r;
                    j0 += jstep;
                    j1 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                }
                it = i6+istep;
                i6 = i5+istep;
                i5 = i4+istep;
                i4 = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i7 = i6+iinc;
        if (i7>=imax) i7 = i7-imax;
	
        /* if factor is 8 */
        if (ifac==8) {
            if (mu==1) {
                c1 = 1.0;
                c2 = P707;
            } else if (mu==3) {
                c1 = -1.0;
                c2 = -P707;
            } else if (mu==5) {
                c1 = 1.0;
                c2 = -P707;
            } else {
                c1 = -1.0;
                c2 = P707;
            }
            c3 = c1*c2;
            for (l=0; l<m; l++) {
                j0 = i0;
                j1 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                j7 = i7;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j0]+z[j4];
                    t1i = z[j0+1]+z[j4+1];
                    t2r = z[j0]-z[j4];
                    t2i = z[j0+1]-z[j4+1];
                    t3r = z[j1]+z[j5];
                    t3i = z[j1+1]+z[j5+1];
                    t4r = z[j1]-z[j5];
                    t4i = z[j1+1]-z[j5+1];
                    t5r = z[j2]+z[j6];
                    t5i = z[j2+1]+z[j6+1];
                    t6r = c1*(z[j2]-z[j6]);
                    t6i = c1*(z[j2+1]-z[j6+1]);
                    t7r = z[j3]+z[j7];
                    t7i = z[j3+1]+z[j7+1];
                    t8r = z[j3]-z[j7];
                    t8i = z[j3+1]-z[j7+1];
                    t9r = t1r+t5r;
                    t9i = t1i+t5i;
                    t10r = t3r+t7r;
                    t10i = t3i+t7i;
                    t11r = c2*(t4r-t8r);
                    t11i = c2*(t4i-t8i);
                    t12r = c3*(t4r+t8r);
                    t12i = c3*(t4i+t8i);
                    y1r = t2r+t11r;
                    y1i = t2i+t11i;
                    y2r = t1r-t5r;
                    y2i = t1i-t5i;
                    y3r = t2r-t11r;
                    y3i = t2i-t11i;
                    y5r = t12r-t6r;
                    y5i = t12i-t6i;
                    y6r = c1*(t3r-t7r);
                    y6i = c1*(t3i-t7i);
                    y7r = t12r+t6r;
                    y7i = t12i+t6i;
                    z[j0] = t9r+t10r;
                    z[j0+1] = t9i+t10i;
                    z[j1] = y1r-y7i;
                    z[j1+1] = y1i+y7r;
                    z[j2] = y2r-y6i;
                    z[j2+1] = y2i+y6r;
                    z[j3] = y3r-y5i;
                    z[j3+1] = y3i+y5r;
                    z[j4] = t9r-t10r;
                    z[j4+1] = t9i-t10i;
                    z[j5] = y3r+y5i;
                    z[j5+1] = y3i-y5r;
                    z[j6] = y2r+y6i;
                    z[j6+1] = y2i-y6r;
                    z[j7] = y1r+y7i;
                    z[j7+1] = y1i-y7r;
                    j0 += jstep;
                    j1 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                    j7 += jstep;
                }
                it = i7+istep;
                i7 = i6+istep;
                i6 = i5+istep;
                i5 = i4+istep;
                i4 = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i8 = i7+iinc;
        if (i8>=imax) i8 = i8-imax;
	
        /* if factor is 9 */
        if (ifac==9) {
            if (mu==1) {
                c1 = P866;
                c2 = P766;
                c3 = P642;
                c4 = P173;
                c5 = P984;
            } else if (mu==2) {
                c1 = -P866;
                c2 = P173;
                c3 = P984;
                c4 = -P939;
                c5 = P342;
            } else if (mu==4) {
                c1 = P866;
                c2 = -P939;
                c3 = P342;
                c4 = P766;
                c5 = -P642;
            } else if (mu==5) {
                c1 = -P866;
                c2 = -P939;
                c3 = -P342;
                c4 = P766;
                c5 = P642;
            } else if (mu==7) {
                c1 = P866;
                c2 = P173;
                c3 = -P984;
                c4 = -P939;
                c5 = -P342;
            } else {
                c1 = -P866;
                c2 = P766;
                c3 = -P642;
                c4 = P173;
                c5 = -P984;
            }
            c6 = c1*c2;
            c7 = c1*c3;
            c8 = c1*c4;
            c9 = c1*c5;
            for (l=0; l<m; l++) {
                j0 = i0;
                j1 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                j7 = i7;
                j8 = i8;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j3]+z[j6];
                    t1i = z[j3+1]+z[j6+1];
                    t2r = z[j0]-0.5*t1r;
                    t2i = z[j0+1]-0.5*t1i;
                    t3r = c1*(z[j3]-z[j6]);
                    t3i = c1*(z[j3+1]-z[j6+1]);
                    t4r = z[j0]+t1r;
                    t4i = z[j0+1]+t1i;
                    t5r = z[j4]+z[j7];
                    t5i = z[j4+1]+z[j7+1];
                    t6r = z[j1]-0.5*t5r;
                    t6i = z[j1+1]-0.5*t5i;
                    t7r = z[j4]-z[j7];
                    t7i = z[j4+1]-z[j7+1];
                    t8r = z[j1]+t5r;
                    t8i = z[j1+1]+t5i;
                    t9r = z[j2]+z[j5];
                    t9i = z[j2+1]+z[j5+1];
                    t10r = z[j8]-0.5*t9r;
                    t10i = z[j8+1]-0.5*t9i;
                    t11r = z[j2]-z[j5];
                    t11i = z[j2+1]-z[j5+1];
                    t12r = z[j8]+t9r;
                    t12i = z[j8+1]+t9i;
                    t13r = t8r+t12r;
                    t13i = t8i+t12i;
                    t14r = t6r+t10r;
                    t14i = t6i+t10i;
                    t15r = t6r-t10r;
                    t15i = t6i-t10i;
                    t16r = t7r+t11r;
                    t16i = t7i+t11i;
                    t17r = t7r-t11r;
                    t17i = t7i-t11i;
                    t18r = c2*t14r-c7*t17r;
                    t18i = c2*t14i-c7*t17i;
                    t19r = c4*t14r+c9*t17r;
                    t19i = c4*t14i+c9*t17i;
                    t20r = c3*t15r+c6*t16r;
                    t20i = c3*t15i+c6*t16i;
                    t21r = c5*t15r-c8*t16r;
                    t21i = c5*t15i-c8*t16i;
                    t22r = t18r+t19r;
                    t22i = t18i+t19i;
                    t23r = t20r-t21r;
                    t23i = t20i-t21i;
                    y1r = t2r+t18r;
                    y1i = t2i+t18i;
                    y2r = t2r+t19r;
                    y2i = t2i+t19i;
                    y3r = t4r-0.5*t13r;
                    y3i = t4i-0.5*t13i;
                    y4r = t2r-t22r;
                    y4i = t2i-t22i;
                    y5r = t3r-t23r;
                    y5i = t3i-t23i;
                    y6r = c1*(t8r-t12r);
                    y6i = c1*(t8i-t12i);
                    y7r = t21r-t3r;
                    y7i = t21i-t3i;
                    y8r = t3r+t20r;
                    y8i = t3i+t20i;
                    z[j0] = t4r+t13r;
                    z[j0+1] = t4i+t13i;
                    z[j1] = y1r-y8i;
                    z[j1+1] = y1i+y8r;
                    z[j2] = y2r-y7i;
                    z[j2+1] = y2i+y7r;
                    z[j3] = y3r-y6i;
                    z[j3+1] = y3i+y6r;
                    z[j4] = y4r-y5i;
                    z[j4+1] = y4i+y5r;
                    z[j5] = y4r+y5i;
                    z[j5+1] = y4i-y5r;
                    z[j6] = y3r+y6i;
                    z[j6+1] = y3i-y6r;
                    z[j7] = y2r+y7i;
                    z[j7+1] = y2i-y7r;
                    z[j8] = y1r+y8i;
                    z[j8+1] = y1i-y8r;
                    j0 += jstep;
                    j1 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                    j7 += jstep;
                    j8 += jstep;
                }
                it = i8+istep;
                i8 = i7+istep;
                i7 = i6+istep;
                i6 = i5+istep;
                i5 = i4+istep;
                i4 = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i9 = i8+iinc;
        if (i9>=imax) i9 = i9-imax;
        i10 = i9+iinc;
        if (i10>=imax) i10 = i10-imax;
	
        /* if factor is 11 */
        if (ifac==11) {
            if (mu==1) {
                c1 = P841;
                c2 = P415;
                c3 = -P142;
                c4 = -P654;
                c5 = -P959;
                c6 = P540;
                c7 = P909;
                c8 = P989;
                c9 = P755;
                c10 = P281;
            } else if (mu==2) {
                c1 = P415;
                c2 = -P654;
                c3 = -P959;
                c4 = -P142;
                c5 = P841;
                c6 = P909;
                c7 = P755;
                c8 = -P281;
                c9 = -P989;
                c10 = -P540;
            } else if (mu==3) {
                c1 = -P142;
                c2 = -P959;
                c3 = P415;
                c4 = P841;
                c5 = -P654;
                c6 = P989;
                c7 = -P281;
                c8 = -P909;
                c9 = P540;
                c10 = P755;
            } else if (mu==4) {
                c1 = -P654;
                c2 = -P142;
                c3 = P841;
                c4 = -P959;
                c5 = P415;
                c6 = P755;
                c7 = -P989;
                c8 = P540;
                c9 = P281;
                c10 = -P909;
            } else if (mu==5) {
                c1 = -P959;
                c2 = P841;
                c3 = -P654;
                c4 = P415;
                c5 = -P142;
                c6 = P281;
                c7 = -P540;
                c8 = P755;
                c9 = -P909;
                c10 = P989;
            } else if (mu==6) {
                c1 = -P959;
                c2 = P841;
                c3 = -P654;
                c4 = P415;
                c5 = -P142;
                c6 = -P281;
                c7 = P540;
                c8 = -P755;
                c9 = P909;
                c10 = -P989;
            } else if (mu==7) {
                c1 = -P654;
                c2 = -P142;
                c3 = P841;
                c4 = -P959;
                c5 = P415;
                c6 = -P755;
                c7 = P989;
                c8 = -P540;
                c9 = -P281;
                c10 = P909;
            } else if (mu==8) {
                c1 = -P142;
                c2 = -P959;
                c3 = P415;
                c4 = P841;
                c5 = -P654;
                c6 = -P989;
                c7 = P281;
                c8 = P909;
                c9 = -P540;
                c10 = -P755;
            } else if (mu==9) {
                c1 = P415;
                c2 = -P654;
                c3 = -P959;
                c4 = -P142;
                c5 = P841;
                c6 = -P909;
                c7 = -P755;
                c8 = P281;
                c9 = P989;
                c10 = P540;
            } else {
                c1 = P841;
                c2 = P415;
                c3 = -P142;
                c4 = -P654;
                c5 = -P959;
                c6 = -P540;
                c7 = -P909;
                c8 = -P989;
                c9 = -P755;
                c10 = -P281;
            }
            for (l=0; l<m; l++) {
                j0 = i0;
                j1 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                j7 = i7;
                j8 = i8;
                j9 = i9;
                j10 = i10;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j1]+z[j10];
                    t1i = z[j1+1]+z[j10+1];
                    t2r = z[j2]+z[j9];
                    t2i = z[j2+1]+z[j9+1];
                    t3r = z[j3]+z[j8];
                    t3i = z[j3+1]+z[j8+1];
                    t4r = z[j4]+z[j7];
                    t4i = z[j4+1]+z[j7+1];
                    t5r = z[j5]+z[j6];
                    t5i = z[j5+1]+z[j6+1];
                    t6r = z[j1]-z[j10];
                    t6i = z[j1+1]-z[j10+1];
                    t7r = z[j2]-z[j9];
                    t7i = z[j2+1]-z[j9+1];
                    t8r = z[j3]-z[j8];
                    t8i = z[j3+1]-z[j8+1];
                    t9r = z[j4]-z[j7];
                    t9i = z[j4+1]-z[j7+1];
                    t10r = z[j5]-z[j6];
                    t10i = z[j5+1]-z[j6+1];
                    t11r = z[j0]-0.5*t5r;
                    t11i = z[j0+1]-0.5*t5i;
                    t12r = t1r-t5r;
                    t12i = t1i-t5i;
                    t13r = t2r-t5r;
                    t13i = t2i-t5i;
                    t14r = t3r-t5r;
                    t14i = t3i-t5i;
                    t15r = t4r-t5r;
                    t15i = t4i-t5i;
                    y1r = t11r+c1*t12r+c2*t13r+c3*t14r+c4*t15r;
                    y1i = t11i+c1*t12i+c2*t13i+c3*t14i+c4*t15i;
                    y2r = t11r+c2*t12r+c4*t13r+c5*t14r+c3*t15r;
                    y2i = t11i+c2*t12i+c4*t13i+c5*t14i+c3*t15i;
                    y3r = t11r+c3*t12r+c5*t13r+c2*t14r+c1*t15r;
                    y3i = t11i+c3*t12i+c5*t13i+c2*t14i+c1*t15i;
                    y4r = t11r+c4*t12r+c3*t13r+c1*t14r+c5*t15r;
                    y4i = t11i+c4*t12i+c3*t13i+c1*t14i+c5*t15i;
                    y5r = t11r+c5*t12r+c1*t13r+c4*t14r+c2*t15r;
                    y5i = t11i+c5*t12i+c1*t13i+c4*t14i+c2*t15i;
                    y6r = c10*t6r-c6*t7r+c9*t8r-c7*t9r+c8*t10r;
                    y6i = c10*t6i-c6*t7i+c9*t8i-c7*t9i+c8*t10i;
                    y7r = c9*t6r-c8*t7r+c6*t8r+c10*t9r-c7*t10r;
                    y7i = c9*t6i-c8*t7i+c6*t8i+c10*t9i-c7*t10i;
                    y8r = c8*t6r-c10*t7r-c7*t8r+c6*t9r+c9*t10r;
                    y8i = c8*t6i-c10*t7i-c7*t8i+c6*t9i+c9*t10i;
                    y9r = c7*t6r+c9*t7r-c10*t8r-c8*t9r-c6*t10r;
                    y9i = c7*t6i+c9*t7i-c10*t8i-c8*t9i-c6*t10i;
                    y10r = c6*t6r+c7*t7r+c8*t8r+c9*t9r+c10*t10r;
                    y10i = c6*t6i+c7*t7i+c8*t8i+c9*t9i+c10*t10i;
                    z[j0] = z[j0]+t1r+t2r+t3r+t4r+t5r;
                    z[j0+1] = z[j0+1]+t1i+t2i+t3i+t4i+t5i;
                    z[j1] = y1r-y10i;
                    z[j1+1] = y1i+y10r;
                    z[j2] = y2r-y9i;
                    z[j2+1] = y2i+y9r;
                    z[j3] = y3r-y8i;
                    z[j3+1] = y3i+y8r;
                    z[j4] = y4r-y7i;
                    z[j4+1] = y4i+y7r;
                    z[j5] = y5r-y6i;
                    z[j5+1] = y5i+y6r;
                    z[j6] = y5r+y6i;
                    z[j6+1] = y5i-y6r;
                    z[j7] = y4r+y7i;
                    z[j7+1] = y4i-y7r;
                    z[j8] = y3r+y8i;
                    z[j8+1] = y3i-y8r;
                    z[j9] = y2r+y9i;
                    z[j9+1] = y2i-y9r;
                    z[j10] = y1r+y10i;
                    z[j10+1] = y1i-y10r;
                    j0 += jstep;
                    j1 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                    j7 += jstep;
                    j8 += jstep;
                    j9 += jstep;
                    j10 += jstep;
                }
                it = i10+istep;
                i10 = i9+istep;
                i9 = i8+istep;
                i8 = i7+istep;
                i7 = i6+istep;
                i6 = i5+istep;
                i5 = i4+istep;
                i4 = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i11 = i10+iinc;
        if (i11>=imax) i11 = i11-imax;
        i12 = i11+iinc;
        if (i12>=imax) i12 = i12-imax;
	
        /* if factor is 13 */
        if (ifac==13) {
            if (mu==1) {
                c1 = P885;
                c2 = P568;
                c3 = P120;
                c4 = -P354;
                c5 = -P748;
                c6 = -P970;
                c7 = P464;
                c8 = P822;
                c9 = P992;
                c10 = P935;
                c11 = P663;
                c12 = P239;
            } else if (mu==2) {
                c1 = P568;
                c2 = -P354;
                c3 = -P970;
                c4 = -P748;
                c5 = P120;
                c6 = P885;
                c7 = P822;
                c8 = P935;
                c9 = P239;
                c10 = -P663;
                c11 = -P992;
                c12 = -P464;
            } else if (mu==3) {
                c1 = P120;
                c2 = -P970;
                c3 = -P354;
                c4 = P885;
                c5 = P568;
                c6 = -P748;
                c7 = P992;
                c8 = P239;
                c9 = -P935;
                c10 = -P464;
                c11 = P822;
                c12 = P663;
            } else if (mu==4) {
                c1 = -P354;
                c2 = -P748;
                c3 = P885;
                c4 = P120;
                c5 = -P970;
                c6 = P568;
                c7 = P935;
                c8 = -P663;
                c9 = -P464;
                c10 = P992;
                c11 = -P239;
                c12 = -P822;
            } else if (mu==5) {
                c1 = -P748;
                c2 = P120;
                c3 = P568;
                c4 = -P970;
                c5 = P885;
                c6 = -P354;
                c7 = P663;
                c8 = -P992;
                c9 = P822;
                c10 = -P239;
                c11 = -P464;
                c12 = P935;
            } else if (mu==6) {
                c1 = -P970;
                c2 = P885;
                c3 = -P748;
                c4 = P568;
                c5 = -P354;
                c6 = P120;
                c7 = P239;
                c8 = -P464;
                c9 = P663;
                c10 = -P822;
                c11 = P935;
                c12 = -P992;
            } else if (mu==7) {
                c1 = -P970;
                c2 = P885;
                c3 = -P748;
                c4 = P568;
                c5 = -P354;
                c6 = P120;
                c7 = -P239;
                c8 = P464;
                c9 = -P663;
                c10 = P822;
                c11 = -P935;
                c12 = P992;
            } else if (mu==8) {
                c1 = -P748;
                c2 = P120;
                c3 = P568;
                c4 = -P970;
                c5 = P885;
                c6 = -P354;
                c7 = -P663;
                c8 = P992;
                c9 = -P822;
                c10 = P239;
                c11 = P464;
                c12 = -P935;
            } else if (mu==9) {
                c1 = -P354;
                c2 = -P748;
                c3 = P885;
                c4 = P120;
                c5 = -P970;
                c6 = P568;
                c7 = -P935;
                c8 = P663;
                c9 = P464;
                c10 = -P992;
                c11 = P239;
                c12 = P822;
            } else if (mu==10) {
                c1 = P120;
                c2 = -P970;
                c3 = -P354;
                c4 = P885;
                c5 = P568;
                c6 = -P748;
                c7 = -P992;
                c8 = -P239;
                c9 = P935;
                c10 = P464;
                c11 = -P822;
                c12 = -P663;
            } else if (mu==11) {
                c1 = P568;
                c2 = -P354;
                c3 = -P970;
                c4 = -P748;
                c5 = P120;
                c6 = P885;
                c7 = -P822;
                c8 = -P935;
                c9 = -P239;
                c10 = P663;
                c11 = P992;
                c12 = P464;
            } else {
                c1 = P885;
                c2 = P568;
                c3 = P120;
                c4 = -P354;
                c5 = -P748;
                c6 = -P970;
                c7 = -P464;
                c8 = -P822;
                c9 = -P992;
                c10 = -P935;
                c11 = -P663;
                c12 = -P239;
            }
            for (l=0; l<m; l++) {
                j0 = i0;
                j1 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                j7 = i7;
                j8 = i8;
                j9 = i9;
                j10 = i10;
                j11 = i11;
                j12 = i12;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j1]+z[j12];
                    t1i = z[j1+1]+z[j12+1];
                    t2r = z[j2]+z[j11];
                    t2i = z[j2+1]+z[j11+1];
                    t3r = z[j3]+z[j10];
                    t3i = z[j3+1]+z[j10+1];
                    t4r = z[j4]+z[j9];
                    t4i = z[j4+1]+z[j9+1];
                    t5r = z[j5]+z[j8];
                    t5i = z[j5+1]+z[j8+1];
                    t6r = z[j6]+z[j7];
                    t6i = z[j6+1]+z[j7+1];
                    t7r = z[j1]-z[j12];
                    t7i = z[j1+1]-z[j12+1];
                    t8r = z[j2]-z[j11];
                    t8i = z[j2+1]-z[j11+1];
                    t9r = z[j3]-z[j10];
                    t9i = z[j3+1]-z[j10+1];
                    t10r = z[j4]-z[j9];
                    t10i = z[j4+1]-z[j9+1];
                    t11r = z[j5]-z[j8];
                    t11i = z[j5+1]-z[j8+1];
                    t12r = z[j6]-z[j7];
                    t12i = z[j6+1]-z[j7+1];
                    t13r = z[j0]-0.5*t6r;
                    t13i = z[j0+1]-0.5*t6i;
                    t14r = t1r-t6r;
                    t14i = t1i-t6i;
                    t15r = t2r-t6r;
                    t15i = t2i-t6i;
                    t16r = t3r-t6r;
                    t16i = t3i-t6i;
                    t17r = t4r-t6r;
                    t17i = t4i-t6i;
                    t18r = t5r-t6r;
                    t18i = t5i-t6i;
                    y1r = t13r+c1*t14r+c2*t15r+c3*t16r+c4*t17r+c5*t18r;
                    y1i = t13i+c1*t14i+c2*t15i+c3*t16i+c4*t17i+c5*t18i;
                    y2r = t13r+c2*t14r+c4*t15r+c6*t16r+c5*t17r+c3*t18r;
                    y2i = t13i+c2*t14i+c4*t15i+c6*t16i+c5*t17i+c3*t18i;
                    y3r = t13r+c3*t14r+c6*t15r+c4*t16r+c1*t17r+c2*t18r;
                    y3i = t13i+c3*t14i+c6*t15i+c4*t16i+c1*t17i+c2*t18i;
                    y4r = t13r+c4*t14r+c5*t15r+c1*t16r+c3*t17r+c6*t18r;
                    y4i = t13i+c4*t14i+c5*t15i+c1*t16i+c3*t17i+c6*t18i;
                    y5r = t13r+c5*t14r+c3*t15r+c2*t16r+c6*t17r+c1*t18r;
                    y5i = t13i+c5*t14i+c3*t15i+c2*t16i+c6*t17i+c1*t18i;
                    y6r = t13r+c6*t14r+c1*t15r+c5*t16r+c2*t17r+c4*t18r;
                    y6i = t13i+c6*t14i+c1*t15i+c5*t16i+c2*t17i+c4*t18i;
                    y7r = c12*t7r-c7*t8r+c11*t9r-c8*t10r+c10*t11r-c9*t12r;
                    y7i = c12*t7i-c7*t8i+c11*t9i-c8*t10i+c10*t11i-c9*t12i;
                    y8r = c11*t7r-c9*t8r+c8*t9r-c12*t10r-c7*t11r+c10*t12r;
                    y8i = c11*t7i-c9*t8i+c8*t9i-c12*t10i-c7*t11i+c10*t12i;
                    y9r = c10*t7r-c11*t8r-c7*t9r+c9*t10r-c12*t11r-c8*t12r;
                    y9i = c10*t7i-c11*t8i-c7*t9i+c9*t10i-c12*t11i-c8*t12i;
                    y10r = c9*t7r+c12*t8r-c10*t9r-c7*t10r+c8*t11r+c11*t12r;
                    y10i = c9*t7i+c12*t8i-c10*t9i-c7*t10i+c8*t11i+c11*t12i;
                    y11r = c8*t7r+c10*t8r+c12*t9r-c11*t10r-c9*t11r-c7*t12r;
                    y11i = c8*t7i+c10*t8i+c12*t9i-c11*t10i-c9*t11i-c7*t12i;
                    y12r = c7*t7r+c8*t8r+c9*t9r+c10*t10r+c11*t11r+c12*t12r;
                    y12i = c7*t7i+c8*t8i+c9*t9i+c10*t10i+c11*t11i+c12*t12i;
                    z[j0] = z[j0]+t1r+t2r+t3r+t4r+t5r+t6r;
                    z[j0+1] = z[j0+1]+t1i+t2i+t3i+t4i+t5i+t6i;
                    z[j1] = y1r-y12i;
                    z[j1+1] = y1i+y12r;
                    z[j2] = y2r-y11i;
                    z[j2+1] = y2i+y11r;
                    z[j3] = y3r-y10i;
                    z[j3+1] = y3i+y10r;
                    z[j4] = y4r-y9i;
                    z[j4+1] = y4i+y9r;
                    z[j5] = y5r-y8i;
                    z[j5+1] = y5i+y8r;
                    z[j6] = y6r-y7i;
                    z[j6+1] = y6i+y7r;
                    z[j7] = y6r+y7i;
                    z[j7+1] = y6i-y7r;
                    z[j8] = y5r+y8i;
                    z[j8+1] = y5i-y8r;
                    z[j9] = y4r+y9i;
                    z[j9+1] = y4i-y9r;
                    z[j10] = y3r+y10i;
                    z[j10+1] = y3i-y10r;
                    z[j11] = y2r+y11i;
                    z[j11+1] = y2i-y11r;
                    z[j12] = y1r+y12i;
                    z[j12+1] = y1i-y12r;
                    j0 += jstep;
                    j1 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                    j7 += jstep;
                    j8 += jstep;
                    j9 += jstep;
                    j10 += jstep;
                    j11 += jstep;
                    j12 += jstep;
                }
                it = i12+istep;
                i12 = i11+istep;
                i11 = i10+istep;
                i10 = i9+istep;
                i9 = i8+istep;
                i8 = i7+istep;
                i7 = i6+istep;
                i6 = i5+istep;
                i5 = i4+istep;
                i4 = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i13 = i12+iinc;
        if (i13>=imax) i13 = i13-imax;
        i14 = i13+iinc;
        if (i14>=imax) i14 = i14-imax;
        i15 = i14+iinc;
        if (i15>=imax) i15 = i15-imax;

        /* if factor is 16 */
        if (ifac==16) {
            if (mu==1) {
                c1 = 1.0;
                c2 = P923;
                c3 = P382;
                c4 = P707;
            } else if (mu==3) {
                c1 = -1.0;
                c2 = P382;
                c3 = P923;
                c4 = -P707;
            } else if (mu==5) {
                c1 = 1.0;
                c2 = -P382;
                c3 = P923;
                c4 = -P707;
            } else if (mu==7) {
                c1 = -1.0;
                c2 = -P923;
                c3 = P382;
                c4 = P707;
            } else if (mu==9) {
                c1 = 1.0;
                c2 = -P923;
                c3 = -P382;
                c4 = P707;
            } else if (mu==11) {
                c1 = -1.0;
                c2 = -P382;
                c3 = -P923;
                c4 = -P707;
            } else if (mu==13) {
                c1 = 1.0;
                c2 = P382;
                c3 = -P923;
                c4 = -P707;
            } else {
                c1 = -1.0;
                c2 = P923;
                c3 = -P382;
                c4 = P707;
            }
            c5 = c1*c4;
            c6 = c1*c3;
            c7 = c1*c2;
            for (l=0; l<m; l++) {
                j0 = i0;
                j1 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                j7 = i7;
                j8 = i8;
                j9 = i9;
                j10 = i10;
                j11 = i11;
                j12 = i12;
                j13 = i13;
                j14 = i14;
                j15 = i15;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j0]+z[j8];
                    t1i = z[j0+1]+z[j8+1];
                    t2r = z[j4]+z[j12];
                    t2i = z[j4+1]+z[j12+1];
                    t3r = z[j0]-z[j8];
                    t3i = z[j0+1]-z[j8+1];
                    t4r = c1*(z[j4]-z[j12]);
                    t4i = c1*(z[j4+1]-z[j12+1]);
                    t5r = t1r+t2r;
                    t5i = t1i+t2i;
                    t6r = t1r-t2r;
                    t6i = t1i-t2i;
                    t7r = z[j1]+z[j9];
                    t7i = z[j1+1]+z[j9+1];
                    t8r = z[j5]+z[j13];
                    t8i = z[j5+1]+z[j13+1];
                    t9r = z[j1]-z[j9];
                    t9i = z[j1+1]-z[j9+1];
                    t10r = z[j5]-z[j13];
                    t10i = z[j5+1]-z[j13+1];
                    t11r = t7r+t8r;
                    t11i = t7i+t8i;
                    t12r = t7r-t8r;
                    t12i = t7i-t8i;
                    t13r = z[j2]+z[j10];
                    t13i = z[j2+1]+z[j10+1];
                    t14r = z[j6]+z[j14];
                    t14i = z[j6+1]+z[j14+1];
                    t15r = z[j2]-z[j10];
                    t15i = z[j2+1]-z[j10+1];
                    t16r = z[j6]-z[j14];
                    t16i = z[j6+1]-z[j14+1];
                    t17r = t13r+t14r;
                    t17i = t13i+t14i;
                    t18r = c4*(t15r-t16r);
                    t18i = c4*(t15i-t16i);
                    t19r = c5*(t15r+t16r);
                    t19i = c5*(t15i+t16i);
                    t20r = c1*(t13r-t14r);
                    t20i = c1*(t13i-t14i);
                    t21r = z[j3]+z[j11];
                    t21i = z[j3+1]+z[j11+1];
                    t22r = z[j7]+z[j15];
                    t22i = z[j7+1]+z[j15+1];
                    t23r = z[j3]-z[j11];
                    t23i = z[j3+1]-z[j11+1];
                    t24r = z[j7]-z[j15];
                    t24i = z[j7+1]-z[j15+1];
                    t25r = t21r+t22r;
                    t25i = t21i+t22i;
                    t26r = t21r-t22r;
                    t26i = t21i-t22i;
                    t27r = t9r+t24r;
                    t27i = t9i+t24i;
                    t28r = t10r+t23r;
                    t28i = t10i+t23i;
                    t29r = t9r-t24r;
                    t29i = t9i-t24i;
                    t30r = t10r-t23r;
                    t30i = t10i-t23i;
                    t31r = t5r+t17r;
                    t31i = t5i+t17i;
                    t32r = t11r+t25r;
                    t32i = t11i+t25i;
                    t33r = t3r+t18r;
                    t33i = t3i+t18i;
                    t34r = c2*t29r-c6*t30r;
                    t34i = c2*t29i-c6*t30i;
                    t35r = t3r-t18r;
                    t35i = t3i-t18i;
                    t36r = c7*t27r-c3*t28r;
                    t36i = c7*t27i-c3*t28i;
                    t37r = t4r+t19r;
                    t37i = t4i+t19i;
                    t38r = c3*t27r+c7*t28r;
                    t38i = c3*t27i+c7*t28i;
                    t39r = t4r-t19r;
                    t39i = t4i-t19i;
                    t40r = c6*t29r+c2*t30r;
                    t40i = c6*t29i+c2*t30i;
                    t41r = c4*(t12r-t26r);
                    t41i = c4*(t12i-t26i);
                    t42r = c5*(t12r+t26r);
                    t42i = c5*(t12i+t26i);
                    y1r = t33r+t34r;
                    y1i = t33i+t34i;
                    y2r = t6r+t41r;
                    y2i = t6i+t41i;
                    y3r = t35r+t40r;
                    y3i = t35i+t40i;
                    y4r = t5r-t17r;
                    y4i = t5i-t17i;
                    y5r = t35r-t40r;
                    y5i = t35i-t40i;
                    y6r = t6r-t41r;
                    y6i = t6i-t41i;
                    y7r = t33r-t34r;
                    y7i = t33i-t34i;
                    y9r = t38r-t37r;
                    y9i = t38i-t37i;
                    y10r = t42r-t20r;
                    y10i = t42i-t20i;
                    y11r = t36r+t39r;
                    y11i = t36i+t39i;
                    y12r = c1*(t11r-t25r);
                    y12i = c1*(t11i-t25i);
                    y13r = t36r-t39r;
                    y13i = t36i-t39i;
                    y14r = t42r+t20r;
                    y14i = t42i+t20i;
                    y15r = t38r+t37r;
                    y15i = t38i+t37i;
                    z[j0] = t31r+t32r;
                    z[j0+1] = t31i+t32i;
                    z[j1] = y1r-y15i;
                    z[j1+1] = y1i+y15r;
                    z[j2] = y2r-y14i;
                    z[j2+1] = y2i+y14r;
                    z[j3] = y3r-y13i;
                    z[j3+1] = y3i+y13r;
                    z[j4] = y4r-y12i;
                    z[j4+1] = y4i+y12r;
                    z[j5] = y5r-y11i;
                    z[j5+1] = y5i+y11r;
                    z[j6] = y6r-y10i;
                    z[j6+1] = y6i+y10r;
                    z[j7] = y7r-y9i;
                    z[j7+1] = y7i+y9r;
                    z[j8] = t31r-t32r;
                    z[j8+1] = t31i-t32i;
                    z[j9] = y7r+y9i;
                    z[j9+1] = y7i-y9r;
                    z[j10] = y6r+y10i;
                    z[j10+1] = y6i-y10r;
                    z[j11] = y5r+y11i;
                    z[j11+1] = y5i-y11r;
                    z[j12] = y4r+y12i;
                    z[j12+1] = y4i-y12r;
                    z[j13] = y3r+y13i;
                    z[j13+1] = y3i-y13r;
                    z[j14] = y2r+y14i;
                    z[j14+1] = y2i-y14r;
                    z[j15] = y1r+y15i;
                    z[j15+1] = y1i-y15r;
                    j0 += jstep;
                    j1 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                    j7 += jstep;
                    j8 += jstep;
                    j9 += jstep;
                    j10 += jstep;
                    j11 += jstep;
                    j12 += jstep;
                    j13 += jstep;
                    j14 += jstep;
                    j15 += jstep;
                }
                it = i15+istep;
                i15 = i14+istep;
                i14 = i13+istep;
                i13 = i12+istep;
                i12 = i11+istep;
                i11 = i10+istep;
                i10 = i9+istep;
                i9 = i8+istep;
                i8 = i7+istep;
                i7 = i6+istep;
                i6 = i5+istep;
                i5 = i4+istep;
                i4 = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
    }
}


void 
pfa2cc (int isign, int idim, int n1, int n2, complex cz[])
{
    int n,nt,k,kt;

    /* determine transform length, number of transforms, and strides */
    if (idim==1) {
        n = n1;
        nt = n2;
        k = 1;
        kt = n1;
    } else {
        n = n2;
        nt = n1;
        k = n1;
        kt = 1;
    }

    /* do multiple complex to complex transforms */
    pfamcc(isign,n,nt,k,kt,cz);
}

void 
pfa2cr (int isign, int idim, int n1, int n2, complex cz[], REAL rz[])
{
    int i1,i2,j,k,it,jt,kt,n,nt,itmul,itinc;
    REAL *z,*temp,tempr,tempi,sumr,sumi,difr,difi;
    REAL wr,wi,wpr,wpi,wtemp,theta;

    /* if transforming dimension 1 */
    if (idim==1) {

        /* copy input to output and fix dc and nyquist */
        z = (REAL*)cz;
        for (i2=0,jt=0,kt=0; i2<n2; i2++,jt+=n1,kt+=(n1+2)) {
            rz[jt+1] = z[kt]-z[kt+n1];
            rz[jt] = z[kt]+z[kt+n1];
            for (i1=2,j=jt+2,k=kt+2; i1<n1; i1++,j++,k++)
                rz[j] = z[k];
        }
        z = rz;

        /* set transform length, number of transforms and strides */
        n = n1;
        nt = n2;
        itmul = 1;
        itinc = n1;

    /* else, if transforming dimension 2 */
    } else {

        /* copy input to output and fix dc and nyquist */
        z = (REAL*)cz;
        for (i2=1; i2<n2/2; i2++) {
            for (i1=0,j=i2*n1*2; i1<n1*2; i1++,j++)
                rz[j] = z[j];
        }
        for (i1=0,j=n1*n2; i1<n1*2; i1+=2,j+=2) {
            rz[i1+1] = z[i1]-z[j];
            rz[i1] = z[i1]+z[j];
        }
        z = rz;

        /* set transform length, number of transforms and strides */
        n = n2;
        nt = n1;
        itmul = n1;
        itinc = 2;
    }

    /* initialize cosine-sine recurrence */
    theta = 2.0*M_PI/(REAL)n;
    if (isign>0) theta = -theta;
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0+wpr;
    wi = wpi;

    /* twiddle transforms simultaneously */
    for (j=2,k=n-2; j<=n/2; j+=2,k-=2) {
        jt = j*itmul;
        kt = k*itmul;
        for (it=0; it<nt; it++) {
            sumr = z[jt]+z[kt];
            sumi = z[jt+1]+z[kt+1];
            difr = z[jt]-z[kt];
            difi = z[jt+1]-z[kt+1];
            tempr = wi*difr-wr*sumi;
            tempi = wi*sumi+wr*difr;
            z[jt] = sumr+tempr;
            z[jt+1] = difi+tempi;
            z[kt] = sumr-tempr;
            z[kt+1] = tempi-difi;
            jt += itinc;
            kt += itinc;
        }
        wtemp = wr;
        wr += wr*wpr-wi*wpi;
        wi += wi*wpr+wtemp*wpi;
    }

    /* if transforming dimension 1 */
    if (idim==1) {

        /* transform as complex elements */
        pfa2cc(isign,1,n1/2,n2,(complex*)z);

    /* else, if transforming dimension 2 */
    } else {

        /* transform as complex elements */
        pfa2cc(isign,2,n1,n2/2,(complex*)z);

        /* unmerge even and odd vectors */
	   temp = (REAL*)malloc(n1*sizeof(REAL));
        for (i2=0; i2<n2; i2+=2) {
            for (i1=0,j=i2*n1+1; i1<n1; i1++,j+=2)
                temp[i1] = z[j];
            for (i1=0,j=i2*n1,k=i2*n1; i1<n1; i1++,j+=2,k++)
                z[k] = z[j];
            for (i1=0,j=(i2+1)*n1; i1<n1; i1++,j++)
                z[j] = temp[i1];
        }
        free(temp);
    }
}

void 
pfa2rc (int isign, int idim, int n1, int n2, REAL rz[], complex cz[])
{
    int i1,i2,j,k,it,jt,kt,n,nt,itmul,itinc;
    REAL *z,*temp,tempr,tempi,sumr,sumi,difr,difi;
    REAL wr,wi,wpr,wpi,wtemp,theta;

    /* copy input to output while scaling */
    z = (REAL*)cz;
    for (i2=0,j=0; i2<n2; i2++)
        for (i1=0; i1<n1; i1++,j++)
            z[j] = 0.5*rz[j];

    /* if transforming dimension 1 */
    if (idim==1) {

        /* transform as complex elements */
        pfa2cc(isign,1,n1/2,n2,cz);

        /* shift rows to make room for nyquist */
        z = (REAL*)cz;
        for (i2=n2-1; i2>0; i2--) {
            jt = i2*n1+n1-1;
            kt = jt+i2*2;
            for (i1=n1-1,j=jt,k=kt; i1>=0; i1--,j--,k--)
                z[k] = z[j];
        }

        /* set transform length, number of transforms and strides */
        n = n1;
        nt = n2;
        itmul = 1;
        itinc = n1+2;

    /* else, if transforming dimension 2 */
    } else {

        /* merge even and odd vectors */
        temp = z+n1*n2;
        for (i2=0; i2<n2; i2+=2) {
            for (i1=0,j=i2*n1; i1<n1; i1++,j++)
                temp[i1] = z[j];
            for (i1=0,j=(i2+1)*n1,k=i2*n1+1; i1<n1; i1++,j++,k+=2)
                z[k] = z[j];
            for (i1=0,j=i2*n1; i1<n1; i1++,j+=2)
                z[j] = temp[i1];
        }

        /* transform as complex elements */
        pfa2cc(isign,2,n1,n2/2,cz);

        /* set transform length, number of transforms and strides */
        n = n2;
        nt = n1;
        itmul = n1;
        itinc = 2;
    }

    /* fix dc and nyquist for each transform */
    for (it=0,j=0,k=n*itmul; it<nt; it++,j+=itinc,k+=itinc) {
        z[k] = 2.0*(z[j]-z[j+1]);
        z[j] = 2.0*(z[j]+z[j+1]);
        z[k+1] = 0.0;
        z[j+1] = 0.0;
    }

    /* initialize cosine-sine recurrence */
    theta = 2.0*M_PI/(REAL)n;
    if (isign<0) theta = -theta;
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0+wpr;
    wi = wpi;

    /* twiddle transforms simultaneously */
    for (j=2,k=n-2; j<=n/2; j+=2,k-=2) {
        jt = j*itmul;
        kt = k*itmul;
        for (it=0; it<nt; it++) {
            sumr = z[jt]+z[kt];
            sumi = z[jt+1]+z[kt+1];
            difr = z[jt]-z[kt];
            difi = z[jt+1]-z[kt+1];
            tempr = wi*difr+wr*sumi;
            tempi = wi*sumi-wr*difr;
            z[jt] = sumr+tempr;
            z[jt+1] = difi+tempi;
            z[kt] = sumr-tempr;
            z[kt+1] = tempi-difi;
            jt += itinc;
            kt += itinc;
        }
        wtemp = wr;
        wr += wr*wpr-wi*wpi;
        wi += wi*wpr+wtemp*wpi;
    }
}

#if defined TEST
/* #include "cwp.h" */
main()
{
	int nmin,n,no,nfft,nffto;
	REAL cpu,total;
	complex c[720720];
	int npfao2 (int nmin, int nmax);
		
	for (n=0; n<720720; ++n)
		c[n].r = c[n].i = 0.0;

	for (nmin=npfa(100); nmin<=10000; nmin=npfa(nmin+1)) {
		n = npfa(nmin);
		for (nfft=0,total=0.0; total<1.0; ++nfft) {
			cpu = cpusec();
			pfacc(1,n,c);
			total += cpusec()-cpu+FLT_EPSILON;
		}
		no = npfao(nmin,2*nmin);
		for (nffto=0,total=0.0; total<1.0; ++nffto) {
			cpu = cpusec();
			pfacc(1,no,c);
			total += cpusec()-cpu+FLT_EPSILON;
		}
		printf("valid n=%d cost=%g optimal n=%d cost=%g\n",
			n,1.0/nfft,no,1.0/nffto);
	}
}
#endif /* TEST */
