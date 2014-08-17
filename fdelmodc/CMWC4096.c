/* http://en.wikipedia.org/wiki/Multiply-with-carry */
#include<stdlib.h>
#include<limits.h>

/* random number generator which can be used as an alternative for drand48() */

/* http://school.anhb.uwa.edu.au/personalpages/kwessen/shared/Marsaglia03.html*/

static unsigned long Q[4096],c=362436; /* choose random initial c<809430660 and */
                                         /* 4096 random 32-bit integers for Q[]   */
void seedCMWC4096(void)
{
	int i;
	for (i=0; i<4096; i++) {
		Q[i] = lrand48();
	}
	return;
}

unsigned long CMWC4096(void)
{
	unsigned long long t, a=18782LL;
	static unsigned long i=4095;
	unsigned long x,r=0xfffffffe;

	i=(i+1)&4095;
	t=a*Q[i]+c;
	c=(t>>32); 
	x=t+c; 
	if(x<c){x++;c++;}
	return(Q[i]=r-x);    
}

/* replace defaults with five random seed values in calling program */
static unsigned long x=123456789,y=362436069,z=521288629,w=88675123,v=886756453;
unsigned long xorshift(void)
{unsigned long t;
 t=(x^(x>>7)); x=y; y=z; z=w; w=v;
 v=(v^(v<<6))^(t^(t<<13)); return (y+y+1)*v;}


double dcmwc4096(void)
{
	double rd;

//	rd = ((double)xorshift())/((double)ULONG_MAX);
	rd = ((double)CMWC4096())/((double)ULONG_MAX);
	return rd;
}

