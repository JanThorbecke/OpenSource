#include<stdlib.h>
#include<stdint.h>
#include <stdbool.h>
#include<math.h>


#define SetBit(A,k)   ( A[((k)/64)] |= (1 <<  ((k)%64)))
#define ClearBit(A,k) ( A[((k)/64)] &= ~(1 << ((k)%64)))
#define TestBit(A,k)  ( A[((k)/64)] & (1 <<   ((k)%64)))

int fill2dfuv(float *v1, float *v2, size_t n1, size_t n2){
/*********************************************************************

	Fills holes in a 2D floating-point unit vector field by nearest
	neighbour interpolating corners into holes.

	AUTHOR:
		Max Holicki
		The Netherlands

	NOTE: Boundary conditions may cause interpolation artefacts.

	ASSUMPTION: Boundaries are filled.

*********************************************************************/
	float tmpf;
	bool *hole;  // Bitarray of already missing values
	size_t *ind; // Indices of holes
	size_t i, index, i1, i2, n, n12, num;

	/***************************/
	/* Determine Filled Values */
	/***************************/
	num=0;
	n12=n1*n2;
	if(n12==0)return(-1);
	hole=(bool*)calloc(n12,sizeof(bool));
	for(i2=0;i2<n2;i2++){
		for(i1=0;i1<n1;i1++){
			if((v1[i2*n1+i1]==0.)&&(v2[i2*n1+i1]==0.)){hole[i2*n1+i1]=1;num++;}
		}
	}
vmess("num=%zu",num);
	/**********************************/
	/* First pass of filling corners. */
	/**********************************/
	index=0;
	ind=(size_t*)malloc(num*sizeof(size_t));
	for(i2=1;i2<n2-1;i2++){
		for(i1=1;i1<n1-1;i1++){
			if(hole[i2*n1+i1]){
				n=8-hole[(i2-1)*n1+i1-1]-hole[i2*n1+i1-1]-hole[(i2+1)*n1+i1-1]-
				    hole[(i2-1)*n1+i1  ]        -         hole[(i2+1)*n1+i1  ]-
				    hole[(i2-1)*n1+i1+1]-hole[i2*n1+i1+1]-hole[(i2+1)*n1+i1+1];
vmess("i1=%zu i2=%zu n=%zu",i1,i2,n);
				if(n>3){
					/* Interpolate */
					v1[i2*n1+i1]=(v1[i2*n1+i1-1-n1]+v1[i2*n1+i1-1]+v1[i2*n1+i1-1+n1]+
					       v1[i2*n1+i1  -n1]    +    v1[i2*n1+i1  +n1]+
					       v1[i2*n1+i1+1-n1]+v1[i2*n1+i1+1]+v1[i2*n1+i1+1+n1])/((float)n);
					v2[i2*n1+i1]=(v2[i2*n1+i1-1-n1]+v2[i2*n1+i1-1]+v2[i2*n1+i1-1+n1]+
					       v2[i2*n1+i1  -n1]    +    v2[i2*n1+i1  +n1]+
					       v2[i2*n1+i1+1-n1]+v2[i2*n1+i1+1]+v2[i2*n1+i1+1+n1])/((float)n);
					/* Normalize */
					tmpf=sqrtf(v1[i2*n1+i1]*v1[i2*n1+i1]+v2[i2*n1+i1]*v2[i2*n1+i1]);
					v1[i2*n1+i1]/=tmpf;v2[i2*n1+i1]/=tmpf;
					/* Clear Bit & Decrement */
					hole[i2*n1+i1]=0;
					num--; // Decrement number of remaining holes
				}else{
					ind[index++]=i2*n1+i1;
				}
			}
		}
	}
//	ind[index]=n12;
vmess("%Zu",index);
	/**********************************/
	/* Remaining Passes To Fill Holes */
	/**********************************/
	num=index; // Placeholder for number of holes remaining
	while(index>0){
		for(i=0;i<num;i++){
			if(hole[ind[i]]){
				n=8-hole[ind[i]-1-n1]-hole[ind[i]-1]-hole[ind[i]-1+n1]-
				    hole[ind[i]  -n1]        -       hole[ind[i]  +n1]-
				    hole[ind[i]+1-n1]-hole[ind[i]+1]-hole[ind[i]+1+n1];
//vmess("n=%zu i1=%Zu i2=%Zu v1=%f v2=%f",n,ind[i]/n1,ind[i]%n1,v1[ind[i]],v2[ind[i]]);
				if(n>3){
					/* Interpolate */
					v1[ind[i]]=(v1[ind[i]-1-n1]+v1[ind[i]-1]+v1[ind[i]-1+n1]+
					            v1[ind[i]  -n1]             +v1[ind[i]  +n1]+
					            v1[ind[i]+1-n1]+v1[ind[i]+1]+v1[ind[i]+1+n1])/((float)n);
					v2[ind[i]]=(v2[ind[i]-1-n1]+v2[ind[i]-1]+v2[ind[i]-1+n1]+
					            v2[ind[i]  -n1]             +v2[ind[i]  +n1]+
					            v2[ind[i]+1-n1]+v2[ind[i]+1]+v2[ind[i]+1+n1])/((float)n);
					/* Normalize */
					tmpf=sqrtf(v1[ind[i]]*v1[ind[i]]+v2[ind[i]]*v2[ind[i]]);
					v1[ind[i]]/=tmpf;v2[ind[i]]/=tmpf;
					/* Clear Bit & Decrement */
					hole[ind[i]]=0;
					index--; // Decrement number of remaining holes
				}
			}
		}
vmess("ind=%d",index);
	}

	/***************/
	/* Free Arrays */
	/***************/
	free(hole);
	free(ind);

	return(0);
}
