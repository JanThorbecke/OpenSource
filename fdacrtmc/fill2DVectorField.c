#include<stdlib.h>
#include<stdint.h>
#include<math.h>

#define SetBit(A,k)   ( A[((k)/32)] |= (1 <<  ((k)%32)))
#define ClearBit(A,k) ( A[((k)/32)] &= ~(1 << ((k)%32)))
#define TestBit(A,k)  ( (A[((k)/32)] & (1 <<   ((k)%32)))!=0)

#pragma optimize("", off) 
int fill2DVectorField(float *v1, float *v2, float p1, float p2, size_t n1, size_t n2){
/*********************************************************************

	Fills holes in a 2D floating-point vector field by nearest
	neighbour interpolating corners into holes.

	AUTHOR:
		Max Holicki
		The Netherlands

	NOTE 1: Boundary conditions may cause interpolation artefacts.
	NOTE 2: p1 & p2 are best computed based on the wavelength of
	        acoustic wavefield as they control how dominant the 
			preferred direction is.

	ASSUMPTION: Boundaries are filled.

*********************************************************************/
	uint32_t *hole; // Bitarray of already missing values
	size_t *indH, *indF; // Indices of holes (filled)
	float *tmp1, *tmp2, *tmpm; //Temporary fields to replace v1 and v2
	float *mag; // Gradient Magnitude
	float tmpf;
	size_t indexH1, indexH2, indexF, i1, i2, n, n12, num;

	/***************************/
	/* Determine Filled Values */
	/***************************/
	num=0;
	n12=n1*n2;
	if(n12==0)return(-1);
	hole=(uint32_t*)calloc((n12-1)/32+1,sizeof(uint32_t));
	mag =(float*)calloc(n12,sizeof(float));
	for(i2=0;i2<n2-0;i2++){
		for(i1=0;i1<n1-0;i1++){
			if((v1[i2*n1+i1]==0.)&&(v2[i2*n1+i1]==0.)){SetBit(hole,i2*n1+i1);num++;}else{
				mag[i2*n1+i1]=sqrtf(v1[i2*n1+i1]*v1[i2*n1+i1]+v2[i2*n1+i1]*v2[i2*n1+i1]);
			}
		}
	}

	/**********************************/
	/* First pass of filling corners. */
	/**********************************/
	indexH1=0;indexF=0;
	indH=(size_t*)malloc(num*sizeof(size_t));
	indF=(size_t*)malloc(num*sizeof(size_t));
	tmp1=(float*) malloc(num*sizeof(float ));
	tmp2=(float*) malloc(num*sizeof(float ));
	tmpm=(float*) malloc(num*sizeof(float ));
	for(i2=1;i2<n2-1;i2++){
		for(i1=1;i1<n1-1;i1++){
			if(TestBit(hole,i2*n1+i1)){ // Check if hole
				// Compute number of filled neighbours
				n=8-(TestBit(hole,i2*n1+i1-1-n1)+TestBit(hole,i2*n1+i1-1)+TestBit(hole,i2*n1+i1-1+n1)+
				     TestBit(hole,i2*n1+i1  -n1)                         +TestBit(hole,i2*n1+i1  +n1)+
				     TestBit(hole,i2*n1+i1+1-n1)+TestBit(hole,i2*n1+i1+1)+TestBit(hole,i2*n1+i1+1+n1));
				if(n>3){ //Interpolate if more than three neighbours
					/* Determine Max Magnitude */
					tmpf=fmaxf(fmaxf(fmaxf(fmax(fmax(fmax(fmaxf(mag[i2*n1+i1-1-n1],mag[i2*n1+i1-1]),
					     mag[i2*n1+i1-1+n1]),mag[i2*n1+i1-n1]),mag[i2*n1+i1+n1]),mag[i2*n1+i1+1-n1]),
						 mag[i2*n1+i1+1]),mag[i2*n1+i1+1+n1]);
					/* Interpolate */
					tmp1[indexF]=(v1[i2*n1+i1-1-n1]+v1[i2*n1+i1-1]+v1[i2*n1+i1-1+n1]+
					              v1[i2*n1+i1  -n1]               +v1[i2*n1+i1  +n1]+
					              v1[i2*n1+i1+1-n1]+v1[i2*n1+i1+1]+v1[i2*n1+i1+1+n1]+
								  ((float)(8-n))*p1*tmpf)*0.125;
					tmp2[indexF]=(v2[i2*n1+i1-1-n1]+v2[i2*n1+i1-1]+v2[i2*n1+i1-1+n1]+
					              v2[i2*n1+i1  -n1]               +v2[i2*n1+i1  +n1]+
					              v2[i2*n1+i1+1-n1]+v2[i2*n1+i1+1]+v2[i2*n1+i1+1+n1]+
								  ((float)(8-n))*p2*tmpf)*0.125;
					/* Compute Magnitude */
					tmpm[indexF] =sqrtf(tmp1[indexF]*tmp1[indexF]+tmp2[indexF]*tmp2[indexF]);
					/* Clear Bit & Decrement */
					indF[indexF++]=i2*n1+i1;
				}else{ //Otherwise store as empty
					indH[indexH1++]=i2*n1+i1;
				}
			}
		}
	}
	indF[indexF]=n12;indH[indexH1]=n12;
	/* Copy over interpolated values */
	indexF=0;
	while(indF[indexF]!=n12){
		v1[indF[indexF]] =tmp1[indexF];
		v2[indF[indexF]] =tmp2[indexF];
		mag[indF[indexF]]=tmpm[indexF];
		ClearBit(hole,indF[indexF]);
		indexF++; num--; // Decrement number of remaining holes
	}

	/**********************************/
	/* Remaining Passes To Fill Holes */
	/**********************************/
	while(num>0){
		indexH1=0;indexH2=0;indexF=0;
		while(indH[indexH1]!=n12){ // Loop over array of holes.
			if(TestBit(hole,indH[indexH1])){ // Check if hole
				// Compute number of filled neighbours
				n=8-(TestBit(hole,indH[indexH1]-1-n1)+TestBit(hole,indH[indexH1]-1)+TestBit(hole,indH[indexH1]-1+n1)+
				     TestBit(hole,indH[indexH1]  -n1)                              +TestBit(hole,indH[indexH1]  +n1)+
				     TestBit(hole,indH[indexH1]+1-n1)+TestBit(hole,indH[indexH1]+1)+TestBit(hole,indH[indexH1]+1+n1));
				if(n>2){ //Interpolate if more than three neighbours
					/* Determine Max Magnitude */
					tmpf=fmaxf(fmaxf(fmaxf(fmax(fmax(fmax(fmaxf(mag[indH[indexH1]-1-n1],mag[indH[indexH1]-1]),
					     mag[indH[indexH1]-1+n1]),mag[indH[indexH1]-n1]),mag[indH[indexH1]+n1]),
						 mag[indH[indexH1]+1-n1]),mag[indH[indexH1]+1]),mag[indH[indexH1]+1+n1]);
					/* Interpolate */
					tmp1[indexF]=(v1[indH[indexH1]-1-n1]+v1[indH[indexH1]-1]+v1[indH[indexH1]-1+n1]+
					              v1[indH[indexH1]  -n1]                    +v1[indH[indexH1]  +n1]+
					              v1[indH[indexH1]+1-n1]+v1[indH[indexH1]+1]+v1[indH[indexH1]+1+n1]+
								  ((float)(8-n))*p1*tmpf)*0.125;
					tmp2[indexF]=(v2[indH[indexH1]-1-n1]+v2[indH[indexH1]-1]+v2[indH[indexH1]-1+n1]+
					              v2[indH[indexH1]  -n1]                    +v2[indH[indexH1]  +n1]+
					              v2[indH[indexH1]+1-n1]+v2[indH[indexH1]+1]+v2[indH[indexH1]+1+n1]+
								  ((float)(8-n))*p2*tmpf)*0.125;
					/* Compute Magnitude */
					tmpm[indexF] =sqrtf(tmp1[indexF]*tmp1[indexF]+tmp2[indexF]*tmp2[indexF]);
					/* Clear Bit & Decrement */
					indF[indexF++]=indH[indexH1];
				}else if(indexH1!=indexH2){ //Otherwise store as empty
					indH[indexH2++]=indH[indexH1];
				}else{indexH2++;}
			}
			indexH1++;
		}
		indF[indexF]=n12;indH[indexH2]=n12;
		/* Copy over interpolated values */
		indexF=0;
		while(indF[indexF]!=n12){
			v1[indF[indexF]] =tmp1[indexF];
			v2[indF[indexF]] =tmp2[indexF];
			mag[indF[indexF]]=tmpm[indexF];
			ClearBit(hole,indF[indexF]);
			indexF++;num--; // Decrement number of remaining holes
		}
	}
	
	/**************************/
	/* Normalize Vector Field */
	/**************************/
	for(i2=0;i2<n2;i2++){
		for(i1=0;i1<n1;i1++){
			if(mag[i2*n1+i1]!=0.){v1[i2*n1+i1]/=mag[i2*n1+i1];v2[i2*n1+i1]/=mag[i2*n1+i1];}
			else{mag[i2*n1+i1]=sqrtf(p1*p1+p2*p2);v1[i2*n1+i1]=p1/mag[i2*n1+i1];v2[i2*n1+i1]=p2/mag[i2*n1+i1];}
			// If amplitudes are too small, we assume preference direction
		}
	}

	/***************/
	/* Free Arrays */
	/***************/
	free(hole);
	free(indH);
	free(indF);
	free(tmp1);
	free(tmp2);
	free(tmpm);

	return(0);
}