#include"fdacrtmc.h"

int mvAvg2d3(size_t n1, size_t n2, float* in, float* out){
/* mvAvg2d3(n1,n2,in,out) computes the 3x3 moving average
 * of 2d input "in" and returns it in output "out". "n1"
 * is the size of the slow dimension while "n2" is the
 * size of the fast dimension. The boundaries are not
 * touched. The array size should be at least 3-by-3,
 * this is not checked. Input & output may NOT be equal.
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float sum[3];
	float c=0.111111111111111111111111111111111111111111111;
	size_t i1,i2,j;

	for(i2=0;i2<n2;i2++)out[i2]=in[i2]; //Copy Data
	for(i1=1;i1<n1-1;i1++){ //Loop over slow dimension
		sum[0]=in[(i1-1)*n2  ]+in[i1*n2  ]+in[(i1+1)*n2  ];
		sum[1]=in[(i1-1)*n2+1]+in[i1*n2+1]+in[(i1+1)*n2+1];
		out[i1*n2]=in[i1*n2];
		for(i2=1,j=0;j<(n2-2)/3;j++){ //Loop over fast dimension - Inner 3-loop unrolled
			sum[2]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]);
			sum[0]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]);
			sum[1]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]);
		}
		// Take care of remaining values
		if(i2<n2-1){
			sum[2]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]);
		}else continue;
		if(i2<n2-1){
			sum[0]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]);
		}
		out[(i1+1)*n2-1]=in[(i1+1)*n2-1];
	}
	for(i2=0;i2<n2;i2++)out[(n1-1)*n2+i2]=in[(n1-1)*n2+i2]; //Copy Data
	return(0);
}

int mvAvg2d3Sgn(size_t n1, size_t n2, float* in, float* out){
/* mvAvg2d3Sgn(n1,n2,in,out) signs the input values according
 * to the sign of the local 3-by-3 average at that point. See
 * mvAvg2d3 for details.
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float sum[3];
	size_t i1,i2,j;

	for(i2=0;i2<n2;i2++)out[i2]=in[i2]; //Copy Data
	for(i1=1;i1<n1-1;i1++){ //Loop over slow dimension
		sum[0]=in[(i1-1)*n2  ]+in[i1*n2  ]+in[(i1+1)*n2  ];
		sum[1]=in[(i1-1)*n2+1]+in[i1*n2+1]+in[(i1+1)*n2+1];
		out[i1*n2]=in[i1*n2];
		for(i2=1,j=0;j<(n2-2)/3;j++){ //Loop over fast dimension - Inner 3-loop unrolled
			sum[2]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]);i2++;
			sum[0]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]);i2++;
			sum[1]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]);i2++;
		}
		// Take care of remaining values
		if(i2<n2-1){
			sum[2]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]);i2++;
		}else continue;
		if(i2<n2-1){
			sum[0]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]);i2++;
		}
		out[(i1+1)*n2-1]=in[(i1+1)*n2-1];
	}
	for(i2=0;i2<n2;i2++)out[(n1-1)*n2+i2]=in[(n1-1)*n2+i2]; //Copy Data
	return(0);
}

int mvAvg2d3Embd(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2){
/* mvAvg2d3Embd(n2,in,out,s1,s2,e1,e2) is the embedded form of mvAvg2d3(n1,n2,in,out). "s1" is
 * the start of the area to average in the slow dimension. "s2" is the start of the area to
 * average in the fast dimension. "e1" is the end of the area to average in the slow dimension.
 * "e2" is the end of the area to average in the fast dimension. Input & output MAY be equal.
 * Requires additional memory. Note: s1,s2>0 & e1<n1 & e2<n2
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float *tmp;
	float sum[3];
	float c=0.111111111111111111111111111111111111111111111;
	size_t i1,i1e,i2,i2e,j;
	size_t n1e=e1-s1;
	size_t n2e=e2-s2;

	tmp=(float*)malloc(n1e*n2e*sizeof(float));
	
	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){ //Loop over slow dimension
		sum[0]=in[(i1-1)*n2+s2-1]+in[i1*n2+s2-1]+in[(i1+1)*n2+s2-1];
		sum[1]=in[(i1-1)*n2+s2  ]+in[i1*n2+s2  ]+in[(i1+1)*n2+s2  ];
		for(i2=s2,i2e=0,j=0;j<n2e/3;j++){ //Loop over fast dimension - Inner 3-loop unrolled
			sum[2]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]);
			sum[0]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]);
			sum[1]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]);
		}
		// Take care of remaining values
		if(i2<e2){
			sum[2]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]);
		}else continue;
		if(i2<e2){
			sum[0]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]);
		}
	}
	//Copy to output
	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){
		for(i2=s2,i2e=0;i2<e2;i2++,i2e++){
			out[i1*n2+i2]=tmp[i1e*n2e+i2e];
		}
	}
	free(tmp);
	return(0);
}

int mvAvg2d3EmbdSgn(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2){
/* mvAvg2d3EmbdSgn(n2,in,out,s1,s2,e1,e2) is the embedded form of mvAvg2d3Sgn(n1,n2,in,out). See
 * mvAvg2d3Embd for details.
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float *tmp;
	float sum[3];
	size_t i1,i1e,i2,i2e,j;
	size_t n1e=e1-s1;
	size_t n2e=e2-s2;

	tmp=(float*)malloc(n1e*n2e*sizeof(float));
	
	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){ //Loop over slow dimension
		sum[0]=in[(i1-1)*n2+s2-1]+in[i1*n2+s2-1]+in[(i1+1)*n2+s2-1];
		sum[1]=in[(i1-1)*n2+s2  ]+in[i1*n2+s2  ]+in[(i1+1)*n2+s2  ];
		for(i2=s2,i2e=0,j=0;j<n2e/3;j++){ //Loop over fast dimension - Inner 3-loop unrolled
			sum[2]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]);
			sum[0]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]);
			sum[1]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]);
		}
		// Take care of remaining values
		if(i2<e2){
			sum[2]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]);
		}else continue;
		if(i2<e2){
			sum[0]=in[(i1-1)*n2+i2+1]+in[i1*n2+i2+1]+in[(i1+1)*n2+i2+1];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]);
		}
	}
	//Copy to output
	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){
		for(i2=s2,i2e=0;i2<e2;i2++,i2e++){
			out[i1*n2+i2]=tmp[i1e*n2e+i2e];
		}
	}
	free(tmp);
	return(0);
}

int mvAvg2d5(size_t n1, size_t n2, float* in, float* out){
/* mvAvg2d5(n1,n2,in,out) computes the 5x5 moving average
 * of 2d input "in" and returns it in output "out". "n1"
 * is the size of the slow dimension while "n2" is the
 * size of the fast dimension. The boundaries are not
 * touched. The array size should be at least 5-by-5,
 * this is not checked. Input & output may NOT be equal.
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float sum[5];
	float c=0.04;
	size_t i1,i2,j;

	for(i2=0;i2<n2;i2++)out[   i2]=in[   i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[n2+i2]=in[n2+i2]; //Copy Data
	for(i1=2;i1<n1-2;i1++){ //Loop over slow dimension
		sum[0]=in[(i1-2)*n2  ]+in[(i1-1)*n2  ]+in[i1*n2  ]+in[(i1+1)*n2  ]+in[(i1+2)*n2  ];
		sum[1]=in[(i1-2)*n2+1]+in[(i1-1)*n2+1]+in[i1*n2+1]+in[(i1+1)*n2+1]+in[(i1+2)*n2+1];
		sum[2]=in[(i1-2)*n2+2]+in[(i1-1)*n2+2]+in[i1*n2+2]+in[(i1+1)*n2+2]+in[(i1+2)*n2+2];
		sum[3]=in[(i1-2)*n2+3]+in[(i1-1)*n2+3]+in[i1*n2+3]+in[(i1+1)*n2+3]+in[(i1+2)*n2+3];
		out[i1*n2]=in[i1*n2];out[i1*n2+1]=in[i1*n2+1];
		for(i2=2,j=0;j<(n2-4)/5;j++){ //Loop over fast dimension - Inner 5-loop unrolled
			sum[4]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
			sum[0]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
			sum[1]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
			sum[2]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
			sum[3]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
		}
		// Take care of remaining values
		if(i2<n2-2){
			sum[4]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
		}else continue;
		if(i2<n2-2){
			sum[0]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
		}else continue;
		if(i2<n2-2){
			sum[1]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
		}else continue;
		if(i2<n2-2){
			sum[2]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
		}
		out[(i1+1)*n2-2]=in[(i1+1)*n2-2];out[(i1+1)*n2-1]=in[(i1+1)*n2-1];
	}
	for(i2=0;i2<n2;i2++)out[(n1-2)*n2+i2]=in[(n1-2)*n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[(n1-1)*n2+i2]=in[(n1-1)*n2+i2]; //Copy Data
	return(0);
}

int mvAvg2d5Sgn(size_t n1, size_t n2, float* in, float* out){
/* mvAvg2d5Sgn(n1,n2,in,out) signs the input values according
 * to the sign of the local 5-by-5 average at that point. See
 * mvAvg2d5 for details.
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float sum[5];
	size_t i1,i2,j;

	for(i2=0;i2<n2;i2++)out[   i2]=in[   i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[n2+i2]=in[n2+i2]; //Copy Data
	for(i1=2;i1<n1-2;i1++){ //Loop over slow dimension
		sum[0]=in[(i1-2)*n2  ]+in[(i1-1)*n2  ]+in[i1*n2  ]+in[(i1+1)*n2  ]+in[(i1+2)*n2  ];
		sum[1]=in[(i1-2)*n2+1]+in[(i1-1)*n2+1]+in[i1*n2+1]+in[(i1+1)*n2+1]+in[(i1+2)*n2+1];
		sum[2]=in[(i1-2)*n2+2]+in[(i1-1)*n2+2]+in[i1*n2+2]+in[(i1+1)*n2+2]+in[(i1+2)*n2+2];
		sum[3]=in[(i1-2)*n2+3]+in[(i1-1)*n2+3]+in[i1*n2+3]+in[(i1+1)*n2+3]+in[(i1+2)*n2+3];
		out[i1*n2]=in[i1*n2];out[i1*n2+1]=in[i1*n2+1];
		for(i2=2,j=0;j<(n2-4)/5;j++){ //Loop over fast dimension - Inner 5-loop unrolled
			sum[4]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);i2++;
			sum[0]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);i2++;
			sum[1]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);i2++;
			sum[2]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);i2++;
			sum[3]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);i2++;
		}
		// Take care of remaining values
		if(i2<n2-2){
			sum[4]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);i2++;
		}else continue;
		if(i2<n2-2){
			sum[0]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);i2++;
		}else continue;
		if(i2<n2-2){
			sum[1]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);i2++;
		}else continue;
		if(i2<n2-2){
			sum[2]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);i2++;
		}
		out[(i1+1)*n2-2]=in[(i1+1)*n2-2];out[(i1+1)*n2-1]=in[(i1+1)*n2-1];
	}
	for(i2=0;i2<n2;i2++)out[(n1-2)*n2+i2]=in[(n1-2)*n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[(n1-1)*n2+i2]=in[(n1-1)*n2+i2]; //Copy Data
	return(0);
}

int mvAvg2d5Embd(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2){
/* mvAvg2d5Embd(n2,in,out,s1,s2,e1,e2) is the embedded form of mvAvg2d5(n1,n2,in,out). "s1" is
 * the start of the area to average in the slow dimension. "s2" is the start of the area to
 * average in the fast dimension. "e1" is the end of the area to average in the slow dimension.
 * "e2" is the end of the area to average in the fast dimension. Input & output MAY be equal.
 * Requires additional memory. Note: s1,s2>1 & e1<n1-1 & e2<n2-1
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float sum[5];
	float *tmp;
	float c=0.04;
	size_t i1,i1e,i2,i2e,j;
	size_t n1e=e1-s1;
	size_t n2e=e2-s2;

	tmp=(float*)malloc(n1e*n2e*sizeof(float));

	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){ //Loop over slow dimension
		sum[0]=in[(i1-2)*n2+s2-2]+in[(i1-1)*n2+s2-2]+in[i1*n2+s2-2]+in[(i1+1)*n2+s2+1]+in[(i1+2)*n2+s2-2];
		sum[1]=in[(i1-2)*n2+s2-1]+in[(i1-1)*n2+s2-1]+in[i1*n2+s2-1]+in[(i1+1)*n2+s2+1]+in[(i1+2)*n2+s2-1];
		sum[2]=in[(i1-2)*n2+s2  ]+in[(i1-1)*n2+s2  ]+in[i1*n2+s2  ]+in[(i1+1)*n2+s2+1]+in[(i1+2)*n2+s2  ];
		sum[3]=in[(i1-2)*n2+s2+1]+in[(i1-1)*n2+s2+1]+in[i1*n2+s2+1]+in[(i1+1)*n2+s2+1]+in[(i1+2)*n2+s2+1];
		for(i2=s2,i2e=0,j=0;j<n2e/5;j++){ //Loop over fast dimension - Inner 5-loop unrolled
			sum[4]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
			sum[0]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
			sum[1]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
			sum[2]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
			sum[3]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
		}
		// Take care of remaining values
		if(i2<e2){
			sum[4]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
		}else continue;
		if(i2<e2){
			sum[0]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
		}else continue;
		if(i2<e2){
			sum[1]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
		}else continue;
		if(i2<e2){
			sum[2]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
		}
	}
	//Copy to output
	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){
		for(i2=s2,i2e=0;i2<e2;i2++,i2e++){
			out[i1*n2+i2]=tmp[i1e*n2e+i2e];
		}
	}
	free(tmp);
	return(0);
}

int mvAvg2d5EmbdSgn(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2){
/* mvAvg2d5Embd(n2,in,out,s1,s2,e1,e2) is the embedded form of mvAvg2d5(n1,n2,in,out). "s1" is
 * the start of the area to average in the slow dimension. "s2" is the start of the area to
 * average in the fast dimension. "e1" is the end of the area to average in the slow dimension.
 * "e2" is the end of the area to average in the fast dimension. Input & output MAY be equal.
 * Requires additional memory. Note: s1,s2>1 & e1<n1-1 & e2<n2-1
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float sum[5];
	float *tmp;
	size_t i1,i1e,i2,i2e,j;
	size_t n1e=e1-s1;
	size_t n2e=e2-s2;
	
	tmp=(float*)malloc(n1e*n2e*sizeof(float));

	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){ //Loop over slow dimension
		sum[0]=in[(i1-2)*n2+s2-2]+in[(i1-1)*n2+s2-2]+in[i1*n2+s2-2]+in[(i1+1)*n2+s2+1]+in[(i1+2)*n2+s2-2];
		sum[1]=in[(i1-2)*n2+s2-1]+in[(i1-1)*n2+s2-1]+in[i1*n2+s2-1]+in[(i1+1)*n2+s2+1]+in[(i1+2)*n2+s2-1];
		sum[2]=in[(i1-2)*n2+s2  ]+in[(i1-1)*n2+s2  ]+in[i1*n2+s2  ]+in[(i1+1)*n2+s2+1]+in[(i1+2)*n2+s2  ];
		sum[3]=in[(i1-2)*n2+s2+1]+in[(i1-1)*n2+s2+1]+in[i1*n2+s2+1]+in[(i1+1)*n2+s2+1]+in[(i1+2)*n2+s2+1];
		for(i2=s2,i2e=0,j=0;j<n2e/5;j++){ //Loop over fast dimension - Inner 5-loop unrolled
			sum[4]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
			sum[0]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
			sum[1]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
			sum[2]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
			sum[3]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
		}
		// Take care of remaining values
		if(i2<e2){
			sum[4]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
		}else continue;
		if(i2<e2){
			sum[0]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
		}else continue;
		if(i2<e2){
			sum[1]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
		}else continue;
		if(i2<e2){
			sum[2]=in[(i1-2)*n2+i2+2]+in[(i1-1)*n2+i2+2]+in[i1*n2+i2+2]+in[(i1+1)*n2+i2+2]+in[(i1+2)*n2+i2+2];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]);
		}
	}
	//Copy to output
	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){
		for(i2=s2,i2e=0;i2<e2;i2++,i2e++){
			out[i1*n2+i2]=tmp[i1e*n2e+i2e];
		}
	}
	free(tmp);
	return(0);
}

int mvAvg2d7(size_t n1, size_t n2, float* in, float* out){
/* mvAvg2d7(n1,n2,in,out) computes the 7x7 moving average
 * of 2d input "in" and returns it in output "out". "n1"
 * is the size of the slow dimension while "n2" is the
 * size of the fast dimension. The boundaries are not
 * touched. The array size should be at least 7-by-7,
 * this is not checked. Input & output may NOT be equal.
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float sum[7];
	float c=1.0/49.0;
	size_t i1,i2,j;

	for(i2=0;i2<n2;i2++)out[     i2]=in[     i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[  n2+i2]=in[  n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[2*n2+i2]=in[2*n2+i2]; //Copy Data
	for(i1=3;i1<n1-3;i1++){ //Loop over slow dimension
		sum[0]=in[(i1-3)*n2  ]+in[(i1-2)*n2  ]+in[(i1-1)*n2  ]+in[i1*n2  ]+in[(i1+1)*n2  ]+in[(i1+2)*n2  ]+in[(i1+3)*n2  ];
		sum[1]=in[(i1-3)*n2+1]+in[(i1-2)*n2+1]+in[(i1-1)*n2+1]+in[i1*n2+1]+in[(i1+1)*n2+1]+in[(i1+2)*n2+1]+in[(i1+3)*n2+1];
		sum[2]=in[(i1-3)*n2+2]+in[(i1-2)*n2+2]+in[(i1-1)*n2+2]+in[i1*n2+2]+in[(i1+1)*n2+2]+in[(i1+2)*n2+2]+in[(i1+3)*n2+2];
		sum[3]=in[(i1-3)*n2+3]+in[(i1-2)*n2+3]+in[(i1-1)*n2+3]+in[i1*n2+3]+in[(i1+1)*n2+3]+in[(i1+2)*n2+3]+in[(i1+3)*n2+3];
		sum[4]=in[(i1-3)*n2+4]+in[(i1-2)*n2+4]+in[(i1-1)*n2+4]+in[i1*n2+4]+in[(i1+1)*n2+4]+in[(i1+2)*n2+4]+in[(i1+3)*n2+4];
		sum[5]=in[(i1-3)*n2+5]+in[(i1-2)*n2+5]+in[(i1-1)*n2+5]+in[i1*n2+5]+in[(i1+1)*n2+5]+in[(i1+2)*n2+5]+in[(i1+3)*n2+5];
		out[i1*n2]=in[i1*n2];out[i1*n2+1]=in[i1*n2+1];out[i1*n2+2]=in[i1*n2+2];
		for(i2=3,j=0;j<(n2-6)/7;j++){ //Loop over fast dimension - Inner 7-loop unrolled
			sum[6]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
			sum[0]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
			sum[1]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
			sum[2]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
			sum[3]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
			sum[4]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
			sum[5]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}
		// Take care of remaining values
		if(i2<n2-3){
			sum[6]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<n2-3){
			sum[0]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<n2-3){
			sum[1]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<n2-3){
			sum[2]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<n2-3){
			sum[3]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<n2-3){
			sum[4]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}
		out[(i1+1)*n2-3]=in[(i1+1)*n2-3];out[(i1+1)*n2-2]=in[(i1+1)*n2-2];out[(i1+1)*n2-1]=in[(i1+1)*n2-1];
	}
	for(i2=0;i2<n2;i2++)out[(n1-3)*n2+i2]=in[(n1-3)*n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[(n1-2)*n2+i2]=in[(n1-2)*n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[(n1-1)*n2+i2]=in[(n1-1)*n2+i2]; //Copy Data
	return(0);
}

int mvAvg2d7Sgn(size_t n1, size_t n2, float* in, float* out){
/* mvAvg2d7Sgn(n1,n2,in,out) signs the input values according
 * to the sign of the local 7-by-7 average at that point. See
 * mvAvg2d7 for details.
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float sum[7];
	size_t i1,i2,j;

	for(i2=0;i2<n2;i2++)out[     i2]=in[     i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[  n2+i2]=in[  n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[2*n2+i2]=in[2*n2+i2]; //Copy Data
	for(i1=3;i1<n1-3;i1++){ //Loop over slow dimension
		sum[0]=in[(i1-3)*n2  ]+in[(i1-2)*n2  ]+in[(i1-1)*n2  ]+in[i1*n2  ]+in[(i1+1)*n2  ]+in[(i1+2)*n2  ]+in[(i1+3)*n2  ];
		sum[1]=in[(i1-3)*n2+1]+in[(i1-2)*n2+1]+in[(i1-1)*n2+1]+in[i1*n2+1]+in[(i1+1)*n2+1]+in[(i1+2)*n2+1]+in[(i1+3)*n2+1];
		sum[2]=in[(i1-3)*n2+2]+in[(i1-2)*n2+2]+in[(i1-1)*n2+2]+in[i1*n2+2]+in[(i1+1)*n2+2]+in[(i1+2)*n2+2]+in[(i1+3)*n2+2];
		sum[3]=in[(i1-3)*n2+3]+in[(i1-2)*n2+3]+in[(i1-1)*n2+3]+in[i1*n2+3]+in[(i1+1)*n2+3]+in[(i1+2)*n2+3]+in[(i1+3)*n2+3];
		sum[4]=in[(i1-3)*n2+4]+in[(i1-2)*n2+4]+in[(i1-1)*n2+4]+in[i1*n2+4]+in[(i1+1)*n2+4]+in[(i1+2)*n2+4]+in[(i1+3)*n2+4];
		sum[5]=in[(i1-3)*n2+5]+in[(i1-2)*n2+5]+in[(i1-1)*n2+5]+in[i1*n2+5]+in[(i1+1)*n2+5]+in[(i1+2)*n2+5]+in[(i1+3)*n2+5];
		out[i1*n2]=in[i1*n2];out[i1*n2+1]=in[i1*n2+1];out[i1*n2+2]=in[i1*n2+2];
		for(i2=3,j=0;j<(n2-6)/5;j++){ //Loop over fast dimension - Inner 7-loop unrolled
			sum[6]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);i2++;
			sum[0]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);i2++;
			sum[1]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);i2++;
			sum[2]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);i2++;
			sum[3]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);i2++;
			sum[4]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);i2++;
			sum[5]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);i2++;
		}
		// Take care of remaining values
		if(i2<n2-3){
			sum[6]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[1]+sum[2]+sum[3]+sum[4]+sum[5]);i2++;
		}else continue;
		if(i2<n2-3){
			sum[0]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[1]+sum[2]+sum[3]+sum[4]+sum[5]);i2++;
		}else continue;
		if(i2<n2-3){
			sum[1]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[1]+sum[2]+sum[3]+sum[4]+sum[5]);i2++;
		}else continue;
		if(i2<n2-3){
			sum[2]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[1]+sum[2]+sum[3]+sum[4]+sum[5]);i2++;
		}else continue;
		if(i2<n2-3){
			sum[3]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[1]+sum[2]+sum[3]+sum[4]+sum[5]);i2++;
		}else continue;
		if(i2<n2-3){
			sum[4]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[1]+sum[2]+sum[3]+sum[4]+sum[5]);i2++;
		}
		out[(i1+1)*n2-3]=in[(i1+1)*n2-3];out[(i1+1)*n2-2]=in[(i1+1)*n2-2];out[(i1+1)*n2-1]=in[(i1+1)*n2-1];
	}
	for(i2=0;i2<n2;i2++)out[(n1-3)*n2+i2]=in[(n1-3)*n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[(n1-2)*n2+i2]=in[(n1-2)*n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[(n1-1)*n2+i2]=in[(n1-1)*n2+i2]; //Copy Data
	return(0);
}

int mvAvg2d7Embd(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2){
/* mvAvg2d7Embd(n2,in,out,s1,s2,e1,e2) is the embedded form of mvAvg2d7(n1,n2,in,out). "s1" is
 * the start of the area to average in the slow dimension. "s2" is the start of the area to
 * average in the fast dimension. "e1" is the end of the area to average in the slow dimension.
 * "e2" is the end of the area to average in the fast dimension. Input & output MAY be equal.
 * Requires additional memory. Note: s1,s2>1 & e1<n1-1 & e2<n2-1
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float sum[7];
	float *tmp;
	float c=1.0/49.0;
	size_t i1,i1e,i2,i2e,j;
	size_t n1e=e1-s1;
	size_t n2e=e2-s2;

	tmp=(float*)malloc(n1e*n2e*sizeof(float));

	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){ //Loop over slow dimension
		sum[0]=in[(i1-3)*n2-3]+in[(i1-2)*n2+s2-3]+in[(i1-1)*n2+s2-3]+in[i1*n2+s2-3]+in[(i1+1)*n2+s2-3]+in[(i1+2)*n2+s2-3]+in[(i1+3)*n2+s2-3];
		sum[1]=in[(i1-3)*n2-2]+in[(i1-2)*n2+s2-2]+in[(i1-1)*n2+s2-2]+in[i1*n2+s2-2]+in[(i1+1)*n2+s2-2]+in[(i1+2)*n2+s2-2]+in[(i1+3)*n2+s2-2];
		sum[2]=in[(i1-3)*n2-1]+in[(i1-2)*n2+s2-1]+in[(i1-1)*n2+s2-1]+in[i1*n2+s2-1]+in[(i1+1)*n2+s2-1]+in[(i1+2)*n2+s2-1]+in[(i1+3)*n2+s2-1];
		sum[3]=in[(i1-3)*n2  ]+in[(i1-2)*n2+s2  ]+in[(i1-1)*n2+s2  ]+in[i1*n2+s2  ]+in[(i1+1)*n2+s2  ]+in[(i1+2)*n2+s2  ]+in[(i1+3)*n2+s2  ];
		sum[4]=in[(i1-3)*n2+1]+in[(i1-2)*n2+s2+1]+in[(i1-1)*n2+s2+1]+in[i1*n2+s2+1]+in[(i1+1)*n2+s2+1]+in[(i1+2)*n2+s2+1]+in[(i1+3)*n2+s2+1];
		sum[5]=in[(i1-3)*n2+2]+in[(i1-2)*n2+s2+2]+in[(i1-1)*n2+s2+2]+in[i1*n2+s2+2]+in[(i1+1)*n2+s2+2]+in[(i1+2)*n2+s2+2]+in[(i1+3)*n2+s2+2];
		for(i2=s2,i2e=0,j=0;j<n2e/7;j++){ //Loop over fast dimension - Inner 7-loop unrolled
			sum[6]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
			sum[0]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
			sum[1]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
			sum[2]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
			sum[3]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
			sum[4]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
			sum[5]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}
		// Take care of remaining values
		if(i2<e2){
			sum[6]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<e2){
			sum[0]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<e2){
			sum[1]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<e2){
			sum[2]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<e2){
			sum[3]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<e2){
			sum[4]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}
	}
	//Copy to output
	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){
		for(i2=s2,i2e=0;i2<e2;i2++,i2e++){
			out[i1*n2+i2]=tmp[i1e*n2e+i2e];
		}
	}
	free(tmp);
	return(0);
}

int mvAvg2d7EmbdSgn(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2){
/* mvAvg2d7Embd(n2,in,out,s1,s2,e1,e2) is the embedded form of mvAvg2d7(n1,n2,in,out). "s1" is
 * the start of the area to average in the slow dimension. "s2" is the start of the area to
 * average in the fast dimension. "e1" is the end of the area to average in the slow dimension.
 * "e2" is the end of the area to average in the fast dimension. Input & output MAY be equal.
 * Requires additional memory. Note: s1,s2>1 & e1<n1-1 & e2<n2-1
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float sum[7];
	float *tmp;
	size_t i1,i1e,i2,i2e,j;
	size_t n1e=e1-s1;
	size_t n2e=e2-s2;
	
	tmp=(float*)malloc(n1e*n2e*sizeof(float));

	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){ //Loop over slow dimension
		sum[0]=in[(i1-3)*n2-3]+in[(i1-2)*n2+s2-3]+in[(i1-1)*n2+s2-3]+in[i1*n2+s2-3]+in[(i1+1)*n2+s2-3]+in[(i1+2)*n2+s2-3]+in[(i1+3)*n2+s2-3];
		sum[1]=in[(i1-3)*n2-2]+in[(i1-2)*n2+s2-2]+in[(i1-1)*n2+s2-2]+in[i1*n2+s2-2]+in[(i1+1)*n2+s2-2]+in[(i1+2)*n2+s2-2]+in[(i1+3)*n2+s2-2];
		sum[2]=in[(i1-3)*n2-1]+in[(i1-2)*n2+s2-1]+in[(i1-1)*n2+s2-1]+in[i1*n2+s2-1]+in[(i1+1)*n2+s2-1]+in[(i1+2)*n2+s2-1]+in[(i1+3)*n2+s2-1];
		sum[3]=in[(i1-3)*n2  ]+in[(i1-2)*n2+s2  ]+in[(i1-1)*n2+s2  ]+in[i1*n2+s2  ]+in[(i1+1)*n2+s2  ]+in[(i1+2)*n2+s2  ]+in[(i1+3)*n2+s2  ];
		sum[4]=in[(i1-3)*n2+1]+in[(i1-2)*n2+s2+1]+in[(i1-1)*n2+s2+1]+in[i1*n2+s2+1]+in[(i1+1)*n2+s2+1]+in[(i1+2)*n2+s2+1]+in[(i1+3)*n2+s2+1];
		sum[5]=in[(i1-3)*n2+2]+in[(i1-2)*n2+s2+2]+in[(i1-1)*n2+s2+2]+in[i1*n2+s2+2]+in[(i1+1)*n2+s2+2]+in[(i1+2)*n2+s2+2]+in[(i1+3)*n2+s2+2];
		for(i2=s2,i2e=0,j=0;j<n2e/7;j++){ //Loop over fast dimension - Inner 7-loop unrolled
			sum[6]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
			sum[0]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[1]+sum[2]+sum[3]+sum[4]+sum[5]);
			sum[1]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[1]+sum[2]+sum[3]+sum[4]+sum[5]);
			sum[2]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[1]+sum[2]+sum[3]+sum[4]+sum[5]);
			sum[3]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[1]+sum[2]+sum[3]+sum[4]+sum[5]);
			sum[4]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[1]+sum[2]+sum[3]+sum[4]+sum[5]);
			sum[5]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[1]+sum[2]+sum[3]+sum[4]+sum[5]);
		}
		// Take care of remaining values
		if(i2<e2){
			sum[6]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<e2){
			sum[0]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<e2){
			sum[1]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<e2){
			sum[2]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<e2){
			sum[3]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<e2){
			sum[4]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}
	}
	//Copy to output
	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){
		for(i2=s2,i2e=0;i2<e2;i2++,i2e++){
			out[i1*n2+i2]=tmp[i1e*n2e+i2e];
		}
	}
	free(tmp);
	return(0);
}

int mvAvg2d9(size_t n1, size_t n2, float* in, float* out){
/* mvAvg2d9(n1,n2,in,out) computes the 9x9 moving average
 * of 2d input "in" and returns it in output "out". "n1"
 * is the size of the slow dimension while "n2" is the
 * size of the fast dimension. The boundaries are not
 * touched. The array size should be at least 9-by-9,
 * this is not checked. Input & output may NOT be equal.
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float sum[9];
	float c=1.0/81.0;
	size_t i1,i2,j;

	for(i2=0;i2<n2;i2++)out[     i2]=in[     i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[  n2+i2]=in[  n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[2*n2+i2]=in[2*n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[3*n2+i2]=in[3*n2+i2]; //Copy Data
	for(i1=4;i1<n1-4;i1++){ //Loop over slow dimension
		sum[0]=in[(i1-4)*n2  ]+in[(i1-3)*n2  ]+in[(i1-2)*n2  ]+in[(i1-1)*n2  ]+in[i1*n2  ]+in[(i1+1)*n2  ]+in[(i1+2)*n2  ]+in[(i1+3)*n2  ]+in[(i1+4)*n2  ];
		sum[1]=in[(i1-4)*n2+1]+in[(i1-3)*n2+1]+in[(i1-2)*n2+1]+in[(i1-1)*n2+1]+in[i1*n2+1]+in[(i1+1)*n2+1]+in[(i1+2)*n2+1]+in[(i1+3)*n2+1]+in[(i1+4)*n2+1];
		sum[2]=in[(i1-4)*n2+2]+in[(i1-3)*n2+2]+in[(i1-2)*n2+2]+in[(i1-1)*n2+2]+in[i1*n2+2]+in[(i1+1)*n2+2]+in[(i1+2)*n2+2]+in[(i1+3)*n2+2]+in[(i1+4)*n2+2];
		sum[3]=in[(i1-4)*n2+3]+in[(i1-3)*n2+3]+in[(i1-2)*n2+3]+in[(i1-1)*n2+3]+in[i1*n2+3]+in[(i1+1)*n2+3]+in[(i1+2)*n2+3]+in[(i1+3)*n2+3]+in[(i1+4)*n2+3];
		sum[4]=in[(i1-4)*n2+4]+in[(i1-3)*n2+4]+in[(i1-2)*n2+4]+in[(i1-1)*n2+4]+in[i1*n2+4]+in[(i1+1)*n2+4]+in[(i1+2)*n2+4]+in[(i1+3)*n2+4]+in[(i1+4)*n2+4];
		sum[5]=in[(i1-4)*n2+5]+in[(i1-3)*n2+5]+in[(i1-2)*n2+5]+in[(i1-1)*n2+5]+in[i1*n2+5]+in[(i1+1)*n2+5]+in[(i1+2)*n2+5]+in[(i1+3)*n2+5]+in[(i1+4)*n2+5];
		sum[6]=in[(i1-4)*n2+6]+in[(i1-3)*n2+6]+in[(i1-2)*n2+6]+in[(i1-1)*n2+6]+in[i1*n2+6]+in[(i1+1)*n2+6]+in[(i1+2)*n2+6]+in[(i1+3)*n2+6]+in[(i1+4)*n2+6];
		sum[7]=in[(i1-4)*n2+7]+in[(i1-3)*n2+7]+in[(i1-2)*n2+7]+in[(i1-1)*n2+7]+in[i1*n2+7]+in[(i1+1)*n2+7]+in[(i1+2)*n2+7]+in[(i1+3)*n2+7]+in[(i1+4)*n2+7];
		out[i1*n2]=in[i1*n2];out[i1*n2+1]=in[i1*n2+1];out[i1*n2+2]=in[i1*n2+2];out[i1*n2+3]=in[i1*n2+3];
		for(i2=4,j=0;j<(n2-8)/9;j++){ //Loop over fast dimension - Inner 9-loop unrolled
			sum[8]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[0]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[1]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[2]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[3]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[4]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[5]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[6]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[7]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}
		// Take care of remaining values
		if(i2<n2-4){
			sum[8]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<n2-4){
			sum[0]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<n2-4){
			sum[1]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<n2-4){
			sum[2]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<n2-4){
			sum[3]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<n2-4){
			sum[4]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<n2-4){
			sum[5]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<n2-4){
			sum[6]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}
		out[(i1+1)*n2-4]=in[(i1+1)*n2-4];out[(i1+1)*n2-3]=in[(i1+1)*n2-3];out[(i1+1)*n2-2]=in[(i1+1)*n2-2];out[(i1+1)*n2-1]=in[(i1+1)*n2-1];
	}
	for(i2=0;i2<n2;i2++)out[(n1-4)*n2+i2]=in[(n1-4)*n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[(n1-3)*n2+i2]=in[(n1-3)*n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[(n1-2)*n2+i2]=in[(n1-2)*n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[(n1-1)*n2+i2]=in[(n1-1)*n2+i2]; //Copy Data
	return(0);
}

int mvAvg2d9Sgn(size_t n1, size_t n2, float* in, float* out){
/* mvAvg2d9Sgn(n1,n2,in,out) signs the input values according
 * to the sign of the local 9-by-9 average at that point. See
 * mvAvg2d7 for details.
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float sum[7];
	size_t i1,i2,j;

	for(i2=0;i2<n2;i2++)out[     i2]=in[     i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[  n2+i2]=in[  n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[2*n2+i2]=in[2*n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[3*n2+i2]=in[3*n2+i2]; //Copy Data
	for(i1=4;i1<n1-2;i1++){ //Loop over slow dimension
		sum[0]=in[(i1-3)*n2  ]+in[(i1-2)*n2  ]+in[(i1-1)*n2  ]+in[i1*n2  ]+in[(i1+1)*n2  ]+in[(i1+2)*n2  ]+in[(i1+3)*n2  ];
		sum[1]=in[(i1-3)*n2+1]+in[(i1-2)*n2+1]+in[(i1-1)*n2+1]+in[i1*n2+1]+in[(i1+1)*n2+1]+in[(i1+2)*n2+1]+in[(i1+3)*n2+1];
		sum[2]=in[(i1-3)*n2+2]+in[(i1-2)*n2+2]+in[(i1-1)*n2+2]+in[i1*n2+2]+in[(i1+1)*n2+2]+in[(i1+2)*n2+2]+in[(i1+3)*n2+2];
		sum[3]=in[(i1-3)*n2+3]+in[(i1-2)*n2+3]+in[(i1-1)*n2+3]+in[i1*n2+3]+in[(i1+1)*n2+3]+in[(i1+2)*n2+3]+in[(i1+3)*n2+3];
		sum[4]=in[(i1-3)*n2+4]+in[(i1-2)*n2+4]+in[(i1-1)*n2+4]+in[i1*n2+4]+in[(i1+1)*n2+4]+in[(i1+2)*n2+4]+in[(i1+3)*n2+4];
		sum[5]=in[(i1-3)*n2+5]+in[(i1-2)*n2+5]+in[(i1-1)*n2+5]+in[i1*n2+5]+in[(i1+1)*n2+5]+in[(i1+2)*n2+5]+in[(i1+3)*n2+5];
		out[i1*n2]=in[i1*n2];out[i1*n2+1]=in[i1*n2+1];out[i1*n2+2]=in[i1*n2+2];out[i1*n2+3]=in[i1*n2+3];
		for(i2=3,j=0;j<(n2-6)/5;j++){ //Loop over fast dimension - Inner 7-loop unrolled
			sum[8]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
			sum[0]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
			sum[1]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
			sum[2]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
			sum[3]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
			sum[4]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
			sum[5]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
			sum[6]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
			sum[7]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
		}
		// Take care of remaining values
		if(i2<n2-4){
			sum[8]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
		}else continue;
		if(i2<n2-4){
			sum[0]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
		}else continue;
		if(i2<n2-4){
			sum[1]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
		}else continue;
		if(i2<n2-4){
			sum[2]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
		}else continue;
		if(i2<n2-4){
			sum[3]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
		}else continue;
		if(i2<n2-4){
			sum[4]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
		}else continue;
		if(i2<n2-4){
			sum[5]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
		}else continue;
		if(i2<n2-4){
			sum[6]==in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			out[i1*n2+i2]=copysignf(in[i1*n2+i2],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);i2++;
		}
		out[(i1+1)*n2-4]=in[(i1+1)*n2-4];out[(i1+1)*n2-3]=in[(i1+1)*n2-3];out[(i1+1)*n2-2]=in[(i1+1)*n2-2];out[(i1+1)*n2-1]=in[(i1+1)*n2-1];
	}
	for(i2=0;i2<n2;i2++)out[(n1-4)*n2+i2]=in[(n1-4)*n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[(n1-3)*n2+i2]=in[(n1-3)*n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[(n1-2)*n2+i2]=in[(n1-2)*n2+i2]; //Copy Data
	for(i2=0;i2<n2;i2++)out[(n1-1)*n2+i2]=in[(n1-1)*n2+i2]; //Copy Data
	return(0);
}

int mvAvg2d9Embd(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2){
/* mvAvg2d9Embd(n2,in,out,s1,s2,e1,e2) is the embedded form of mvAvg2d7(n1,n2,in,out). "s1" is
 * the start of the area to average in the slow dimension. "s2" is the start of the area to
 * average in the fast dimension. "e1" is the end of the area to average in the slow dimension.
 * "e2" is the end of the area to average in the fast dimension. Input & output MAY be equal.
 * Requires additional memory. Note: s1,s2>1 & e1<n1-1 & e2<n2-1
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float sum[9];
	float *tmp;
	float c=1.0/81.0;
	size_t i1,i1e,i2,i2e,j;
	size_t n1e=e1-s1;
	size_t n2e=e2-s2;

	tmp=(float*)malloc(n1e*n2e*sizeof(float));

	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){ //Loop over slow dimension
		sum[0]=in[(i1-4)*n2-4]+in[(i1-3)*n2+s2-4]+in[(i1-2)*n2+s2-4]+in[(i1-1)*n2+s2-3]+in[i1*n2+s2-3]+in[(i1+1)*n2+s2-3]+in[(i1+2)*n2+s2-3]+in[(i1+3)*n2+s2-3]+in[(i1+4)*n2+s2-3];
		sum[1]=in[(i1-4)*n2-3]+in[(i1-3)*n2+s2-3]+in[(i1-2)*n2+s2-3]+in[(i1-1)*n2+s2-3]+in[i1*n2+s2-3]+in[(i1+1)*n2+s2-3]+in[(i1+2)*n2+s2-3]+in[(i1+3)*n2+s2-3]+in[(i1+4)*n2+s2-3];
		sum[2]=in[(i1-4)*n2-2]+in[(i1-3)*n2+s2-2]+in[(i1-2)*n2+s2-2]+in[(i1-1)*n2+s2-2]+in[i1*n2+s2-2]+in[(i1+1)*n2+s2-2]+in[(i1+2)*n2+s2-2]+in[(i1+3)*n2+s2-2]+in[(i1+4)*n2+s2-2];
		sum[3]=in[(i1-4)*n2-1]+in[(i1-3)*n2+s2-1]+in[(i1-2)*n2+s2-1]+in[(i1-1)*n2+s2-1]+in[i1*n2+s2-1]+in[(i1+1)*n2+s2-1]+in[(i1+2)*n2+s2-1]+in[(i1+3)*n2+s2-1]+in[(i1+4)*n2+s2-1];
		sum[4]=in[(i1-4)*n2  ]+in[(i1-3)*n2+s2  ]+in[(i1-2)*n2+s2  ]+in[(i1-1)*n2+s2  ]+in[i1*n2+s2  ]+in[(i1+1)*n2+s2  ]+in[(i1+2)*n2+s2  ]+in[(i1+3)*n2+s2  ]+in[(i1+4)*n2+s2  ];
		sum[5]=in[(i1-4)*n2+1]+in[(i1-3)*n2+s2+1]+in[(i1-2)*n2+s2+1]+in[(i1-1)*n2+s2+1]+in[i1*n2+s2+1]+in[(i1+1)*n2+s2+1]+in[(i1+2)*n2+s2+1]+in[(i1+3)*n2+s2+1]+in[(i1+4)*n2+s2+1];
		sum[6]=in[(i1-4)*n2+2]+in[(i1-3)*n2+s2+2]+in[(i1-2)*n2+s2+2]+in[(i1-1)*n2+s2+2]+in[i1*n2+s2+2]+in[(i1+1)*n2+s2+2]+in[(i1+2)*n2+s2+2]+in[(i1+3)*n2+s2+2]+in[(i1+4)*n2+s2+2];
		sum[7]=in[(i1-4)*n2+3]+in[(i1-3)*n2+s2+3]+in[(i1-2)*n2+s2+3]+in[(i1-1)*n2+s2+3]+in[i1*n2+s2+3]+in[(i1+1)*n2+s2+3]+in[(i1+2)*n2+s2+3]+in[(i1+3)*n2+s2+3]+in[(i1+4)*n2+s2+3];
		for(i2=s2,i2e=0,j=0;j<n2e/9;j++){ //Loop over fast dimension - Inner 9-loop unrolled
			sum[8]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[0]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[1]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[2]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[3]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[4]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[5]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[6]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[7]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}
		// Take care of remaining values
		if(i2<e2){
			sum[8]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<e2){
			sum[0]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<e2){
			sum[1]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<e2){
			sum[2]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<e2){
			sum[3]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<e2){
			sum[4]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<e2){
			sum[5]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}else continue;
		if(i2<e2){
			sum[6]=in[(i1-3)*n2+i2+3]+in[(i1-2)*n2+i2+3]+in[(i1-1)*n2+i2+3]+in[i1*n2+i2+3]+in[(i1+1)*n2+i2+3]+in[(i1+2)*n2+i2+3]+in[(i1+3)*n2+i2+3];i2++;
			tmp[i1e*n2e+i2e++]=c*(sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]);
		}
	}
	//Copy to output
	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){
		for(i2=s2,i2e=0;i2<e2;i2++,i2e++){
			out[i1*n2+i2]=tmp[i1e*n2e+i2e];
		}
	}
	free(tmp);
	return(0);
}

int mvAvg2d9EmbdSgn(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2){
/* mvAvg2d9Embd(n2,in,out,s1,s2,e1,e2) is the embedded form of mvAvg2d7(n1,n2,in,out). "s1" is
 * the start of the area to average in the slow dimension. "s2" is the start of the area to
 * average in the fast dimension. "e1" is the end of the area to average in the slow dimension.
 * "e2" is the end of the area to average in the fast dimension. Input & output MAY be equal.
 * Requires additional memory. Note: s1,s2>1 & e1<n1-1 & e2<n2-1
 *
 * AUTHOR:
 *      Max Holicki
 *      Delft University of Technology
 *      21.11.2017
 */
	float sum[9];
	float *tmp;
	size_t i1,i1e,i2,i2e,j;
	size_t n1e=e1-s1;
	size_t n2e=e2-s2;
	
	tmp=(float*)malloc(n1e*n2e*sizeof(float));

	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){ //Loop over slow dimension
		sum[0]=in[(i1-4)*n2-4]+in[(i1-3)*n2+s2-4]+in[(i1-2)*n2+s2-4]+in[(i1-1)*n2+s2-3]+in[i1*n2+s2-3]+in[(i1+1)*n2+s2-3]+in[(i1+2)*n2+s2-3]+in[(i1+3)*n2+s2-3]+in[(i1+4)*n2+s2-3];
		sum[1]=in[(i1-4)*n2-3]+in[(i1-3)*n2+s2-3]+in[(i1-2)*n2+s2-3]+in[(i1-1)*n2+s2-3]+in[i1*n2+s2-3]+in[(i1+1)*n2+s2-3]+in[(i1+2)*n2+s2-3]+in[(i1+3)*n2+s2-3]+in[(i1+4)*n2+s2-3];
		sum[2]=in[(i1-4)*n2-2]+in[(i1-3)*n2+s2-2]+in[(i1-2)*n2+s2-2]+in[(i1-1)*n2+s2-2]+in[i1*n2+s2-2]+in[(i1+1)*n2+s2-2]+in[(i1+2)*n2+s2-2]+in[(i1+3)*n2+s2-2]+in[(i1+4)*n2+s2-2];
		sum[3]=in[(i1-4)*n2-1]+in[(i1-3)*n2+s2-1]+in[(i1-2)*n2+s2-1]+in[(i1-1)*n2+s2-1]+in[i1*n2+s2-1]+in[(i1+1)*n2+s2-1]+in[(i1+2)*n2+s2-1]+in[(i1+3)*n2+s2-1]+in[(i1+4)*n2+s2-1];
		sum[4]=in[(i1-4)*n2  ]+in[(i1-3)*n2+s2  ]+in[(i1-2)*n2+s2  ]+in[(i1-1)*n2+s2  ]+in[i1*n2+s2  ]+in[(i1+1)*n2+s2  ]+in[(i1+2)*n2+s2  ]+in[(i1+3)*n2+s2  ]+in[(i1+4)*n2+s2  ];
		sum[5]=in[(i1-4)*n2+1]+in[(i1-3)*n2+s2+1]+in[(i1-2)*n2+s2+1]+in[(i1-1)*n2+s2+1]+in[i1*n2+s2+1]+in[(i1+1)*n2+s2+1]+in[(i1+2)*n2+s2+1]+in[(i1+3)*n2+s2+1]+in[(i1+4)*n2+s2+1];
		sum[6]=in[(i1-4)*n2+2]+in[(i1-3)*n2+s2+2]+in[(i1-2)*n2+s2+2]+in[(i1-1)*n2+s2+2]+in[i1*n2+s2+2]+in[(i1+1)*n2+s2+2]+in[(i1+2)*n2+s2+2]+in[(i1+3)*n2+s2+2]+in[(i1+4)*n2+s2+2];
		sum[7]=in[(i1-4)*n2+3]+in[(i1-3)*n2+s2+3]+in[(i1-2)*n2+s2+3]+in[(i1-1)*n2+s2+3]+in[i1*n2+s2+3]+in[(i1+1)*n2+s2+3]+in[(i1+2)*n2+s2+3]+in[(i1+3)*n2+s2+3]+in[(i1+4)*n2+s2+3];
		for(i2=s2,i2e=0,j=0;j<n2e/9;j++){ //Loop over fast dimension - Inner 5-loop unrolled
			sum[8]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[0]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[1]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[2]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[3]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[4]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[5]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[6]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
			sum[7]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}
		// Take care of remaining values
		if(i2<e2){
			sum[8]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<e2){
			sum[0]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<e2){
			sum[1]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<e2){
			sum[2]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<e2){
			sum[3]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<e2){
			sum[4]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<e2){
			sum[5]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}else continue;
		if(i2<e2){
			sum[6]=in[(i1-4)*n2+i2+4]+in[(i1-3)*n2+i2+4]+in[(i1-2)*n2+i2+4]+in[(i1-1)*n2+i2+4]+in[i1*n2+i2+4]+in[(i1+1)*n2+i2+4]+in[(i1+2)*n2+i2+4]+in[(i1+3)*n2+i2+4]+in[(i1+4)*n2+i2+4];
			tmp[i1e*n2e+i2e++]=copysignf(in[i1*n2+i2++],sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]);
		}
	}
	//Copy to output
	for(i1=s1,i1e=0;i1<e1;i1++,i1e++){
		for(i2=s2,i2e=0;i2<e2;i2++,i2e++){
			out[i1*n2+i2]=tmp[i1e*n2e+i2e];
		}
	}
	free(tmp);
	return(0);
}
