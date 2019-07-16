#include<complex.h>
#include<fftw3.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<errno.h>
#include"par.h"
#include"fdacrtmc.h"

static char WisdomStem[1024]=WISDOMDIR;

/* Initialize Multithreaded FFTW */
int initFFTW(){
	char *base;
	int i;
	int n,m=0;

	/* Update Wisdom Stem */
	i=omp_get_max_threads();
	n=i;while(n!=0){n/=10;++m;}
	n=strlen(WisdomStem);
	if(2+n+m>1024)return(-1);
	sprintf(&WisdomStem[n],"n%d_",i);

	/* Initialize Multithreaded FFTW */
	if(!fftw_init_threads())i=1;else i=0;
        fftw_plan_with_nthreads(omp_get_max_threads());
	return(i);
}

/* Get Wisdom Functions */
int getWisdom_1d_r2r(size_t n1){
	// Load fftw wisdom from file
	char *FileName;
	FILE *fp;
	size_t n,m=0;

	// Count number of digits
	n=n1;while(n!=0){n/=10;++m;}

	// Create Filename
	FileName=(char *)malloc((strlen(WisdomStem)+9+m)*sizeof(char));
	memcpy(FileName,WisdomStem,strlen(WisdomStem));
	sprintf(&FileName[strlen(WisdomStem)],"fft_r2r-%zu",n1);

	// Load Wisdom
	if((fp=fopen(FileName,"r"))==NULL){
		if(errno==ENOENT)return(1);else{
			free(FileName);
			return(-1);
		}
	}else{
		fftw_import_wisdom_from_file(fp);
		fclose(fp);
		free(FileName);
		return(0);
	}
}
int getWisdom_1d_r2c(size_t n1){
	// Load fftw wisdom from file
	char *FileName;
	FILE *fp;
	size_t n,m=0;

	// Count number of digits
	n=n1;while(n!=0){n/=10;++m;}

	// Create Filename
	FileName=(char *)malloc((strlen(WisdomStem)+9+m)*sizeof(char));
	memcpy(FileName,WisdomStem,strlen(WisdomStem));
	sprintf(&FileName[strlen(WisdomStem)],"fft_r2c-%zu",n1);

	// Load Wisdom
	if((fp=fopen(FileName,"r"))==NULL){
		if(errno==ENOENT)return(1);else{
			free(FileName);
			return(-1);
		}
	}else{
		fftw_import_wisdom_from_file(fp);
		fclose(fp);
		free(FileName);
		return(0);
	}
}
int getWisdom_1d_c2c(size_t n1){
	// Load fftw wisdom from file
	char *FileName;
	FILE *fp;
	size_t n,m=0;

	// Count number of digits
	n=n1;while(n!=0){n/=10;++m;}

	// Create Filename
	FileName=(char *)malloc((strlen(WisdomStem)+9+m)*sizeof(char));
	memcpy(FileName,WisdomStem,strlen(WisdomStem));
	sprintf(&FileName[strlen(WisdomStem)],"fft_c2c-%zu",n1);

	// Load Wisdom
	if((fp=fopen(FileName,"r"))==NULL){
		if(errno==ENOENT)return(1);else{
			free(FileName);
			return(-1);
		}
	}else{
		fftw_import_wisdom_from_file(fp);
		fclose(fp);
		free(FileName);
		return(0);
	}
}
int getWisdom_2d_r2c(size_t n1,size_t n2){
	// Load fftw wisdom from file
	char *FileName;
	FILE *fp;
	size_t n,m1=0,m2=0;

	// Count number of digits
	n=n1;while(n!=0){n/=10;++m1;}
	n=n2;while(n!=0){n/=10;++m2;}

	// Create Filename
	FileName=(char *)malloc((strlen(WisdomStem)+10+m1+m2)*sizeof(char));
	memcpy(FileName,WisdomStem,strlen(WisdomStem));
	sprintf(&FileName[strlen(WisdomStem)],"fft_r2c-%zu_%zu",n1,n2);

	// Load Wisdom
	if((fp=fopen(FileName,"r"))==NULL){
		if(errno==ENOENT)return(1);else{
			free(FileName);
			return(-1);
		}
	}else{
		fftw_import_wisdom_from_file(fp);
		fclose(fp);
		free(FileName);
		return(0);
	}
}
int getWisdom_2d_c2c(size_t n1,size_t n2){
	// Load fftw wisdom from file
	char *FileName;
	FILE *fp;
	size_t n,m1=0,m2=0;

	// Count number of digits
	n=n1;while(n!=0){n/=10;++m1;}
	n=n2;while(n!=0){n/=10;++m2;}

	// Create Filename
	FileName=(char *)malloc((strlen(WisdomStem)+10+m1+m2)*sizeof(char));
	memcpy(FileName,WisdomStem,strlen(WisdomStem));
	sprintf(&FileName[strlen(WisdomStem)],"fft_c2c-%zu_%zu",n1,n2);

	// Load Wisdom
	if((fp=fopen(FileName,"r"))==NULL){
		if(errno==ENOENT)return(1);else{
			free(FileName);
			return(-1);
		}
	}else{
		fftw_import_wisdom_from_file(fp);
		fclose(fp);
		free(FileName);
		return(0);
	}
}
int getWisdom_2d_c2c_1(size_t n1,size_t n2){
	// Load fftw wisdom from file
	char *FileName;
	FILE *fp;
	size_t n,m1=0,m2=0;

	// Count number of digits
	n=n1;while(n!=0){n/=10;++m1;}
	n=n2;while(n!=0){n/=10;++m2;}

	// Create Filename
	FileName=(char *)malloc((strlen(WisdomStem)+12+m1+m2)*sizeof(char));
	memcpy(FileName,WisdomStem,strlen(WisdomStem));
	sprintf(&FileName[strlen(WisdomStem)],"fft_c2c_1-%zu_%zu",n1,n2);

	// Load Wisdom
	if((fp=fopen(FileName,"r"))==NULL){
		if(errno==ENOENT)return(1);else{
			free(FileName);
			return(-1);
		}
	}else{
		fftw_import_wisdom_from_file(fp);
		fclose(fp);
		free(FileName);
		return(0);
	}
}
int getWisdom_2d_c2c_2(size_t n1,size_t n2){
	// Load fftw wisdom from file
	char *FileName;
	FILE *fp;
	size_t n,m1=0,m2=0;

	// Count number of digits
	n=n1;while(n!=0){n/=10;++m1;}
	n=n2;while(n!=0){n/=10;++m2;}

	// Create Filename
	FileName=(char *)malloc((strlen(WisdomStem)+12+m1+m2)*sizeof(char));
	memcpy(FileName,WisdomStem,strlen(WisdomStem));
	sprintf(&FileName[strlen(WisdomStem)],"fft_c2c_2-%zu_%zu",n1,n2);

	// Load Wisdom
	if((fp=fopen(FileName,"r"))==NULL){
		if(errno==ENOENT)return(1);else{
			free(FileName);
			return(-1);
		}
	}else{
		fftw_import_wisdom_from_file(fp);
		fclose(fp);
		free(FileName);
		return(0);
	}
}

/* Write Wisdom Functions */
int writeStuff(){
	char *FileName;
	FILE *fp;

	// Create Filename
	FileName=(char *)malloc((strlen(WisdomStem)+3)*sizeof(char));
	memcpy(FileName,WisdomStem,strlen(WisdomStem));
	FileName[strlen(WisdomStem)  ]='t';
	FileName[strlen(WisdomStem)+1]='m';
	FileName[strlen(WisdomStem)+2]='p';
	FileName[strlen(WisdomStem)+3]='\0';
	
	// Write Garbage
	if((fp=fopen(FileName,"w"))==NULL){
		free(FileName);
		return(-1);
	}
	fftw_export_wisdom_to_file(fp);
	fclose(fp);
	free(FileName);

	// Forget Wisdom
	fftw_forget_wisdom();
	return(0);
}
int writeWisdom_1d_r2r(size_t n1){
	// Write fftw wisdom to file
	char *FileName;
	FILE *fp;
	size_t n,m=0;

	// Count number of digits
	n=n1;while(n!=0){n/=10;++m;}
	// Create Filename
	FileName=(char *)malloc((strlen(WisdomStem)+9+m)*sizeof(char));
	memcpy(FileName,WisdomStem,strlen(WisdomStem));
	sprintf(&FileName[strlen(WisdomStem)],"fft_r2r-%zu",n1);
	// Write Wisdom
	if((fp=fopen(FileName,"w"))==NULL){
		free(FileName);
		return(-1);
	}else{
		fftw_export_wisdom_to_file(fp);
		fclose(fp);
		free(FileName);
		return(0);
	}
}
int writeWisdom_1d_r2c(size_t n1){
	// Write fftw wisdom to file
	char *FileName;
	FILE *fp;
	size_t n,m=0;

	// Count number of digits
	n=n1;while(n!=0){n/=10;++m;}
	// Create Filename
	FileName=(char *)malloc((strlen(WisdomStem)+9+m)*sizeof(char));
	memcpy(FileName,WisdomStem,strlen(WisdomStem));
	sprintf(&FileName[strlen(WisdomStem)],"fft_r2c-%zu",n1);
	// Write Wisdom
	if((fp=fopen(FileName,"w"))==NULL){
		free(FileName);
		return(-1);
	}else{
		fftw_export_wisdom_to_file(fp);
		fclose(fp);
		free(FileName);
		return(0);
	}
}
int writeWisdom_1d_c2c(size_t n1){
	// Write fftw wisdom to file
	char *FileName;
	FILE *fp;
	size_t n,m=0;

	// Count number of digits
	n=n1;while(n!=0){n/=10;++m;}
	// Create Filename
	FileName=(char *)malloc((strlen(WisdomStem)+9+m)*sizeof(char));
	memcpy(FileName,WisdomStem,strlen(WisdomStem));
	sprintf(&FileName[strlen(WisdomStem)],"fft_c2c-%zu",n1);
	// Write Wisdom
	if((fp=fopen(FileName,"w"))==NULL){
		free(FileName);
		return(-1);
	}else{
		fftw_export_wisdom_to_file(fp);
		fclose(fp);
		free(FileName);
		return(0);
	}
}
int writeWisdom_2d_r2c(size_t n1,size_t n2){
	// Write fftw wisdom to file
	char *FileName;
	FILE *fp;
	size_t n,m1=0,m2=0;

	// Count number of digits
	n=n1;while(n!=0){n/=10;++m1;}
	n=n2;while(n!=0){n/=10;++m2;}
	// Create Filename
	FileName=(char *)malloc((strlen(WisdomStem)+10+m1+m2)*sizeof(char));
	memcpy(FileName,WisdomStem,strlen(WisdomStem));
	sprintf(&FileName[strlen(WisdomStem)],"fft_r2c-%zu_%zu",n1,n2);;
	// Write Wisdom
	if((fp=fopen(FileName,"w"))==NULL){
		verr("Could not write 2d real-complex FFTw wisdom to file: %s",FileName);
		free(FileName);
		return(-1);
	}else{
		fftw_export_wisdom_to_file(fp);
		fclose(fp);
		free(FileName);
		return(0);
	}
}
int writeWisdom_2d_c2c(size_t n1,size_t n2){
	// Write fftw wisdom to file
	char *FileName;
	FILE *fp;
	size_t n,m1=0,m2=0;

	// Count number of digits
	n=n1;while(n!=0){n/=10;++m1;}
	n=n2;while(n!=0){n/=10;++m2;}
	// Create Filename
	FileName=(char *)malloc((strlen(WisdomStem)+10+m1+m2)*sizeof(char));
	memcpy(FileName,WisdomStem,strlen(WisdomStem));
	sprintf(&FileName[strlen(WisdomStem)],"fft_c2c-%zu_%zu",n1,n2);
	// Write Wisdom
	if((fp=fopen(FileName,"w"))==NULL){
		free(FileName);
		return(-1);
	}else{
		fftw_export_wisdom_to_file(fp);
		fclose(fp);
		free(FileName);
		return(0);
	}
}
int writeWisdom_2d_c2c_1(size_t n1,size_t n2){
//TODO: Maybe wrong
	// Write fftw wisdom to file
	char *FileName;
	FILE *fp;
	size_t n,m1=0,m2=0;

	// Count number of digits
	n=n1;while(n!=0){n/=10;++m1;}
	n=n2;while(n!=0){n/=10;++m2;}
	// Create Filename
	FileName=(char *)malloc((strlen(WisdomStem)+12+m1+m2)*sizeof(char));
	memcpy(FileName,WisdomStem,strlen(WisdomStem));
	sprintf(&FileName[strlen(WisdomStem)],"fft_c2c_1-%zu_%zu",n1,n2);
	// Write Wisdom
	if((fp=fopen(FileName,"w"))==NULL){
		free(FileName);
		return(-1);
	}else{
		fftw_export_wisdom_to_file(fp);
		fclose(fp);
		free(FileName);
		return(0);
	}
}
int writeWisdom_2d_c2c_2(size_t n1,size_t n2){
	// Write fftw wisdom to file
	char *FileName;
	FILE *fp;
	size_t n,m1=0,m2=0;

	// Count number of digits
	n=n1;while(n!=0){n/=10;++m1;}
	n=n2;while(n!=0){n/=10;++m2;}
	// Create Filename
	FileName=(char *)malloc((strlen(WisdomStem)+12+m1+m2)*sizeof(char));
	memcpy(FileName,WisdomStem,strlen(WisdomStem));
	sprintf(&FileName[strlen(WisdomStem)],"fft_c2c_2-%zu_%zu",n1,n2);
	// Write Wisdom
	if((fp=fopen(FileName,"w"))==NULL){
		free(FileName);
		return(-1);
	}else{
		fftw_export_wisdom_to_file(fp);
		fclose(fp);
		free(FileName);
		return(0);
	}
}

/* Plan FFT */
int PlanFFT_1d_r2r(double *in, size_t n1, fftw_plan *fftp, fftw_plan *ifftp){
	double *pDouble, *out;

	// Allocate Output Array
	out=(double *)fftw_malloc(n1*sizeof(fftw_complex));

	// Test Forward Transform
	if(getWisdom_1d_r2r(n1)) return(-1);
	*fftp=fftw_plan_r2r_1d((int)n1,in,out,FFTW_R2HC,FFTW_WISDOM_ONLY);
	if(!(*fftp)){fftw_free(out);return(-1);}

	// Allocate Temporary Input & Copy Memory
	pDouble=(double *)fftw_malloc(n1*sizeof(double));
	memcpy(pDouble,in,n1*sizeof(double));

	// Test Inverse Transform
	fftw_execute_r2r(*fftp,pDouble,out);
	*ifftp=fftw_plan_r2r_1d((int)n1,out,pDouble,FFTW_HC2R,FFTW_WISDOM_ONLY);

	// Free Arrays
	fftw_free(pDouble);fftw_free(out);
	if(!(*ifftp)) return(-1);else return(0);
}
int PlanFFT_1d_r2c(double *in, size_t n1, fftw_plan *fftp, fftw_plan *ifftp){
	double *pDouble;
	fftw_complex *out;

	// Allocate Output Array
	out=(fftw_complex *)fftw_malloc((n1/2+1)*sizeof(fftw_complex));

	// Test Forward Transform
	if(getWisdom_1d_r2c(n1)) return(-1);
	*fftp=fftw_plan_dft_r2c_1d((int)n1,in,out,FFTW_WISDOM_ONLY);
	if(!(*fftp)){fftw_free(out);return(-1);}

	// Allocate Temporary Input & Copy Memory
	pDouble=(double *)fftw_malloc(n1*sizeof(double));
	memcpy(pDouble,in,n1*sizeof(double));

	// Test Inverse Transform
	fftw_execute_dft_r2c(*fftp,pDouble,out);
	*ifftp=fftw_plan_dft_c2r_1d((int)n1,out,pDouble,FFTW_WISDOM_ONLY);

	// Free Arrays
	fftw_free(pDouble);fftw_free(out);
	if(!(*ifftp)) return(-1);else return(0);
}
int PlanFFT_1d_c2c(fftw_complex *in, size_t n1, fftw_plan *fftp, fftw_plan *ifftp){
	fftw_complex *pComplex, *out;

	// Allocate Output Array
	out=(fftw_complex *)fftw_malloc(n1*sizeof(fftw_complex));

	// Test Forward Transform
	if(getWisdom_1d_c2c(n1)) return(-1);
	*fftp=fftw_plan_dft_1d((int)n1,in,out,FFTW_FORWARD,FFTW_WISDOM_ONLY);
	if(!(*fftp)){fftw_free(out);return(-1);}

	// Allocate Temporary Input & Copy Memory
	pComplex=(fftw_complex *)fftw_malloc(n1*sizeof(fftw_complex));
	memcpy(pComplex,in,n1*sizeof(fftw_complex));

	// Test Inverse Transform
	fftw_execute_dft(*fftp,pComplex,out);
	*ifftp=fftw_plan_dft_1d((int)n1,out,pComplex,FFTW_BACKWARD,FFTW_WISDOM_ONLY);

	// Free Arrays
	fftw_free(pComplex);fftw_free(out);
	if(!(*ifftp)) return(-1);else return(0);
}
int PlanFFT_2d_r2c(double *in, size_t n1, size_t n2, fftw_plan *fftp, fftw_plan *ifftp){
	double *pDouble;
	fftw_complex *out;
	size_t n12;

	// Compute Array Size
	n12=n1*n2;

	// Allocate Output Array
	out=(fftw_complex *)fftw_malloc((n2/2+1)*n1*sizeof(fftw_complex));

	// Test Forward Transform
	if(getWisdom_2d_r2c(n1,n2)) return(-1);
	*fftp=fftw_plan_dft_r2c_2d((int)n1,(int)n2,in,out,FFTW_WISDOM_ONLY);
	if(!(*fftp)){fftw_free(out);return(-1);}

	// Allocate Temporary Input & Copy Memory
	pDouble=(double *)fftw_malloc(n12*sizeof(double));
	memcpy(pDouble,in,n12*sizeof(double));

	// Test Inverse Transform
	fftw_execute_dft_r2c(*fftp,pDouble,out);
	*ifftp=fftw_plan_dft_c2r_2d((int)n1,(int)n2,out,pDouble,FFTW_WISDOM_ONLY);

	// Free Arrays
	fftw_free(pDouble);fftw_free(out);
	if(!(*ifftp)) return(-1);else return(0);
}
int PlanFFT_2d_c2c(fftw_complex *in, size_t n1, size_t n2, fftw_plan *fftp, fftw_plan *ifftp){
	fftw_complex *pComplex, *out;
	size_t n12;

	// Compute Array Size
	n12=n1*n2;

	// Allocate Output Array
	out=(fftw_complex *)fftw_malloc(n12*sizeof(fftw_complex));

	// Test Forward Transform
	if(getWisdom_2d_c2c(n1,n2)) return(-1);
	*fftp=fftw_plan_dft_2d((int)n1,(int)n2,in,out,FFTW_FORWARD,FFTW_WISDOM_ONLY);
	if(!(*fftp)){fftw_free(out);return(-1);}

	// Allocate Temporary Input & Copy Memory
	pComplex=(fftw_complex *)fftw_malloc(n12*sizeof(fftw_complex));
	memcpy(pComplex,in,n12*sizeof(fftw_complex));

	// Test Inverse Transform
	fftw_execute_dft(*fftp,pComplex,out);
	*ifftp=fftw_plan_dft_2d((int)n1,(int)n2,out,pComplex,FFTW_BACKWARD,FFTW_WISDOM_ONLY);

	// Free Arrays
	fftw_free(pComplex);fftw_free(out);
	if(!(*ifftp)) return(-1);else return(0);
}
int PlanFFT_2d_c2c_1(fftw_complex *in, size_t n1, size_t n2, fftw_plan *fftp, fftw_plan *ifftp){
//TODO: Maybe wrong
	fftw_complex *pComplex, *out;
	int *n;
	size_t n12;

	// Compute Array Size
	n12=n1*n2;

	// Set Transform Length
	n=(int*)malloc(sizeof(int));n[0]=(int)n1;

	// Allocate Output Array
	out=(fftw_complex *)fftw_malloc(n12*sizeof(fftw_complex));

	// Test Forward Transform
	if(getWisdom_2d_c2c_1(n1,n2)) return(-1);
	*fftp=fftw_plan_many_dft(1,n,(int)n2,in,n,(int)n2,1,out,n,(int)n2,1,FFTW_FORWARD,FFTW_WISDOM_ONLY);
	if(!(*fftp)){fftw_free(out);return(-1);}

	// Allocate Temporary Input & Copy Memory
	pComplex=(fftw_complex *)fftw_malloc(n12*sizeof(fftw_complex));
	memcpy(pComplex,in,n12*sizeof(fftw_complex));

	// Test Inverse Transform
	fftw_execute_dft(*fftp,pComplex,out);
	*ifftp=fftw_plan_many_dft(1,n,(int)n2,out,n,(int)n2,1,pComplex,n,(int)n2,1,FFTW_BACKWARD,FFTW_WISDOM_ONLY);

	// Free Arrays
	fftw_free(pComplex);fftw_free(out);
	if(!(*ifftp)) return(-1);else return(0);
}
int PlanFFT_2d_c2c_2(fftw_complex *in, size_t n1, size_t n2, fftw_plan *fftp, fftw_plan *ifftp){
	fftw_complex *pComplex, *out;
	int *n;
	size_t n12;

	// Compute Array Size
	n12=n1*n2;

	// Set Transform Length
	n=(int*)malloc(sizeof(int));n[0]=(int)n2;

	// Allocate Output Array
	out=(fftw_complex *)fftw_malloc(n12*sizeof(fftw_complex));

	// Test Forward Transform
	if(getWisdom_2d_c2c_2(n1,n2)) return(-1);
	*fftp=fftw_plan_many_dft(1,n,(int)n1,in,n,1,(int)n2,out,n,1,(int)n2,FFTW_FORWARD,FFTW_WISDOM_ONLY);
	if(!(*fftp)){fftw_free(out);return(-1);}

	// Allocate Temporary Input & Copy Memory
	pComplex=(fftw_complex *)fftw_malloc(n12*sizeof(fftw_complex));
	memcpy(pComplex,in,n12*sizeof(fftw_complex));

	// Test Inverse Transform
	fftw_execute_dft(*fftp,pComplex,out);
	*ifftp=fftw_plan_many_dft(1,n,(int)n1,out,n,1,(int)n2,pComplex,n,1,(int)n2,FFTW_BACKWARD,FFTW_WISDOM_ONLY);

	// Free Arrays
	fftw_free(pComplex);fftw_free(out);
	if(!(*ifftp)) return(-1);else return(0);
}

/* Train FFT */
int TrainGarbage(){
        fftw_complex *in, *out;
        fftw_plan fftp, ifftp;
        size_t i;

        // Initialize Arrays
	in =(fftw_complex*)fftw_malloc(sizeof(fftw_complex));
	out=(fftw_complex*)fftw_malloc(sizeof(fftw_complex));

	// Seed Input Array
	*in=((double)rand()/(double)RAND_MAX)+I*((double)rand()/(double)RAND_MAX); //Initialize Training Value

	// Generate Wisdom
	fftp=fftw_plan_dft_1d((int)1,in,out,FFTW_FORWARD,FFTW_EXHAUSTIVE);
	fftw_execute_dft(fftp,in,out);
	ifftp=fftw_plan_dft_1d((int)1,out,in,FFTW_BACKWARD,FFTW_EXHAUSTIVE);
	// Write Wisdom to Disk
	writeStuff();
//	writeWisdom_1d_c2c(n1);

	// Destroy Plans & Free Wisdom
	fftw_destroy_plan(fftp);
	fftw_destroy_plan(ifftp);
	fftw_cleanup();

	// Free Arrays
	fftw_free(in);fftw_free(out);
	return(0);
}
int TrainFFTw_1d_r2r(size_t n1){
	double *in, *out;
	fftw_plan fftp, ifftp;
	size_t i;

	// Initialize Arrays
	in =(double*)fftw_malloc(n1*sizeof(double));
	out=(double*)fftw_malloc(n1*sizeof(double));

	// Seed Input Array
	for(i=0;i<n1;i++)in[i]=((double)rand()/(double)RAND_MAX); //Initialize Training Matrix

	// Generate Wisdom
	fftp=fftw_plan_r2r_1d((int)n1,in,out,FFTW_R2HC,FFTW_EXHAUSTIVE);
	fftw_execute_r2r(fftp,in,out);
	ifftp=fftw_plan_r2r_1d((int)n1,out,in,FFTW_HC2R,FFTW_EXHAUSTIVE);

	// Write Wisdom to Disk
	writeWisdom_1d_r2r(n1);

	// Destroy Plans & Free Wisdom
	fftw_destroy_plan(fftp);
	fftw_destroy_plan(ifftp);
	fftw_cleanup();

	// Free Arrays
	fftw_free(in);fftw_free(out);
	return(0);
}
int TrainFFTw_1d_r2c(size_t n1){
	double *in;
	fftw_complex *out;
	fftw_plan fftp, ifftp;
	size_t i;

	// Initialize Arrays
	in =(double*)fftw_malloc(n1*sizeof(double));
	out=(fftw_complex*)fftw_malloc(n1*sizeof(fftw_complex));

	// Seed Input Array
	for(i=0;i<n1;i++)in[i]=((double)rand()/(double)RAND_MAX); //Initialize Training Matrix

	// Generate Wisdom
	fftp=fftw_plan_dft_r2c_1d((int)n1,in,out,FFTW_EXHAUSTIVE);
	fftw_execute_dft_r2c(fftp,in,out);
	ifftp=fftw_plan_dft_c2r_1d((int)n1,out,in,FFTW_EXHAUSTIVE);

	// Write Wisdom to Disk
	writeWisdom_1d_r2c(n1);

	// Destroy Plans & Free Wisdom
	fftw_destroy_plan(fftp);
	fftw_destroy_plan(ifftp);
	fftw_cleanup();

	// Free Arrays
	fftw_free(in);fftw_free(out);
	return(0);
}
int TrainFFTw_1d_c2c(size_t n1){
	fftw_complex *in, *out;
	fftw_plan fftp, ifftp;
	size_t i;

	// Initialize Arrays
	in =(fftw_complex*)fftw_malloc(n1*sizeof(fftw_complex));
	out=(fftw_complex*)fftw_malloc(n1*sizeof(fftw_complex));

	// Seed Input Array
	for(i=0;i<n1;i++)in[i]=((double)rand()/(double)RAND_MAX)+I*((double)rand()/(double)RAND_MAX); //Initialize Training Matrix

	// Generate Wisdom
	fftp=fftw_plan_dft_1d((int)n1,in,out,FFTW_FORWARD,FFTW_EXHAUSTIVE);
	fftw_execute_dft(fftp,in,out);
	ifftp=fftw_plan_dft_1d((int)n1,out,in,FFTW_BACKWARD,FFTW_EXHAUSTIVE);

	// Write Wisdom to Disk
	writeWisdom_1d_c2c(n1);

	// Destroy Plans & Free Wisdom
	fftw_destroy_plan(fftp);
	fftw_destroy_plan(ifftp);
	fftw_cleanup();

	// Free Arrays
	fftw_free(in);fftw_free(out);
	return(0);
}
int TrainFFTw_2d_r2c(size_t n1, size_t n2){
	double *in;
	fftw_complex *out;
	fftw_plan fftp, ifftp;
	size_t i, n12;

	// Precompute Constants
	n12=n1*n2;

	// Initialize Arrays
	in =(double*)fftw_malloc(n12*sizeof(double));
	out=(fftw_complex*)fftw_malloc((n2/2+1)*n1*sizeof(fftw_complex));

	// Seed Input Array
	for(i=0;i<n12;i++)in[i]=((double)rand()/(double)RAND_MAX)+I*((double)rand()/(double)RAND_MAX); //Initialize Training Matrix

	// Generate Wisdom
	fftp=fftw_plan_dft_r2c_2d((int)n1,(int)n2,in,out,FFTW_EXHAUSTIVE);
	fftw_execute_dft_r2c(fftp,in,out);
	ifftp=fftw_plan_dft_c2r_2d((int)n1,(int)n2,out,in,FFTW_EXHAUSTIVE);

	// Write Wisdom to Disk
	writeWisdom_2d_r2c(n1,n2);

	// Destroy Plans & Free Wisdom
	fftw_destroy_plan(fftp);
	fftw_destroy_plan(ifftp);
	fftw_cleanup();

	// Free Arrays
	fftw_free(in);fftw_free(out);
	return(0);
}
int TrainFFTw_2d_c2c(size_t n1, size_t n2){
	fftw_complex *in, *out;
	fftw_plan fftp, ifftp;
	size_t i, n12;

	// Precompute Constants
	n12=n1*n2;

	// Initialize Arrays
	in =(fftw_complex*)fftw_malloc(n12*sizeof(fftw_complex));
	out=(fftw_complex*)fftw_malloc(n12*sizeof(fftw_complex));

	// Seed Input Array
	for(i=0;i<n12;i++)in[i]=((double)rand()/(double)RAND_MAX)+I*((double)rand()/(double)RAND_MAX); //Initialize Training Matrix

	// Generate Wisdom
	fftp=fftw_plan_dft_2d((int)n1,(int)n2,in,out,FFTW_FORWARD,FFTW_EXHAUSTIVE);
	fftw_execute_dft(fftp,in,out);
	ifftp=fftw_plan_dft_2d((int)n1,(int)n2,out,in,FFTW_BACKWARD,FFTW_EXHAUSTIVE);

	// Write Wisdom to Disk
	i=(size_t)writeWisdom_2d_c2c(n1,n2);

	// Destroy Plans & Free Wisdom
	fftw_destroy_plan(fftp);
	fftw_destroy_plan(ifftp);
	fftw_cleanup();

	// Free Arrays
	fftw_free(in);fftw_free(out);
	return((int)i);
}
int TrainFFTw_2d_c2c_1(size_t n1, size_t n2){
//TODO: Maybe wrong
	fftw_complex *in, *out;
	fftw_plan fftp, ifftp;
	int *n;
	size_t i, n12;

	// Precompute Constants
	n12=n1*n2;

	// Set Transform Length
	n=(int*)malloc(sizeof(int));n[0]=(int)n1;

	// Initialize Arrays
	in =(fftw_complex*)fftw_malloc(n12*sizeof(fftw_complex));
	out=(fftw_complex*)fftw_malloc(n12*sizeof(fftw_complex));

	// Seed Input Array
	for(i=0;i<n12;i++)in[i]=((double)rand()/(double)RAND_MAX)+I*((double)rand()/(double)RAND_MAX); //Initialize Training Matrix

	// Generate Wisdom
	fftp=fftw_plan_many_dft(1,n,(int)n2,in,n,(int)n2,1,out,n,(int)n2,1,FFTW_FORWARD,FFTW_EXHAUSTIVE);
	fftw_execute_dft(fftp,in,out);
	ifftp=fftw_plan_many_dft(1,n,(int)n2,in,n,(int)n2,1,out,n,(int)n2,1,FFTW_BACKWARD,FFTW_EXHAUSTIVE);

	// Write Wisdom to Disk
	writeWisdom_2d_c2c_1(n1,n2);

	// Destroy Plans & Free Wisdom
	fftw_destroy_plan(fftp);
	fftw_destroy_plan(ifftp);
	fftw_cleanup();

	// Free Arrays
	free(n);fftw_free(in);fftw_free(out);
	return(0);
}
int TrainFFTw_2d_c2c_2(size_t n1, size_t n2){
	fftw_complex *in, *out;
	fftw_plan fftp, ifftp;
	int *n;
	size_t i, n12;

	// Precompute Constants
	n12=n1*n2;

	// Set Transform Length
	n=(int*)malloc(sizeof(int));n[0]=(int)n2;

	// Initialize Arrays
	in =(fftw_complex*)fftw_malloc(n12*sizeof(fftw_complex));
	out=(fftw_complex*)fftw_malloc(n12*sizeof(fftw_complex));

	// Seed Input Array
	for(i=0;i<n12;i++)in[i]=((double)rand()/(double)RAND_MAX)+I*((double)rand()/(double)RAND_MAX); //Initialize Training Matrix

	// Generate Wisdom
	fftp=fftw_plan_many_dft(1,n,(int)n1,in,n,1,(int)n2,out,n,1,(int)n2,FFTW_FORWARD,FFTW_EXHAUSTIVE);
	fftw_execute_dft(fftp,in,out);
	ifftp=fftw_plan_many_dft(1,n,(int)n1,out,n,1,(int)n2,in,n,1,(int)n2,FFTW_BACKWARD,FFTW_EXHAUSTIVE);

	// Write Wisdom to Disk
	writeWisdom_2d_c2c_2(n1,n2);

	// Free Arrays
	free(n);
	fftw_free(in);
	fftw_free(out);

	// Destroy Plans & Free Wisdom
	fftw_destroy_plan(fftp);
	fftw_destroy_plan(ifftp);
	fftw_cleanup();

	return(0);
}

/* Initialize FFTw Plans */
int initializeFFTwPlans(fftPlansPar* fftPlans){
	/* Initialize Allocated FFTw Plans */
	// Generic Plans
	fftPlans->fft_1d_r2r   =NULL; //One  dimensional forward real         to half-complex FFT plan
	fftPlans->ifft_1d_r2r  =NULL; //One  dimensional inverse half-complex to real         FFT plan
	fftPlans->fft_1d_r2c   =NULL; //One  dimensional forward real         to      complex FFT plan
	fftPlans->ifft_1d_c2r  =NULL; //One  dimensional inverse      complex to real         FFT plan
	fftPlans->fft_1d_c2c   =NULL; //One  dimensional forward      complex to      complex FFT plan
	fftPlans->ifft_1d_c2c  =NULL; //One  dimensional inverse      complex to      complex FFT plan
	fftPlans->fft_2d_r2c   =NULL; //Two  dimensional forward real         to      complex FFT plan
	fftPlans->ifft_2d_c2r  =NULL; //Two  dimensional inverse real         to      complex FFT plan
	fftPlans->fft_2d_c2c   =NULL; //Two  dimensional forward      complex to      complex FFT plan
	fftPlans->ifft_2d_c2c  =NULL; //Two  dimensional inverse      complex to      complex FFT plan
	fftPlans->fft_2d_r2c_1 =NULL; //Slow dimension   forward real         to      complex FFT plan
	fftPlans->ifft_2d_c2r_1=NULL; //Slow dimension   inverse real         to      complex FFT plan
	fftPlans->fft_2d_c2c_1 =NULL; //Slow dimension   forward      complex to      complex FFT plan
	fftPlans->ifft_2d_c2c_1=NULL; //Slow dimension   inverse      complex to      complex FFT plan
	fftPlans->fft_2d_r2c_2 =NULL; //Fast dimension   forward real         to      complex FFT plan
	fftPlans->ifft_2d_c2r_2=NULL; //Fast dimension   inverse real         to      complex FFT plan
	fftPlans->fft_2d_c2c_2 =NULL; //Fast dimension   forward      complex to      complex FFT plan
	fftPlans->ifft_2d_c2c_2=NULL; //Fast dimension   inverse      complex to      complex FFT plan
	//////// Special Plans ////////
	// Up-Down Plane-Wave Wavefield Decomposition
	fftPlans->fft_2d_r2c_TZ    =NULL; // (t,z)  to (w,kz)   Double Fourier Transform
	fftPlans->ifft_2d_c2r_WKz  =NULL; // (w,kz) to (t,z)    Double Fourier Transform
	fftPlans->fft_2d_c2c_2_WZ  =NULL; // (w,z)   to (w,kz)  Single Fourier Transform
	fftPlans->ifft_2d_c2c_2_WKz=NULL; // (w,kz)  to (w,z)   Single Fourier Transform
	fftPlans->fft_2d_c2c_2_TZ  =NULL; // (t,z)  to (t,kz)   Single Fourier Transform
	fftPlans->ifft_2d_c2c_2_TKz=NULL; // (t,kz) to (t,z)    Single Fourier Transform
	// 2D Wavenumber Transform
	fftPlans->fft_2d_r2c_ZX    =NULL; // (z,x)   to (kz,kx) Double Fourier Transform
	fftPlans->ifft_2d_c2r_KzKx =NULL; // (kz,kx) to (z,x)   Double Fourier Transform
	// Hilbert Plans
	fftPlans->fft_1d_c2c_x  =NULL; // Forward Horizontal Hilbert Transform Plans
	fftPlans->ifft_1d_c2c_Kx=NULL; // Inverse Horizontal Hilbert Transform Plans
	fftPlans->fft_1d_c2c_z  =NULL; // Forward Vertical   Hilbert Transform Plans
	fftPlans->ifft_1d_c2c_Kz=NULL; // Inverse Vertical   Hilbert Transform Plans
	return(0);
}

/* Destroy FFTw Plans */
int destroyFFTwPlans(fftPlansPar* fftPlans){
/* Destroys Allocated FFTw Plans */
	// Generic Plans
	if(fftPlans->fft_1d_r2r   ){fftw_destroy_plan(fftPlans->fft_1d_r2r   );fftPlans->fft_1d_r2r   =NULL;} //One  dimensional forward real         to half-complex FFT plan
	if(fftPlans->ifft_1d_r2r  ){fftw_destroy_plan(fftPlans->ifft_1d_r2r  );fftPlans->ifft_1d_r2r  =NULL;} //One  dimensional inverse half-complex to real         FFT plan
	if(fftPlans->fft_1d_r2c   ){fftw_destroy_plan(fftPlans->fft_1d_r2c   );fftPlans->fft_1d_r2c   =NULL;} //One  dimensional forward real         to      complex FFT plan
	if(fftPlans->ifft_1d_c2r  ){fftw_destroy_plan(fftPlans->ifft_1d_c2r  );fftPlans->ifft_1d_c2r  =NULL;} //One  dimensional inverse      complex to real         FFT plan
	if(fftPlans->fft_1d_c2c   ){fftw_destroy_plan(fftPlans->fft_1d_c2c   );fftPlans->fft_1d_c2c   =NULL;} //One  dimensional forward      complex to      complex FFT plan
	if(fftPlans->ifft_1d_c2c  ){fftw_destroy_plan(fftPlans->ifft_1d_c2c  );fftPlans->ifft_1d_c2c  =NULL;} //One  dimensional inverse      complex to      complex FFT plan
	if(fftPlans->fft_2d_r2c   ){fftw_destroy_plan(fftPlans->fft_2d_r2c   );fftPlans->fft_2d_r2c   =NULL;} //Two  dimensional forward real         to      complex FFT plan
	if(fftPlans->ifft_2d_c2r  ){fftw_destroy_plan(fftPlans->ifft_2d_c2r  );fftPlans->ifft_2d_c2r  =NULL;} //Two  dimensional inverse real         to      complex FFT plan
	if(fftPlans->fft_2d_c2c   ){fftw_destroy_plan(fftPlans->fft_2d_c2c   );fftPlans->fft_2d_c2c   =NULL;} //Two  dimensional forward      complex to      complex FFT plan
	if(fftPlans->ifft_2d_c2c  ){fftw_destroy_plan(fftPlans->ifft_2d_c2c  );fftPlans->ifft_2d_c2c  =NULL;} //Two  dimensional inverse      complex to      complex FFT plan
	if(fftPlans->fft_2d_r2c_1 ){fftw_destroy_plan(fftPlans->fft_2d_r2c_1 );fftPlans->fft_2d_r2c_1 =NULL;} //Slow dimension   forward real         to      complex FFT plan
	if(fftPlans->ifft_2d_c2r_1){fftw_destroy_plan(fftPlans->ifft_2d_c2r_1);fftPlans->ifft_2d_c2r_1=NULL;} //Slow dimension   inverse real         to      complex FFT plan
	if(fftPlans->fft_2d_c2c_1 ){fftw_destroy_plan(fftPlans->fft_2d_c2c_1 );fftPlans->fft_2d_c2c_1 =NULL;} //Slow dimension   forward      complex to      complex FFT plan
	if(fftPlans->ifft_2d_c2c_1){fftw_destroy_plan(fftPlans->ifft_2d_c2c_1);fftPlans->ifft_2d_c2c_1=NULL;} //Slow dimension   inverse      complex to      complex FFT plan
	if(fftPlans->fft_2d_r2c_2 ){fftw_destroy_plan(fftPlans->fft_2d_r2c_2 );fftPlans->fft_2d_r2c_2 =NULL;} //Fast dimension   forward real         to      complex FFT plan
	if(fftPlans->ifft_2d_c2r_2){fftw_destroy_plan(fftPlans->ifft_2d_c2r_2);fftPlans->ifft_2d_c2r_2=NULL;} //Fast dimension   inverse real         to      complex FFT plan
	if(fftPlans->fft_2d_c2c_2 ){fftw_destroy_plan(fftPlans->fft_2d_c2c_2 );fftPlans->fft_2d_c2c_2 =NULL;} //Fast dimension   forward      complex to      complex FFT plan
	if(fftPlans->ifft_2d_c2c_2){fftw_destroy_plan(fftPlans->ifft_2d_c2c_2);fftPlans->ifft_2d_c2c_2=NULL;} //Fast dimension   inverse      complex to      complex FFT plan
	//////// Special Plans ////////
	// Up-Down Plane-Wave Wavefield Decomposition
	if(fftPlans->fft_2d_r2c_TZ    ){fftw_destroy_plan(fftPlans->fft_2d_r2c_TZ    );fftPlans->fft_2d_r2c_TZ    =NULL;} // (t,z)   to (w,kz)  Double Fourier Transform
	if(fftPlans->ifft_2d_c2r_WKz  ){fftw_destroy_plan(fftPlans->ifft_2d_c2r_WKz  );fftPlans->ifft_2d_c2r_WKz  =NULL;} // (w,kz)  to (t,z)   Double Fourier Transform
	if(fftPlans->fft_2d_c2c_2_WZ  ){fftw_destroy_plan(fftPlans->fft_2d_c2c_2_WZ  );fftPlans->fft_2d_c2c_2_WZ  =NULL;} // (w,z)   to (w,kz)  Single Fourier Transform
	if(fftPlans->ifft_2d_c2c_2_WKz){fftw_destroy_plan(fftPlans->ifft_2d_c2c_2_WKz);fftPlans->ifft_2d_c2c_2_WKz=NULL;} // (w,kz)  to (w,z)   Single Fourier Transform
	if(fftPlans->fft_2d_c2c_2_TZ  ){fftw_destroy_plan(fftPlans->fft_2d_c2c_2_TZ  );fftPlans->fft_2d_c2c_2_TZ  =NULL;} // (t,z)   to (t,kz)  Single Fourier Transform
	if(fftPlans->ifft_2d_c2c_2_TKz){fftw_destroy_plan(fftPlans->ifft_2d_c2c_2_TKz);fftPlans->ifft_2d_c2c_2_TKz=NULL;} // (t,kz)  to (t,z)   Single Fourier Transform
	// 2D Wavenumber Transform
	if(fftPlans->fft_2d_r2c_ZX    ){fftw_destroy_plan(fftPlans->fft_2d_r2c_ZX    );fftPlans->fft_2d_r2c_ZX    =NULL;} // (z,x)   to (kz,kx) Double Fourier Transform
	if(fftPlans->ifft_2d_c2r_KzKx ){fftw_destroy_plan(fftPlans->ifft_2d_c2r_KzKx );fftPlans->ifft_2d_c2r_KzKx =NULL;} // (kz,kx) to (z,x)   Double Fourier Transform
	// Hilbert Plans
	if(fftPlans->fft_1d_c2c_x  ){fftw_destroy_plan(fftPlans->fft_1d_c2c_x  );fftPlans->fft_1d_c2c_x  =NULL;} // Forward Horizontal Hilbert Transform Plans
	if(fftPlans->ifft_1d_c2c_Kx){fftw_destroy_plan(fftPlans->ifft_1d_c2c_Kx);fftPlans->ifft_1d_c2c_Kx=NULL;} // Inverse Horizontal Hilbert Transform Plans
	if(fftPlans->fft_1d_c2c_z  ){fftw_destroy_plan(fftPlans->fft_1d_c2c_z  );fftPlans->fft_1d_c2c_z  =NULL;} // Forward Vertical   Hilbert Transform Plans
	if(fftPlans->ifft_1d_c2c_Kz){fftw_destroy_plan(fftPlans->ifft_1d_c2c_Kz);fftPlans->ifft_1d_c2c_Kz=NULL;} // Inverse Vertical   Hilbert Transform Plans
	return(0);
}

/* Global Planning Functions */
int CreateUDPlaneWaveImagingFFTPlans(fftPlansPar *fftPlans,size_t nt,size_t nz){
	double *inD;
	fftw_complex *inC;
	size_t i, ntz;

	// Precompute Constants
	ntz=nt*nz;

	// Allocate Arrays
	inD=(double*)fftw_malloc(ntz*sizeof(double));
	for(i=0;i<ntz;i++)inD[i]=((double)rand()/(double)RAND_MAX); //Initialize Training Matrix

	// Plan 2D Real to Complex FFT
	if(PlanFFT_2d_r2c(inD,nt,nz,&(fftPlans->fft_2d_r2c_TZ),&(fftPlans->ifft_2d_c2r_WKz))){ //Check if Plan can Be Created
		destroyFFTwPlans(fftPlans);
		fftw_cleanup();
		if(TrainFFTw_2d_r2c(nt,nz))verr("Could not save FFTw wisdom, WisdomStem not set correctly?");
		if(PlanFFT_2d_r2c(inD,nt,nz,&(fftPlans->fft_2d_r2c_TZ),&(fftPlans->ifft_2d_c2r_WKz)))verr("Could not load generated FFTw wisdom!");
	}

	// Free Arrays
	fftw_free(inD);
	return(0);
}
int Create2DWavenumberTransformPlan(fftPlansPar *fftPlans,size_t nz,size_t nx){
	double *inD;
	size_t i, nxz;

	// Allocate Arrays
	nxz=nx*nz;
	inD=(double*)fftw_malloc(nxz*sizeof(double));
	for(i=0;i<nxz;i++)inD[i]=((double)rand()/(double)RAND_MAX); //Initialize Training Matrix

	// Plan 2D Real to Complex FFT
	if(PlanFFT_2d_r2c(inD,nx,nz,&(fftPlans->fft_2d_r2c_ZX),&(fftPlans->ifft_2d_c2r_KzKx))){ //Check if Plan can Be Created
		destroyFFTwPlans(fftPlans);
		fftw_cleanup();
		if(TrainFFTw_2d_r2c(nx,nz))verr("Could not save FFTw wisdom, WisdomStem not set correctly?");
		if(PlanFFT_2d_r2c(inD,nx,nz,&(fftPlans->fft_2d_r2c_ZX),&(fftPlans->ifft_2d_c2r_KzKx)))verr("Could not load generated FFTw wisdom!");
	}
	// Free Arrays
	free(inD);
	return(0);
}
int Create1DWavenumberFrequencyTransformPlans(fftPlansPar *fftPlans,size_t nx,size_t nz,size_t nt){
	fftw_complex *in;
	size_t i, n;

	// Find largest dimension
	if(nt<nx)n=nx;else n=nt;if(n<nz)n=nz;

	/*******************/
	/* Plan Transforms */
	/*******************/
	in=(fftw_complex*)fftw_malloc(n*sizeof(fftw_complex));
	for(i=0;i<nt;i++)in[i]=((double)rand()/(double)RAND_MAX); //Initialize Training Matrix
	// Plan 1D Complex to Complex FFTs
	if(PlanFFT_1d_c2c(in,nx,&(fftPlans->fft_1d_c2c_x),&(fftPlans->ifft_1d_c2c_Kx))){ //Check if Plan can Be Created
		destroyFFTwPlans(fftPlans);
		fftw_cleanup();
		if(TrainFFTw_1d_c2c(nx))verr("Could not save FFTw wisdom, WisdomStem not set correctly?");
		if(PlanFFT_1d_c2c(in,nx,&(fftPlans->fft_1d_c2c_x),&(fftPlans->ifft_1d_c2c_Kx)))verr("Could not load generated FFTw wisdom!");
	}
	if(PlanFFT_1d_c2c(in,nz,&(fftPlans->fft_1d_c2c_z),&(fftPlans->ifft_1d_c2c_Kz))){ //Check if Plan can Be Created
		destroyFFTwPlans(fftPlans);
		fftw_cleanup();
		if(TrainFFTw_1d_c2c(nz))verr("Could not save FFTw wisdom, WisdomStem not set correctly?");
		if(PlanFFT_1d_c2c(in,nx,&(fftPlans->fft_1d_c2c_x),&(fftPlans->ifft_1d_c2c_Kx)))verr("Could not load generated FFTw wisdom!");
		if(PlanFFT_1d_c2c(in,nz,&(fftPlans->fft_1d_c2c_z),&(fftPlans->ifft_1d_c2c_Kz)))verr("Could not load generated FFTw wisdom!");
	}
	if(PlanFFT_1d_c2c(in,nt,&(fftPlans->fft_1d_c2c_t),&(fftPlans->ifft_1d_c2c_W))){ //Check if Plan can Be Created
		destroyFFTwPlans(fftPlans);
		fftw_cleanup();
		if(TrainFFTw_1d_c2c(nt))verr("Could not save FFTw wisdom, WisdomStem not set correctly?");
		if(PlanFFT_1d_c2c(in,nx,&(fftPlans->fft_1d_c2c_x),&(fftPlans->ifft_1d_c2c_Kx)))verr("Could not load generated FFTw wisdom!");
		if(PlanFFT_1d_c2c(in,nz,&(fftPlans->fft_1d_c2c_z),&(fftPlans->ifft_1d_c2c_Kz)))verr("Could not load generated FFTw wisdom!");
		if(PlanFFT_1d_c2c(in,nt,&(fftPlans->fft_1d_c2c_t),&(fftPlans->ifft_1d_c2c_W )))verr("Could not load generated FFTw wisdom!");
	}
	
	free(in);
	return(0);
}
int Create1DWavenumberTransformPlans(fftPlansPar *fftPlans,size_t nx,size_t nz){
	fftw_complex *in;
	size_t i, n;

	// Find largest dimension
	if(nz<nx)n=nx;else n=nz;

	/*******************/
	/* Plan Transforms */
	/*******************/
	in=(fftw_complex*)fftw_malloc(n*sizeof(fftw_complex));
	for(i=0;i<n;i++)in[i]=((double)rand()/(double)RAND_MAX); //Initialize Training Matrix
	// TODO: I have to plan at least once before loading wisdom.
	// Don't know why I need the following line to properly read stuffa
	// Hence I plan a 1 element FFT and forget the wisdom.
	TrainGarbage();

	// Plan 1D Complex to Complex FFTs
	if(PlanFFT_1d_c2c(in,nx,&(fftPlans->fft_1d_c2c_x),&(fftPlans->ifft_1d_c2c_Kx))){
		destroyFFTwPlans(fftPlans);
		fftw_forget_wisdom();
		fftw_cleanup();
		if(TrainFFTw_1d_c2c(nx))verr("Could not save FFTw wisdom, WisdomStem not set correctly?");
		if(PlanFFT_1d_c2c(in,nx,&(fftPlans->fft_1d_c2c_x),&(fftPlans->ifft_1d_c2c_Kx)))verr("Could not load generated FFTw wisdom!");
	}
	if(PlanFFT_1d_c2c(in,nz,&(fftPlans->fft_1d_c2c_z),&(fftPlans->ifft_1d_c2c_Kz))){ //Check if Plan can Be Created
		destroyFFTwPlans(fftPlans);
		fftw_forget_wisdom();
		fftw_cleanup();
		if(TrainFFTw_1d_c2c(nz))verr("Could not save FFTw wisdom, WisdomStem not set correctly?");
		if(PlanFFT_1d_c2c(in,nx,&(fftPlans->fft_1d_c2c_x),&(fftPlans->ifft_1d_c2c_Kx)))verr("Could not load generated FFTw wisdom!");
		if(PlanFFT_1d_c2c(in,nz,&(fftPlans->fft_1d_c2c_z),&(fftPlans->ifft_1d_c2c_Kz)))verr("Could not load generated FFTw wisdom!");
	}

	fftw_free(in);
	return(0);
}

/* Wisdom Path Related Functions */
int setWisdomPath(char *path){
	if(strlen(path)+1>sizeof(WisdomStem)){
		vmess("Error not able to change wisdom stem %s",WisdomStem);
		vmess("to new wisdom stem %s",path);
		verr(" Path to long, maximum path length is %zu!",sizeof(WisdomStem));
		return(-1);
	}
	strcpy(WisdomStem,path);
	return(0);
}
void printWisdomPath(){
	vmess("Wisdom Stem: %s",WisdomStem);
}
