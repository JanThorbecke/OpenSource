#include <stdio.h>


/**
* prints an example parameter file for makemod
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/
void plotexample()
{
	printf("\nfile_base=example.su\n");
	printf("cp0=1200\n");
	printf("cs0=600\n");
	printf("ro0=1000\n");
	printf("sizex=2000\n");
	printf("sizez=2500\n");
	printf("dx=5\n");
	printf("dz=10\n\n");
	printf("gradt=1\n");
	printf("writeint=1\n");
	printf("verbose=1\n\n");

	printf("intt=def\n");
	printf("grad=150\n");
	printf("poly=2\n");
	printf("cp=1800,2500,1400,2500,3000\n");
	printf("cs=1200,1800,1000,1800,1200\n");
	printf("ro=1000\n");
	printf("x=0,500,1000,1500,2000\n");
	printf("z=100,400,150,500,100\n\n");

	printf("intt=sin\n");
	printf("var=50,100\n");
	printf("grad=0\n");
	printf("poly=0\n");
	printf("cp=2400\n");
	printf("cs=1900\n");
	printf("ro=1400\n");
	printf("x=0,500,1000,1500,2000\n");
	printf("z=800,800,900,1000,900\n\n");

	printf("intt=rough\n");
	printf("var=50,1.4,10\n");
	printf("grad=50\n");
	printf("poly=1\n");
	printf("cp=3000,3300\n");
	printf("cs=2200,2300\n");
	printf("ro=1800\n");
	printf("x=0,500,1000,1500,2000\n");
	printf("z=1400,1800,1800,1300,1500\n\n");

	printf("intt=fract\n");
	printf("var=60,30,1.6,1,1.8,10\n");
	printf("grad=0\n");
	printf("poly=0\n");
	printf("cp=3800\n");
	printf("cs=2600\n");
	printf("ro=2100\n");
	printf("x=0,2000\n");
	printf("z=2100,2100\n\n");

	printf("intt=diffr\n");
	printf("cp=100\n");
	printf("cs=100\n");
	printf("ro=100\n");
	printf("x=500,1500\n");
	printf("z=1300,1800\n\n");

	return;
}
