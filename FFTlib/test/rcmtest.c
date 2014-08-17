#include <genfft.h>
#include <time.h>

main () {

	int l,j,i,n,sign, isign, ldr, ldc;
	int N, Nmax=8193, Nlot=150, Nitcc; 
	float diff, sumi, sumr, scl;
	double t0, t1, t, bias, k;
	float *data, *c_data;
	complex *cdata;
	char *machine=getenv("MACHINE");

	t = 0.0;
	t0 = wallclock_time();
	t1 = wallclock_time();
	bias = t1-t0;
	fprintf(stderr,"bias wallclock_time = %f\n",bias);

	data   = (float *) malloc (Nlot*Nmax*sizeof(float));
	c_data = (float *) malloc (Nlot*Nmax*sizeof(float));
	cdata = (complex *) calloc (Nlot*((Nmax+2)/2),sizeof(complex));

	N = 16;
	k = 5.0;
	sign = 1;
	isign = -1;
	while (N <= Nmax) {

		/* Initialize the data */

		ldr = N;
		ldc = (N+2)/2;
		for (j=0;j<Nlot;j++) {
			for (i=0;i<N;i++) {
				c_data[j*ldr+i]   = (float)-0.1+0.5*(N/2-i);
//				c_data[j*ldr+i]   = 0.0;
			}
//			c_data[j*ldr]   = 1.0;
		}
		t = 0.0;

        /* FFT */
		for (i=0; i<100; i++) {
			for (j=0;j<ldr*Nlot;j++) data[j] = c_data[j];
			t0 = wallclock_time();
			rcmfft(&data[0], &cdata[0], N, Nlot, ldr, ldc, sign);
			crmfft(&cdata[0], &data[0], N, Nlot, ldc, ldr, isign);
			t1 = wallclock_time();
			t += t1-t0;
		}

         /* Compare the data */
		scl = 1.0/(float)N;
		for (j=0; j<Nlot; j++) {
			for (i=0;i<N;i++) {
				l = j*ldr+i;
	/*
            	fprintf(stderr,"%s: l = %d (%d,%d) data = %f C-data = %f\n", machine, l, j, i, data[l]*scl, c_data[l]);
	*/
				if (c_data[l] != 0.0) {
					diff = fabs((data[l]*scl - c_data[l]) / c_data[l]);
				}
				else {
					diff = 0.0;
				}
				if (diff >= 5.0e-3) {
					fprintf(stderr,"Bad values at %i, diff=%12.2e.\n",l,diff);
					fprintf(stderr,"data[i] = %12.2e c_data[i] = %12.2e.\n",data[l], c_data[l]);
					break;
				}
/*
*/
			}
		}

		fprintf(stderr,"%s: N = %d*%d wallclock_time = %f\n",machine,Nlot,N,t-bias);

	/* find next N valid for pfacc */

		N = pow(2.0,k);
		k += 1.0;

   	} 
}

