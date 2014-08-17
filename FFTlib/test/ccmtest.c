#include <genfft.h>
#include <time.h>

main () {

	int j,i,n,sign, isign;
	int N, Nmax=4097, Nitcc, Nlot=513, ld1; 
	float diff, sumi, sumr, scl;
	double t0, t1, t, bias, k;
	complex *data, *c_data;
	char *machine=getenv("MACHINE");

	t = 0.0;
	t0 = wallclock_time();
	t1 = wallclock_time();
	bias = t1-t0;
	fprintf(stderr,"bias wallclock_time = %f\n",bias);

	data   = (complex *) malloc (Nlot*Nmax*sizeof(complex));
	c_data = (complex *) malloc (Nlot*Nmax*sizeof(complex));

	N = 256;
	k = log(N)/log(2)+1;
	sign = 1;
	isign = -1;
	while (N <= Nmax) {

		ld1 = N;
		/* Initialize the data */

		for (j=0;j<250;j++) {
			for (i=0;i<N;i++) {
				c_data[j*ld1+i].r = (float)-0.1+0.5*(N/2-i);
				c_data[j*ld1+i].i = (float)0.3+j*0.01;
			}
			for (i=N;i<ld1;i++) {
				c_data[j*ld1+i].r = (float)0.0;
				c_data[j*ld1+i].i = (float)0.0;
			}
		}
		t = 0.0;

        /* FFT */
		for (i=0; i<10; i++) {
			for (j=0;j<ld1*250;j++) data[j] = c_data[j];
			t0 = wallclock_time();
			ccmfft(data, N, 250, ld1, sign);
			ccmfft(data, N, 250, ld1, isign);
			t1 = wallclock_time();
			t += t1-t0;
		}

         /* Compare the data */

		scl = 1.0/(float)N;
		for (i=0; i<ld1*250; i++) {
/*
            fprintf(stderr,"%s: i = %d data = %f C-data = %f\n", machine,i, data[i].r*scl, c_data[i].r);
            fprintf(stderr,"%s: i = %d data = %f C-data = %f\n", machine,  i, data[i].i*scl, c_data[i].i);
*/
			if (c_data[i].r != 0.0 && c_data[i].i != 0.0) {
				diff = fabs((data[i].r*scl - c_data[i].r) / c_data[i].r);
				diff += fabs((data[i].i*scl - c_data[i].i) / c_data[i].i);
			}
			else diff = 0.0;

			if (diff >= 5.0e-3) {
				fprintf(stderr,"Bad values at %i, diff=%12.2e.\n",i,diff);
				fprintf(stderr,"data[i].r = %12.2e.\n",data[i].r);
				fprintf(stderr,"c_data[i].r = %12.2e.\n",c_data[i].r);
				exit(1);
			}
		}

		fprintf(stderr,"%s: N = %d*250 wallclock_time = %f\n",machine,N,t-bias);

		/* Initialize the data */

		for (j=0;j<Nlot;j++) {
			for (i=0;i<N;i++) {
				c_data[j*ld1+i].r = (float)-0.1+0.5*(N/2-i);
				c_data[j*ld1+i].i = (float)0.3+j*0.01;
			}
			for (i=N;i<ld1;i++) {
				c_data[j*ld1+i].r = (float)0.0;
				c_data[j*ld1+i].i = (float)0.0;
			}
		}
		t = 0.0;

        /* FFT */
		for (i=0; i<10; i++) {
			for (j=0;j<ld1*Nlot;j++) data[j] = c_data[j];
			t0 = wallclock_time();
			ccmfft(data, N, Nlot, ld1, sign);
			ccmfft(data, N, Nlot, ld1, isign);
			t1 = wallclock_time();
			t += t1-t0;
		}

         /* Compare the data */

		scl = 1.0/(float)N;
		for (i=0; i<ld1*Nlot; i++) {
/*
            fprintf(stderr,"%s: i = %d data = %f C-data = %f\n", machine, i, data[i].r*scl, c_data[i].r);
            fprintf(stderr,"%s: i = %d data = %f C-data = %f\n", machine, i, data[i].i*scl, c_data[i].i);
*/
			if (c_data[i].r != 0.0 && c_data[i].i != 0.0) {
				diff = fabs((data[i].r*scl - c_data[i].r) / c_data[i].r);
				diff += fabs((data[i].i*scl - c_data[i].i) / c_data[i].i);
			}
			else diff=0.0;

			if (diff >= 5.0e-3) {
				fprintf(stderr,"Bad values at %i, diff=%12.2e.\n",i,diff);
				fprintf(stderr,"data[i].r = %12.2e.\n",data[i].r);
				fprintf(stderr,"c_data[i].r = %12.2e.\n",c_data[i].r);
				exit(1);
			}
		}

		fprintf(stderr,"%s: N = %d*%d wallclock_time = %f\n",machine,N,Nlot,t-bias);
	/* find next N valid for pfacc */

		N += 1;
		N = pow(2.0,k);
		k += 1.0;

   	} 
}

