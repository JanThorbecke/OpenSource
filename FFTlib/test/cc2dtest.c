#include <genfft.h>
#include <time.h>

main () {

	int j,i,n,sign, isign;
	int N, Nmax=1024, Nitcc, Nlot=1024, ld1;
	float diff, sumi, sumr, scl;
	double t0, t1, t, bias, k;
	complex *data, *c_data;
	char machine[128];

	t = 0.0;
	t0 = wallclock_time();
	t1 = wallclock_time();
	bias = t1-t0;
	fprintf(stderr,"bias wallclock_time = %f\n",bias);

	data   = (complex *) malloc (Nlot*Nmax*sizeof(complex));
	c_data = (complex *) malloc (Nlot*Nmax*sizeof(complex));

	N = 64;
	k = log(N)/log(2)+1;
	sign = 1;
	isign = -1;
	while (N <= Nmax) {

		ld1 = N;
		/* Initialize the data */

		for (j=0;j<N;j++) {
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
			for (j=0;j<ld1*N;j++) data[j] = c_data[j];
			t0 = wallclock_time();
			cc2dfft(data, N, N, ld1, sign);
			cc2dfft(data, N, N, ld1, isign);
			t1 = wallclock_time();
			t += t1-t0;
		}

         /* Compare the data */

		scl = 1.0/(float)(N*N);
		for (i=0; i<ld1*N; i++) {
/*
            fprintf(stderr,"%s: i = %d data = %12.6e C-data = %12.6e\n", getenv("MACHINE"), i, data[i].r*scl, c_data[i].r);
            fprintf(stderr,"%s: i = %d data = %12.6e C-data = %12.6e\n", getenv("MACHINE"), i, data[i].i*scl, c_data[i].i);
*/
			if (c_data[i].r != 0.0 && c_data[i].i != 0.0) {
				diff = fabs((data[i].r*scl - c_data[i].r) / c_data[i].r);
				diff += fabs((data[i].i*scl - c_data[i].i) / c_data[i].i);
			}
			else diff = 0.0;

			if (diff >= 5.0e-3) {
				fprintf(stderr,"Bad values at %i, diff=%12.2e.\n",i,diff);
				fprintf(stderr,"data[i].r   = %12.2e.\n",data[i].r);
				fprintf(stderr,"c_data[i].r = %12.2e.\n",c_data[i].r);
				exit(1);
			}
		}

		fprintf(stderr,"%s: N = %d*%d wallclock_time = %f\n",getenv("MACHINE"),N,N,t-bias);

	/* find next N */

		N = pow(2.0,k);
		k += 1.0;
	}

/******************** NEXT PART ***************************/

	N = 128;
	Nlot = 125;
	while (Nlot <= 2000) {

		ld1 = N;
		/* Initialize the data */

		for (j=0;j<Nlot;j++) {
			for (i=0;i<N;i++) {
				c_data[j*ld1+i].r = (float)-0.1+0.5*(N/2-i);
				c_data[j*ld1+i].i = (float)0.3+j*0.01;
//				c_data[j*ld1+i].r = 0.0;
//				c_data[j*ld1+i].i = 0.0;
			}
//			c_data[j*ld1].r = 1.0;
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
			cc2dfft(data, N, Nlot, ld1, sign);
			cc2dfft(data, N, Nlot, ld1, isign);
			t1 = wallclock_time();
			t += t1-t0;
		}

         /* Compare the data */

		scl = 1.0/(float)(N*Nlot);
		for (i=0; i<ld1*Nlot; i++) {
/*
            fprintf(stderr,"%s: i = %d data = %12.6e C-data = %12.6e\n", getenv("MACHINE"), i, data[i].r*scl, c_data[i].r);
            fprintf(stderr,"%s: i = %d data = %12.6e C-data = %12.6e\n", getenv("MACHINE"), i, data[i].i*scl, c_data[i].i);
*/
			if (c_data[i].r != 0.0 && c_data[i].i != 0.0) {
				diff = fabs((data[i].r*scl - c_data[i].r) / c_data[i].r);
				diff += fabs((data[i].i*scl - c_data[i].i) / c_data[i].i);
			}
			else diff=0.0;

			if (diff >= 5.0e-3) {
				fprintf(stderr,"Bad values at %i, diff=%12.2e.\n",i,diff);
				fprintf(stderr,"data[i]   = %12.6e %12.6e.\n",data[i].r, data[i].i);
				fprintf(stderr,"c_data[i] = %12.6e %12.6e.\n",c_data[i].r, c_data[i].i);
				exit(1);
			}
		}

		fprintf(stderr,"%s: N = %d*%d wallclock_time = %f\n",getenv("MACHINE"),N,Nlot,t-bias);
	/* find next N valid for pfacc */

		Nlot *= 2;

   	} 
}

