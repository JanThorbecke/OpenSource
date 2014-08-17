#include <genfft.h>
#include <time.h>

main () {

	int j,i,n,sign, isign;
	int N, Nmax=8192, Nitcc; 
	float diff, sumi, sumr, scl;
	double t0, t1, t, bias, k;
	complex *data, *c_data;
	char *machine=getenv("MACHINE");

	t = 0.0;
	t0 = wallclock_time();
	t1 = wallclock_time();
	bias = t1-t0;
	fprintf(stderr,"bias wallclock_time = %f\n",bias);

	data   = (complex *) malloc (Nmax*sizeof(complex));
	c_data = (complex *) malloc (Nmax*sizeof(complex));

	N = 8;
	k = 3.0;
	sign = 1;
	isign = -1;
	while (N <= Nmax) {

		/* Initialize the data */

		for (i=0;i<N;i++) {
			c_data[i].r   = (float)-0.1+0.5*(N/2-i);
			c_data[i].i   = (float)0.3;
		}
		t = 0.0;

        /* FFT */
		for (i=0; i<2500; i++) {
			for (j=0;j<N;j++) data[j] = c_data[j];
			t0 = wallclock_time();
			cc1fft(data, N, sign);
			cc1fft(data, N, isign);
			t1 = wallclock_time();
			t += t1-t0;
		}

/*
		if (NINT(pow(2.0, (double)NINT(log((double)N)/log(2.0)))) != N) {
		if (npfa(N) == N) fprintf(stderr,"SU_fft ");
		else fprintf(stderr,"dft ");
		}
		else fprintf(stderr,"Mayer ");
*/


         /* Compare the data */
		scl = 1.0/(float)N;
		for (i=0; i<N; i++) {
/*
            fprintf(stderr,"%s: i = %d data = %f C-data = %f\n", machine, i, data[i].r*scl, c_data[i].r);
            fprintf(stderr,"%s: i = %d data = %f C-data = %f\n", machine, i, data[i].i*scl, c_data[i].i);
*/
			diff = fabs((data[i].r*scl - c_data[i].r) );
			diff += fabs((data[i].i*scl - c_data[i].i) );
			diff = diff/(c_data[i].r+c_data[i].i);
			if (diff >= 5.0e-3) {
				fprintf(stderr,"Bad values at %i, diff=%12.2e.\n",i,diff);
				fprintf(stderr,"data[i].r = %12.6e %12.6e.\n",data[i].r*scl, data[i].i*scl);
				fprintf(stderr,"c_data[i].r = %12.6e %12.6e.\n",c_data[i].r, c_data[i].i);
				exit(1);
			}
		}

		fprintf(stderr,"%s: N = %d wallclock_time = %f\n",machine,N,t-bias);

	/* find next N valid for pfacc */

		N += 1;
		N = pow(2.0,k);
		k += 1.0;

   	} 
}

