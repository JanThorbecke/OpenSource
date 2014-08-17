#include <genfft.h>
#include <time.h>

main () {

	int j,i,n,sign, isign;
	int N, Nmax=8192, Nitcc; 
	float diff, sumi, sumr, scl;
	double t0, t1, t, bias;
	float k;
	float *data, *c_data;
	complex *cdata;
	char *machine=getenv("MACHINE");

	t = 0.0;
	t0 = wallclock_time();
	t1 = wallclock_time();
	bias = t1-t0;
	fprintf(stderr,"bias wallclock_time = %f\n",bias);

	data   = (float *) malloc (Nmax*sizeof(float));
	c_data = (float *) malloc (Nmax*sizeof(float));
	cdata = (complex *) malloc ((Nmax+2)/2*sizeof(complex));

	N = 16;
	k = 5.0;
	sign = 1;
	isign = -1;
	while (N <= Nmax) {

		/* Initialize the data */

		for (i=0;i<N;i++) {
			c_data[i]   = (float)-0.1+0.5*(N/2-i);
//			c_data[i]   = 0.0;
		}
//		c_data[0]   = 1.0;
		t = 0.0;

        /* FFT */
		for (i=0; i<2500; i++) {
			for (j=0;j<N;j++) data[j] = c_data[j];
			t0 = wallclock_time();
			rc1fft(data, cdata, N, sign);
			cr1fft(cdata, data, N, isign);
			t1 = wallclock_time();
			t += t1-t0;
		}

         /* Compare the data */
		scl = 1.0/(float)N;
		for (i=0; i<N; i++) {
/*
            fprintf(stderr,"%s: i = %d data = %f C-data = %f\n", machine, i, data[i].r*scl, c_data[i].r);
            fprintf(stderr,"%s: i = %d data = %f C-data = %f\n", machine, i, data[i].i*scl, c_data[i].i);
*/
			if (c_data[i] != 0.0) {
				diff = fabs((data[i]*scl - c_data[i]) / c_data[i]);
			}
			else {
				diff = 0.0;
			}
			if (diff >= 5.0e-3) {
				fprintf(stderr,"Bad values at %i, diff=%12.6e.\n",i,diff);
				fprintf(stderr,"data[i] = %12.6e.\n",data[i]);
				fprintf(stderr,"c_data[i] = %12.6e.\n",c_data[i]);
				break;
			}
		}

		fprintf(stderr,"%s: N = %d wallclock_time = %f\n",machine,N,t-bias);

	/* find next N valid for pfacc */

		N += 1;
		N = pow(2.0,k);
		k += 1.0;

   	} 
}

