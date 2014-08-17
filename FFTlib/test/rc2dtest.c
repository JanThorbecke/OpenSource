#include <genfft.h>
#include <time.h>

main () {

	int l,j,i,n,sign, isign, ldr, ldc;
	int N, Nmax=2048, Nlot=150, Nitcc; 
	float diff, sumi, sumr, scl;
	double t0, t1, t, bias, k;
	float *data, *c_data;
	complex *cdata;
	char machine[128];

	t = 0.0;
	t0 = wallclock_time();
	t1 = wallclock_time();
	bias = t1-t0;
	fprintf(stderr,"bias wallclock_time = %f\n",bias);

	data   = (float *) malloc (Nlot*Nmax*sizeof(float));
	c_data = (float *) malloc (Nlot*Nmax*sizeof(float));
	cdata = (complex *) malloc (Nlot*((Nmax+2)/2)*sizeof(complex));

	N = 16;
	k = log(N)/log(2)+1;
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
//			c_data[j*ldr+0]   = 1.0;
		}
		t = 0.0;

        /* FFT */
		for (i=0; i<10; i++) {
			for (j=0;j<ldr*Nlot;j++) data[j] = c_data[j];
			t0 = wallclock_time();
			rc2dfft(data, cdata, N, Nlot, ldr, ldc, sign);
			cr2dfft(cdata, data, N, Nlot, ldc, ldr, isign);
			t1 = wallclock_time();
			t += t1-t0;
		}

         /* Compare the data */
		scl = 1.0/(float)(N*Nlot);
		for (j=0; j<Nlot; j++) {
			for (l=0; l<N; l++) {
				i = j*ldr+l;
	/*
            	fprintf(stderr,"%s: i = %d data = %f C-data = %f\n", getenv("MACHINE"), i, data[i]*scl, c_data[i]);
	*/
				if (c_data[i] != 0.0) {
					diff = fabs((data[i]*scl - c_data[i]) / c_data[i]);
				}
				else {
					diff = 0.0;
				}
/*
				if (diff >= 5.0e-3) {
					fprintf(stderr,"Bad values at %i, diff=%12.2e.\n",i,diff);
					fprintf(stderr,"data[i] = %12.2e c_data[i] = %12.2e.\n",data[i], c_data[i]);
					break;
				}
*/
			}
		}

		fprintf(stderr,"%s: N = %d*%d wallclock_time = %f\n",getenv("MACHINE"),Nlot,N,t-bias);

	/* find next N valid for pfacc */

		N = pow(2.0,k);
		k += 1.0;

   	} 

/******************** NEXT PART ***************************/

    N = 128;
    Nlot = 125;
    while (Nlot <= 2000) {

        ldr = N;
		ldc = (N+2)/2;
        /* Initialize the data */

        for (j=0;j<Nlot;j++) {
            for (i=0;i<N;i++) {
				c_data[j*ldr+i]   = (float)-0.1+0.5*(N/2-i);
            }
        }
        t = 0.0;

        /* FFT */
        for (i=0; i<10; i++) {
            for (j=0;j<ldr*Nlot;j++) data[j] = c_data[j];
            t0 = wallclock_time();
			rc2dfft(data, cdata, N, Nlot, ldr, ldc, sign);
			cr2dfft(cdata, data, N, Nlot, ldc, ldr, isign);
            t1 = wallclock_time();
            t += t1-t0;
        }

         /* Compare the data */

        scl = 1.0/(float)(N*Nlot);
		for (j=0; j<Nlot; j++) {
			for (l=0; l<N; l++) {
				i = j*ldr+l;
	/*
            	fprintf(stderr,"%s: i = %d data = %f C-data = %f\n", getenv("MACHINE"), i, data[i].r*scl, c_data[i].r);
            	fprintf(stderr,"%s: i = %d data = %f C-data = %f\n", getenv("MACHINE"), i, data[i].i*scl, c_data[i].i);
	*/
				if (c_data[i] != 0.0) {
					diff = fabs((data[i]*scl - c_data[i]) / c_data[i]);
				}
				else {
					diff = 0.0;
				}
				if (diff >= 5.0e-3) {
					fprintf(stderr,"Bad values at %i, diff=%12.2e.\n",i,diff);
					fprintf(stderr,"data[i] = %12.2e c_data[i] = %12.2e.\n",data[i], c_data[i]);
					break;
				}
			}
			if (diff >= 5.0e-3) break;
        }

        fprintf(stderr,"%s: N = %d*%d wallclock_time = %f\n",getenv("MACHINE"),N,Nlot,t-bias);
    /* find next N valid for pfacc */

        Nlot *= 2;

    }

}

