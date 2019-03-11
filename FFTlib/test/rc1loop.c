#include <genfft.h>
#include <time.h>

void crdft(complex *cdata, REAL *rdata, int n, int sign);
void rcdft(REAL *rdata, complex *cdata, int n, int sign);


int main () {

	int j,i,n,k,sign, isign;
	int N, Nmax=600, Nitcc; 
	float diff, sumi, sumr, scl;
	double t0, t1, t, bias;
	float *wtimes;
	float *data, *c_data;
	complex *cdata;
	char *machine=getenv("MACHINE");

	t = 0.0;
	t0 = wallclock_time();
	t1 = wallclock_time();
	bias = t1-t0;
	fprintf(stderr,"bias wallclock_time = %f\n",bias);

	wtimes = (float *) malloc ((Nmax-N+1)*sizeof(float));
	cdata = (complex *) malloc ((Nmax+2)/2*sizeof(complex));

	N = 500;
	sign = 1;
	isign = -1;
#pragma omp parallel 
#pragma omp shared(c_data, N, sign, isign, t, Nmax)
#pragma omp private(i, data, cdata, j, t0, t1)
	data   = (float *) malloc (Nmax*sizeof(float));
	c_data = (float *) malloc (Nmax*sizeof(float));
	while (N <= Nmax) {

		/* Initialize the data */

#pragma omp master
{
		for (i=0;i<N;i++) {
			c_data[i]   = (float)-0.1+0.5*(N/2-i);
		}
		t = 0.0;
}

        /* FFT */
#pragma omp for
		for (i=0; i<100; i++) {
			for (j=0;j<N;j++) data[j] = c_data[j];
			t0 = wallclock_time();
//            kiss_fftr_cfg st = kiss_fftr_alloc( N ,1 ,0,0);
//            kiss_fftri( st ,(complex*)buf,(float*)bufout );
//            kiss_fftr( st ,(float*)buf,(complex*)bufout );
//            free(st);

			rc1fft(data, cdata, N, sign);
			cr1fft(cdata, data, N, isign);
			//rcdft(data, cdata, N, sign);
			//crdft(cdata, data, N, isign);
			t1 = wallclock_time();
#pragma omp critical
{
			t += t1-t0;
}
		}

         /* Compare the data */
		scl = 1.0/(float)N;
		for (i=0; i<N; i++) {
/*
            fprintf(stderr,"%s: i = %d data = %f C-data = %f diff = %f\n", machine, i, data[i]*scl, c_data[i], c_data[i]-data[i]*scl);
*/
			if (c_data[i] != 0.0) {
				diff = fabs((data[i]*scl - c_data[i]) / c_data[i]);
			}
			else {
				diff = 0.0;
			}
			if (diff >= 5.0e-3) {
				fprintf(stderr,"Bad values at %i, diff=%12.6e.\n",i,diff);
				fprintf(stderr,"data[i] = %12.6e.\n",data[i]*scl);
				fprintf(stderr,"c_data[i] = %12.6e.\n",c_data[i]);
				break;
			}
		}
//exit(1);
		fprintf(stderr,"%s: N = %d wallclock_time = %f\n",machine,N,t-bias);
		wtimes[k++] = t-bias;

	/* find next N valid for pfacc */

		N += 1;

   	} 
	return 0;
}

