#include <time.h>
#include <sys/time.h>
#include <stdio.h>

double wallclock_time(void)
{
	struct timeval s_val;
	static struct timeval b_val;
	double time;
	static int base=0;

	gettimeofday(&s_val,0);

	if (!base) {
		b_val = s_val;
		base = 1;
		return 0.0;
	}

	time = (double)(s_val.tv_sec-b_val.tv_sec) + 
		   (double)(1e-6*((double)s_val.tv_usec-(double)b_val.tv_usec));

	return (double)time;
}

double wallclock_time_(void)
{
	return (double)wallclock_time();
}

