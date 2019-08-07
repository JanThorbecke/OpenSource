unsigned int gcd_uint(unsigned int u, unsigned int v){
	unsigned int shift, t;

	if (u == v) return u;
	if (u == 0) return v;
	if (v == 0) return u;

	/* Let shift both by K, where K is the greatest power of 2 dividing both u and v. */
	for(shift = 0;((u|v)&1)==0;++shift){u>>=1;v >>= 1;}
	while((u&1)==0) u>>=1;

	do{
		while((v&1)==0) v>>=1;
		if(u>v){t=v;v=u;u=t;}
	}while (v != 0);

	return u << shift;
}