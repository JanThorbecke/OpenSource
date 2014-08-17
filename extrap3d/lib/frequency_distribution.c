#ifdef MPI
#include <mpi.h>
#endif
#include <assert.h>
#include <stdio.h>

int frequency_distribution(int nw_low, int nw, int npes, int pe, int maxlw, int *freq_index, int type )
{
	int my_n_local_w=0;
	int min_n_local_w, n_local_w;
	int leftover_w;
	int index=0, ipe, dpe, iw, i;

	min_n_local_w  = nw/npes;
	leftover_w     = nw % npes;

	assert( 0 <= pe && pe < npes );
	assert( 0 < nw_low );
	assert( 0 <= maxlw && maxlw <= nw );


	if (type == 0) { /* simple blocks */
		for( ipe=0; ipe<npes; ipe++ ) {
			n_local_w = min_n_local_w + ( ipe<leftover_w ? 1 : 0 );
			if( pe == ipe ) {
				my_n_local_w = n_local_w;
				for( i=0; i<my_n_local_w; i++ ) {
					freq_index[i] = nw_low+index+i;
				}
			}
			index += n_local_w;
		}
	}
	else if (type == 1) { /* round robin */
		ipe=0;
		for(iw=nw_low; iw<nw+nw_low; iw++ ) {
			if( ipe==pe ) freq_index[my_n_local_w++] = iw;
			ipe = (ipe+1)%npes;
		}
	}
	else if (type == 2)  { /* min_high */
		ipe=0;
		dpe=1;
		for( iw=nw-1+nw_low; iw>= nw_low; iw-- ) {
			if( ipe==pe ) freq_index[my_n_local_w++] = iw;
			ipe += dpe;
			if( ipe== npes || ipe==-1 ) { dpe*=-1; ipe+=dpe; }
			assert( 0 <= ipe && ipe < npes );
		}
	}
	assert( min_n_local_w <= my_n_local_w  &&  my_n_local_w <= maxlw );
	assert( leftover_w > 0  ||  min_n_local_w == my_n_local_w );

	return my_n_local_w;
}

#ifdef MPI
int check_frequency_distribution( int nw_low, int nw, int npes, int pe,
				  int maxlw, int *freq_index, int nlw )
{
	/* setup the array to check that each frequency is on 1 and only 1 processor. */
	int *frequency_found, nw_on_pe, ipe, iw, f_index, valid, nw_on_ipe;

	frequency_found = (int *)malloc( nw*sizeof(int) );

	// count the number of times a frequency is found on a processor.
	for( ipe=0; ipe<npes; ipe++ ) {
		nw_on_ipe = (pe==ipe) ? nlw : -1;
		MPI_Bcast( &nw_on_ipe, 1, MPI_INT, ipe, MPI_COMM_WORLD );
		assert( nw_on_ipe >= 0 && nw_on_ipe<= maxlw );
		if( nw_on_ipe > 0 ) {
			for( iw=0; iw<nw_on_ipe; iw++ ) {
				f_index = (pe==ipe) ? freq_index[iw] : -1;
         		MPI_Bcast( &f_index, 1, MPI_INT, ipe, MPI_COMM_WORLD );
	 			assert( nw_low<=f_index && f_index<nw+nw_low );
	 			frequency_found[f_index-nw_low]++;
      		}
		}
    }

	// check the results.
	valid = 0;
	for( iw=0; iw<nw; iw++ ) {
		valid = valid && frequency_found[iw]==1;
//		assert( frequency_found[iw] == 1 );
		fprintf(stderr,"frequency_found[%d] = %d\n", iw, frequency_found[iw]);
	}

	free(frequency_found);
	return valid;
}
#endif
