#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"fdacrtmc.h"
#include"zfp.h"

static zfp_field* field;
static zfp_stream* zfp;
static void* buffer;
size_t bufsize;
static bitstream* stream;

static float* farr; /* Static pointer to float array with same size as migration snapshots */

int initializeCompressionField(migPar *mig){

	//Initialize Snapshot Storage Before & After (De)Compression
	farr=(float*)malloc(mig->sizem*sizeof(float)); //Allocate Snapshot Buffer

	//Initialize (De)Compression Operator
	field=zfp_field_2d(farr,zfp_type_float,(uint) mig->nz,(uint) mig->nx); //Initialize Compression Field
	zfp = zfp_stream_open(NULL);                                           //Compressed stream and parameters
//	zfp_stream_set_accuracy( zfp,0.000001);                                //Set tolerance for fixed-accuracy  mode
//	zfp_stream_set_accuracy( zfp,0.001   );                                //Set tolerance for fixed-accuracy  mode
	zfp_stream_set_precision(zfp,5       );                                //Set tolerance for fixed-precision mode
	bufsize=zfp_stream_maximum_size(zfp,field);                            //Determine conservative buffer size
	buffer=malloc(bufsize);                                                //Compressed stream storage
	stream=stream_open(buffer, bufsize);                                   //Bit stream for compression
	zfp_stream_set_bit_stream(zfp, stream);                                //Associate with compressed stream

	return(0);
}

int deallocateCompressionField(){
	zfp_field_free(field);
	zfp_stream_close(zfp);
	stream_close(stream);
	free(buffer);
	free(farr);
	return(0);
}

float* getFArrPointer(){
	return(farr);
}

size_t compressSnapshot(void** out){
	size_t zfpsize;

	//Compress
	zfp_stream_rewind(zfp);                     //Rewind stream
	zfpsize=zfp_compress(zfp,field);            //Compress Data
//	zfp_stream_flush(zfp);                      //Is this necessary

	//Copy to output
	*out=malloc(zfpsize+sizeof(size_t));        //Allocate space for compressed data.
	*((size_t*)*out)=zfpsize;                   //Store compressed size as first value.
	memcpy(*out+sizeof(size_t),buffer,zfpsize); //Then store the rest of the compressed buffer.

	return(zfpsize);
}

size_t decompresSnapshot(void* in){
	size_t zfpsize;

	//Determine Size
	zfpsize=*((size_t*)in);

	//Copy Bitstream
	memcpy(buffer,in+sizeof(size_t),zfpsize);

	//Decompess to floating point array farr
	zfp_stream_rewind(zfp);
	return(zfp_decompress(zfp,field));
}
