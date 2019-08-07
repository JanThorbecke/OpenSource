inline void transpose_scalar_block(float *A, float *B, const int lda, const int ldb, const int block_size){
	size_t i,j;

#pragma omp parallel for
    for(i=0;i<block_size;i++){
        for(j=0;j<block_size;j++){
            B[j*ldb+i]=A[i*lda+j];
        }
    }
}

inline void transpose_block(float *A, float *B, const int n, const int m, const int lda, const int ldb, const int block_size){
	size_t i,j;

#pragma omp parallel for
    for(i=0;i<n;i+=block_size){
        for(j=0;j<m;j+=block_size){
            transpose_scalar_block(&A[i*lda+j],&B[j*ldb+i],lda,ldb,block_size);
        }
    }
}