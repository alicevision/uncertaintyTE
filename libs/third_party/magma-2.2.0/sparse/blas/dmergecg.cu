/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @generated from sparse/blas/zmergecg.cu, normal z -> d, Sun Nov 20 20:20:43 2016
       @author Hartwig Anzt

*/
#include "magmasparse_internal.h"

#define BLOCK_SIZE 512

#define PRECISION_d


// These routines merge multiple kernels from dmergecg into one
// for a description see 
// "Reformulated Conjugate Gradient for the Energy-Aware 
// Solution of Linear Systems on GPUs (ICPP '13)

// accelerated reduction for one vector
__global__ void
magma_dcgreduce_kernel_spmv1( 
    int Gs,
    int n, 
    double * vtmp,
    double * vtmp2 )
{
    extern __shared__ double temp[];    
    int Idx = threadIdx.x;
    int blockSize = 128;
    int gridSize = blockSize  * 2 * gridDim.x; 
    temp[Idx] = MAGMA_D_MAKE( 0.0, 0.0);
    int i = blockIdx.x * ( blockSize * 2 ) + Idx;   
    while (i < Gs ) {
        temp[ Idx  ] += vtmp[ i ]; 
        temp[ Idx  ] += ( i + blockSize < Gs ) ? vtmp[ i + blockSize ] 
                                                : MAGMA_D_MAKE( 0.0, 0.0); 
        i += gridSize;
    }
    __syncthreads();
    if ( Idx < 64 ) {
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ) {
            temp[ Idx ] += temp[ Idx + 32 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 8  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 4  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 2  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 1  ]; __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ) {
            volatile double *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ) {
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    if ( Idx == 0 ) {
        vtmp2[ blockIdx.x ] = temp[ 0 ];
    }
}


// accelerated reduction for two vectors
__global__ void
magma_dcgreduce_kernel_spmv2( 
    int Gs,
    int n, 
    double * vtmp,
    double * vtmp2 )
{
    extern __shared__ double temp[];     
    int Idx = threadIdx.x;
    int blockSize = 128;
    int gridSize = blockSize  * 2 * gridDim.x; 
    int j;

    for( j=0; j<2; j++){
        int i = blockIdx.x * ( blockSize * 2 ) + Idx;   
        temp[Idx+j*(blockSize)] = MAGMA_D_ZERO;
        while (i < Gs ) {
            temp[ Idx+j*(blockSize)  ] += vtmp[ i+j*n ]; 
            temp[ Idx+j*(blockSize)  ] += 
                ( i + (blockSize) < Gs ) ? vtmp[ i+j*n + (blockSize) ] 
                                                : MAGMA_D_ZERO;
            i += gridSize;
        }
    }
    __syncthreads();
    if ( Idx < 64 ){
        for( j=0; j<2; j++){
            temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 64 ];
        }
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            for( j=0; j<2; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 32 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 16 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 8 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 4 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 2 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*(blockSize) ] += temp[ Idx+j*(blockSize) + 1 ];
                __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile double *temp2 = temp;
            for( j=0; j<2; j++){
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 32 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 16 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 8 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 4 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 2 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 1 ];
            }
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            for( j=0; j<2; j++){
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 32 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 16 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 8 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 4 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 2 ];
                temp2[ Idx+j*(blockSize) ] += temp2[ Idx+j*(blockSize) + 1 ];
            }
        }
    #endif
    if ( Idx == 0 ){
        for( j=0; j<2; j++){
            vtmp2[ blockIdx.x+j*n ] = temp[ j*(blockSize) ];
        }
    }
}



// computes the SpMV using CSR and the first step of the reduction
__global__ void
magma_dcgmerge_spmvcsr_kernel(  
    int n,
    double * dval, 
    magma_index_t * drowptr, 
    magma_index_t * dcolind,
    double * d,
    double * z,
    double * vtmp )
{
    extern __shared__ double temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;
    int j;

    temp[ Idx ] = MAGMA_D_MAKE( 0.0, 0.0);

    if( i<n ) {
        double dot = MAGMA_D_ZERO;
        int start = drowptr[ i ];
        int end = drowptr[ i+1 ];
        for( j=start; j<end; j++)
            dot += dval[ j ] * d[ dcolind[j] ];
        z[ i ] =  dot;
        temp[ Idx ] =  d[ i ] * dot;
    }

    __syncthreads();
    if ( Idx < 128 ) {
        temp[ Idx ] += temp[ Idx + 128 ];
    }
    __syncthreads();
    if ( Idx < 64 ) {
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ) {
            temp[ Idx ] += temp[ Idx + 32 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 8  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 4  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 2  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 1  ]; __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ) {
            volatile double *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ) {
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif

    if ( Idx == 0 ) {
            vtmp[ blockIdx.x ] = temp[ 0 ];
    }
}

// computes the SpMV using ELL and the first step of the reduction
__global__ void
magma_dcgmerge_spmvell_kernel(  
    int n,
    int num_cols_per_row,
    double * dval, 
    magma_index_t * dcolind,
    double * d,
    double * z,
    double * vtmp )
{
    extern __shared__ double temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;

    temp[ Idx ] = MAGMA_D_MAKE( 0.0, 0.0);

    if(i < n ) {
        double dot = MAGMA_D_MAKE(0.0, 0.0);
        for ( int k = 0; k < num_cols_per_row; k++ ) {
            int col = dcolind [ n * k + i ];
            double val = dval [ n * k + i ];
            if( val != 0)
                dot += val * d[ col ];
        }
        z[ i ] =  dot;
        temp[ Idx ] = d[ i ] * dot;
    }

    __syncthreads();
    if ( Idx < 128 ) {
        temp[ Idx ] += temp[ Idx + 128 ];
    }
    __syncthreads();
    if ( Idx < 64 ) {
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ) {
            temp[ Idx ] += temp[ Idx + 32 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 8  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 4  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 2  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 1  ]; __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ) {
            volatile double *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ) {
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif

    if ( Idx == 0 ) {
            vtmp[ blockIdx.x ] = temp[ 0 ];
    }
}


// computes the SpMV using ELLPACK and the first step of the reduction
__global__ void
magma_dcgmerge_spmvellpack_kernel(  
    int n,
    int num_cols_per_row,
    double * dval, 
    magma_index_t * dcolind,
    double * d,
    double * z,
    double * vtmp )
{
    extern __shared__ double temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;

    temp[ Idx ] = MAGMA_D_MAKE( 0.0, 0.0);

    if(i < n ) {
        double dot = MAGMA_D_MAKE(0.0, 0.0);
        for ( int k = 0; k < num_cols_per_row; k++ ) {
            int col = dcolind [ num_cols_per_row * i + k ];
            double val = dval [ num_cols_per_row * i + k ];
            if( val != 0)
                dot += val * d[ col ];
        }
        z[ i ] =  dot;
        temp[ Idx ] = d[ i ] * dot;
    }

    __syncthreads();
    if ( Idx < 128 ) {
        temp[ Idx ] += temp[ Idx + 128 ];
    }
    __syncthreads();
    if ( Idx < 64 ) {
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ) {
            temp[ Idx ] += temp[ Idx + 32 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 8  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 4  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 2  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 1  ]; __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ) {
            volatile double *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ) {
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif

    if ( Idx == 0 ) {
            vtmp[ blockIdx.x ] = temp[ 0 ];
    }
}


// computes the SpMV using SELL alignment 1 and the first step of the reduction
__global__ void
magma_dcgmerge_spmvell_kernelb1(  
    int n,
    int blocksize,
    double * dval, 
    magma_index_t * dcolind,
    magma_index_t * drowptr,
    double * d,
    double * z,
    double * vtmp )
{
    extern __shared__ double temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;

    temp[ Idx ] = MAGMA_D_MAKE( 0.0, 0.0);
    
    int idx = threadIdx.x;      // local row
    int bdx = blockIdx.x; // global block index
    int row = bdx * 256 + idx;  // global row index
    // int lblocksize = ( row + blocksize < num_rows) ? blocksize : ( num_rows - blocksize * (row/blocksize) );
    int lrow = threadIdx.x%blocksize; // local row;
    
    if( row < n ) {
        int offset = drowptr[ row/blocksize ];
        int border = (drowptr[ row/blocksize+1 ]-offset)/blocksize;
    
        double dot = MAGMA_D_MAKE(0.0, 0.0);
        for ( int n = 0; n < border; n++) { 
            int col = dcolind [ offset+ blocksize * n + lrow ];
            double val = dval[ offset+ blocksize * n + lrow ];
            dot = dot + val * d [ col ];
        }
        z[ i ] = dot;
        temp[ Idx ] = d[ i ] * dot;
    }
    
/*
    if(i < n ) {
        int offset = drowptr[ blockIdx.x ];
        int border = (drowptr[ blockIdx.x+1 ]-offset)/blocksize;
        double dot = MAGMA_D_MAKE(0.0, 0.0);
        for ( int k = 0; k < border; k++){ 
            int col = dcolind [ offset+ blocksize * k + threadIdx.x ];
            double val = dval[offset+ blocksize * k + threadIdx.x];
            if( val != 0){
                  dot += val*d[col];
            }
        }
        
        
        //double dot = MAGMA_D_MAKE(0.0, 0.0);
        //for ( int k = 0; k < num_cols_per_row; k++ ) {
        //    int col = dcolind [ n * k + i ];
        //    double val = dval [ n * k + i ];
        //    if( val != 0)
        //        dot += val * d[ col ];
        //}
        z[ i ] =  dot;
        temp[ Idx ] = d[ i ] * dot;
    }*/

    __syncthreads();
    if ( Idx < 128 ) {
        temp[ Idx ] += temp[ Idx + 128 ];
    }
    __syncthreads();
    if ( Idx < 64 ) {
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ) {
            temp[ Idx ] += temp[ Idx + 32 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 8  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 4  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 2  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 1  ]; __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ) {
            volatile double *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ) {
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif

    if ( Idx == 0 ) {
            vtmp[ blockIdx.x ] = temp[ 0 ];
    }
}


// computes the SpMV using ELLRT 8 threads per row
__global__ void
magma_dcgmerge_spmvellpackrt_kernel_8(  
    int n,
    double * dval, 
    magma_index_t * dcolind,
    magma_index_t * drowlength,
    double * d,
    double * z,
    double * vtmp,
    magma_int_t T, 
    magma_int_t alignment  )
{
    int idx = blockIdx.y * gridDim.x * blockDim.x + 
              blockDim.x * blockIdx.x + threadIdx.x;  // global thread index
    int idb = threadIdx.x;   // local thread index
    int idp = idb%T;  // number of threads assigned to one row
    int i = idx/T;  // row index
    
    extern __shared__ double shared[];

    if(i < n ) {
        double dot = MAGMA_D_MAKE(0.0, 0.0);
        int max_ = magma_ceildiv( drowlength[i], T );  
            // number of elements each thread handles

        for ( int k = 0; k < max_; k++ ) {
            // original code in paper (not working for me)
            //double val = dval[ k*(T*alignment)+(i*T)+idp ];  
            //int col = dcolind [ k*(T*alignment)+(i*T)+idp ];    

            // new code (working for me)        
            double val = dval[ k*(T)+(i*alignment)+idp ];
            int col = dcolind [ k*(T)+(i*alignment)+idp ];

            dot += val * d[ col ];
        }
        shared[idb]  = dot;
        if( idp < 4 ) {
            shared[idb]+=shared[idb+4];
            if( idp < 2 ) shared[idb]+=shared[idb+2];
            if( idp == 0 ) {
                z[i] = (shared[idb]+shared[idb+1]);
            }
        }
    }
}

// computes the SpMV using ELLRT 8 threads per row
__global__ void
magma_dcgmerge_spmvellpackrt_kernel_16(  
    int n,
    double * dval, 
    magma_index_t * dcolind,
    magma_index_t * drowlength,
    double * d,
    double * z,
    double * vtmp,
    magma_int_t T, 
    magma_int_t alignment  )
{
    int idx = blockIdx.y * gridDim.x * blockDim.x + 
              blockDim.x * blockIdx.x + threadIdx.x;  // global thread index
    int idb = threadIdx.x;   // local thread index
    int idp = idb%T;  // number of threads assigned to one row
    int i = idx/T;  // row index
    
    extern __shared__ double shared[];

    if(i < n ) {
        double dot = MAGMA_D_MAKE(0.0, 0.0);
        int max_ = magma_ceildiv( drowlength[i], T );  
            // number of elements each thread handles

        for ( int k = 0; k < max_; k++ ) {
            // original code in paper (not working for me)
            //double val = dval[ k*(T*alignment)+(i*T)+idp ];  
            //int col = dcolind [ k*(T*alignment)+(i*T)+idp ];    

            // new code (working for me)        
            double val = dval[ k*(T)+(i*alignment)+idp ];
            int col = dcolind [ k*(T)+(i*alignment)+idp ];

            dot += val * d[ col ];
        }
        shared[idb]  = dot;
        if( idp < 8 ) {
            shared[idb]+=shared[idb+8];
            if( idp < 4 ) shared[idb]+=shared[idb+4];
            if( idp < 2 ) shared[idb]+=shared[idb+2];
            if( idp == 0 ) {
                z[i] = (shared[idb]+shared[idb+1]);
            }
        }
    }
}

// computes the SpMV using ELLRT 8 threads per row
__global__ void
magma_dcgmerge_spmvellpackrt_kernel_32(  
    int n,
    double * dval, 
    magma_index_t * dcolind,
    magma_index_t * drowlength,
    double * d,
    double * z,
    double * vtmp,
    magma_int_t T, 
    magma_int_t alignment  )
{
    int idx = blockIdx.y * gridDim.x * blockDim.x + 
              blockDim.x * blockIdx.x + threadIdx.x;  // global thread index
    int idb = threadIdx.x;   // local thread index
    int idp = idb%T;  // number of threads assigned to one row
    int i = idx/T;  // row index
    
    extern __shared__ double shared[];

    if(i < n ) {
        double dot = MAGMA_D_MAKE(0.0, 0.0);
        int max_ = magma_ceildiv( drowlength[i], T );  
            // number of elements each thread handles

        for ( int k = 0; k < max_; k++ ) {
            // original code in paper (not working for me)
            //double val = dval[ k*(T*alignment)+(i*T)+idp ];  
            //int col = dcolind [ k*(T*alignment)+(i*T)+idp ];    

            // new code (working for me)        
            double val = dval[ k*(T)+(i*alignment)+idp ];
            int col = dcolind [ k*(T)+(i*alignment)+idp ];

            dot += val * d[ col ];
        }
        shared[idb]  = dot;
        if( idp < 16 ) {
            shared[idb]+=shared[idb+16];
            if( idp < 8 ) shared[idb]+=shared[idb+8];
            if( idp < 4 ) shared[idb]+=shared[idb+4];
            if( idp < 2 ) shared[idb]+=shared[idb+2];
            if( idp == 0 ) {
                z[i] = (shared[idb]+shared[idb+1]);
            }
        }
    }
}


// additional kernel necessary to compute first reduction step
__global__ void
magma_dcgmerge_spmvellpackrt_kernel2(  
    int n,
    double * z,
    double * d,
    double * vtmp2 )
{
    extern __shared__ double temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;
    

    temp[ Idx ] = ( i < n ) ? z[i]*d[i] : MAGMA_D_MAKE(0.0, 0.0);
    __syncthreads();
    if ( Idx < 128 ) {
        temp[ Idx ] += temp[ Idx + 128 ];
    }
    __syncthreads();
    if ( Idx < 64 ) {
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ) {
            temp[ Idx ] += temp[ Idx + 32 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 8  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 4  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 2  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 1  ]; __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ) {
            volatile double *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ) {
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif

    if ( Idx == 0 ) {
            vtmp2[ blockIdx.x ] = temp[ 0 ];
    }
}



// computes the SpMV using SELLC
__global__ void
magma_dcgmerge_spmvsellc_kernel(   
    int num_rows, 
    int blocksize,
    double * dval, 
    magma_index_t * dcolind,
    magma_index_t * drowptr,
    double * d,
    double * z,
    double * vtmp)
{
    extern __shared__ double temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;
    int offset = drowptr[ blockIdx.x ];
    int border = (drowptr[ blockIdx.x+1 ]-offset)/blocksize;

    temp[ Idx ] = MAGMA_D_MAKE( 0.0, 0.0);


    if(i < num_rows ) {
        double dot = MAGMA_D_MAKE(0.0, 0.0);
        for ( int n = 0; n < border; n ++) {
            int col = dcolind [offset+ blocksize * n + Idx ];
            double val = dval[offset+ blocksize * n + Idx];
            if( val != 0) {
                  dot=dot+val*d[col];
            }
        }
        z[ i ] =  dot;
        temp[ Idx ] = d[ i ] * dot;
    }
    __syncthreads();
    if ( Idx < 128 ) {
        temp[ Idx ] += temp[ Idx + 128 ];
    }
    __syncthreads();
    if ( Idx < 64 ) {
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ) {
            temp[ Idx ] += temp[ Idx + 32 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 8  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 4  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 2  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 1  ]; __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ) {
            volatile double *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ) {
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif

    if ( Idx == 0 ) {
            vtmp[ blockIdx.x ] = temp[ 0 ];
    }
}


// SELLP SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void
magma_dcgmerge_spmvsellpt_kernel_8( 
    int num_rows, 
    int blocksize,
    int T,
    double * dval, 
    magma_index_t * dcolind,
    magma_index_t * drowptr,
    double * d,
    double * z)
{
   // T threads assigned to each row
    int idx = threadIdx.y;      // thread in row
    int idy = threadIdx.x;      // local row
    int ldx = idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index

    extern __shared__ double shared[];

    if(row < num_rows ) {
        double dot = MAGMA_D_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles
        for ( int k = 0; k < max_; k++ ) {
            double val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    dcolind[ offset + ldx + block*k ];
            dot += val * d[ col ];
        }
        shared[ldx]  = dot;

        __syncthreads();
        if( idx < 4 ) {
            shared[ldx]+=shared[ldx+blocksize*4];              
            __syncthreads();
            if( idx < 2 ) shared[ldx]+=shared[ldx+blocksize*2];   
            __syncthreads();
            if( idx == 0 ) {
                z[row] = 
                (shared[ldx]+shared[ldx+blocksize*1]);
            }
        }
    }
}
// SELLP SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void
magma_dcgmerge_spmvsellpt_kernel_16( 
    int num_rows, 
    int blocksize,
    int T,
    double * dval, 
    magma_index_t * dcolind,
    magma_index_t * drowptr,
    double * d,
    double * z)
{
   // T threads assigned to each row
    int idx = threadIdx.y;      // thread in row
    int idy = threadIdx.x;      // local row
    int ldx = idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index

    extern __shared__ double shared[];

    if(row < num_rows ) {
        double dot = MAGMA_D_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles
        for ( int k = 0; k < max_; k++ ) {
            double val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    dcolind[ offset + ldx + block*k ];
            dot += val * d[ col ];
        }
        shared[ldx]  = dot;

        __syncthreads();
        if( idx < 8 ) {
            shared[ldx]+=shared[ldx+blocksize*8];              
            __syncthreads();
            if( idx < 4 ) shared[ldx]+=shared[ldx+blocksize*4];   
            __syncthreads();
            if( idx < 2 ) shared[ldx]+=shared[ldx+blocksize*2];   
            __syncthreads();
            if( idx == 0 ) {
                z[row] = 
                (shared[ldx]+shared[ldx+blocksize*1]);
            }
        }
    }
}


// SELLP SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void
magma_dcgmerge_spmvsellpt_kernel_32( 
    int num_rows, 
    int blocksize,
    int T,
    double * dval, 
    magma_index_t * dcolind,
    magma_index_t * drowptr,
    double * d,
    double * z)
{
   // T threads assigned to each row
    int idx = threadIdx.y;      // thread in row
    int idy = threadIdx.x;      // local row
    int ldx = idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index

    extern __shared__ double shared[];

    if(row < num_rows ) {
        double dot = MAGMA_D_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles
        for ( int k = 0; k < max_; k++ ) {
            double val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    dcolind[ offset + ldx + block*k ];
            dot += val * d[ col ];
        }
        shared[ldx]  = dot;

        __syncthreads();
        if( idx < 16 ) {
            shared[ldx]+=shared[ldx+blocksize*16];              
            __syncthreads();
            if( idx < 8 ) shared[ldx]+=shared[ldx+blocksize*8];  
            __syncthreads();
            if( idx < 4 ) shared[ldx]+=shared[ldx+blocksize*4];   
            __syncthreads();
            if( idx < 2 ) shared[ldx]+=shared[ldx+blocksize*2];   
            __syncthreads();
            if( idx == 0 ) {
                z[row] = 
                (shared[ldx]+shared[ldx+blocksize*1]);
            }
        }
    }
}


// kernel to handle scalars
__global__ void // rho = beta/tmp; gamma = beta;
magma_dcg_rhokernel(  
    double * skp ) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i==0 ) {
        double tmp = skp[1];
        skp[3] = tmp/skp[4];
        skp[2] = tmp;
    }
}

/**
    Purpose
    -------

    Merges the first SpmV using different formats with the dot product 
    and the computation of rho

    Arguments
    ---------

    @param[in]
    A           magma_d_matrix
                input matrix 

    @param[in]
    d1          magmaDouble_ptr 
                temporary vector

    @param[in]
    d2          magmaDouble_ptr 
                temporary vector

    @param[in]
    dd          magmaDouble_ptr 
                input vector d

    @param[out]
    dz          magmaDouble_ptr 
                input vector z

    @param[out]
    skp         magmaDouble_ptr 
                array for parameters ( skp[3]=rho )

    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_dgegpuk
    ********************************************************************/

extern "C" magma_int_t
magma_dcgmerge_spmv1(
    magma_d_matrix A,
    magmaDouble_ptr d1,
    magmaDouble_ptr d2,
    magmaDouble_ptr dd,
    magmaDouble_ptr dz,
    magmaDouble_ptr skp,
    magma_queue_t queue )
{
    int local_block_size=256;
    dim3 Bs( local_block_size );
    dim3 Gs( magma_ceildiv( A.num_rows, local_block_size ) );
    dim3 Gs_next;
    int Ms =  local_block_size * sizeof( double ); 
    magmaDouble_ptr aux1 = d1, aux2 = d2;
    int b = 1;        

    if ( A.storage_type == Magma_CSR )
        magma_dcgmerge_spmvcsr_kernel<<< Gs, Bs, Ms, queue->cuda_stream() >>>
        ( A.num_rows, A.dval, A.drow, A.dcol, dd, dz, d1 );
    else if ( A.storage_type == Magma_ELLPACKT )
        magma_dcgmerge_spmvellpack_kernel<<< Gs, Bs, Ms, queue->cuda_stream() >>>
        ( A.num_rows, A.max_nnz_row, A.dval, A.dcol, dd, dz, d1 );
    else if ( A.storage_type == Magma_ELL )
        magma_dcgmerge_spmvell_kernel<<< Gs, Bs, Ms, queue->cuda_stream() >>>
        ( A.num_rows, A.max_nnz_row, A.dval, A.dcol, dd, dz, d1 );
    else if ( A.storage_type == Magma_CUCSR ) {
        cusparseHandle_t cusparseHandle = 0;
        cusparseMatDescr_t descr = 0;
        double c_one = MAGMA_D_ONE;
        double c_zero = MAGMA_D_ZERO;
        cusparseCreate( &cusparseHandle );
        cusparseSetStream( cusparseHandle, queue->cuda_stream() );
        cusparseCreateMatDescr( &descr );
        cusparseSetMatType( descr, CUSPARSE_MATRIX_TYPE_GENERAL );
        cusparseSetMatIndexBase( descr, CUSPARSE_INDEX_BASE_ZERO );
        cusparseDcsrmv( cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE,
        A.num_rows, A.num_cols, A.nnz, &c_one, descr,
        A.dval, A.drow, A.dcol, dd, &c_zero, dz );
        cusparseDestroyMatDescr( descr );
        cusparseDestroy( cusparseHandle );
        cusparseHandle = 0;
        descr = 0;
        magma_dcgmerge_spmvellpackrt_kernel2<<< Gs, Bs, Ms, queue->cuda_stream() >>>
                      ( A.num_rows, dz, dd, d1 );
    }
    else if ( A.storage_type == Magma_SELLP && A.alignment == 1 ) {
            magma_dcgmerge_spmvell_kernelb1<<< Gs, Bs, Ms, queue->cuda_stream() >>>
            ( A.num_rows, A.blocksize, 
                A.dval, A.dcol, A.drow, dd, dz, d1 );
    }
    else if ( A.storage_type == Magma_SELLP && A.alignment > 1) {
            int num_threadssellp = A.blocksize*A.alignment;
            magma_int_t arch = magma_getdevice_arch();
            if ( arch < 200 && num_threadssellp > 256 )
                printf("error: too much shared memory requested.\n");

            dim3 block( A.blocksize, A.alignment, 1);
            int dimgrid1 = int( sqrt( double( A.numblocks )));
            int dimgrid2 = magma_ceildiv( A.numblocks, dimgrid1 );

            dim3 gridsellp( dimgrid1, dimgrid2, 1);
            int Mssellp = num_threadssellp * sizeof( double );

            if ( A.alignment == 8)
                magma_dcgmerge_spmvsellpt_kernel_8
                <<< gridsellp, block, Mssellp, queue->cuda_stream() >>>
                ( A.num_rows, A.blocksize, A.alignment, 
                    A.dval, A.dcol, A.drow, dd, dz);

            else if ( A.alignment == 16)
                magma_dcgmerge_spmvsellpt_kernel_16
                <<< gridsellp, block, Mssellp, queue->cuda_stream() >>>
                ( A.num_rows, A.blocksize, A.alignment, 
                    A.dval, A.dcol, A.drow, dd, dz);

            else if ( A.alignment == 32)
                magma_dcgmerge_spmvsellpt_kernel_32
                <<< gridsellp, block, Mssellp, queue->cuda_stream() >>>
                ( A.num_rows, A.blocksize, A.alignment, 
                    A.dval, A.dcol, A.drow, dd, dz);

            else
                printf("error: alignment not supported.\n");

        // in case of using SELLP, we can't efficiently merge the 
        // dot product and the first reduction loop into the SpMV kernel
        // as the SpMV grid would result in low occupancy.
        magma_dcgmerge_spmvellpackrt_kernel2<<< Gs, Bs, Ms, queue->cuda_stream() >>>
                              ( A.num_rows, dz, dd, d1 );
    }
    else if ( A.storage_type == Magma_ELLRT ) {
        // in case of using ELLRT, we need a different grid, assigning
        // threads_per_row processors to each row
        // the block size is num_threads
        // fixed values


    int num_blocks = magma_ceildiv( A.num_rows, A.blocksize );

    int num_threads = A.alignment*A.blocksize;

    int real_row_length = magma_roundup( A.max_nnz_row, A.alignment );

    magma_int_t arch = magma_getdevice_arch();
    if ( arch < 200 && num_threads > 256 )
        printf("error: too much shared memory requested.\n");

    int dimgrid1 = int( sqrt( double( num_blocks )));
    int dimgrid2 = magma_ceildiv( num_blocks, dimgrid1 );
    dim3 gridellrt( dimgrid1, dimgrid2, 1);

    int Mellrt = A.alignment * A.blocksize * sizeof( double );
    // printf("launch kernel: %dx%d %d %d\n", grid.x, grid.y, num_threads , Ms);

    if ( A.alignment == 32 ) {
        magma_dcgmerge_spmvellpackrt_kernel_32
                <<< gridellrt, num_threads , Mellrt, queue->cuda_stream() >>>
                 ( A.num_rows, A.dval, A.dcol, A.drow, dd, dz, d1, 
                                                 A.alignment, real_row_length );
    }
    else if ( A.alignment == 16 ) {
        magma_dcgmerge_spmvellpackrt_kernel_16
                <<< gridellrt, num_threads , Mellrt, queue->cuda_stream() >>>
                 ( A.num_rows, A.dval, A.dcol, A.drow, dd, dz, d1, 
                                                 A.alignment, real_row_length );
    }
    else if ( A.alignment == 8 ) {
        magma_dcgmerge_spmvellpackrt_kernel_8
                <<< gridellrt, num_threads , Mellrt, queue->cuda_stream() >>>
                 ( A.num_rows, A.dval, A.dcol, A.drow, dd, dz, d1, 
                                                 A.alignment, real_row_length );
    }
    else {
        printf("error: alignment %d not supported.\n", int(A.alignment) );
        return MAGMA_ERR_NOT_SUPPORTED;
    }
        // in case of using ELLRT, we can't efficiently merge the 
        // dot product and the first reduction loop into the SpMV kernel
        // as the SpMV grid would result in low occupancy.

        magma_dcgmerge_spmvellpackrt_kernel2<<< Gs, Bs, Ms, queue->cuda_stream() >>>
                              ( A.num_rows, dz, dd, d1 );
    }

    while( Gs.x > 1 ) {
        Gs_next.x = magma_ceildiv( Gs.x, Bs.x );
        if ( Gs_next.x == 1 ) Gs_next.x = 2;
        magma_dcgreduce_kernel_spmv1<<< Gs_next.x/2, Bs.x/2, Ms/2, queue->cuda_stream()>>> 
                                        ( Gs.x,  A.num_rows, aux1, aux2 );
        Gs_next.x = Gs_next.x /2;
        Gs.x = Gs_next.x;
        b = 1 - b;
        if ( b ) { aux1 = d1; aux2 = d2; }
        else   { aux2 = d1; aux1 = d2; }
    }


    magma_dcopyvector( 1, aux1, 1, skp+4, 1, queue );
    dim3 Bs2( 2 );
    dim3 Gs2( 1 );
    magma_dcg_rhokernel<<< Gs2, Bs2, 0, queue->cuda_stream()>>>( skp );

    return MAGMA_SUCCESS;
}


/* -------------------------------------------------------------------------- */

// updates x and r and computes the first part of the dot product r*r
__global__ void
magma_dcgmerge_xrbeta_kernel(  
    int n, 
    double * x, 
    double * r,
    double * d,
    double * z,
    double * skp,
    double * vtmp )
{
    extern __shared__ double temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;

    double rho = skp[3];
    double mrho = MAGMA_D_MAKE( -1.0, 0.0)*rho;

    temp[ Idx ] = MAGMA_D_MAKE( 0.0, 0.0);

    if( i<n ) {
        x[i] += rho * d[i];
        r[i] += mrho * z[i];
        temp[ Idx ] = r[i] * r[i];
    }
    __syncthreads();
    if ( Idx < 128 ) {
        temp[ Idx ] += temp[ Idx + 128 ];
    }
    __syncthreads();
    if ( Idx < 64 ) {
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ) {
            temp[ Idx ] += temp[ Idx + 32 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 8  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 4  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 2  ]; __syncthreads();
            temp[ Idx ] += temp[ Idx + 1  ]; __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ) {
            volatile double *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ) {
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif

    if ( Idx == 0 ) {
            vtmp[ blockIdx.x ] = temp[ 0 ];
    }
}

// kernel to handle scalars
__global__ void //alpha = beta / gamma
magma_dcg_alphabetakernel(  
    double * skp )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i==0 ) {
        double tmp1 = skp[1];
        skp[0] =  tmp1/skp[2];
        //printf("beta=%e\n", MAGMA_D_REAL(tmp1));
    }
}

// update search Krylov vector d
__global__ void
magma_dcg_d_kernel(  
    int n, 
    double * skp,
    double * r,
    double * d )
{
    int i   = blockIdx.x * blockDim.x + threadIdx.x;

    double alpha = skp[0];

    if( i<n ) {
        d[i] = r[i] + alpha * d[i];
    }
}



/**
    Purpose
    -------

    Merges the update of r and x with the dot product and performs then 
    the update for the Krylov vector d

    Arguments
    ---------

    @param[in]
    n           int
                dimension n

    @param[in]
    d1          magmaDouble_ptr 
                temporary vector

    @param[in]
    d2          magmaDouble_ptr 
                temporary vector

    @param[in,out]
    dx          magmaDouble_ptr
                input vector x

    @param[in,out]
    dr          magmaDouble_ptr 
                input/output vector r

    @param[in]
    dd          magmaDouble_ptr 
                input vector d

    @param[in]
    dz          magmaDouble_ptr 
                input vector z
    @param[in]
    skp         magmaDouble_ptr 
                array for parameters

    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_dsygpuk
    ********************************************************************/

extern "C" magma_int_t
magma_dcgmerge_xrbeta(
    magma_int_t n,
    magmaDouble_ptr d1,
    magmaDouble_ptr d2,
    magmaDouble_ptr dx,
    magmaDouble_ptr dr,
    magmaDouble_ptr dd,
    magmaDouble_ptr dz, 
    magmaDouble_ptr skp,
    magma_queue_t queue )
{
    int local_block_size=256;
    dim3 Bs( local_block_size );
    dim3 Gs( magma_ceildiv( n, local_block_size ) );
    dim3 Gs_next;
    int Ms =  2*local_block_size * sizeof( double ); 
    magmaDouble_ptr aux1 = d1, aux2 = d2;
    int b = 1;        
    magma_dcgmerge_xrbeta_kernel<<< Gs, Bs, Ms, queue->cuda_stream()>>>
                                    ( n, dx, dr, dd, dz, skp, d1);  



    while( Gs.x > 1 ) {
        Gs_next.x = magma_ceildiv( Gs.x, Bs.x );
        if ( Gs_next.x == 1 ) Gs_next.x = 2;
        magma_dcgreduce_kernel_spmv1<<< Gs_next.x/2, Bs.x/2, Ms/2, queue->cuda_stream()>>> 
                                    ( Gs.x, n, aux1, aux2 );
        Gs_next.x = Gs_next.x /2;
        Gs.x = Gs_next.x;
        b = 1 - b;
        if ( b ) { aux1 = d1; aux2 = d2; }
        else   { aux2 = d1; aux1 = d2; }
    }


    magma_dcopyvector( 1, aux1, 1, skp+1, 1, queue );
    dim3 Bs2( 2 );
    dim3 Gs2( 1 );
    magma_dcg_alphabetakernel<<< Gs2, Bs2, 0, queue->cuda_stream()>>>( skp );

    dim3 Bs3( local_block_size );
    dim3 Gs3( magma_ceildiv( n, local_block_size ) );
    magma_dcg_d_kernel<<< Gs3, Bs3, 0, queue->cuda_stream()>>>( n, skp, dr, dd );  

    return MAGMA_SUCCESS;
}

/* -------------------------------------------------------------------------- */

// updates x and r
__global__ void
magma_dpcgmerge_xrbeta_kernel(  
    int n, 
    double * x, 
    double * r,
    double * d,
    double * z,
    double * skp )
{
    extern __shared__ double temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;

    double rho = skp[3];
    double mrho = MAGMA_D_MAKE( -1.0, 0.0)*rho;

    if( i<n ) {
        x[i] += rho * d[i];
        r[i] += mrho * z[i];
    }
}


// dot product for multiple vectors
__global__ void
magma_dmddot_one_kernel_1( 
    int n, 
    double * v0,
    double * w0,
    double * vtmp)
{
    extern __shared__ double temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;
    int j;

    // 1 vectors v(i)/w(i)
    
    temp[ Idx ]                 = ( i < n ) ?
                v0[ i ] * w0[ i ] : MAGMA_D_ZERO;
    temp[ Idx + blockDim.x ]    = ( i < n ) ?
                v0[ i ] * v0[ i ] : MAGMA_D_ZERO;
    
    __syncthreads();
    if ( Idx < 128 ){
        for( j=0; j<2; j++){
            temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 128 ];
        }
    }
    __syncthreads();
    if ( Idx < 64 ){
        for( j=0; j<2; j++){
            temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 64 ];
        }
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 32 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 16 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 8 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 4 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 2 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 1 ];
                __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile double *temp2 = temp;
            for( j=0; j<2; j++){
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 32 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 16 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 8 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 4 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 2 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 1 ];
            }
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            for( j=0; j<2; j++){
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 32 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 16 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 8 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 4 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 2 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 1 ];
            }
        }
    #endif  
    
    if ( Idx == 0 ){
            vtmp[ blockIdx.x ] = temp[ 0 ];
            vtmp[ blockIdx.x+n ] = temp[ blockDim.x ];
    }
}

/**
    Purpose
    -------

    Merges the update of r and x with the dot product and performs then 
    the update for the Krylov vector d

    Arguments
    ---------

    @param[in]
    n           int
                dimension n

    @param[in,out]
    dx          magmaDouble_ptr
                input vector x

    @param[in,out]
    dr          magmaDouble_ptr 
                input/output vector r

    @param[in]
    dd          magmaDouble_ptr 
                input vector d

    @param[in]
    dz          magmaDouble_ptr 
                input vector z
    @param[in]
    skp         magmaDouble_ptr 
                array for parameters

    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_dsygpuk
    ********************************************************************/

extern "C" magma_int_t
magma_dpcgmerge_xrbeta1(
    magma_int_t n,
    magmaDouble_ptr dx,
    magmaDouble_ptr dr,
    magmaDouble_ptr dd,
    magmaDouble_ptr dz, 
    magmaDouble_ptr skp,
    magma_queue_t queue )
{
    int local_block_size=256;
    dim3 Bs( local_block_size );
    dim3 Gs( magma_ceildiv( n, local_block_size ) );
    magma_dpcgmerge_xrbeta_kernel<<< Gs, Bs, 0, queue->cuda_stream()>>>
                                    ( n, dx, dr, dd, dz, skp );  
                                    
    return MAGMA_SUCCESS;
}

/* -------------------------------------------------------------------------- */


/**
    Purpose
    -------

    Merges the update of r and x with the dot product and performs then 
    the update for the Krylov vector d

    Arguments
    ---------

    @param[in]
    n           int
                dimension n

    @param[in]
    d1          magmaDouble_ptr 
                temporary vector

    @param[in]
    d2          magmaDouble_ptr 
                temporary vector

    @param[in]
    dh          magmaDouble_ptr
                input vector x

    @param[in]
    dr          magmaDouble_ptr 
                input/output vector r
                
    @param[in]
    dd          magmaDouble_ptr 
                input/output vector d

    @param[in]
    skp         magmaDouble_ptr 
                array for parameters

    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_dsygpuk
    ********************************************************************/

extern "C" magma_int_t
magma_dpcgmerge_xrbeta2(
    magma_int_t n,
    magmaDouble_ptr d1,
    magmaDouble_ptr d2,
    magmaDouble_ptr dh,
    magmaDouble_ptr dr, 
    magmaDouble_ptr dd, 
    magmaDouble_ptr skp,
    magma_queue_t queue )
{
    int local_block_size=256;
    dim3 Bs( local_block_size );
    dim3 Gs( magma_ceildiv( n, local_block_size ) );
    dim3 Gs_next;
    int Ms =  4*local_block_size * sizeof( double ); 
    magmaDouble_ptr aux1 = d1, aux2 = d2;
    int b = 1;        
                                    
    magma_dmddot_one_kernel_1<<< Gs, Bs, Ms, queue->cuda_stream()>>>
                                    ( n, dr, dh, d1);  

    while( Gs.x > 1 ) {
        Gs_next.x = magma_ceildiv( Gs.x, Bs.x );
        if ( Gs_next.x == 1 ) Gs_next.x = 2;
        magma_dcgreduce_kernel_spmv2<<< Gs_next.x/2, Bs.x/2, Ms/2, queue->cuda_stream()>>> 
                                    ( Gs.x, n, aux1, aux2 );
        Gs_next.x = Gs_next.x /2;
        Gs.x = Gs_next.x;
        b = 1 - b;
        if ( b ) { aux1 = d1; aux2 = d2; }
        else   { aux2 = d1; aux1 = d2; }
    }


    magma_dcopyvector( 1, aux1, 1, skp+1, 1, queue );
    magma_dcopyvector( 1, aux1+n, 1, skp+6, 1, queue );
    dim3 Bs2( 2 );
    dim3 Gs2( 1 );
    magma_dcg_alphabetakernel<<< Gs2, Bs2, 0, queue->cuda_stream()>>>( skp );

    dim3 Bs3( local_block_size );
    dim3 Gs3( magma_ceildiv( n, local_block_size ) );
    magma_dcg_d_kernel<<< Gs3, Bs3, 0, queue->cuda_stream()>>>( n, skp, dh, dd );  

    return MAGMA_SUCCESS;
}

/* -------------------------------------------------------------------------- */



// updates x and r
__global__ void
magma_djcgmerge_xrbeta_kernel(  
    int n, 
    double * diag, 
    double * x,     
    double * r,
    double * d,
    double * z,
    double * h,
    double * vtmp,
    double * skp )
{
    extern __shared__ double temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;
    int j;

    double rho = skp[3];
    double mrho = MAGMA_D_MAKE( -1.0, 0.0)*rho;

    if( i<n ) {
        x[i] += rho * d[i];
        r[i] += mrho * z[i];
        h[i] = r[i] * diag[i];
    }
    __syncthreads();
    temp[ Idx ]                 = ( i < n ) ?
                h[ i ] * r[ i ] : MAGMA_D_ZERO;
    temp[ Idx + blockDim.x ]    = ( i < n ) ?
                r[ i ] * r[ i ] : MAGMA_D_ZERO;
    
    __syncthreads();
    if ( Idx < 128 ){
        for( j=0; j<2; j++){
            temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 128 ];
        }
    }
    __syncthreads();
    if ( Idx < 64 ){
        for( j=0; j<2; j++){
            temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 64 ];
        }
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 32 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 16 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 8 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 4 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 2 ];
                __syncthreads();
            for( j=0; j<2; j++)
                temp[ Idx+j*blockDim.x ] += temp[ Idx+j*blockDim.x + 1 ];
                __syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile double *temp2 = temp;
            for( j=0; j<2; j++){
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 32 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 16 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 8 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 4 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 2 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 1 ];
            }
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            for( j=0; j<2; j++){
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 32 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 16 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 8 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 4 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 2 ];
                temp2[ Idx+j*blockDim.x ] += temp2[ Idx+j*blockDim.x + 1 ];
            }
        }
    #endif  
    
    if ( Idx == 0 ){
            vtmp[ blockIdx.x ] = temp[ 0 ];
            vtmp[ blockIdx.x+n ] = temp[ blockDim.x ];
    }
}


/**
    Purpose
    -------

    Merges the update of r and x with the dot product and performs then 
    the update for the Krylov vector d

    Arguments
    ---------

    @param[in]
    n           int
                dimension n

    @param[in]
    d1          magmaDouble_ptr 
                temporary vector

    @param[in]
    d2          magmaDouble_ptr 
                temporary vector
                
    @param[in]
    diag        magmaDouble_ptr 
                inverse diagonal (Jacobi preconditioner)

    @param[in]
    dx          magmaDouble_ptr
                iteration vector x

    @param[in]
    dr          magmaDouble_ptr 
                input/output vector r
                
    @param[in]
    dd          magmaDouble_ptr
                input vector d

                
    @param[in]
    dz          magmaDouble_ptr
                input vector z
                
    @param[in]
    dh          magmaDouble_ptr
                input vector h

    @param[in]
    skp         magmaDouble_ptr 
                array for parameters

    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_dsygpuk
    ********************************************************************/

extern "C" magma_int_t
magma_djcgmerge_xrbeta(
    magma_int_t n,
    magmaDouble_ptr d1,
    magmaDouble_ptr d2,
    magmaDouble_ptr diag,
    magmaDouble_ptr dx,
    magmaDouble_ptr dr,
    magmaDouble_ptr dd,
    magmaDouble_ptr dz,
    magmaDouble_ptr dh, 
    magmaDouble_ptr skp,
    magma_queue_t queue )
{
    int local_block_size=256;
    dim3 Bs( local_block_size );
    dim3 Gs( magma_ceildiv( n, local_block_size ) );
    dim3 Gs_next;
    int Ms =  4*local_block_size * sizeof( double ); 
    magmaDouble_ptr aux1 = d1, aux2 = d2;
    int b = 1;    
                                    
    magma_djcgmerge_xrbeta_kernel<<< Gs, Bs, Ms, queue->cuda_stream() >>>
                                    ( n, diag, dx, dr, dd, dz, dh, d1, skp );  
                                    
    while( Gs.x > 1 ) {
        Gs_next.x = magma_ceildiv( Gs.x, Bs.x );
        if ( Gs_next.x == 1 ) Gs_next.x = 2;
        magma_dcgreduce_kernel_spmv2<<< Gs_next.x/2, Bs.x/2, Ms/2, queue->cuda_stream() >>> 
                                    ( Gs.x, n, aux1, aux2 );
        Gs_next.x = Gs_next.x /2;
        Gs.x = Gs_next.x;
        b = 1 - b;
        if ( b ) { aux1 = d1; aux2 = d2; }
        else   { aux2 = d1; aux1 = d2; }
    }


    magma_dcopyvector( 1, aux1, 1, skp+1, 1, queue );
    magma_dcopyvector( 1, aux1+n, 1, skp+6, 1, queue );
    dim3 Bs2( 2 );
    dim3 Gs2( 1 );
    magma_dcg_alphabetakernel<<< Gs2, Bs2, 0, queue->cuda_stream()>>>( skp );

    dim3 Bs3( local_block_size );
    dim3 Gs3( magma_ceildiv( n, local_block_size ) );
    magma_dcg_d_kernel<<< Gs3, Bs3, 0, queue->cuda_stream()>>>( n, skp, dh, dd );  

    return MAGMA_SUCCESS;
}

/* -------------------------------------------------------------------------- */
