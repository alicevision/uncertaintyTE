/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @generated from magmablas/zgeadd2.cu, normal z -> s, Sun Nov 20 20:20:27 2016
       @author Mark Gates
*/
#include "magma_internal.h"

#define BLK_X 64
#define BLK_Y 32

/*
    Divides matrix into ceil( m/BLK_X ) x ceil( n/BLK_Y ) blocks.
    Each block has BLK_X threads.
    Each thread loops across one row, updating BLK_Y entries.

    Code similar to slaset.
*/
__global__
void sgeadd2_full(
    int m, int n,
    float alpha,
    const float *dA, int ldda,
    float beta,
    float       *dB, int lddb )
{
    int ind = blockIdx.x*BLK_X + threadIdx.x;
    int iby = blockIdx.y*BLK_Y;
    /* check if full block-column */
    bool full = (iby + BLK_Y <= n);
    /* do only rows inside matrix */
    if ( ind < m ) {
        dA += ind + iby*ldda;
        dB += ind + iby*lddb;
        if ( full ) {
            // full block-column
            #pragma unroll
            for( int j=0; j < BLK_Y; ++j ) {
                dB[j*lddb] = alpha*dA[j*ldda] + beta*dB[j*lddb];
            }
        }
        else {
            // partial block-column
            for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
                dB[j*lddb] = alpha*dA[j*ldda] + beta*dB[j*lddb];
            }
        }
    }
}


/***************************************************************************//**
    Purpose
    -------
    ZGEADD2 adds two matrices, dB = alpha*dA + beta*dB.
    
    @see ZGEADD for dB = alpha*dA + dB, lacking beta.
    
    Arguments
    ---------
    
    @param[in]
    m       INTEGER
            The number of rows of the matrix dA.  M >= 0.
    
    @param[in]
    n       INTEGER
            The number of columns of the matrix dA.  N >= 0.
    
    @param[in]
    alpha   REAL
            The scalar alpha.
            
    @param[in]
    dA      REAL array, dimension (LDDA,N)
            The m by n matrix dA.
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).
    
    @param[in]
    beta    REAL
            The scalar beta.
            
    @param[in,out]
    dB      REAL array, dimension (LDDB,N)
            The m by n matrix dB.
    
    @param[in]
    lddb    INTEGER
            The leading dimension of the array dB.  LDDB >= max(1,M).
    
    @param[in]
    queue   magma_queue_t
            Queue to execute in.

    @ingroup magma_geadd
*******************************************************************************/
extern "C" void
magmablas_sgeadd2(
    magma_int_t m, magma_int_t n,
    float alpha,
    magmaFloat_const_ptr dA, magma_int_t ldda,
    float beta,
    magmaFloat_ptr       dB, magma_int_t lddb,
    magma_queue_t queue )
{
    magma_int_t info = 0;
    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( ldda < max(1,m))
        info = -5;
    else if ( lddb < max(1,m))
        info = -7;
    
    if ( info != 0 ) {
        magma_xerbla( __func__, -(info) );
        return;
    }
    
    if ( m == 0 || n == 0 )
        return;
    
    dim3 threads( BLK_X, 1 );
    dim3 grid( magma_ceildiv( m, BLK_X ), magma_ceildiv( n, BLK_Y ) );
    
    sgeadd2_full<<< grid, threads, 0, queue->cuda_stream() >>>
        ( m, n, alpha, dA, ldda, beta, dB, lddb );
}
