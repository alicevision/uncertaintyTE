/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @generated from magmablas/zlarfbx.cu, normal z -> s, Sun Nov 20 20:20:29 2016

*/
#include "magma_internal.h"
#include "commonblas_s.h"
#include "magma_templates.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE 512


/******************************************************************************/
extern "C"
__global__ void 
magma_sgemv_kernel1(int m, const float * __restrict__ V, int ldv, 
                    const float * __restrict__ c, 
                    float *dwork)
{
    const int i = threadIdx.x;
    const float *dV = V + (blockIdx.x) * ldv;

    __shared__ float sum[ BLOCK_SIZE ];
    float lsum;

    /*  lsum := v**H * C  */
    lsum = MAGMA_S_ZERO;
    for (int j = i; j < m; j += BLOCK_SIZE)
       lsum += MAGMA_S_MUL( MAGMA_S_CONJ( dV[j] ), c[j] );
    
    sum[i] = lsum;
    magma_sum_reduce< BLOCK_SIZE >( i, sum );

    __syncthreads();
    if (i == 0)
       dwork [blockIdx.x] = sum[0];
}

/******************************************************************************/
/*
    Call 
        magma_sgemv_kernel3<<< n, BLOCK_SIZE, 0, queue->cuda_stream() >>>(m, V, ldv, c, dwork, tau)
    to compute
        SGEMV( "Conjugate transpose", m, n, -tau[0], V, ldv, c, 1, zero, dwork, 1)
        and to set c[0] to 1.
    i.e., 
        work = -tau[0] V**H c
*/
extern "C"
__global__ void
magma_sgemv_kernel3(int m, const float * __restrict__ V, int ldv, float *c,
                    float *dwork, float *tau)
{
    const int i = threadIdx.x;
    const float *dV = V + (blockIdx.x) * ldv;

    __shared__ float sum[ BLOCK_SIZE ];
    float lsum;

    if (i == 0)
       c[0] = MAGMA_S_ONE;           

    /*  lsum := v**H * C  */
    lsum = MAGMA_S_ZERO;
    for (int j = i; j < m; j += BLOCK_SIZE)
       lsum += MAGMA_S_MUL( MAGMA_S_CONJ( dV[j] ), c[j] );

    sum[i] = lsum;
    magma_sum_reduce< BLOCK_SIZE >( i, sum );

    __syncthreads();
    if (i == 0)
       dwork [blockIdx.x] = -tau[0]*sum[0];
}


/******************************************************************************/
extern "C"
__global__ void
magma_sgemv_kernel2(int m, int n, const float * __restrict__ V, int ldv, 
                    const float * __restrict__ x, float *c)
{
    const int i = threadIdx.x;
    const int j = i + BLOCK_SIZE * blockIdx.x;
    float lsum;

    V += j;

    lsum = MAGMA_S_ZERO;
    if (j < m) {
        for (int k=0; k < n; k++)
            lsum += MAGMA_S_MUL( V[k*ldv], x[k]);
        
        c[j] -= lsum;
    }
}


/******************************************************************************/
/*
    Apply a real block reflector H to a real vector C from the left
    (i.e., C = H C). H is represented in the form
          H = I - V T V**H
    where T is the real k-by-k upper triangular matrix in the 
    representation of the block reflector, and V is a real block of
    k elementary reflectors. 
*/
extern "C" void
magma_slarfbx_gpu(
    magma_int_t m, magma_int_t k,
    magmaFloat_ptr V,  magma_int_t ldv,
    magmaFloat_ptr dT, magma_int_t ldt,
    magmaFloat_ptr c,
    magmaFloat_ptr dwork,
    magma_queue_t queue )
{
    /* dwork = V**H c     */
    magma_sgemv_kernel1
        <<< k, BLOCK_SIZE, 0, queue->cuda_stream() >>>
        (m, V, ldv, c, dwork); 

    /* dwork = T**H dwork */
    magma_strmv_tkernel
        <<< k, k, 0, queue->cuda_stream() >>>
        ( dT, ldt, dwork, dwork+k);
 
    /* c = c - V dwork    */
    dim3  blocks3( magma_ceildiv( m, BLOCK_SIZE ) );
    dim3 threads3( BLOCK_SIZE );     
    magma_sgemv_kernel2
        <<< blocks3, threads3, 0, queue->cuda_stream() >>>
        ( m, k, V, ldv, dwork+k, c);
}
