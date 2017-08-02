/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016
       
       @author Azzam Haidar
       @author Tingxing Dong
       @author Ahmad Abdelfattah
       
       @generated from src/zgetri_outofplace_batched.cpp, normal z -> s, Sun Nov 20 20:20:26 2016
*/
#include "magma_internal.h"
#include "batched_kernel_param.h"

/***************************************************************************//**
    Purpose
    -------
    SGETRI computes the inverse of a matrix using the LU factorization
    computed by SGETRF. This method inverts U and then computes inv(A) by
    solving the system inv(A)*L = inv(U) for inv(A).
    
    Note that it is generally both faster and more accurate to use SGESV,
    or SGETRF and SGETRS, to solve the system AX = B, rather than inverting
    the matrix and multiplying to form X = inv(A)*B. Only in special
    instances should an explicit inverse be computed with this routine.

    Arguments
    ---------
    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    dA_array Array of pointers, dimension (batchCount).
            Each is a REAL array on the GPU, dimension (LDDA,N)
            On entry, the factors L and U from the factorization
            A = P*L*U as computed by SGETRF_GPU.
            On exit, if INFO = 0, the inverse of the original matrix A.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[in]
    dipiv_array Array of pointers, dimension (batchCount), for corresponding matrices.
            Each is an INTEGER array, dimension (N)
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    @param[out]
    dinvA_array Array of pointers, dimension (batchCount).
            Each is a REAL array on the GPU, dimension (LDDIA,N)
            It contains the inverse of the matrix
  
    @param[in]
    lddia   INTEGER
            The leading dimension of the array invA_array.  LDDIA >= max(1,N).

    @param[out]
    info_array  Array of INTEGERs, dimension (batchCount), for corresponding matrices.
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.
      -     > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                  has been completed, but the factor U is exactly
                  singular, and division by zero will occur if it is used
                  to solve a system of equations.

    @param[in]
    batchCount  INTEGER
                The number of matrices to operate on.

    @param[in]
    queue   magma_queue_t
            Queue to execute in.
                  
    @ingroup magma_getri_batched
*******************************************************************************/
extern "C" magma_int_t
magma_sgetri_outofplace_batched( magma_int_t n, 
                  float **dA_array, magma_int_t ldda,
                  magma_int_t **dipiv_array, 
                  float **dinvA_array, magma_int_t lddia,
                  magma_int_t *info_array,
                  magma_int_t batchCount, magma_queue_t queue)
{
    /* Local variables */
    magma_int_t info = 0;
    if (n < 0)
        info = -1;
    else if (ldda < max(1,n))
        info = -3;
    else if (lddia < max(1,n))
        info = -6;

    if (info != 0) {
        magma_xerbla( __func__, -(info) );
        return info;
    }

    /* Quick return if possible */
    if ( n == 0 )
        return info;

    magma_int_t ib, j;
    magma_int_t nb = 256; //256; // BATRF_NB;

    float **dA_displ   = NULL;
    float **dW0_displ  = NULL;
    float **dW1_displ  = NULL;
    float **dW2_displ  = NULL;
    float **dW3_displ  = NULL;
    float **dW4_displ  = NULL;
    float **dinvdiagA_array = NULL;
    float **dwork_array = NULL;
    float **dW5_displ   = NULL;
    magma_malloc((void**)&dA_displ,   batchCount * sizeof(*dA_displ));
    magma_malloc((void**)&dW0_displ,  batchCount * sizeof(*dW0_displ));
    magma_malloc((void**)&dW1_displ,  batchCount * sizeof(*dW1_displ));
    magma_malloc((void**)&dW2_displ,  batchCount * sizeof(*dW2_displ));
    magma_malloc((void**)&dW3_displ,  batchCount * sizeof(*dW3_displ));
    magma_malloc((void**)&dW4_displ,  batchCount * sizeof(*dW4_displ));
    magma_malloc((void**)&dinvdiagA_array, batchCount * sizeof(*dinvdiagA_array));
    magma_malloc((void**)&dwork_array, batchCount * sizeof(*dwork_array));
    magma_malloc((void**)&dW5_displ,  batchCount * sizeof(*dW5_displ));

    float* dinvdiagA;
    float* dwork; // dinvdiagA and dwork are workspace in strsm
    magma_int_t invdiagA_msize = magma_roundup( n, BATRI_NB )*BATRI_NB;
    magma_int_t dwork_msize = n*nb;
    magma_smalloc( &dinvdiagA, invdiagA_msize * batchCount);
    magma_smalloc( &dwork, dwork_msize * batchCount );
    /* check allocation */
    if ( dA_displ  == NULL || dW1_displ == NULL || dW2_displ       == NULL || dW3_displ   == NULL || 
         dW4_displ == NULL || dW5_displ  == NULL || dinvdiagA_array == NULL || dwork_array == NULL || 
         dinvdiagA == NULL || dwork     == NULL ) {
        magma_free(dA_displ);
        magma_free(dW1_displ);
        magma_free(dW2_displ);
        magma_free(dW3_displ);
        magma_free(dW4_displ);
        magma_free(dW5_displ);
        magma_free(dinvdiagA_array);
        magma_free(dwork_array);
        magma_free(dinvdiagA);
        magma_free( dwork );
        info = MAGMA_ERR_DEVICE_ALLOC;
        return info;
    }

    magmablas_slaset( MagmaFull, invdiagA_msize, batchCount, MAGMA_S_ZERO, MAGMA_S_ZERO, dinvdiagA, invdiagA_msize, queue );
    magmablas_slaset( MagmaFull, dwork_msize, batchCount, MAGMA_S_ZERO, MAGMA_S_ZERO, dwork, dwork_msize, queue );
    magma_sset_pointer( dwork_array, dwork, n, 0, 0, dwork_msize, batchCount, queue );
    magma_sset_pointer( dinvdiagA_array, dinvdiagA, STRTRI_BATCHED_NB, 0, 0, invdiagA_msize, batchCount, queue );

    magma_sdisplace_pointers(dA_displ, dA_array, ldda, 0, 0, batchCount, queue);
    // set dinvdiagA to identity
    magmablas_slaset_batched( MagmaFull, n, n, MAGMA_S_ZERO, MAGMA_S_ONE, dinvA_array, lddia, batchCount, queue );

    for (j = 0; j < n; j += nb) {
        ib = min(nb, n-j);
        // dinvdiagA * Piv' = I * U^-1 * L^-1 = U^-1 * L^-1 * I
        // Azzam : optimization can be done:
        //          2- compute invdiagL invdiagU only one time

        //printf(" @ step %d calling solve 1\n",j);
        // solve dwork = L^-1 * I
        magmablas_slaset_batched( MagmaFull, j, ib, MAGMA_S_ZERO, MAGMA_S_ZERO, dwork_array, n, batchCount, queue );
        magma_sdisplace_pointers(dW5_displ, dwork_array, n, j, 0, batchCount, queue);
        magma_sdisplace_pointers(dW0_displ, dinvA_array, lddia, j, j, batchCount, queue);
        magma_sdisplace_pointers(dA_displ, dA_array, ldda, j, j, batchCount, queue);
        
        magmablas_strsm_outofplace_batched( MagmaLeft, MagmaLower, MagmaNoTrans, MagmaUnit, 1,
                n-j, ib,
                MAGMA_S_ONE,
                dA_displ,       ldda, // dA
                dW0_displ,   lddia, // dB
                dW5_displ,        n, // dX //output
                dinvdiagA_array,  invdiagA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue );
        
        //printf(" @ step %d calling solve 2\n",j);
        // solve dinvdiagA = U^-1 * dwork
        magma_sdisplace_pointers(dW5_displ, dwork_array, n, 0, 0, batchCount, queue);
        magma_sdisplace_pointers(dW0_displ, dinvA_array, lddia, 0, j, batchCount, queue);
        magma_sdisplace_pointers(dA_displ, dA_array, ldda, 0, 0, batchCount, queue);
        magmablas_strsm_outofplace_batched( MagmaLeft, MagmaUpper, MagmaNoTrans, MagmaNonUnit, 1,
                n, ib,
                MAGMA_S_ONE,
                dA_displ,       ldda, // dA
                dW5_displ,        n, // dB 
                dW0_displ,   lddia, // dX //output
                dinvdiagA_array,  invdiagA_msize, 
                dW1_displ,   dW2_displ, 
                dW3_displ,   dW4_displ,
                1, batchCount, queue );
    }

    // Apply column interchanges
    magma_slaswp_columnserial_batched( n, dinvA_array, lddia, max(1,n-1), 1, dipiv_array, batchCount, queue );

    magma_queue_sync(queue);

    magma_free(dA_displ);
    magma_free(dW0_displ);
    magma_free(dW1_displ);
    magma_free(dW2_displ);
    magma_free(dW3_displ);
    magma_free(dW4_displ);
    magma_free(dW5_displ);
    magma_free(dinvdiagA_array);
    magma_free(dwork_array);
    magma_free(dinvdiagA);
    magma_free( dwork );
    
    return info;
}
