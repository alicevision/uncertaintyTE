/*
   -- MAGMA (version 2.2.0) --
   Univ. of Tennessee, Knoxville
   Univ. of California, Berkeley
   Univ. of Colorado, Denver
   @date November 2016

   @author Azzam Haidar
   @author Tingxing Dong

   @generated from src/zgetrf_panel_batched.cpp, normal z -> c, Sun Nov 20 20:20:26 2016
*/
#include "magma_internal.h"

/***************************************************************************//**
    Purpose
    -------
    This is an internal routine that might have many assumption.
    Documentation is not fully completed

    CGETRF_PANEL computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.

    The factorization has the form
        A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 3 BLAS version of the algorithm.

    This is a batched version that factors batchCount M-by-N matrices in parallel.
    dA, ipiv, and info become arrays with one entry per matrix.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of each matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of each matrix A.  N >= 0.

    @param[in]
    min_recpnb   INTEGER.
                 Internal use. The recursive nb

    @param[in,out]
    dA_array    Array of pointers, dimension (batchCount).
            Each is a COMPLEX array on the GPU, dimension (LDDA,N).
            On entry, each pointer is an M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    @param[in]
    ldda    INTEGER
            The leading dimension of each array A.  LDDA >= max(1,M).

    @param[out]
    dipiv_array  Array of pointers, dimension (batchCount), for corresponding matrices.
            Each is an INTEGER array, dimension (min(M,N))
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    @param[out]
    dpivinfo_array  Array of pointers, dimension (batchCount), for internal use.

    @param[in,out]
    dX_array       Array of pointers, dimension (batchCount).
             Each is a COMPLEX array X of dimension ( lddx, n ).
             On entry, should be set to 0
             On exit, the solution matrix X

    @param[in]
    dX_length    INTEGER.
                 The size of each workspace matrix dX

    @param[in,out]
    dinvA_array    Array of pointers, dimension (batchCount).
            Each is a COMPLEX array dinvA, a workspace on device.
            If side == MagmaLeft,  dinvA must be of size >= ceil(m/CTRTRI_BATCHED_NB)*CTRTRI_BATCHED_NB*CTRTRI_BATCHED_NB,
            If side == MagmaRight, dinvA must be of size >= ceil(n/CTRTRI_BATCHED_NB)*CTRTRI_BATCHED_NB*CTRTRI_BATCHED_NB,

    @param[in]
    dinvA_length    INTEGER
                   The size of each workspace matrix dinvA
    @param[in]
    dW1_displ  Workspace array of pointers, for internal use.

    @param[in]
    dW2_displ  Workspace array of pointers, for internal use.

    @param[in]
    dW3_displ  Workspace array of pointers, for internal use.

    @param[in]
    dW4_displ  Workspace array of pointers, for internal use.

    @param[in]
    dW5_displ  Workspace array of pointers, for internal use.

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
    gbstep  INTEGER
            internal use.

    @param[in]
    batchCount  INTEGER
                The number of matrices to operate on.

    @param[in]
    queue   magma_queue_t
            Queue to execute in.

    @ingroup magma_getrf_batched
*******************************************************************************/
extern "C" magma_int_t
magma_cgetrf_recpanel_batched(
    magma_int_t m, magma_int_t n, magma_int_t min_recpnb,    
    magmaFloatComplex** dA_array,    magma_int_t ldda,
    magma_int_t** dipiv_array, magma_int_t** dpivinfo_array,
    magmaFloatComplex** dX_array,    magma_int_t dX_length,
    magmaFloatComplex** dinvA_array, magma_int_t dinvA_length,
    magmaFloatComplex** dW1_displ, magmaFloatComplex** dW2_displ,  
    magmaFloatComplex** dW3_displ, magmaFloatComplex** dW4_displ,
    magmaFloatComplex** dW5_displ,
    magma_int_t *info_array, magma_int_t gbstep,  
    magma_int_t batchCount,  magma_queue_t queue)
{
    //magma_int_t DEBUG = 3;
    // Quick return if possible
    if (m == 0 || n == 0) {
        return 0;
    }


    magmaFloatComplex **dA_displ  = NULL;
    magma_malloc((void**)&dA_displ,   batchCount * sizeof(*dA_displ));
    magma_int_t **dipiv_displ = NULL;
    magma_malloc((void**)&dipiv_displ, batchCount * sizeof(*dipiv_displ));
    
    magma_int_t panel_nb = n;
    if (panel_nb <= min_recpnb) {
        //if (DEBUG > 0)printf("calling bottom panel recursive with m=%d nb=%d\n",m,n);
        //  panel factorization
        //magma_cdisplace_pointers(dA_displ, dA_array, ldda, 0, 0, batchCount);
        magma_cgetf2_batched(m, panel_nb,
                           dA_array, ldda,
                           dW1_displ, dW2_displ, dW3_displ,
                           dipiv_array, info_array, gbstep, batchCount, queue);
    }
    else {
        // split A over two [A A2]
        // panel on A1, update on A2 then panel on A1    
        magma_int_t n1 = n/2;
        magma_int_t n2 = n-n1;
        magma_int_t m1 = m;
        magma_int_t m2 = m-n1;
        magma_int_t p1 = 0;
        magma_int_t p2 = n1;
        // panel on A1
        //if (DEBUG > 0)printf("calling recursive panel on A1 with m=%d nb=%d min_recpnb %d\n",m1,n1,min_recpnb);
        magma_cdisplace_pointers(dA_displ, dA_array, ldda, p1, p1, batchCount, queue); 
        magma_idisplace_pointers(dipiv_displ, dipiv_array, 1, p1, 0, batchCount, queue);
        magma_cgetrf_recpanel_batched(
                           m1, n1, min_recpnb,
                           dA_displ, ldda,
                           dipiv_displ, dpivinfo_array,
                           dX_array, dX_length,
                           dinvA_array, dinvA_length,
                           dW1_displ, dW2_displ,
                           dW3_displ, dW4_displ, dW5_displ,
                           info_array, gbstep, batchCount, queue);

        // update A2
        //if (DEBUG > 0)printf("calling TRSM  with             m=%d n=%d\n",m1,n2);
        
        // setup pivinfo 
        setup_pivinfo_batched(dpivinfo_array, dipiv_displ, m1, n1, batchCount, queue);
        magma_cdisplace_pointers(dW5_displ, dA_array, ldda, p1, p2, batchCount, queue); 
        magma_claswp_rowparallel_batched( n2, dW5_displ, ldda,
                           dX_array, n1,
                           0, n1,
                           dpivinfo_array, batchCount, queue );
        magmablas_ctrsm_outofplace_batched( MagmaLeft, MagmaLower, MagmaNoTrans, MagmaUnit, 1,
                              n1, n2,
                              MAGMA_C_ONE,
                              dA_displ,    ldda, // dA
                              dX_array,  n1, // dB
                              dW5_displ,   ldda, // dX
                              dinvA_array, dinvA_length,
                              dW1_displ,   dW2_displ, 
                              dW3_displ,   dW4_displ,
                              0, batchCount, queue );

        magma_cdisplace_pointers(dW1_displ, dA_array, ldda, p2, 0, batchCount, queue); 
        magma_cdisplace_pointers(dA_displ, dA_array, ldda, p2, p2, batchCount, queue); 

        //if (DEBUG > 0)printf("calling update A2(%d,%d) -= A(%d,%d)*A(%d,%d)  with             m=%d n=%d k=%d ldda %d\n",p2,p2,p2,0,p1,p2,m2,n2,n1,ldda);

        magma_cgemm_batched( MagmaNoTrans, MagmaNoTrans, m2, n2, n1, 
                             MAGMA_C_NEG_ONE, dW1_displ, ldda, 
                             dW5_displ, ldda, 
                             MAGMA_C_ONE,  dA_displ, ldda, 
                             batchCount, queue );
        // panel on A2
        //if (DEBUG > 0)printf("calling recursive panel on A2 with m=%d nb=%d min_recpnb %d\n",m2,n2,min_recpnb);
        magma_idisplace_pointers(dipiv_displ, dipiv_array, 1, p2, 0, batchCount, queue);
        magma_cgetrf_recpanel_batched(
                           m2, n2, min_recpnb,
                           dA_displ, ldda,
                           dipiv_displ, dpivinfo_array,
                           dX_array, dX_length,
                           dinvA_array, dinvA_length,
                           dW1_displ, dW2_displ,
                           dW3_displ, dW4_displ, dW5_displ,
                           info_array, gbstep+p2, batchCount, queue);

        // setup pivinfo
        setup_pivinfo_batched(dpivinfo_array, dipiv_displ, m2, n2, batchCount, queue);
        adjust_ipiv_batched(dipiv_displ, n2, n1, batchCount, queue);
        
        magma_cdisplace_pointers(dW1_displ, dA_array, ldda, p2, 0, batchCount, queue); // no need since it is above
        magma_claswp_rowparallel_batched( n1, dW1_displ, ldda,
                           dW1_displ, ldda,
                           n1, n,
                           dpivinfo_array, batchCount, queue );
    }

    magma_free(dA_displ);
    magma_free(dipiv_displ);
    return 0;
}
