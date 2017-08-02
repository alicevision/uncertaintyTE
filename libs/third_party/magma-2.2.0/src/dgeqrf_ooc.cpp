/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @generated from src/zgeqrf_ooc.cpp, normal z -> d, Sun Nov 20 20:20:22 2016

*/
#include <cuda_runtime.h>

#include "magma_internal.h"

/***************************************************************************//**
    Purpose
    -------
    DGEQRF_OOC computes a QR factorization of a DOUBLE PRECISION M-by-N matrix A:
    A = Q * R. This version does not require work space on the GPU
    passed as input. GPU memory is allocated in the routine.
    This is an out-of-core (ooc) version that is similar to magma_dgeqrf but
    the difference is that this version can use a GPU even if the matrix
    does not fit into the GPU memory at once.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in,out]
    A       DOUBLE PRECISION array, dimension (LDA,N)
            On entry, the M-by-N matrix A.
            On exit, the elements on and above the diagonal of the array
            contain the min(M,N)-by-N upper trapezoidal matrix R (R is
            upper triangular if m >= n); the elements below the diagonal,
            with the array TAU, represent the orthogonal matrix Q as a
            product of min(m,n) elementary reflectors (see Further
            Details).
    \n
            Higher performance is achieved if A is in pinned memory, e.g.
            allocated using magma_malloc_pinned.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    @param[out]
    tau     DOUBLE PRECISION array, dimension (min(M,N))
            The scalar factors of the elementary reflectors (see Further
            Details).

    @param[out]
    work    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.
    \n
            Higher performance is achieved if WORK is in pinned memory, e.g.
            allocated using magma_malloc_pinned.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK.  LWORK >= N*NB,
            where NB can be obtained through magma_get_dgeqrf_nb( M, N ).
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.

    Further Details
    ---------------
    The matrix Q is represented as a product of elementary reflectors

        Q = H(1) H(2) . . . H(k), where k = min(m,n).

    Each H(i) has the form

        H(i) = I - tau * v * v'

    where tau is a real scalar, and v is a real vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
    and tau in TAU(i).

    @ingroup magma_geqrf
*******************************************************************************/
extern "C" magma_int_t
magma_dgeqrf_ooc(
    magma_int_t m, magma_int_t n,
    double *A,    magma_int_t lda, double *tau,
    double *work, magma_int_t lwork,
    magma_int_t *info )
{
    #define  A(i_,j_) ( A + (i_) + (j_)*lda )
    #define dA(i_,j_) (dA + (i_) + (j_)*ldda)

    /* Constants */
    const double c_one = MAGMA_D_ONE;
    
    /* Local variables */
    magmaDouble_ptr dA, dwork;
    magma_int_t i, ib, IB, j, min_mn, lddwork, ldda, rows;

    magma_int_t nb = magma_get_dgeqrf_nb( m, n );

    magma_int_t lwkopt = n * nb;
    work[0] = magma_dmake_lwork( lwkopt );
    bool lquery = (lwork == -1);
    *info = 0;
    if (m < 0) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (lda < max(1,m)) {
        *info = -4;
    } else if (lwork < max(1,n) && ! lquery) {
        *info = -7;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
        return *info;
    }

    /* Check how much memory do we have */
    size_t freeMem, totalMem;
    cudaMemGetInfo( &freeMem, &totalMem );
    freeMem /= sizeof(double);
    
    magma_int_t NB = (magma_int_t)(0.8*freeMem/m);
    NB = (NB / nb) * nb;

    if (NB >= n)
        return magma_dgeqrf(m, n, A, lda, tau, work, lwork, info);

    min_mn = min(m,n);
    if (min_mn == 0) {
        work[0] = c_one;
        return *info;
    }

    lddwork = magma_roundup( NB, 32 ) + nb;
    ldda    = magma_roundup( m, 32 );

    if (MAGMA_SUCCESS != magma_dmalloc( &dA, (NB + nb)*ldda + nb*lddwork )) {
        *info = MAGMA_ERR_DEVICE_ALLOC;
        return *info;
    }

    magma_queue_t queues[2];
    magma_device_t cdev;
    magma_getdevice( &cdev );
    magma_queue_create( cdev, &queues[0] );
    magma_queue_create( cdev, &queues[1] );

    magmaDouble_ptr ptr = dA + ldda*NB;
    dwork = dA + ldda*(NB + nb);

    /* start the main loop over the blocks that fit in the GPU memory */
    for (i=0; i < n; i += NB) {
        IB = min( n-i, NB );
        //printf("Processing %5lld columns -- %5lld to %5lld ...\n", (long long) IB, (long long) i, (long long)(i+IB) );

        /* 1. Copy the next part of the matrix to the GPU */
        magma_dsetmatrix_async( m, IB,
                                A(0,i),  lda,
                                dA(0,0), ldda, queues[0] );
        magma_queue_sync( queues[0] );

        /* 2. Update it with the previous transformations */
        for (j=0; j < min(i,min_mn); j += nb) {
            ib = min( min_mn-j, nb );

            /* Get a panel in ptr.                                           */
            //   1. Form the triangular factor of the block reflector
            //   2. Send it to the GPU.
            //   3. Put 0s in the upper triangular part of V.
            //   4. Send V to the GPU in ptr.
            //   5. Update the matrix.
            //   6. Restore the upper part of V.
            rows = m-j;
            lapackf77_dlarft( MagmaForwardStr, MagmaColumnwiseStr,
                              &rows, &ib, A(j,j), &lda, tau+j, work, &ib);
            magma_dsetmatrix_async( ib, ib,
                                    work,  ib,
                                    dwork, lddwork, queues[1] );

            magma_dpanel_to_q( MagmaUpper, ib, A(j,j), lda, work+ib*ib );
            magma_dsetmatrix_async( rows, ib,
                                    A(j,j), lda,
                                    ptr,    rows, queues[1] );
            magma_queue_sync( queues[1] );

            magma_dlarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
                              rows, IB, ib,
                              ptr, rows, dwork,    lddwork,
                              dA(j, 0), ldda, dwork+ib, lddwork, queues[1] );

            magma_dq_to_panel( MagmaUpper, ib, A(j,j), lda, work+ib*ib );
        }

        /* 3. Do a QR on the current part */
        if (i < min_mn)
            magma_dgeqrf2_gpu( m-i, IB, dA(i,0), ldda, tau+i, info );

        /* 4. Copy the current part back to the CPU */
        magma_dgetmatrix_async( m, IB,
                                dA(0,0), ldda,
                                A(0,i),  lda, queues[0] );
    }

    magma_queue_sync( queues[0] );

    magma_queue_destroy( queues[0] );
    magma_queue_destroy( queues[1] );
    magma_free( dA );
    
    return *info;
} /* magma_dgeqrf_ooc */
