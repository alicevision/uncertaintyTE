/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @generated from src/zgelqf_gpu.cpp, normal z -> d, Sun Nov 20 20:20:21 2016

*/
#include "magma_internal.h"

#define REAL

/***************************************************************************//**
    Purpose
    -------
    DGELQF computes an LQ factorization of a DOUBLE PRECISION M-by-N matrix dA:
    dA = L * Q.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in,out]
    dA      DOUBLE PRECISION array on the GPU, dimension (LDDA,N)
            On entry, the M-by-N matrix dA.
            On exit, the elements on and below the diagonal of the array
            contain the m-by-min(m,n) lower trapezoidal matrix L (L is
            lower triangular if m <= n); the elements above the diagonal,
            with the array TAU, represent the orthogonal matrix Q as a
            product of elementary reflectors (see Further Details).

    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).

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
            The dimension of the array WORK.  LWORK >= max(1,M).
            For optimum performance LWORK >= M*NB, where NB is the
            optimal blocksize.
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

       Q = H(k) . . . H(2) H(1), where k = min(m,n).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a real scalar, and v is a real vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
    and tau in TAU(i).

    @ingroup magma_gelqf
*******************************************************************************/
extern "C" magma_int_t
magma_dgelqf_gpu(
    magma_int_t m, magma_int_t n,
    magmaDouble_ptr dA, magma_int_t ldda,
    double *tau,
    double *work, magma_int_t lwork,
    magma_int_t *info)
{
    /* Constants */
    const double c_one = MAGMA_D_ONE;
    const magma_int_t ione = 1;
    MAGMA_UNUSED( ione );  // used only for real

    /* Local variables */
    magmaDouble_ptr dAT=NULL;
    magma_int_t min_mn, maxm, maxn, nb;
    magma_int_t iinfo;

    *info = 0;
    nb = magma_get_dgelqf_nb( m, n );
    min_mn = min( m, n );

    work[0] = magma_dmake_lwork( m*nb );
    bool lquery = (lwork == -1);
    if (m < 0) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,m)) {
        *info = -4;
    } else if (lwork < max(1,m) && ! lquery) {
        *info = -7;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
        return *info;
    }

    /*  Quick return if possible */
    if (min_mn == 0) {
        work[0] = c_one;
        return *info;
    }

    maxm = magma_roundup( m, 32 );
    maxn = magma_roundup( n, 32 );

    magma_int_t lddat = maxn;

    magma_queue_t queue;
    magma_device_t cdev;
    magma_getdevice( &cdev );
    magma_queue_create( cdev, &queue );
    
    if ( m == n ) {
        dAT = dA;
        lddat = ldda;
        magmablas_dtranspose_inplace( m, dAT, ldda, queue );
    }
    else {
        if (MAGMA_SUCCESS != magma_dmalloc( &dAT, maxm*maxn ) ) {
            *info = MAGMA_ERR_DEVICE_ALLOC;
            goto cleanup;
        }
        
        magmablas_dtranspose( m, n, dA, ldda, dAT, lddat, queue );
    }
    
    magma_dgeqrf2_gpu( n, m, dAT, lddat, tau, &iinfo );
    assert( iinfo >= 0 );
    if ( iinfo > 0 ) {
        *info = iinfo;
    }
    
    // conjugate tau
    #ifdef COMPLEX
    lapackf77_dlacgv( &min_mn, tau, &ione );
    #endif
    
    if ( m == n ) {
        magmablas_dtranspose_inplace( m, dAT, lddat, queue );
    }
    else {
        magmablas_dtranspose( n, m, dAT, lddat, dA, ldda, queue );
        magma_free( dAT );
    }

cleanup:
    magma_queue_destroy( queue );
    
    return *info;
} /* magma_dgelqf_gpu */
