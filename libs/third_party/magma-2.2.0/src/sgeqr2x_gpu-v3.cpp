/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016
       
       @author Stan Tomov

       @generated from src/zgeqr2x_gpu-v3.cpp, normal z -> s, Sun Nov 20 20:20:21 2016

*/
#include "magma_internal.h"

/******************************************************************************/
// TODO: how does this differ from larfb_gpu?
extern "C" magma_int_t
magma_slarfb2_gpu(
    magma_int_t m, magma_int_t n, magma_int_t k,
    magmaFloat_const_ptr dV,    magma_int_t lddv,
    magmaFloat_const_ptr dT,    magma_int_t lddt,
    magmaFloat_ptr       dC,    magma_int_t lddc,
    magmaFloat_ptr       dwork, magma_int_t ldwork,
    magma_queue_t queue )
{
    float c_zero    = MAGMA_S_ZERO;
    float c_one     = MAGMA_S_ONE;
    float c_neg_one = MAGMA_S_NEG_ONE;

    if (m <= 0 || n <= 0)
        return MAGMA_SUCCESS;

    // W = C^H V
    magma_sgemm( MagmaConjTrans, MagmaNoTrans,
    //magmablas_sgemm_reduce(
                           n, k, m,
                           c_one,  dC,    lddc,
                                   dV,    lddv,
                           c_zero, dwork, ldwork, queue );

    // W = W T^H = C^H V T^H
    magma_strmm( MagmaRight, MagmaUpper, MagmaNoTrans, MagmaNonUnit,
                 n, k,
                 c_one, dT,    lddt,
                        dwork, ldwork, queue );

    // C = C - V W^H = C - V T V^H C = (I - V T V^H) C = H C
    magma_sgemm( MagmaNoTrans, MagmaConjTrans,
                 m, n, k,
                 c_neg_one, dV,    lddv,
                            dwork, ldwork,
                 c_one,     dC,    lddc, queue );
    
    return MAGMA_SUCCESS;
}


/***************************************************************************//**
    Purpose
    -------
    SGEQR2 computes a QR factorization of a real m by n matrix A:
    A = Q * R.

    This expert routine requires two more arguments than the standard
    sgeqr2, namely, dT and ddA, explained below. The storage for A is
    also not as in the LAPACK's sgeqr2 routine (see below).

    The first is used to output the triangular
    n x n factor T of the block reflector used in the factorization.
    The second holds the diagonal nxn blocks of A, i.e., the diagonal
    submatrices of R. This routine implements the left looking QR.

    This version adds internal blocking.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in,out]
    dA      REAL array, dimension (LDDA,N)
            On entry, the m by n matrix A.
            On exit, the orthogonal matrix Q as a
            product of elementary reflectors (see Further Details).
    \n
            the elements on and above the diagonal of the array
            contain the min(m,n) by n upper trapezoidal matrix R (R is
            upper triangular if m >= n); the elements below the diagonal,
            with the array TAU, represent the orthogonal matrix Q as a
            product of elementary reflectors (see Further Details).

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

    @param[out]
    dtau    REAL array, dimension (min(M,N))
            The scalar factors of the elementary reflectors (see Further
            Details).

    @param[out]
    dT      REAL array, dimension N x N.
            Stores the triangular N x N factor T of the block reflector
            used in the factorization. The lower triangular part is 0.

    @param[out]
    ddA     REAL array, dimension N x N.
            Stores the elements of the upper N x N diagonal block of A.
            LAPACK stores this array in A. There are 0s below the diagonal.

    @param
    dwork   (workspace) REAL array, dimension (3 N)

    @param[out]
    info    INTEGER
      -     = 0: successful exit
      -     < 0: if INFO = -i, the i-th argument had an illegal value

    Further Details
    ---------------
    The matrix Q is represented as a product of elementary reflectors

       Q = H(1) H(2) . . . H(k), where k = min(m,n).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a real scalar, and v is a real vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
    and tau in TAU(i).

    @ingroup magma_geqr2
*******************************************************************************/
extern "C" magma_int_t
magma_sgeqr2x3_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloat_ptr dA, magma_int_t ldda,
    magmaFloat_ptr dtau,
    magmaFloat_ptr dT,
    magmaFloat_ptr ddA,
    magmaFloat_ptr        dwork,
    magma_int_t *info)
{
    #define dA(i_,j_) (dA + (i_) + (j_)*ldda)
    #define BLOCK_SIZE 32

    magma_int_t b, i, min_mn;

    magmaFloat_ptr dnorm = dwork;
    magmaFloat_ptr dwork2 = (magmaFloat_ptr)(dwork + 2*n);

    *info = 0;
    if (m < 0) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,m)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_xerbla( __func__, -(*info) );
        return *info;
    }

    magma_queue_t queue;
    magma_device_t cdev;
    magma_getdevice( &cdev );
    magma_queue_create( cdev, &queue );

    /* Compute the norms of the trailing columns */
    min_mn = min(m,n);
    // magmablas_snrm2_cols( m, min_mn, dA(0,0), ldda, dnorm, queue );

    for (b=0; b < min_mn; b += BLOCK_SIZE) {
        for (i = b; i < min(min_mn, b+BLOCK_SIZE); ++i) {
            /*   Apply H' to A(:,i) from the left */
            if ( i-b > 0)
                magma_slarfbx_gpu( m-b, i-b, dA(b, b), ldda,
                                  dT+b+b*min_mn, min_mn, dA(b, i), dwork2, queue );

            /*   Adjust the dnorm[i] to hold the norm of A(i:m,i) */
            //if ( i > 0 )
            //    magmablas_snrm2_adjust( i, dnorm+i, dA(0, i), queue );
            magmablas_snrm2_cols( m-i, 1, dA(i,i), ldda, dnorm+i, queue );
            
            /*  Generate elementary reflector H(i) to annihilate A(i+1:m,i)
                1. 1 is not yet put on the diagonal of A
                2. Elements above the diagonal are copied in ddA and
                   the ones in A are set to zero
                3. update T */
            magma_slarfgtx_gpu(m-i, dA(i, i), dA(min(i+1,m), i), dtau+i,
                               dnorm+i, ddA + i + i*(n), i,
                               dA(i,0), ldda,  dT, min_mn, dwork2, queue);
        }
        
        /* Apply the transformations to the trailing matrix. */
        //magma_slarfb2_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
        magma_slarfb2_gpu(
                           m-b, min_mn-i, BLOCK_SIZE,
                           dA(b, b), ldda, dT+b+b*min_mn, min_mn,
                           dA(b, i), ldda, dwork2, min_mn-i, queue );
    }

    magma_queue_destroy( queue );

    return *info;
} /* magma_sgeqr2 */
