/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @author Tingxing Dong
       @author Azzam Haidar

       @generated from magmablas/ztrsv.cu, normal z -> c, Sun Nov 20 20:20:30 2016
*/

#include "magma_internal.h"
#include "magma_templates.h"


#define PRECISION_c

#define NB 256  //NB is the 1st level blocking in recursive blocking, NUM_THREADS is the 2ed level, NB=256, NUM_THREADS=64 is optimal for batched

#define NUM_THREADS 128 //64 //128

#define BLOCK_SIZE_N 128
#define DIM_X_N 128
#define DIM_Y_N 1

#define BLOCK_SIZE_T 32
#define DIM_X_T 16
#define DIM_Y_T 8

#include "ctrsv_template_device.cuh"

#define A(i, j)  (A + (i) + (j)*lda)   // A(i, j) means at i row, j column

extern __shared__ magmaFloatComplex shared_data[];


/******************************************************************************/
template< const int BLOCK_SIZE, const int DIM_X, const int DIM_Y,
          const int TILE_SIZE, const int flag, const magma_uplo_t uplo,
          const magma_trans_t trans, const magma_diag_t diag >
__global__ void
ctrsv_notrans_kernel_outplace(
    int n,
    const magmaFloatComplex * __restrict__ A, int lda,
    magmaFloatComplex *b, int incb,
    magmaFloatComplex *x)
{
    ctrsv_notrans_device< BLOCK_SIZE, DIM_X, DIM_Y, TILE_SIZE, flag, uplo, trans, diag >( n, A, lda, b, incb, x);
}


/******************************************************************************/
template< const int BLOCK_SIZE, const int DIM_X, const int DIM_Y,
          const int TILE_SIZE, const int flag, const magma_uplo_t uplo,
          const magma_trans_t trans, const magma_diag_t diag >
__global__ void
ctrsv_trans_kernel_outplace(
    int n,
    const magmaFloatComplex * __restrict__ A, int lda,
    magmaFloatComplex *b, int incb,
    magmaFloatComplex *x)
{
    ctrsv_trans_device< BLOCK_SIZE, DIM_X, DIM_Y, TILE_SIZE, flag, uplo, trans, diag >( n, A, lda, b, incb, x);
}


/******************************************************************************/
extern "C" void
magmablas_ctrsv_outofplace(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t n,
    magmaFloatComplex_const_ptr A, magma_int_t lda,
    magmaFloatComplex_ptr       b, magma_int_t incb,
    magmaFloatComplex_ptr       x,
    magma_queue_t queue,
    magma_int_t flag=0)
{
    /* Check arguments */
    magma_int_t info = 0;
    if ( uplo != MagmaUpper && uplo != MagmaLower ) {
        info = -1;
    } else if ( trans != MagmaNoTrans && trans != MagmaTrans && trans != MagmaConjTrans ) {
        info = -2;
    } else if ( diag != MagmaUnit && diag != MagmaNonUnit ) {
        info = -3;
    } else if (n < 0) {
        info = -5;
    } else if (lda < max(1,n)) {
        info = -8;
    }

    if (info != 0) {
        magma_xerbla( __func__, -(info) );
        return;
    }

    // quick return if possible.
    if (n == 0)
        return;

    dim3 threads( NUM_THREADS );
    dim3 blocks( 1, 1, 1 );
    size_t shmem = n * sizeof(magmaFloatComplex);

    if (trans == MagmaNoTrans)
    {
        if (uplo == MagmaUpper)
        {
            if (diag == MagmaNonUnit)
            {
                if (flag == 0) {
                    ctrsv_notrans_kernel_outplace< BLOCK_SIZE_N, DIM_X_N, DIM_Y_N, MagmaBigTileSize, 0, MagmaUpper, MagmaNoTrans, MagmaNonUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
                else {
                    ctrsv_notrans_kernel_outplace< BLOCK_SIZE_N, DIM_X_N, DIM_Y_N, MagmaBigTileSize, 1, MagmaUpper, MagmaNoTrans, MagmaNonUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
            }
            else if (diag == MagmaUnit)
            {
                if (flag == 0) {
                    ctrsv_notrans_kernel_outplace< BLOCK_SIZE_N, DIM_X_N, DIM_Y_N, MagmaBigTileSize, 0, MagmaUpper, MagmaNoTrans, MagmaUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
                else {
                    ctrsv_notrans_kernel_outplace< BLOCK_SIZE_N, DIM_X_N, DIM_Y_N, MagmaBigTileSize, 1, MagmaUpper, MagmaNoTrans, MagmaUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
            }
        }
        else //Lower
        {
            if (diag == MagmaNonUnit)
            {
                if (flag == 0)
                {
                    ctrsv_notrans_kernel_outplace< BLOCK_SIZE_N, DIM_X_N, DIM_Y_N, MagmaBigTileSize, 0, MagmaLower, MagmaNoTrans, MagmaNonUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
                else {
                    ctrsv_notrans_kernel_outplace< BLOCK_SIZE_N, DIM_X_N, DIM_Y_N, MagmaBigTileSize, 1, MagmaLower, MagmaNoTrans, MagmaNonUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
            }
            else if (diag == MagmaUnit)
            {
                if (flag == 0)
                {
                    ctrsv_notrans_kernel_outplace< BLOCK_SIZE_N, DIM_X_N, DIM_Y_N, MagmaBigTileSize, 0, MagmaLower, MagmaNoTrans, MagmaUnit>
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
                else {
                    ctrsv_notrans_kernel_outplace< BLOCK_SIZE_N, DIM_X_N, DIM_Y_N, MagmaBigTileSize, 1, MagmaLower, MagmaNoTrans, MagmaUnit>
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
            }
        }
    }
    else if (trans == MagmaTrans)
    {
        if (uplo == MagmaUpper)
        {
            if (diag == MagmaNonUnit) {
                if (flag == 0)
                {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T, MagmaBigTileSize, 0, MagmaUpper, MagmaTrans, MagmaNonUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
                else {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T, MagmaBigTileSize, 1, MagmaUpper, MagmaTrans, MagmaNonUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
            }
            else if (diag == MagmaUnit) {
                if (flag == 0)
                {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T, MagmaBigTileSize, 0, MagmaUpper, MagmaTrans, MagmaUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
                else {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T, MagmaBigTileSize, 1, MagmaUpper, MagmaTrans, MagmaUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
            }
        }
        else
        {
            if (diag == MagmaNonUnit) {
                if (flag == 0)
                {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T,MagmaBigTileSize, 0, MagmaLower, MagmaTrans, MagmaNonUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
                else {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T, MagmaBigTileSize, 1, MagmaLower, MagmaTrans, MagmaNonUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
            }
            else if (diag == MagmaUnit) {
                if (flag == 0)
                {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T,MagmaBigTileSize, 0, MagmaLower, MagmaTrans, MagmaUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
                else {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T, MagmaBigTileSize, 1, MagmaLower, MagmaTrans, MagmaUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
            }
        }
    }
    else if (trans == MagmaConjTrans)
    {
        if (uplo == MagmaUpper)
        {
            if (diag == MagmaNonUnit) {
                if (flag == 0)
                {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T, MagmaBigTileSize, 0, MagmaUpper, MagmaConjTrans, MagmaNonUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
                else {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T, MagmaBigTileSize, 1, MagmaUpper, MagmaConjTrans, MagmaNonUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
            }
            else if (diag == MagmaUnit) {
                if (flag == 0)
                {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T, MagmaBigTileSize, 0, MagmaUpper, MagmaConjTrans, MagmaUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
                else {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T, MagmaBigTileSize, 1, MagmaUpper, MagmaConjTrans, MagmaUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
            }
        }
        else
        {
            if (diag == MagmaNonUnit) {
                if (flag == 0)
                {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T,MagmaBigTileSize, 0, MagmaLower, MagmaConjTrans, MagmaNonUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
                else {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T, MagmaBigTileSize, 1, MagmaLower, MagmaConjTrans, MagmaNonUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
            }
            else if (diag == MagmaUnit) {
                if (flag == 0)
                {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T,MagmaBigTileSize, 0, MagmaLower, MagmaConjTrans, MagmaUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
                else {
                    ctrsv_trans_kernel_outplace< BLOCK_SIZE_T, DIM_X_T, DIM_Y_T, MagmaBigTileSize, 1, MagmaLower, MagmaConjTrans, MagmaUnit >
                        <<< blocks, threads, shmem, queue->cuda_stream() >>>
                        (n, A, lda, b, incb, x);
                }
            }
        }
    }
}


/******************************************************************************/
/*
    README: flag decides if the ctrsv_outplace see an updated x or not. 0: No; other: Yes
    In recursive, flag must be nonzero except the 1st call
*/
extern "C" void
magmablas_ctrsv_recursive_outofplace(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t n,
    magmaFloatComplex_const_ptr A, magma_int_t lda,
    magmaFloatComplex_ptr       b, magma_int_t incb,
    magmaFloatComplex_ptr       x,
    magma_queue_t queue)
{
    /* Check arguments */
    magma_int_t info = 0;
    if ( uplo != MagmaUpper && uplo != MagmaLower ) {
        info = -1;
    } else if ( trans != MagmaNoTrans && trans != MagmaTrans && trans != MagmaConjTrans ) {
        info = -2;
    } else if ( diag != MagmaUnit && diag != MagmaNonUnit ) {
        info = -3;
    } else if (n < 0) {
        info = -5;
    } else if (lda < max(1,n)) {
        info = -8;
    }

    if (info != 0) {
        magma_xerbla( __func__, -(info) );
        return;
    }

    // quick return if possible.
    if (n == 0)
        return;

    //Init x with zero
    //magmablas_claset( MagmaFull, n, incb, MAGMA_C_ZERO, MAGMA_C_ZERO, x, n, queue );

    magma_int_t col = n;

    if (trans == MagmaNoTrans)
    {
        for (magma_int_t i=0; i < n; i+= NB)
        {
            magma_int_t jb = min(NB, n-i);

            if (uplo == MagmaUpper)
            {
                col -= jb;
                //assume x_array contains zero elements, magmablas_cgemv will cause slow down
                magma_cgemv( MagmaNoTrans, jb, i, MAGMA_C_ONE, A(col, col+jb), lda,
                             x+col+jb, 1, MAGMA_C_ONE, x+col, 1, queue );
            }
            else
            {
                col = i;
                magma_cgemv( MagmaNoTrans, jb, i, MAGMA_C_ONE, A(col, 0), lda,
                             x, 1, MAGMA_C_ONE, x+col, 1, queue );
            }

            magmablas_ctrsv_outofplace( uplo, trans, diag, jb, A(col, col), lda, b+col, incb, x+col, queue, i );
        }
    }
    else
    {
        for (magma_int_t i=0; i < n; i += NB)
        {
            magma_int_t jb = min(NB, n-i);

            if (uplo == MagmaLower)
            {
                col -= jb;

                magma_cgemv( MagmaConjTrans, i, jb, MAGMA_C_ONE, A(col+jb, col), lda, x+col+jb, 1, MAGMA_C_ONE, x+col, 1, queue );
            }
            else
            {
                col = i;
                
                magma_cgemv( MagmaConjTrans, i, jb, MAGMA_C_ONE, A(0, col), lda, x, 1, MAGMA_C_ONE, x+col, 1, queue );
            }
     
            magmablas_ctrsv_outofplace( uplo, trans, diag, jb, A(col, col), lda, b+col, incb, x+col, queue, i );
        }
    }
}


/***************************************************************************//**
    Purpose
    -------
    ctrsv solves one of the matrix equations on gpu

        op(A)*x = B,   or
        x*op(A) = B,

    where alpha is a scalar, X and B are vectors, A is a unit, or
    non-unit, upper or lower triangular matrix and op(A) is one of

        op(A) = A,    or
        op(A) = A^T,  or
        op(A) = A^H.

    The vector x is overwritten on b.

    Arguments
    ----------
    @param[in]
    uplo    magma_uplo_t.
            On entry, uplo specifies whether the matrix A is an upper or
            lower triangular matrix as follows:
      -     = MagmaUpper:  A is an upper triangular matrix.
      -     = MagmaLower:  A is a  lower triangular matrix.

    @param[in]
    trans  magma_trans_t.
            On entry, trans specifies the form of op(A) to be used in
            the matrix multiplication as follows:
      -     = MagmaNoTrans:    op(A) = A.
      -     = MagmaTrans:      op(A) = A^T.
      -     = MagmaConjTrans:  op(A) = A^H.

    @param[in]
    diag    magma_diag_t.
            On entry, diag specifies whether or not A is unit triangular
            as follows:
      -     = MagmaUnit:     A is assumed to be unit triangular.
      -     = MagmaNonUnit:  A is not assumed to be unit triangular.

    @param[in]
    n       INTEGER.
            On entry, n N specifies the order of the matrix A. n >= 0.

    @param[in]
    dA      COMPLEX array of dimension ( lda, n )
            Before entry with uplo = MagmaUpper, the leading n by n
            upper triangular part of the array A must contain the upper
            triangular matrix and the strictly lower triangular part of
            A is not referenced.
            Before entry with uplo = MagmaLower, the leading n by n
            lower triangular part of the array A must contain the lower
            triangular matrix and the strictly upper triangular part of
            A is not referenced.
            Note that when diag = MagmaUnit, the diagonal elements of
            A are not referenced either, but are assumed to be unity.

    @param[in]
    ldda    INTEGER.
            On entry, lda specifies the first dimension of A.
            lda >= max( 1, n ).

    @param[in]
    db      COMPLEX array of dimension  n
            On exit, b is overwritten with the solution vector X.

    @param[in]
    incb    INTEGER.
            On entry,  incb specifies the increment for the elements of
            b. incb must not be zero.
            Unchanged on exit.

    @param[in]
    queue   magma_queue_t
            Queue to execute in.

    @ingroup magma_trsv
*******************************************************************************/
extern "C" void
magmablas_ctrsv(
    magma_uplo_t uplo, magma_trans_t trans, magma_diag_t diag,
    magma_int_t n,
    magmaFloatComplex_const_ptr dA, magma_int_t ldda,
    magmaFloatComplex_ptr db, magma_int_t incb,
    magma_queue_t queue)
{
    magma_int_t size_x = n * incb;

    magmaFloatComplex_ptr dx=NULL;

    magma_cmalloc( &dx, size_x );

    magmablas_claset( MagmaFull, n, 1, MAGMA_C_ZERO, MAGMA_C_ZERO, dx, n, queue );

    magmablas_ctrsv_recursive_outofplace( uplo, trans, diag, n, dA, ldda, db, incb, dx, queue );

    magmablas_clacpy( MagmaFull, n, 1, dx, n, db, n, queue );

    magma_free( dx );
}
