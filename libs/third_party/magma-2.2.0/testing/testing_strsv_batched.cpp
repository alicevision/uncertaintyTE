/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @generated from testing/testing_ztrsv_batched.cpp, normal z -> s, Thu Jul 27 17:29:49 2017
       @author Tingxing Dong

*/
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <cuda.h>  // for CUDA_VERSION

// includes, project
#include "flops.h"
#include "magma_v2.h"
#include "magma_lapack.h"
#include "testings.h"

#if defined(_OPENMP)
#include <omp.h>
#include "../control/magma_threadsetting.h"  // internal header
#endif


#define h_A(i,j,s) (h_A + (i) + (j)*lda + (s)*lda*Ak)


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing strsm_batched
*/
int main( int argc, char** argv)
{
    TESTING_CHECK( magma_init() );
    magma_print_environment();

    real_Double_t   gflops, magma_perf, magma_time=0, cublas_perf=0, cublas_time=0, cpu_perf=0, cpu_time=0;
    float          magma_error, cublas_error, normx, normr, normA, work[1];
    magma_int_t i, j, s, N, info;
    magma_int_t Ak;
    magma_int_t sizeA, sizeB;
    magma_int_t lda, ldda;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    magma_int_t *ipiv;

    float c_zero = MAGMA_S_ZERO;
    
    float *h_A, *h_b, *h_bcublas, *h_bmagma, *h_blapack, *h_x;
    float *d_A, *d_b;
    float **d_A_array = NULL;
    float **d_b_array = NULL;
    
    float **dwork_array = NULL;

    float c_neg_one = MAGMA_S_NEG_ONE;
    float alpha = MAGMA_S_ONE;
    int status = 0;
    magma_int_t batchCount;

    magma_opts opts( MagmaOptsBatched );
    opts.parse_opts( argc, argv );
    batchCount = opts.batchcount;
    
    float tol = opts.tolerance * lapackf77_slamch("E");

    printf("%% uplo = %s, transA = %s, diag = %s \n",
           lapack_uplo_const(opts.uplo),
           lapack_trans_const(opts.transA), lapack_diag_const(opts.diag) );
    printf("%% BatchCount   N   MAGMA Gflop/s (ms)   CUBLAS Gflop/s (ms)    CPU Gflop/s (ms)      MAGMA   CUBLAS error\n");
    printf("%%========================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.msize[itest];

            gflops = FLOPS_STRSM(opts.side, N, 1) / 1e9 * batchCount;

            lda = Ak = N;
            
            ldda = magma_roundup( lda, opts.align );  // multiple of 32 by default

            sizeA = lda*Ak*batchCount;
            sizeB = N*batchCount;

            TESTING_CHECK( magma_smalloc_cpu( &h_A,       sizeA  ));
            TESTING_CHECK( magma_smalloc_cpu( &h_b,       sizeB   ));
            TESTING_CHECK( magma_smalloc_cpu( &h_x,       sizeB   ));
            TESTING_CHECK( magma_smalloc_cpu( &h_blapack, sizeB   ));
            TESTING_CHECK( magma_smalloc_cpu( &h_bcublas, sizeB   ));
            TESTING_CHECK( magma_smalloc_cpu( &h_bmagma,  sizeB   ));
            TESTING_CHECK( magma_imalloc_cpu( &ipiv,      Ak      ));
            
            TESTING_CHECK( magma_smalloc( &d_A,       ldda*Ak*batchCount ));
            TESTING_CHECK( magma_smalloc( &d_b,       N*batchCount  ));
            
            TESTING_CHECK( magma_malloc( (void**) &d_A_array,   batchCount * sizeof(float*) ));
            TESTING_CHECK( magma_malloc( (void**) &d_b_array,   batchCount * sizeof(float*) ));
            TESTING_CHECK( magma_malloc( (void**) &dwork_array, batchCount * sizeof(float*) ));


            magmaFloat_ptr dwork=NULL; // invA and work are workspace in strsm
            magma_int_t dwork_batchSize = N;
            TESTING_CHECK( magma_smalloc( &dwork, dwork_batchSize * batchCount ));
    
            magma_sset_pointer( dwork_array, dwork, N, 0, 0, dwork_batchSize, batchCount, opts.queue );

            memset( h_bmagma, 0, batchCount*N*sizeof(float) );
            magmablas_slaset( MagmaFull, N, batchCount, c_zero, c_zero, dwork, N, opts.queue );

            /* Initialize the matrices */
            /* Factor A into LU to get well-conditioned triangular matrix.
             * Copy L to U, since L seems okay when used with non-unit diagonal
             * (i.e., from U), while U fails when used with unit diagonal. */
            lapackf77_slarnv( &ione, ISEED, &sizeA, h_A );

            for (s=0; s < batchCount; s++) {
                lapackf77_sgetrf( &Ak, &Ak, h_A + s * lda * Ak, &lda, ipiv, &info );
                for( j = 0; j < Ak; ++j ) {
                    for( i = 0; i < j; ++i ) {
                        *h_A(i,j,s) = *h_A(j,i,s);
                    }
                }
            }

            lapackf77_slarnv( &ione, ISEED, &sizeB, h_b );
            memcpy( h_blapack, h_b, sizeB*sizeof(float) );
            
            /* =====================================================================
               Performs operation using MAGMABLAS
               =================================================================== */
            magma_ssetmatrix( Ak, Ak*batchCount, h_A, lda, d_A, ldda, opts.queue );
            magma_ssetmatrix( N,  batchCount, h_b, N, d_b, N, opts.queue );

            magma_sset_pointer( d_A_array, d_A, ldda, 0, 0, ldda*Ak, batchCount, opts.queue );
            magma_sset_pointer( d_b_array, d_b, N, 0, 0, N, batchCount, opts.queue );
            magma_sset_pointer( dwork_array, dwork, N, 0, 0, N, batchCount, opts.queue );

            magma_time = magma_sync_wtime( opts.queue );

            magmablas_strsv_work_batched(opts.uplo, opts.transA, opts.diag, 
                                    N, d_A_array, ldda,
                                    d_b_array, 1, dwork_array, batchCount, opts.queue); 

            magma_time = magma_sync_wtime( opts.queue ) - magma_time;
            magma_perf = gflops / magma_time;
            magma_sgetmatrix( N, batchCount, dwork, N, h_bmagma, N, opts.queue );


            /* =====================================================================
               Performs operation using CUBLAS
               =================================================================== */
            magma_ssetmatrix( N, batchCount, h_b, N, d_b, N, opts.queue );
            magma_sset_pointer( d_b_array, d_b, N, 0, 0, N, batchCount, opts.queue );

            // CUBLAS version <= 6.0 has float **            dA_array, no cast needed.
            // CUBLAS version    6.5 has float const**       dA_array, requiring cast.
            // Correctly, it should be   float const* const* dA_array, to avoid requiring cast.
            #if CUDA_VERSION >= 6050
                cublas_time = magma_sync_wtime( opts.queue );
                cublasStrsmBatched(
                    opts.handle, cublas_side_const(MagmaLeft), cublas_uplo_const(opts.uplo),
                    cublas_trans_const(opts.transA), cublas_diag_const(opts.diag),
                    int(N), 1, &alpha,
                    (const float**) d_A_array, int(ldda),
                    d_b_array, int(N), int(batchCount) );
                cublas_time = magma_sync_wtime( opts.queue ) - cublas_time;
                cublas_perf = gflops / cublas_time;
            #else
                MAGMA_UNUSED( alpha );
            #endif

            magma_sgetmatrix( N, batchCount, d_b, N, h_bcublas, N, opts.queue );
            
            /* =====================================================================
               Performs operation using CPU BLAS
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_wtime();
                #if !defined (BATCHED_DISABLE_PARCPU) && defined(_OPENMP)
                magma_int_t nthreads = magma_get_lapack_numthreads();
                magma_set_lapack_numthreads(1);
                magma_set_omp_numthreads(nthreads);
                #pragma omp parallel for schedule(dynamic)
                #endif
                for (s=0; s < batchCount; s++) {
                    blasf77_strsv(
                        lapack_uplo_const(opts.uplo),
                        lapack_trans_const(opts.transA), lapack_diag_const(opts.diag),
                        &N, 
                        h_A + s * lda * Ak, &lda,
                        h_blapack + s * N,  &ione);
                }
                #if !defined (BATCHED_DISABLE_PARCPU) && defined(_OPENMP)
                    magma_set_lapack_numthreads(nthreads);
                #endif
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            // ||b - Ax|| / (||A||*||x||)
            magma_error  = 0;
            cublas_error = 0;
            for (s=0; s < batchCount; s++) {
                // error for CUBLAS
                normA = lapackf77_slange( "F", &N, &N, h_A + s * lda * Ak, &lda, work );
                float err;

                #if CUDA_VERSION >= 6050
                normx = lapackf77_slange( "F", &N, &ione, h_bcublas + s * N, &ione, work );
                blasf77_strmv( lapack_uplo_const(opts.uplo), lapack_trans_const(opts.transA), lapack_diag_const(opts.diag),
                               &N,
                               h_A + s * lda * Ak, &lda,
                               h_bcublas + s * N, &ione );

                blasf77_saxpy( &N, &c_neg_one, h_b + s * N, &ione, h_bcublas + s * N, &ione );
                normr = lapackf77_slange( "F", &N, &ione, h_bcublas + s * N, &N, work );
                err = normr / (normA*normx);
                
                if ( isnan(err) || isinf(err) ) {
                    printf("error for matrix %lld cublas_error = %7.2f where normr=%7.2f normx=%7.2f and normA=%7.2f\n", 
                            (long long) s, err, normr, normx, normA);
                    cublas_error = err;
                    break;
                }
                cublas_error = max( err, cublas_error );
                #endif

                // error for MAGMA
                normx = lapackf77_slange( "F", &N, &ione, h_bmagma + s * N, &ione, work );
                blasf77_strmv( lapack_uplo_const(opts.uplo), lapack_trans_const(opts.transA), lapack_diag_const(opts.diag),
                               &N,
                               h_A + s * lda * Ak, &lda,
                               h_bmagma  + s * N, &ione );

                blasf77_saxpy( &N, &c_neg_one, h_b + s * N, &ione, h_bmagma + s * N, &ione );
                normr = lapackf77_slange( "F", &N, &ione, h_bmagma + s * N, &N, work );
                err = normr / (normA*normx);

                if ( isnan(err) || isinf(err) ) {
                    printf("error for matrix %lld magma_error = %7.2f where normr=%7.2f normx=%7.2f and normA=%7.2f\n", 
                            (long long) s, err, normr, normx, normA);
                    magma_error = err;
                    break;
                }
                magma_error = max( err, magma_error );
            }
            bool okay = (magma_error < tol && cublas_error < tol);
            status += ! okay;

            if ( opts.lapack ) {
                printf("%10lld %5lld    %7.2f (%7.2f)     %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %8.2e   %s\n",
                        (long long) batchCount, (long long) N,
                        magma_perf,  1000.*magma_time,
                        cublas_perf, 1000.*cublas_time,
                        cpu_perf,    1000.*cpu_time,
                        magma_error, cublas_error,
                        (okay ? "ok" : "failed"));
            }
            else {
                printf("%10lld %5lld    %7.2f (%7.2f)     %7.2f (%7.2f)     ---   (  ---  )   %8.2e   %8.2e   %s\n",
                        (long long) batchCount, (long long) N,
                        magma_perf,  1000.*magma_time,
                        cublas_perf, 1000.*cublas_time,
                        magma_error, cublas_error,
                        (okay ? "ok" : "failed"));
            }
            
            magma_free_cpu( h_A );
            magma_free_cpu( h_b );
            magma_free_cpu( h_x );
            magma_free_cpu( h_blapack );
            magma_free_cpu( h_bcublas );
            magma_free_cpu( h_bmagma  );
            magma_free_cpu( ipiv );
            
            magma_free( d_A );
            magma_free( d_b );
            magma_free( d_A_array );
            magma_free( d_b_array );

            magma_free( dwork );
            magma_free( dwork_array );
            
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    opts.cleanup();
    TESTING_CHECK( magma_finalize() );
    return status;
}
