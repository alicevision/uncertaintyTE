/*
   -- MAGMA (version 2.2.0) --
   Univ. of Tennessee, Knoxville
   Univ. of California, Berkeley
   Univ. of Colorado, Denver
   @date November 2016

   @author Mark gates
   @author Azzam Haidar
   @author Tingxing Dong

   @generated from testing/testing_zposv_batched.cpp, normal z -> c, Thu Jul 27 17:29:59 2017
 */
// includes, system
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_v2.h"
#include "magma_lapack.h"
#include "testings.h"

#if defined(_OPENMP)
#include <omp.h>
#include "../control/magma_threadsetting.h"  // internal header
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cposv_batched
*/
int main(int argc, char **argv)
{
    TESTING_CHECK( magma_init() );
    magma_print_environment();

    real_Double_t   gflops, cpu_perf, cpu_time, gpu_perf, gpu_time;
    float          error, Rnorm, Anorm, Xnorm, *work;
    magmaFloatComplex c_one     = MAGMA_C_ONE;
    magmaFloatComplex c_neg_one = MAGMA_C_NEG_ONE;
    magmaFloatComplex *h_A, *h_B, *h_X;
    magmaFloatComplex_ptr d_A, d_B;
    magma_int_t *cpu_info;
    magma_int_t *dinfo_array;
    magma_int_t N, nrhs, lda, ldb, ldda, lddb, info, sizeA, sizeB;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    int status = 0;
    magma_int_t batchCount;

    magmaFloatComplex **dA_array = NULL;
    magmaFloatComplex **dB_array = NULL;

    magma_opts opts( MagmaOptsBatched );
    opts.parse_opts( argc, argv );
    
    float tol = opts.tolerance * lapackf77_slamch("E");
    magma_queue_t queue = opts.queue;

    nrhs = opts.nrhs;
    batchCount = opts.batchcount;

    printf("%% uplo = %s\n", lapack_uplo_const(opts.uplo) );
    printf("%% BatchCount   N  NRHS   CPU Gflop/s (sec)   GPU Gflop/s (sec)   ||B - AX|| / N*||A||*||X||\n");
    printf("%%==========================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda    = N;
            ldb    = lda;
            ldda   = magma_roundup( N, opts.align );  // multiple of 32 by default
            lddb   = ldda;
            gflops = ( FLOPS_CPOTRF( N) + FLOPS_CPOTRS( N, nrhs ) ) / 1e9 * batchCount;
            
            sizeA = lda*N*batchCount;
            sizeB = ldb*nrhs*batchCount;

            TESTING_CHECK( magma_cmalloc_cpu( &h_A, sizeA ));
            TESTING_CHECK( magma_cmalloc_cpu( &h_B, sizeB ));
            TESTING_CHECK( magma_cmalloc_cpu( &h_X, sizeB ));
            TESTING_CHECK( magma_smalloc_cpu( &work, N ));
            TESTING_CHECK( magma_imalloc_cpu( &cpu_info, batchCount ));

            TESTING_CHECK( magma_cmalloc( &d_A, ldda*N*batchCount    ));
            TESTING_CHECK( magma_cmalloc( &d_B, lddb*nrhs*batchCount ));
            TESTING_CHECK( magma_imalloc( &dinfo_array, batchCount ));

            TESTING_CHECK( magma_malloc( (void**) &dA_array, batchCount * sizeof(magmaFloatComplex*) ));
            TESTING_CHECK( magma_malloc( (void**) &dB_array, batchCount * sizeof(magmaFloatComplex*) ));

            /* Initialize the matrices */
            lapackf77_clarnv( &ione, ISEED, &sizeA, h_A );
            lapackf77_clarnv( &ione, ISEED, &sizeB, h_B );

            for (int i=0; i < batchCount; i++)
            {
                magma_cmake_hpd( N, h_A + i * lda * N, lda ); // need modification
            }

            magma_csetmatrix( N, N*batchCount,    h_A, lda, d_A, ldda, opts.queue );
            magma_csetmatrix( N, nrhs*batchCount, h_B, ldb, d_B, lddb, opts.queue );
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            magma_cset_pointer( dA_array, d_A, ldda, 0, 0, ldda*N, batchCount, queue );
            magma_cset_pointer( dB_array, d_B, lddb, 0, 0, lddb*nrhs, batchCount, queue );

            gpu_time = magma_sync_wtime( opts.queue );
            info = magma_cposv_batched(opts.uplo, N, nrhs, dA_array, ldda, dB_array, lddb, dinfo_array, batchCount, queue);
            gpu_time = magma_sync_wtime( opts.queue ) - gpu_time;
            gpu_perf = gflops / gpu_time;
            // check correctness of results throught "dinfo_magma" and correctness of argument throught "info"
            magma_getvector( batchCount, sizeof(magma_int_t), dinfo_array, 1, cpu_info, 1, opts.queue );
            for (int i=0; i < batchCount; i++)
            {
                if (cpu_info[i] != 0 ) {
                    printf("magma_cposv_batched matrix %lld returned internal error %lld\n",
                            (long long) i, (long long) cpu_info[i] );
                }
            }
            if (info != 0) {
                printf("magma_cposv_batched returned argument error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
            }
            
            //=====================================================================
            // Residual
            //=====================================================================
            magma_cgetmatrix( N, nrhs*batchCount, d_B, lddb, h_X, ldb, opts.queue );

            error = 0;
            for (magma_int_t s=0; s < batchCount; s++)
            {
                Anorm = lapackf77_clange("I", &N, &N,    h_A + s * lda * N, &lda, work);
                Xnorm = lapackf77_clange("I", &N, &nrhs, h_X + s * ldb * nrhs, &ldb, work);
            
                blasf77_cgemm( MagmaNoTransStr, MagmaNoTransStr, &N, &nrhs, &N,
                           &c_one,     h_A + s * lda * N, &lda,
                                       h_X + s * ldb * nrhs, &ldb,
                           &c_neg_one, h_B + s * ldb * nrhs, &ldb);
            
                Rnorm = lapackf77_clange("I", &N, &nrhs, h_B + s * ldb * nrhs, &ldb, work);
                float err = Rnorm/(N*Anorm*Xnorm);
                
                if ( isnan(err) || isinf(err) ) {
                    error = err;
                    break;
                }
                error = max( err, error );
            }
            bool okay = (error < tol);
            status += ! okay;

            /* ====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_wtime();
                // #define BATCHED_DISABLE_PARCPU
                #if !defined (BATCHED_DISABLE_PARCPU) && defined(_OPENMP)
                magma_int_t nthreads = magma_get_lapack_numthreads();
                magma_set_lapack_numthreads(1);
                magma_set_omp_numthreads(nthreads);
                #pragma omp parallel for schedule(dynamic)
                #endif
                for (magma_int_t s=0; s < batchCount; s++)
                {
                    magma_int_t locinfo;
                    lapackf77_cposv( lapack_uplo_const(opts.uplo), &N, &nrhs, h_A + s * lda * N, &lda, h_B + s * ldb * nrhs, &ldb, &locinfo );
                    if (locinfo != 0) {
                        printf("lapackf77_cposv matrix %lld returned error %lld: %s.\n", 
                               (long long) s, (long long) locinfo, magma_strerror( locinfo ));
                    }
                }
                #if !defined (BATCHED_DISABLE_PARCPU) && defined(_OPENMP)
                    magma_set_lapack_numthreads(nthreads);
                #endif
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                
                printf( "%10lld %5lld %5lld   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
                        (long long) batchCount, (long long) N, (long long) nrhs, cpu_perf, cpu_time, gpu_perf, gpu_time,
                        error, (okay ? "ok" : "failed"));
            }
            else {
                printf( "%10lld %5lld %5lld     ---   (  ---  )   %7.2f (%7.2f)   %8.2e   %s\n",
                        (long long) batchCount, (long long) N, (long long) nrhs, gpu_perf, gpu_time,
                        error, (okay ? "ok" : "failed"));
            }
            
            magma_free_cpu( h_A );
            magma_free_cpu( h_B );
            magma_free_cpu( h_X );
            magma_free_cpu( work );
            magma_free_cpu( cpu_info );
            
            magma_free( d_A );
            magma_free( d_B );

            magma_free( dinfo_array );

            magma_free( dA_array );
            magma_free( dB_array );

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
