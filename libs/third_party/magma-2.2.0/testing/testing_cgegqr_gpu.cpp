/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @generated from testing/testing_zgegqr_gpu.cpp, normal z -> c, Thu Jul 27 17:28:34 2017
       @author Stan Tomov

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_v2.h"
#include "magma_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing cgegqr
*/
int main( int argc, char** argv)
{
    TESTING_CHECK( magma_init() );
    magma_print_environment();

    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float           error, e1, e2, e3, e4, e5, *work;
    magmaFloatComplex c_neg_one = MAGMA_C_NEG_ONE;
    magmaFloatComplex c_one     = MAGMA_C_ONE;
    magmaFloatComplex c_zero    = MAGMA_C_ZERO;
    magmaFloatComplex *h_A, *h_R, *tau, *dtau, *h_work, *h_rwork, tmp[1], unused[1];

    magmaFloatComplex_ptr d_A, dwork;
    magma_int_t M, N, n2, lda, ldda, lwork, info, min_mn;
    magma_int_t ione     = 1, ldwork;
    magma_int_t ISEED[4] = {0,0,0,1};
    int status = 0;

    magma_opts opts;
    opts.parse_opts( argc, argv );
    opts.lapack |= opts.check;  // check (-c) implies lapack (-l)
    
    // versions 1...4 are valid
    if (opts.version < 1 || opts.version > 4) {
        printf("Unknown version %lld; exiting\n", (long long) opts.version );
        return -1;
    }
    
    float tol = 10. * opts.tolerance * lapackf77_slamch("E");
    
    printf("%% version %lld\n", (long long) opts.version );
    printf("%% M     N     CPU Gflop/s (ms)    GPU Gflop/s (ms)      ||I-Q'Q||_F / M     ||I-Q'Q||_I / M    ||A-Q R||_I\n");
    printf("%%                                                       MAGMA  /  LAPACK    MAGMA  /  LAPACK\n");
    printf("%%=========================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];

            if (N > 128) {
                printf("%5lld %5lld   skipping because cgegqr requires N <= 128\n",
                        (long long) M, (long long) N);
                continue;
            }
            if (M < N) {
                printf("%5lld %5lld   skipping because cgegqr requires M >= N\n",
                        (long long) M, (long long) N);
                continue;
            }

            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N;
            ldda   = magma_roundup( M, opts.align );  // multiple of 32 by default
            gflops = FLOPS_CGEQRF( M, N ) / 1e9 +  FLOPS_CUNGQR( M, N, N ) / 1e9;
            
            // query for workspace size
            lwork = -1;
            lapackf77_cgeqrf( &M, &N, unused, &M, unused, tmp, &lwork, &info );
            lwork = (magma_int_t)MAGMA_C_REAL( tmp[0] );
            lwork = max(lwork, 3*N*N);
            
            ldwork = N*N;
            if (opts.version == 2) {
                ldwork = 3*N*N + min_mn + 2;
            }

            TESTING_CHECK( magma_cmalloc_pinned( &tau,    min_mn ));
            TESTING_CHECK( magma_cmalloc_pinned( &h_work, lwork  ));
            TESTING_CHECK( magma_cmalloc_pinned( &h_rwork, lwork  ));

            TESTING_CHECK( magma_cmalloc_cpu( &h_A,   n2     ));
            TESTING_CHECK( magma_cmalloc_cpu( &h_R,   n2     ));
            TESTING_CHECK( magma_smalloc_cpu( &work,  M      ));
            
            TESTING_CHECK( magma_cmalloc( &d_A,   ldda*N ));
            TESTING_CHECK( magma_cmalloc( &dtau,  min_mn ));
            TESTING_CHECK( magma_cmalloc( &dwork, ldwork ));

            /* Initialize the matrix */
            lapackf77_clarnv( &ione, ISEED, &n2, h_A );

            lapackf77_clacpy( MagmaFullStr, &M, &N, h_A, &lda, h_R, &lda );
            magma_csetmatrix( M, N, h_R, lda, d_A, ldda, opts.queue );
            
            // warmup
            if ( opts.warmup ) {
                magma_cgegqr_gpu( 1, M, N, d_A, ldda, dwork, h_work, &info );
                magma_csetmatrix( M, N, h_R, lda, d_A, ldda, opts.queue );
            }
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            gpu_time = magma_sync_wtime( opts.queue );
            magma_cgegqr_gpu( opts.version, M, N, d_A, ldda, dwork, h_rwork, &info );
            gpu_time = magma_sync_wtime( opts.queue ) - gpu_time;
            gpu_perf = gflops / gpu_time;
            if (info != 0) {
                printf("magma_cgegqr returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
            }

            magma_cgetmatrix( M, N, d_A, ldda, h_R, lda, opts.queue );

            // Regenerate R
            // blasf77_cgemm("t", "n", &N, &N, &M, &c_one, h_R, &lda, h_A, &lda, &c_zero, h_rwork, &N);
            // magma_cprint(N, N, h_work, N);

            blasf77_ctrmm("r", "u", "n", "n", &M, &N, &c_one, h_rwork, &N, h_R, &lda);
            blasf77_caxpy( &n2, &c_neg_one, h_A, &ione, h_R, &ione );
            e5 = lapackf77_clange("i", &M, &N, h_R, &lda, work) /
                 lapackf77_clange("i", &M, &N, h_A, &lda, work);
            magma_cgetmatrix( M, N, d_A, ldda, h_R, lda, opts.queue );
 
            if ( opts.lapack ) {
                /* =====================================================================
                   Performs operation using LAPACK
                   =================================================================== */
                cpu_time = magma_wtime();

                /* Orthogonalize on the CPU */
                lapackf77_cgeqrf(&M, &N, h_A, &lda, tau, h_work, &lwork, &info);
                lapackf77_cungqr(&M, &N, &N, h_A, &lda, tau, h_work, &lwork, &info );

                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0) {
                    printf("lapackf77_cungqr returned error %lld: %s.\n",
                           (long long) info, magma_strerror( info ));
                }
                
                /* =====================================================================
                   Check the result compared to LAPACK
                   =================================================================== */
                blasf77_cgemm("c", "n", &N, &N, &M, &c_one, h_R, &lda, h_R, &lda, &c_zero, h_work, &N);
                for (magma_int_t ii = 0; ii < N*N; ii += N+1 ) {
                    h_work[ii] = MAGMA_C_SUB(h_work[ii], c_one);
                }
                e1 = lapackf77_clange("f", &N, &N, h_work, &N, work) / N;
                e3 = lapackf77_clange("i", &N, &N, h_work, &N, work) / N;

                blasf77_cgemm("c", "n", &N, &N, &M, &c_one, h_A, &lda, h_A, &lda, &c_zero, h_work, &N);
                for (magma_int_t ii = 0; ii < N*N; ii += N+1 ) {
                    h_work[ii] = MAGMA_C_SUB(h_work[ii], c_one);
                }
                e2 = lapackf77_clange("f", &N, &N, h_work, &N, work) / N;
                e4 = lapackf77_clange("i", &N, &N, h_work, &N, work) / N;

                if (opts.version != 4)
                    error = e1;
                else
                    error = e1 / (10.*max(M,N));

                printf("%5lld %5lld   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e / %8.2e   %8.2e / %8.2e   %8.2e  %s\n",
                       (long long) M, (long long) N, cpu_perf, 1000.*cpu_time, gpu_perf, 1000.*gpu_time,
                       e1, e2, e3, e4, e5,
                       (error < tol ? "ok" : "failed"));
                status += ! (error < tol);
            }
            else {
                printf("%5lld %5lld     ---   (  ---  )   %7.2f (%7.2f)     ---  \n",
                       (long long) M, (long long) N, gpu_perf, 1000.*gpu_time );
            }
            
            magma_free_pinned( tau    );
            magma_free_pinned( h_work );
            magma_free_pinned( h_rwork );
           
            magma_free_cpu( h_A  );
            magma_free_cpu( h_R  );
            magma_free_cpu( work );

            magma_free( d_A   );
            magma_free( dtau  );
            magma_free( dwork );

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
