/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @generated from testing/testing_zgeqp3_gpu.cpp, normal z -> d, Thu Jul 27 17:28:39 2017

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

#define REAL

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dgeqp3_gpu
*/
int main( int argc, char** argv)
{
    TESTING_CHECK( magma_init() );
    magma_print_environment();
    
    real_Double_t    gflops, gpu_perf, gpu_time, cpu_perf=0, cpu_time=0;
    double *h_A, *h_R, *tau, *h_work;
    magmaDouble_ptr d_A, dtau, d_work;
    magma_int_t *jpvt;
    magma_int_t M, N, K, n2, lda, lwork, j, info, min_mn, nb;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    int status = 0;
    
    magma_opts opts;
    opts.parse_opts( argc, argv );

    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    printf("%% M     N     CPU Gflop/s (sec)   GPU Gflop/s (sec)   ||A*P - Q*R||_F\n");
    printf("%%====================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            K = opts.ksize[itest];
            if ( M < K || N < K || K <= 0 ) {
                printf( "%5lld %5lld %5lld   skipping because dgeqp3 requires M >= K, N >= K, K(the rank) >= 0\n",
                        (long long) M, (long long) N, (long long) K );
                continue;
            }
 
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N;
            nb     = magma_get_dgeqp3_nb( M, N );
            gflops = FLOPS_DGEQRF( M, N ) / 1e9;
            
            lwork = ( N+1 )*nb;
            #ifdef REAL
            lwork += 2*N;
            #endif
            if ( opts.check )
                lwork = max( lwork, M*N + N );
            
            #ifdef COMPLEX
            double *rwork;
            magmaDouble_ptr drwork;
            TESTING_CHECK( magma_dmalloc( &drwork, 2*N + (N+1)*nb ));
            TESTING_CHECK( magma_dmalloc_cpu( &rwork,  2*N ));
            #endif
            TESTING_CHECK( magma_imalloc_cpu( &jpvt,   N      ));
            TESTING_CHECK( magma_dmalloc_cpu( &tau,    min_mn ));
            TESTING_CHECK( magma_dmalloc_cpu( &h_A,    n2     ));
            
            TESTING_CHECK( magma_dmalloc_pinned( &h_R,    n2     ));
            TESTING_CHECK( magma_dmalloc_pinned( &h_work, lwork  ));
            
            TESTING_CHECK( magma_dmalloc( &dtau,   min_mn ));
            TESTING_CHECK( magma_dmalloc( &d_A,    lda*N  ));
            TESTING_CHECK( magma_dmalloc( &d_work, lwork  ));
            
            /* Initialize the matrix */
            lapackf77_dlarnv( &ione, ISEED, &n2, h_R );

            /* Make h_A of rank K */
            double alpha = MAGMA_D_MAKE(  1., 0. );
            double beta  = MAGMA_D_MAKE(  0., 0. );
            blasf77_dgemm("N", "N", &M, &N, &K, &alpha, h_R, &lda, h_R, &lda,
                          &beta, h_A, &lda);

            lapackf77_dlacpy( MagmaFullStr, &M, &N, h_A, &lda, h_R, &lda );
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                for( j = 0; j < N; j++)
                    jpvt[j] = 0;
                
                cpu_time = magma_wtime();
                lapackf77_dgeqp3( &M, &N, h_R, &lda, jpvt, tau, h_work, &lwork,
                                  #ifdef COMPLEX
                                  rwork,
                                  #endif
                                  &info );
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0) {
                    printf("lapack_dgeqp3 returned error %lld: %s.\n",
                           (long long) info, magma_strerror( info ));
                }
            }
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            lapackf77_dlacpy( MagmaFullStr, &M, &N, h_A, &lda, h_R, &lda );
            for( j = 0; j < N; j++)
                jpvt[j] = 0;
            
            /* copy A to gpu */
            magma_dsetmatrix( M, N, h_R, lda, d_A, lda, opts.queue );

            /* call gpu-interface */
            gpu_time = magma_wtime();
            magma_dgeqp3_gpu( M, N, d_A, lda, jpvt, dtau, d_work, lwork,
                              #ifdef COMPLEX
                              drwork,
                              #endif
                              &info );
            gpu_time = magma_wtime() - gpu_time;
            
            /* copy outputs to cpu */
            magma_dgetmatrix( M, N, d_A, lda, h_R, lda, opts.queue );
            magma_dgetvector( min_mn, dtau, 1, tau, 1, opts.queue );
            
            gpu_perf = gflops / gpu_time;
            if (info != 0) {
                printf("magma_dgeqp3 returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            if ( opts.lapack ) {
                printf("%5lld %5lld   %7.2f (%7.2f)   %7.2f (%7.2f)",
                       (long long) M, (long long) N, cpu_perf, cpu_time, gpu_perf, gpu_time );
            }
            else {
                printf("%5lld %5lld     ---   (  ---  )   %7.2f (%7.2f)",
                       (long long) M, (long long) N, gpu_perf, gpu_time );
            }
            if ( opts.check ) {
                double error, ulp;
                ulp = lapackf77_dlamch( "P" );
                
                // Compute norm( A*P - Q*R )
                error = lapackf77_dqpt01( &M, &N, &min_mn, h_A, h_R, &lda,
                                          tau, jpvt, h_work, &lwork );
                error *= ulp;
                printf("   %8.2e   %s\n", error, (error < tol ? "ok" : "failed"));
                status += ! (error < tol);
            }
            else {
                printf("     ---  \n");
            }
            
            #ifdef COMPLEX
            magma_free_cpu( rwork  );
            magma_free( drwork );
            #endif
            magma_free_cpu( jpvt   );
            magma_free_cpu( tau    );
            magma_free_cpu( h_A    );
            
            magma_free_pinned( h_R    );
            magma_free_pinned( h_work );
            
            magma_free( dtau   );
            magma_free( d_A    );
            magma_free( d_work );
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
