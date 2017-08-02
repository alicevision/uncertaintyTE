/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @generated from testing/testing_zlacpy.cpp, normal z -> s, Thu Jul 27 17:27:43 2017
       @author Mark Gates
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "magma_v2.h"
#include "magma_lapack.h"
#include "testings.h"

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing slacpy
*/
int main( int argc, char** argv)
{
    TESTING_CHECK( magma_init() );
    magma_print_environment();

    real_Double_t    gbytes, gpu_perf, gpu_time, cpu_perf, cpu_time;
    float           error, work[1];
    float  c_neg_one = MAGMA_S_NEG_ONE;
    float *h_A, *h_B, *h_R;
    magmaFloat_ptr d_A, d_B;
    magma_int_t M, N, size, lda, ldb, ldda, lddb;
    magma_int_t ione     = 1;
    int status = 0;
    
    magma_opts opts;
    opts.parse_opts( argc, argv );

    magma_uplo_t uplo[] = { MagmaLower, MagmaUpper, MagmaFull };
    
    printf("%% uplo    M     N   CPU GByte/s (ms)    GPU GByte/s (ms)    check\n");
    printf("%%================================================================\n");
    for( int iuplo = 0; iuplo < 3; ++iuplo ) {
      for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            lda    = M;
            ldb    = lda;
            ldda   = magma_roundup( M, opts.align );  // multiple of 32 by default
            lddb   = ldda;
            size   = lda*N;
            if ( uplo[iuplo] == MagmaLower ) {
                // load & save lower trapezoid (with diagonal)
                if ( M > N ) {
                    gbytes = 2. * sizeof(float) * (1.*M*N - 0.5*N*(N-1)) / 1e9;
                } else {
                    gbytes = 2. * sizeof(float) * 0.5*M*(M+1) / 1e9;
                }
            }
            else if ( uplo[iuplo] == MagmaUpper ) {
                // load & save upper trapezoid (with diagonal)
                if ( N > M ) {
                    gbytes = 2. * sizeof(float) * (1.*M*N - 0.5*M*(M-1)) / 1e9;
                } else {
                    gbytes = 2. * sizeof(float) * 0.5*N*(N+1) / 1e9;
                }
            }
            else {
                // load & save entire matrix
                gbytes = 2. * sizeof(float) * 1.*M*N / 1e9;
            }
    
            TESTING_CHECK( magma_smalloc_cpu( &h_A, size   ));
            TESTING_CHECK( magma_smalloc_cpu( &h_B, size   ));
            TESTING_CHECK( magma_smalloc_cpu( &h_R, size   ));
            
            TESTING_CHECK( magma_smalloc( &d_A, ldda*N ));
            TESTING_CHECK( magma_smalloc( &d_B, lddb*N ));
            
            /* Initialize the matrix */
            for( int j = 0; j < N; ++j ) {
                for( int i = 0; i < M; ++i ) {
                    h_A[i + j*lda] = MAGMA_S_MAKE( i + j/10000., j );
                    h_B[i + j*ldb] = MAGMA_S_MAKE( i - j/10000. + 10000., j );
                }
            }
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            magma_ssetmatrix( M, N, h_A, lda, d_A, ldda, opts.queue );
            magma_ssetmatrix( M, N, h_B, ldb, d_B, lddb, opts.queue );
            
            gpu_time = magma_sync_wtime( opts.queue );
            //magmablas_slacpy( uplo[iuplo], M-2, N-2, d_A+1+ldda, ldda, d_B+1+lddb, lddb, opts.queue );  // inset by 1 row & col
            magmablas_slacpy( uplo[iuplo], M, N, d_A, ldda, d_B, lddb, opts.queue );
            gpu_time = magma_sync_wtime( opts.queue ) - gpu_time;
            gpu_perf = gbytes / gpu_time;
            
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            cpu_time = magma_wtime();
            //magma_int_t M2 = M-2;  // inset by 1 row & col
            //magma_int_t N2 = N-2;
            //lapackf77_slacpy( lapack_uplo_const(uplo[iuplo]), &M2, &N2, h_A+1+lda, &lda, h_B+1+ldb, &ldb );
            lapackf77_slacpy( lapack_uplo_const(uplo[iuplo]), &M, &N, h_A, &lda, h_B, &ldb );
            cpu_time = magma_wtime() - cpu_time;
            cpu_perf = gbytes / cpu_time;
            
            if ( opts.verbose ) {
                printf( "A= " );  magma_sprint(     M, N, h_A, lda );
                printf( "B= " );  magma_sprint(     M, N, h_B, ldb );
                printf( "dA=" );  magma_sprint_gpu( M, N, d_A, ldda, opts.queue );
                printf( "dB=" );  magma_sprint_gpu( M, N, d_B, lddb, opts.queue );
            }
            
            /* =====================================================================
               Check the result
               =================================================================== */
            magma_sgetmatrix( M, N, d_B, lddb, h_R, lda, opts.queue );
            
            blasf77_saxpy(&size, &c_neg_one, h_B, &ione, h_R, &ione);
            error = lapackf77_slange("f", &M, &N, h_R, &lda, work);

            printf("%5s %5lld %5lld   %7.2f (%7.2f)   %7.2f (%7.2f)   %s\n",
                   lapack_uplo_const(uplo[iuplo]), (long long) M, (long long) N,
                   cpu_perf, cpu_time*1000., gpu_perf, gpu_time*1000.,
                   (error == 0. ? "ok" : "failed") );
            status += ! (error == 0.);
            
            magma_free_cpu( h_A );
            magma_free_cpu( h_B );
            magma_free_cpu( h_R );
            
            magma_free( d_A );
            magma_free( d_B );
            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
      }
      printf( "\n" );
    }

    opts.cleanup();
    TESTING_CHECK( magma_finalize() );
    return status;
}
