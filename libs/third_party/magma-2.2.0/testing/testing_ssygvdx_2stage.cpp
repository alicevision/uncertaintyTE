/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @author Raffaele Solca
       @author Azzam Haidar
       @author Mark Gates

       @generated from testing/testing_zhegvdx_2stage.cpp, normal z -> s, Thu Jul 27 17:29:20 2017

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

#include "../control/magma_threadsetting.h"  // internal header

#define REAL

/* ////////////////////////////////////////////////////////////////////////////
   -- Testing ssygvdx
*/
int main( int argc, char** argv)
{
    TESTING_CHECK( magma_init() );
    magma_print_environment();

    /* Constants */
    const float c_zero    = MAGMA_S_ZERO;
    const float c_one     = MAGMA_S_ONE;
    const float c_neg_one = MAGMA_S_NEG_ONE;
    const magma_int_t ione = 1;
    
    /* Local variables */
    real_Double_t gpu_time;

    float *h_A, *h_R, *h_B, *h_S, *h_work;

    #ifdef COMPLEX
    float *rwork;
    magma_int_t lrwork;
    #endif

    float *w1, *w2, result[2]={0,0};
    magma_int_t *iwork;
    magma_int_t N, n2, info, lda, lwork, liwork;
    magma_int_t ISEED[4] = {0,0,0,1};
    int status = 0;

    magma_opts opts;
    opts.parse_opts( argc, argv );
    
    float tol    = opts.tolerance * lapackf77_slamch("E");
    float tolulp = opts.tolerance * lapackf77_slamch("P");

    magma_range_t range = MagmaRangeAll;
    if (opts.fraction != 1)
        range = MagmaRangeI;

    // pass ngpu = -1 to test multi-GPU code using 1 gpu
    magma_int_t abs_ngpu = abs( opts.ngpu );
    
    printf("%% itype = %lld, jobz = %s, range = %s, uplo = %s, fraction = %6.4f, ngpu = %lld\n",
           (long long) opts.itype, lapack_vec_const(opts.jobz), lapack_range_const(range), lapack_uplo_const(opts.uplo),
           opts.fraction, (long long) abs_ngpu);

    if (opts.itype == 1) {
        printf("%%   N     M   GPU Time (sec)   |AZ-BZD|   |D - D_magma|\n");
    }                                                   
    else if (opts.itype == 2) {                      
        printf("%%   N     M   GPU Time (sec)   |ABZ-ZD|   |D - D_magma|\n");
    }                                                   
    else if (opts.itype == 3) {                      
        printf("%%   N     M   GPU Time (sec)   |BAZ-ZD|   |D - D_magma|\n");
    }                                     
        printf("%%======================================================\n");
    magma_int_t threads = magma_get_parallel_numthreads();
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            N = opts.nsize[itest];
            lda = N;
            n2  = lda*N;
            
            // TODO: test vl-vu range
            magma_int_t m1 = 0;
            float vl = 0;
            float vu = 0;
            magma_int_t il = 0;
            magma_int_t iu = 0;
            if (opts.fraction == 0) {
                il = max( 1, magma_int_t(0.1*N) );
                iu = max( 1, magma_int_t(0.3*N) );
            }
            else {
                il = 1;
                iu = max( 1, magma_int_t(opts.fraction*N) );
            }

            magma_ssyevdx_getworksize(N, threads, (opts.jobz == MagmaVec),
                                     &lwork,
                                     #ifdef COMPLEX
                                     &lrwork,
                                     #endif
                                     &liwork);
            /* Allocate host memory for the matrix */
            TESTING_CHECK( magma_smalloc_cpu( &h_A,    n2 ));
            TESTING_CHECK( magma_smalloc_cpu( &h_B,    n2 ));
            TESTING_CHECK( magma_smalloc_cpu( &w1,     N ));
            TESTING_CHECK( magma_smalloc_cpu( &w2,     N ));
            TESTING_CHECK( magma_imalloc_cpu( &iwork,  liwork ));
            
            TESTING_CHECK( magma_smalloc_pinned( &h_R,    n2 ));
            TESTING_CHECK( magma_smalloc_pinned( &h_S,    n2 ));
            TESTING_CHECK( magma_smalloc_pinned( &h_work, max( lwork, N*N ) ));  // check needs N*N
            #ifdef COMPLEX
            TESTING_CHECK( magma_smalloc_pinned( &rwork,  lrwork ));
            #endif

            /* Initialize the matrix */
            lapackf77_slarnv( &ione, ISEED, &n2, h_A );
            lapackf77_slarnv( &ione, ISEED, &n2, h_B );
            magma_smake_hpd( N, h_B, lda );
            magma_smake_symmetric( N, h_A, lda );

            lapackf77_slacpy( MagmaFullStr, &N, &N, h_A, &lda, h_R, &lda );
            lapackf77_slacpy( MagmaFullStr, &N, &N, h_B, &lda, h_S, &lda );

            // ===================================================================
            // Performs operation using MAGMA
            // ===================================================================
            gpu_time = magma_wtime();
            if (opts.ngpu == 1) {
                magma_ssygvdx_2stage( opts.itype, opts.jobz, range, opts.uplo,
                                      N, h_R, lda, h_S, lda, vl, vu, il, iu, &m1, w1,
                                      h_work, lwork,
                                      #ifdef COMPLEX
                                      rwork, lrwork,
                                      #endif
                                      iwork, liwork,
                                      &info );
            }
            else {
                magma_ssygvdx_2stage_m( abs_ngpu, opts.itype, opts.jobz, range, opts.uplo,
                                        N, h_R, lda, h_S, lda, vl, vu, il, iu, &m1, w1,
                                        h_work, lwork,
                                        #ifdef COMPLEX
                                        rwork, lrwork,
                                        #endif
                                        iwork, liwork,
                                        &info );
            }
            gpu_time = magma_wtime() - gpu_time;
            if (info != 0) {
                printf("magma_ssygvdx_2stage returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
            }
            
            if ( opts.check ) {
                /* =====================================================================
                   Check the results following the LAPACK's [zc]hegvdx routine.
                   A x = lambda B x is solved
                   and the following 3 tests computed:
                   (1)    | A Z - B Z D | / ( |A| |Z| N )  (itype = 1)
                          | A B Z - Z D | / ( |A| |Z| N )  (itype = 2)
                          | B A Z - Z D | / ( |A| |Z| N )  (itype = 3)
                   (2)    | D(with V, magma) - D(w/o V, lapack) | / | D |
                   =================================================================== */
                #ifdef REAL
                float *rwork = h_work + N*N;
                #endif
                
                if ( opts.jobz != MagmaNoVec ) {
                    result[0] = 1.;
                    result[0] /= safe_lapackf77_slansy("1", lapack_uplo_const(opts.uplo), &N, h_A, &lda, rwork);
                    result[0] /= lapackf77_slange("1", &N, &m1, h_R, &lda, rwork);
                    
                    if (opts.itype == 1) {
                        blasf77_ssymm("L", lapack_uplo_const(opts.uplo), &N, &m1, &c_one, h_A, &lda, h_R, &lda, &c_zero, h_work, &N);
                        for (int i=0; i < m1; ++i)
                            blasf77_sscal(&N, &w1[i], &h_R[i*N], &ione);
                        blasf77_ssymm("L", lapack_uplo_const(opts.uplo), &N, &m1, &c_neg_one, h_B, &lda, h_R, &lda, &c_one, h_work, &N);
                        result[0] *= lapackf77_slange("1", &N, &m1, h_work, &N, rwork)/N;
                    }
                    else if (opts.itype == 2) {
                        blasf77_ssymm("L", lapack_uplo_const(opts.uplo), &N, &m1, &c_one, h_B, &lda, h_R, &lda, &c_zero, h_work, &N);
                        for (int i=0; i < m1; ++i)
                            blasf77_sscal(&N, &w1[i], &h_R[i*N], &ione);
                        blasf77_ssymm("L", lapack_uplo_const(opts.uplo), &N, &m1, &c_one, h_A, &lda, h_work, &N, &c_neg_one, h_R, &lda);
                        result[0] *= lapackf77_slange("1", &N, &m1, h_R, &lda, rwork)/N;
                    }
                    else if (opts.itype == 3) {
                        blasf77_ssymm("L", lapack_uplo_const(opts.uplo), &N, &m1, &c_one, h_A, &lda, h_R, &lda, &c_zero, h_work, &N);
                        for (int i=0; i < m1; ++i)
                            blasf77_sscal(&N, &w1[i], &h_R[i*N], &ione);
                        blasf77_ssymm("L", lapack_uplo_const(opts.uplo), &N, &m1, &c_one, h_B, &lda, h_work, &N, &c_neg_one, h_R, &lda);
                        result[0] *= lapackf77_slange("1", &N, &m1, h_R, &lda, rwork)/N;
                    }
                }
                
                lapackf77_slacpy( MagmaFullStr, &N, &N, h_A, &lda, h_R, &lda );
                lapackf77_slacpy( MagmaFullStr, &N, &N, h_B, &lda, h_S, &lda );
                
                lapackf77_ssygvd( &opts.itype, "N", lapack_uplo_const(opts.uplo), &N,
                                  h_R, &lda, h_S, &lda, w2,
                                  h_work, &lwork,
                                  #ifdef COMPLEX
                                  rwork, &lrwork,
                                  #endif
                                  iwork, &liwork,
                                  &info );
                if (info != 0) {
                    printf("lapackf77_ssygvd returned error %lld: %s.\n",
                           (long long) info, magma_strerror( info ));
                }
                
                float maxw=0, diff=0;
                for (int j=0; j < m1; j++) {
                    maxw = max(maxw, fabs(w1[j]));
                    maxw = max(maxw, fabs(w2[j]));
                    diff = max(diff, fabs(w1[j] - w2[j]));
                }
                result[1] = diff / (m1*maxw);
            }
            
            /* =====================================================================
               Print execution time
               =================================================================== */
            printf("%5lld %5lld   %9.4f     ",
                   (long long) N, (long long) m1, gpu_time);
            if ( opts.check ) {
                bool okay = (result[1] < tolulp);
                if ( opts.jobz != MagmaNoVec ) {
                    okay = okay && (result[0] < tol);
                    printf("   %8.2e", result[0] );
                }
                else {
                    printf("     ---   ");
                }
                printf("        %8.2e  %s\n", result[1], (okay ? "ok" : "failed"));
                status += ! okay;
            }
            else {
                printf("     ---\n");
            }
            
            magma_free_cpu( h_A   );
            magma_free_cpu( h_B   );
            magma_free_cpu( w1    );
            magma_free_cpu( w2    );
            magma_free_cpu( iwork );
            
            magma_free_pinned( h_R );
            magma_free_pinned( h_S );
            magma_free_pinned( h_work );
            #ifdef COMPLEX
            magma_free_pinned( rwork );
            #endif
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
