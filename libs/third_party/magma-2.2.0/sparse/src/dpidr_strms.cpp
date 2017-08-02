/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @author Hartwig Anzt
       @author Eduardo Ponce
       @author Moritz Kreutzer

       @generated from sparse/src/zpidr_strms.cpp, normal z -> d, Sun Nov 20 20:20:46 2016
*/

#include "magmasparse_internal.h"

#define RTOLERANCE     lapackf77_dlamch( "E" )
#define ATOLERANCE     lapackf77_dlamch( "E" )


/**
    Purpose
    -------

    Solves a system of linear equations
       A * X = B
    where A is a real symmetric N-by-N positive definite matrix A.
    This is a GPU implementation of the preconditioned Induced Dimension
    Reduction method applying kernel fusion and kernel overlap.

    Arguments
    ---------

    @param[in]
    A           magma_d_matrix
                input matrix A

    @param[in]
    b           magma_d_matrix
                RHS b

    @param[in,out]
    x           magma_d_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_d_solver_par*
                solver parameters

    @param[in]
    precond_par magma_d_preconditioner*
                preconditioner

    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_dposv
    ********************************************************************/


extern "C" magma_int_t
magma_dpidr_strms(
    magma_d_matrix A, magma_d_matrix b, magma_d_matrix *x,
    magma_d_solver_par *solver_par,
    magma_d_preconditioner *precond_par,
    magma_queue_t queue )
{
    magma_int_t info = MAGMA_NOTCONVERGED;

    // prepare solver feedback
    solver_par->solver = Magma_PIDRMERGE;
    solver_par->numiter = 0;
    solver_par->spmv_count = 0;
    solver_par->init_res = 0.0;
    solver_par->final_res = 0.0;
    solver_par->iter_res = 0.0;
    solver_par->runtime = 0.0;

    // constants
    const double c_zero = MAGMA_D_ZERO;
    const double c_one = MAGMA_D_ONE;
    const double c_n_one = MAGMA_D_NEG_ONE;

    // internal user options
    const magma_int_t smoothing = 1;   // 0 = disable, 1 = enable
    const double angle = 0.7;          // [0-1]

    // local variables
    magma_int_t iseed[4] = {0, 0, 0, 1};
    magma_int_t dof;
    magma_int_t s;
    magma_int_t distr;
    magma_int_t k, i, sk;
    magma_int_t innerflag;
    magma_int_t ldd;
    magma_int_t q;
    double residual;
    double nrm;
    double nrmb;
    double nrmr;
    double nrmt;
    double rho;
    double om;
    double gamma;

    // matrices and vectors
    magma_d_matrix dxs = {Magma_CSR};
    magma_d_matrix dr = {Magma_CSR}, drs = {Magma_CSR};
    magma_d_matrix dP = {Magma_CSR}, dP1 = {Magma_CSR};
    magma_d_matrix dG = {Magma_CSR}, dGcol = {Magma_CSR};
    magma_d_matrix dU = {Magma_CSR};
    magma_d_matrix dM = {Magma_CSR};
    magma_d_matrix df = {Magma_CSR};
    magma_d_matrix dt = {Magma_CSR}, dtt = {Magma_CSR};
    magma_d_matrix dc = {Magma_CSR};
    magma_d_matrix dv = {Magma_CSR};
    magma_d_matrix dlu = {Magma_CSR};
    magma_d_matrix dskp = {Magma_CSR};
    magma_d_matrix dalpha = {Magma_CSR};
    magma_d_matrix dbeta = {Magma_CSR};
    double *hMdiag = NULL;
    double *hskp = NULL;
    double *halpha = NULL;
    double *hbeta = NULL;
    double *d1 = NULL, *d2 = NULL;
    
    // queue variables
    const magma_int_t nqueues = 3;     // number of queues
    magma_queue_t queues[nqueues];    

    // chronometry
    real_Double_t tempo1, tempo2;

    // create additional queues
    queues[0] = queue;
    for ( q = 1; q < nqueues; q++ ) {
        magma_queue_create( queue->device(), &(queues[q]) );
    }

    // initial s space
    // TODO: add option for 's' (shadow space number)
    // Hack: uses '--restart' option as the shadow space number.
    //       This is not a good idea because the default value of restart option is used to detect
    //       if the user provided a custom restart. This means that if the default restart value
    //       is changed then the code will think it was the user (unless the default value is
    //       also updated in the 'if' statement below.
    s = 1;
    if ( solver_par->restart != 50 ) {
        if ( solver_par->restart > A.num_cols ) {
            s = A.num_cols;
        } else {
            s = solver_par->restart;
        }
    }
    solver_par->restart = s;

    // set max iterations
    solver_par->maxiter = min( 2 * A.num_cols, solver_par->maxiter );

    // check if matrix A is square
    if ( A.num_rows != A.num_cols ) {
        //printf("Matrix A is not square.\n");
        info = MAGMA_ERR_NOT_SUPPORTED;
        goto cleanup;
    }

    // |b|
    nrmb = magma_dnrm2( b.num_rows, b.dval, 1, queue );
    if ( nrmb == 0.0 ) {
        magma_dscal( x->num_rows, MAGMA_D_ZERO, x->dval, 1, queue );
        info = MAGMA_SUCCESS;
        goto cleanup;
    }

    // t = 0
    // make t twice as large to contain both, dt and dr
    ldd = magma_roundup( b.num_rows, 32 );
    CHECK( magma_dvinit( &dt, Magma_DEV, ldd, 2, c_zero, queue ));
    dt.num_rows = b.num_rows;
    dt.num_cols = 1;
    dt.nnz = dt.num_rows;

    // redirect the dr.dval to the second part of dt
    CHECK( magma_dvinit( &dr, Magma_DEV, b.num_rows, 1, c_zero, queue ));
    magma_free( dr.dval );
    dr.dval = dt.dval + ldd;

    // r = b - A x
    CHECK( magma_dresidualvec( A, b, *x, &dr, &nrmr, queue ));
    
    // |r|
    solver_par->init_res = nrmr;
    solver_par->final_res = solver_par->init_res;
    solver_par->iter_res = solver_par->init_res;
    if ( solver_par->verbose > 0 ) {
        solver_par->res_vec[0] = (real_Double_t)nrmr;
    }

    // check if initial is guess good enough
    if ( nrmr <= solver_par->atol ||
        nrmr/nrmb <= solver_par->rtol ) {
        info = MAGMA_SUCCESS;
        goto cleanup;
    }

    // P = randn(n, s)
    // P = ortho(P)
//---------------------------------------
    // P = 0.0
    CHECK( magma_dvinit( &dP, Magma_CPU, A.num_cols, s, c_zero, queue ));

    // P = randn(n, s)
    distr = 3;        // 1 = unif (0,1), 2 = unif (-1,1), 3 = normal (0,1) 
    dof = dP.num_rows * dP.num_cols;
    lapackf77_dlarnv( &distr, iseed, &dof, dP.val );

    // transfer P to device
    CHECK( magma_dmtransfer( dP, &dP1, Magma_CPU, Magma_DEV, queue ));
    magma_dmfree( &dP, queue );

    // P = ortho(P1)
    if ( dP1.num_cols > 1 ) {
        // P = magma_dqr(P1), QR factorization
        CHECK( magma_dqr( dP1.num_rows, dP1.num_cols, dP1, dP1.ld, &dP, NULL, queue ));
    } else {
        // P = P1 / |P1|
        nrm = magma_dnrm2( dof, dP1.dval, 1, queue );
        nrm = 1.0 / nrm;
        magma_dscal( dof, nrm, dP1.dval, 1, queue );
        CHECK( magma_dmtransfer( dP1, &dP, Magma_DEV, Magma_DEV, queue ));
    }
    magma_dmfree( &dP1, queue );
//---------------------------------------

    // allocate memory for the scalar products
    CHECK( magma_dmalloc_pinned( &hskp, 5 ));
    CHECK( magma_dvinit( &dskp, Magma_DEV, 4, 1, c_zero, queue ));

    CHECK( magma_dmalloc_pinned( &halpha, s ));
    CHECK( magma_dvinit( &dalpha, Magma_DEV, s, 1, c_zero, queue ));

    CHECK( magma_dmalloc_pinned( &hbeta, s ));
    CHECK( magma_dvinit( &dbeta, Magma_DEV, s, 1, c_zero, queue ));
    
    // workspace for merged dot product
    CHECK( magma_dmalloc( &d1, max(2, s) * b.num_rows ));
    CHECK( magma_dmalloc( &d2, max(2, s) * b.num_rows ));

    // smoothing enabled
    if ( smoothing > 0 ) {
        // set smoothing solution vector
        CHECK( magma_dmtransfer( *x, &dxs, Magma_DEV, Magma_DEV, queue ));

        // tt = 0
        // make tt twice as large to contain both, dtt and drs
        ldd = magma_roundup( b.num_rows, 32 );
        CHECK( magma_dvinit( &dtt, Magma_DEV, ldd, 2, c_zero, queue ));
        dtt.num_rows = dr.num_rows;
        dtt.num_cols = 1;
        dtt.nnz = dtt.num_rows;

        // redirect the drs.dval to the second part of dtt
        CHECK( magma_dvinit( &drs, Magma_DEV, dr.num_rows, 1, c_zero, queue ));
        magma_free( drs.dval );
        drs.dval = dtt.dval + ldd;

        // set smoothing residual vector
        magma_dcopyvector( dr.num_rows, dr.dval, 1, drs.dval, 1, queue );
    }

    // G(n,s) = 0
    if ( s > 1 ) {
        ldd = magma_roundup( A.num_rows, 32 );
        CHECK( magma_dvinit( &dG, Magma_DEV, ldd, s, c_zero, queue ));
        dG.num_rows = A.num_rows;
    } else {
        CHECK( magma_dvinit( &dG, Magma_DEV, A.num_rows, s, c_zero, queue ));
    }

    // dGcol represents a single column of dG, array pointer is set inside loop
    CHECK( magma_dvinit( &dGcol, Magma_DEV, dG.num_rows, 1, c_zero, queue ));
    magma_free( dGcol.dval );

    // U(n,s) = 0
    if ( s > 1 ) {
        ldd = magma_roundup( A.num_cols, 32 );
        CHECK( magma_dvinit( &dU, Magma_DEV, ldd, s, c_zero, queue ));
        dU.num_rows = A.num_cols;
    } else {
        CHECK( magma_dvinit( &dU, Magma_DEV, A.num_cols, s, c_zero, queue ));
    }

    // M(s,s) = I
    CHECK( magma_dvinit( &dM, Magma_DEV, s, s, c_zero, queue ));
    CHECK( magma_dmalloc_pinned( &hMdiag, s ));
    magmablas_dlaset( MagmaFull, dM.num_rows, dM.num_cols, c_zero, c_one, dM.dval, dM.ld, queue );

    // f = 0
    CHECK( magma_dvinit( &df, Magma_DEV, dP.num_cols, 1, c_zero, queue ));

    // c = 0
    CHECK( magma_dvinit( &dc, Magma_DEV, dM.num_cols, 1, c_zero, queue ));

    // v = r
    CHECK( magma_dmtransfer( dr, &dv, Magma_DEV, Magma_DEV, queue ));

    // lu = 0
    CHECK( magma_dvinit( &dlu, Magma_DEV, dr.num_rows, 1, c_zero, queue ));

    //--------------START TIME---------------
    // chronometry
    tempo1 = magma_sync_wtime( queue );
    if ( solver_par->verbose > 0 ) {
        solver_par->timing[0] = 0.0;
    }

    om = MAGMA_D_ONE;
    gamma = MAGMA_D_ZERO;
    innerflag = 0;

    // start iteration
    do
    {
        solver_par->numiter++;

        // new RHS for small systems
        // f = P' r
        // Q1
        magma_dgemvmdot_shfl( dP.num_rows, dP.num_cols, dP.dval, dr.dval, d1, d2, df.dval, queues[1] );

        // skp[4] = f(k)
        // Q1
        magma_dgetvector_async( 1, df.dval, 1, &hskp[4], 1, queues[1] );

        // c(k:s) = f(k:s)
        // Q1
        magma_dcopyvector_async( s, df.dval, 1, dc.dval, 1, queues[1] );

        // c(k:s) = M(k:s,k:s) \ f(k:s)
        // Q1
        magma_dtrsv( MagmaLower, MagmaNoTrans, MagmaNonUnit, s, dM.dval, dM.ld, dc.dval, 1, queues[1] );

        // shadow space loop
        for ( k = 0; k < s; ++k ) {
            sk = s - k;
            dGcol.dval = dG.dval + k * dG.ld;

            // v = r - G(:,k:s) c(k:s)
            // Q1
            magmablas_dgemv( MagmaNoTrans, dG.num_rows, sk, c_n_one, dGcol.dval, dG.ld, &dc.dval[k], 1, c_one, dv.dval, 1, queues[1] );

            // preconditioning operation 
            // v = L \ v;
            // v = U \ v;
            // Q1
            CHECK( magma_d_applyprecond_left( MagmaNoTrans, A, dv, &dlu, precond_par, queues[1] )); 
            CHECK( magma_d_applyprecond_right( MagmaNoTrans, A, dlu, &dv, precond_par, queues[1] )); 

            // sync Q0 --> U(:,k) = U(:,k) - U(:,1:k) * alpha(1:k)
            magma_queue_sync( queues[0] );

            // U(:,k) = om * v + U(:,k:s) c(k:s)
            // Q1
            magmablas_dgemv( MagmaNoTrans, dU.num_rows, sk, c_one, &dU.dval[k*dU.ld], dU.ld, &dc.dval[k], 1, om, dv.dval, 1, queues[1] );

            // G(:,k) = A U(:,k)
            // Q1
            CHECK( magma_d_spmv( c_one, A, dv, c_zero, dGcol, queues[1] ));
            solver_par->spmv_count++;

            // bi-orthogonalize the new basis vectors
            for ( i = 0; i < k; ++i ) {
                // alpha = P(:,i)' G(:,k)
                // Q1
                halpha[i] = magma_ddot( dP.num_rows, &dP.dval[i*dP.ld], 1, dGcol.dval, 1, queues[1] );
                // implicit sync Q1 --> alpha = P(:,i)' G(:,k) 

                // alpha = alpha / M(i,i)
                halpha[i] = halpha[i] / hMdiag[i];
                    
                // G(:,k) = G(:,k) - alpha * G(:,i)
                // Q1
                magma_daxpy( dG.num_rows, -halpha[i], &dG.dval[i*dG.ld], 1, dGcol.dval, 1, queues[1] );
            }

            // sync Q1 --> compute new G, skp[4] = f(k
            magma_queue_sync( queues[1] );

            // new column of M = P'G, first k-1 entries are zero
            // M(k:s,k) = P(:,k:s)' G(:,k)
            // Q2
            magma_dgemvmdot_shfl( dP.num_rows, sk, &dP.dval[k*dP.ld], dGcol.dval, d1, d2, &dM.dval[k*dM.ld+k], queues[2] );

            // U(:,k) = v
            // Q0
            magma_dcopyvector_async( dU.num_rows, dv.dval, 1, &dU.dval[k*dU.ld], 1, queues[0] );

            // non-first s iteration
            if ( k > 0 ) {
                // alpha = dalpha
                // Q0
                magma_dsetvector_async( k, halpha, 1, dalpha.dval, 1, queues[0] );

                // U update outside of loop using GEMV
                // U(:,k) = U(:,k) - U(:,1:k) * alpha(1:k)
                // Q0
                magmablas_dgemv( MagmaNoTrans, dU.num_rows, k, c_n_one, dU.dval, dU.ld, dalpha.dval, 1, c_one, &dU.dval[k*dU.ld], 1, queues[0] );
            }

            // Mdiag(k) = M(k,k)
            // Q2
            magma_dgetvector( 1, &dM.dval[k*dM.ld+k], 1, &hMdiag[k], 1, queues[2] );
            // implicit sync Q2 --> Mdiag(k) = M(k,k)

            // check M(k,k) == 0
            if ( MAGMA_D_EQUAL(hMdiag[k], MAGMA_D_ZERO) ) {
                innerflag = 1;
                info = MAGMA_DIVERGENCE;
                break;
            }

            // beta = f(k) / M(k,k)
            hbeta[k] = hskp[4] / hMdiag[k];

            // check for nan
            if ( magma_d_isnan( hbeta[k] ) || magma_d_isinf( hbeta[k] )) {
                innerflag = 1;
                info = MAGMA_DIVERGENCE;
                break;
            }

            // non-last s iteration 
            if ( (k + 1) < s ) {
                // f(k+1:s) = f(k+1:s) - beta * M(k+1:s,k)
                // Q1
                magma_daxpy( sk-1, -hbeta[k], &dM.dval[k*dM.ld+(k+1)], 1, &df.dval[k+1], 1, queues[1] );

                // c(k+1:s) = f(k+1:s)
                // Q1
                magma_dcopyvector_async( sk-1, &df.dval[k+1], 1, &dc.dval[k+1], 1, queues[1] );

                // c(k+1:s) = M(k+1:s,k+1:s) \ f(k+1:s)
                // Q1
                magma_dtrsv( MagmaLower, MagmaNoTrans, MagmaNonUnit, sk-1, &dM.dval[(k+1)*dM.ld+(k+1)], dM.ld, &dc.dval[k+1], 1, queues[1] );

                // skp[4] = f(k+1)
                // Q1
                magma_dgetvector_async( 1, &df.dval[k+1], 1, &hskp[4], 1, queues[1] );
            }

            // r = r - beta * G(:,k)
            // Q2
            magma_daxpy( dr.num_rows, -hbeta[k], dGcol.dval, 1, dr.dval, 1, queues[2] );

            // smoothing disabled
            if ( smoothing <= 0 ) {
                // |r|
                // Q2
                nrmr = magma_dnrm2( dr.num_rows, dr.dval, 1, queues[2] );           
                // implicit sync Q2 --> |r|

                // v = r
                // Q1
                magma_dcopyvector_async( dr.num_rows, dr.dval, 1, dv.dval, 1, queues[1] );

            // smoothing enabled
            } else {
                // x = x + beta * U(:,k)
                // Q0
                magma_daxpy( x->num_rows, hbeta[k], &dU.dval[k*dU.ld], 1, x->dval, 1, queues[0] );

                // smoothing operation
//---------------------------------------
                // t = rs - r
                // Q2
                magma_didr_smoothing_1( drs.num_rows, drs.num_cols, drs.dval, dr.dval, dtt.dval, queues[2] );

                // t't
                // t'rs
                // Q2
                CHECK( magma_dgemvmdot_shfl( dt.ld, 2, dtt.dval, dtt.dval, d1, d2, &dskp.dval[2], queues[2] ));

                // skp[2-3] = dskp[2-3]
                // Q2
                magma_dgetvector( 2, &dskp.dval[2], 1, &hskp[2], 1, queues[2] );
                // implicit sync Q2 --> skp = dskp

                // gamma = (t' * rs) / (t' * t)
                gamma = hskp[3] / hskp[2];
                
                // xs = xs - gamma * (xs - x) 
                // Q0
                magma_didr_smoothing_2( dxs.num_rows, dxs.num_cols, -gamma, x->dval, dxs.dval, queues[0] );

                // v = r
                // Q1
                magma_dcopyvector_async( dr.num_rows, dr.dval, 1, dv.dval, 1, queues[1] );

                // rs = rs - gamma * t 
                // Q2
                magma_daxpy( drs.num_rows, -gamma, dtt.dval, 1, drs.dval, 1, queues[2] );

                // |rs|
                // Q2
                nrmr = magma_dnrm2( drs.num_rows, drs.dval, 1, queues[2] );       
                // implicit sync Q2 --> |r|
//---------------------------------------
            }

            // store current timing and residual
            if ( solver_par->verbose > 0 ) {
                tempo2 = magma_sync_wtime( queue );
                if ( (solver_par->numiter) % solver_par->verbose == 0 ) {
                    solver_par->res_vec[(solver_par->numiter) / solver_par->verbose]
                            = (real_Double_t)nrmr;
                    solver_par->timing[(solver_par->numiter) / solver_par->verbose]
                            = (real_Double_t)tempo2 - tempo1;
                }
            }

            // check convergence or iteration limit
            if ( nrmr <= solver_par->atol ||
                nrmr/nrmb <= solver_par->rtol ) { 
                s = k + 1; // for the x-update outside the loop
                innerflag = 2;
                info = MAGMA_SUCCESS;
                break;
            }
        }

        // smoothing disabled
        if ( smoothing <= 0 && innerflag != 1 ) {
            // dbeta(1:s) = beta(1:s)
            // Q0
            magma_dsetvector_async( s, hbeta, 1, dbeta.dval, 1, queues[0] );

            // x = x + U(:,1:s) * beta(1:s)
            // Q0
            magmablas_dgemv( MagmaNoTrans, dU.num_rows, s, c_one, dU.dval, dU.ld, dbeta.dval, 1, c_one, x->dval, 1, queues[0] );
        }

        // check convergence or iteration limit or invalid result of inner loop
        if ( innerflag > 0 ) {
            break;
        }

        // preconditioning operation 
        // v = L \ v;
        // v = U \ v;
        // Q2
        CHECK( magma_d_applyprecond_left( MagmaNoTrans, A, dv, &dlu, precond_par, queues[2] )); 
        CHECK( magma_d_applyprecond_right( MagmaNoTrans, A, dlu, &dv, precond_par, queues[2] )); 

        // t = A v
        // Q2
        CHECK( magma_d_spmv( c_one, A, dv, c_zero, dt, queues[2] ));
        solver_par->spmv_count++;

        // computation of a new omega
//---------------------------------------
        // t't
        // t'r
        // Q2
        CHECK( magma_dgemvmdot_shfl( dt.ld, 2, dt.dval, dt.dval, d1, d2, dskp.dval, queues[2] ));

        // skp[0-2] = dskp[0-2]
        // Q2
        magma_dgetvector( 2, dskp.dval, 1, hskp, 1, queues[2] );
        // implicit sync Q2 --> skp = dskp

        // |t|
        nrmt = magma_dsqrt( MAGMA_D_REAL(hskp[0]) );

        // rho = abs((t' * r) / (|t| * |r|))
        rho = MAGMA_D_ABS( MAGMA_D_REAL(hskp[1]) / (nrmt * nrmr) );

        // om = (t' * r) / (|t| * |t|)
        om = hskp[1] / hskp[0]; 
        if ( rho < angle ) {
            om = (om * angle) / rho;
        }
//---------------------------------------
        if ( MAGMA_D_EQUAL(om, MAGMA_D_ZERO) ) {
            info = MAGMA_DIVERGENCE;
            break;
        }

        // sync Q1 --> v = r
        magma_queue_sync( queues[1] );

        // r = r - om * t
        // Q2
        magma_daxpy( dr.num_rows, -om, dt.dval, 1, dr.dval, 1, queues[2] );

        // x = x + om * v
        // Q0
        magma_daxpy( x->num_rows, om, dv.dval, 1, x->dval, 1, queues[0] );

        // smoothing disabled
        if ( smoothing <= 0 ) {
            // |r|
            // Q2
            nrmr = magma_dnrm2( dr.num_rows, dr.dval, 1, queues[2] );           
            // implicit sync Q2 --> |r|

            // v = r
            // Q1
            magma_dcopyvector_async( dr.num_rows, dr.dval, 1, dv.dval, 1, queues[1] );

        // smoothing enabled
        } else {
            // smoothing operation
//---------------------------------------
            // t = rs - r
            // Q2
            magma_didr_smoothing_1( drs.num_rows, drs.num_cols, drs.dval, dr.dval, dtt.dval, queues[2] );

            // t't
            // t'rs
            // Q2
            CHECK( magma_dgemvmdot_shfl( dt.ld, 2, dtt.dval, dtt.dval, d1, d2, &dskp.dval[2], queues[2] ));

            // skp[2-3] = dskp[2-3]
            // Q2
            magma_dgetvector( 2, &dskp.dval[2], 1, &hskp[2], 1, queues[2] );
            // implicit sync Q2 --> skp = dskp

            // gamma = (t' * rs) / (t' * t)
            gamma = hskp[3] / hskp[2];

            // xs = xs - gamma * (xs - x) 
            // Q0
            magma_didr_smoothing_2( dxs.num_rows, dxs.num_cols, -gamma, x->dval, dxs.dval, queues[0] );

            // v = r
            // Q1
            magma_dcopyvector_async( dr.num_rows, dr.dval, 1, dv.dval, 1, queues[1] );

            // rs = rs - gamma * (rs - r) 
            // Q2
            magma_daxpy( drs.num_rows, -gamma, dtt.dval, 1, drs.dval, 1, queues[2] );

            // |rs|
            // Q2
            nrmr = magma_dnrm2( drs.num_rows, drs.dval, 1, queues[2] );           
            // implicit sync Q2 --> |r|
//---------------------------------------
        }

        // store current timing and residual
        if ( solver_par->verbose > 0 ) {
            tempo2 = magma_sync_wtime( queue );
            magma_queue_sync( queue );
            if ( (solver_par->numiter) % solver_par->verbose == 0 ) {
                solver_par->res_vec[(solver_par->numiter) / solver_par->verbose]
                        = (real_Double_t)nrmr;
                solver_par->timing[(solver_par->numiter) / solver_par->verbose]
                        = (real_Double_t)tempo2 - tempo1;
            }
        }

        // check convergence or iteration limit
        if ( nrmr <= solver_par->atol ||
            nrmr/nrmb <= solver_par->rtol ) { 
            info = MAGMA_SUCCESS;
            break;
        }
    }
    while ( solver_par->numiter + 1 <= solver_par->maxiter );

    // sync all queues
    for ( q = 0; q < nqueues; q++ ) {
        magma_queue_sync( queues[q] );
    }

    // smoothing enabled
    if ( smoothing > 0 ) {
        // x = xs
        magma_dcopyvector_async( x->num_rows, dxs.dval, 1, x->dval, 1, queue );

        // r = rs
        magma_dcopyvector_async( dr.num_rows, drs.dval, 1, dr.dval, 1, queue );
    }

    // get last iteration timing
    tempo2 = magma_sync_wtime( queue );
    magma_queue_sync( queue );
    solver_par->runtime = (real_Double_t)tempo2 - tempo1;
//--------------STOP TIME----------------

    // get final stats
    solver_par->iter_res = nrmr;
    CHECK( magma_dresidualvec( A, b, *x, &dr, &residual, queue ));
    solver_par->final_res = residual;

    // set solver conclusion
    if ( info != MAGMA_SUCCESS && info != MAGMA_DIVERGENCE ) {
        if ( solver_par->init_res > solver_par->final_res ) {
            info = MAGMA_SLOW_CONVERGENCE;
        }
    }


cleanup:
    // free resources
    // sync all queues, destory additional queues
    magma_queue_sync( queues[0] );
    for ( q = 1; q < nqueues; q++ ) {
        magma_queue_sync( queues[q] );
        magma_queue_destroy( queues[q] );
    }

    // smoothing enabled
    if ( smoothing > 0 ) {
        drs.dval = NULL;  // needed because its pointer is redirected to dtt
        magma_dmfree( &dxs, queue );
        magma_dmfree( &drs, queue ); 
        magma_dmfree( &dtt, queue );
    }
    dr.dval = NULL;       // needed because its pointer is redirected to dt
    dGcol.dval = NULL;    // needed because its pointer is redirected to dG
    magma_dmfree( &dr, queue );
    magma_dmfree( &dP, queue );
    magma_dmfree( &dP1, queue );
    magma_dmfree( &dG, queue );
    magma_dmfree( &dGcol, queue );
    magma_dmfree( &dU, queue );
    magma_dmfree( &dM, queue );
    magma_dmfree( &df, queue );
    magma_dmfree( &dt, queue );
    magma_dmfree( &dc, queue );
    magma_dmfree( &dv, queue );
    magma_dmfree( &dlu, queue );
    magma_dmfree( &dskp, queue );
    magma_dmfree( &dalpha, queue );
    magma_dmfree( &dbeta, queue );
    magma_free_pinned( hMdiag );
    magma_free_pinned( hskp );
    magma_free_pinned( halpha );
    magma_free_pinned( hbeta );
    magma_free( d1 );
    magma_free( d2 );

    solver_par->info = info;
    return info;
    /* magma_dpidr_strms */
}
