/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @generated from sparse/src/zbicgstab.cpp, normal z -> d, Sun Nov 20 20:20:44 2016
       @author Hartwig Anzt

*/
#include "magmasparse_internal.h"

#define RTOLERANCE     lapackf77_dlamch( "E" )
#define ATOLERANCE     lapackf77_dlamch( "E" )


/**
    Purpose
    -------

    Solves a system of linear equations
       A * X = B
    where A is a general matrix.
    This is a GPU implementation of the Biconjugate Gradient Stabilized method.

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
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_dgesv
    ********************************************************************/

extern "C" magma_int_t
magma_dbicgstab(
    magma_d_matrix A, magma_d_matrix b, magma_d_matrix *x,
    magma_d_solver_par *solver_par,
    magma_queue_t queue )
{
    magma_int_t info = MAGMA_NOTCONVERGED;
    
    // prepare solver feedback
    solver_par->solver = Magma_BICGSTAB;
    solver_par->numiter = 0;
    solver_par->spmv_count = 0;

    // some useful variables
    double c_zero = MAGMA_D_ZERO;
    double c_one  = MAGMA_D_ONE;
    double c_neg_one = MAGMA_D_NEG_ONE;
    
    magma_int_t dofs = A.num_rows * b.num_cols;

    // workspace
    magma_d_matrix r={Magma_CSR}, rr={Magma_CSR}, p={Magma_CSR}, v={Magma_CSR}, s={Magma_CSR}, t={Magma_CSR};
    CHECK( magma_dvinit( &r, Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_dvinit( &rr,Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_dvinit( &p, Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_dvinit( &v, Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_dvinit( &s, Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_dvinit( &t, Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));

    
    // solver variables
    double alpha, beta, omega, rho_old, rho_new;
    double betanom, nom0, r0, res, nomb;
    //double den;

    // solver setup
    CHECK(  magma_dresidualvec( A, b, *x, &r, &nom0, queue));
    magma_dcopy( dofs, r.dval, 1, rr.dval, 1, queue );                  // rr = r
    betanom = nom0;
    rho_new = magma_ddot( dofs, r.dval, 1, r.dval, 1, queue );             // rho=<rr,r>
    rho_old = omega = alpha = MAGMA_D_MAKE( 1.0, 0. );
    solver_par->init_res = nom0;

    CHECK( magma_d_spmv( c_one, A, r, c_zero, v, queue ));              // z = A r
    //den = MAGMA_D_REAL( magma_ddot( dofs, v.dval, 1, r.dval, 1), queue ); // den = z' * r

    nomb = magma_dnrm2( dofs, b.dval, 1, queue );
    if ( nomb == 0.0 ){
        nomb=1.0;
    }       
    if ( (r0 = nomb * solver_par->rtol) < ATOLERANCE ){
        r0 = ATOLERANCE;
    }
    
    solver_par->final_res = solver_par->init_res;
    solver_par->iter_res = solver_par->init_res;
    if ( solver_par->verbose > 0 ) {
        solver_par->res_vec[0] = nom0;
        solver_par->timing[0] = 0.0;
    }
    if ( nomb < r0 ) {
        info = MAGMA_SUCCESS;
        goto cleanup;
    }

    //Chronometry
    real_Double_t tempo1, tempo2;
    tempo1 = magma_sync_wtime( queue );


    solver_par->numiter = 0;
    solver_par->spmv_count = 0;
    // start iteration
    do
    {
        solver_par->numiter++;

        rho_new = magma_ddot( dofs, rr.dval, 1, r.dval, 1, queue );  // rho=<rr,r>
        beta = rho_new/rho_old * alpha/omega;   // beta=rho/rho_old *alpha/omega
        magma_dscal( dofs, beta, p.dval, 1, queue );                 // p = beta*p
        magma_daxpy( dofs, c_neg_one * omega * beta, v.dval, 1 , p.dval, 1, queue );
                                                        // p = p-omega*beta*v
        magma_daxpy( dofs, c_one, r.dval, 1, p.dval, 1, queue );      // p = p+r
        CHECK( magma_d_spmv( c_one, A, p, c_zero, v, queue ));      // v = Ap
        solver_par->spmv_count++;
        alpha = rho_new / magma_ddot( dofs, rr.dval, 1, v.dval, 1, queue );
        magma_dcopy( dofs, r.dval, 1 , s.dval, 1, queue );            // s=r
        magma_daxpy( dofs, c_neg_one * alpha, v.dval, 1 , s.dval, 1, queue ); // s=s-alpha*v

        CHECK( magma_d_spmv( c_one, A, s, c_zero, t, queue ));       // t=As
        solver_par->spmv_count++;
        omega = magma_ddot( dofs, t.dval, 1, s.dval, 1, queue )   // omega = <s,t>/<t,t>
                   / magma_ddot( dofs, t.dval, 1, t.dval, 1, queue );

        magma_daxpy( dofs, alpha, p.dval, 1 , x->dval, 1, queue );     // x=x+alpha*p
        magma_daxpy( dofs, omega, s.dval, 1 , x->dval, 1, queue );     // x=x+omega*s

        magma_dcopy( dofs, s.dval, 1 , r.dval, 1, queue );             // r=s
        magma_daxpy( dofs, c_neg_one * omega, t.dval, 1 , r.dval, 1, queue ); // r=r-omega*t
        res = betanom = magma_dnrm2( dofs, r.dval, 1, queue );

        rho_old = rho_new;                                    // rho_old=rho

        if ( solver_par->verbose > 0 ) {
            tempo2 = magma_sync_wtime( queue );
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) res;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }

        if ( res/nomb <= solver_par->rtol || res <= solver_par->atol ){
            break;
        }
    }
    while ( solver_par->numiter+1 <= solver_par->maxiter );
    
    tempo2 = magma_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    double residual;
    CHECK(  magma_dresidualvec( A, b, *x, &r, &residual, queue));
    solver_par->iter_res = res;
    solver_par->final_res = residual;

    if ( solver_par->numiter < solver_par->maxiter ) {
        info = MAGMA_SUCCESS;
    } else if ( solver_par->init_res > solver_par->final_res ) {
        if ( solver_par->verbose > 0 ) {
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) betanom;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }
        info = MAGMA_SLOW_CONVERGENCE;
        if( solver_par->iter_res < solver_par->rtol*solver_par->init_res ||
            solver_par->iter_res < solver_par->atol ) {
            info = MAGMA_SUCCESS;
        }
    }
    else {
        if ( solver_par->verbose > 0 ) {
            if ( (solver_par->numiter)%solver_par->verbose==0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) betanom;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }
        info = MAGMA_DIVERGENCE;
    }
    
cleanup:
    magma_dmfree(&r, queue );
    magma_dmfree(&rr, queue );
    magma_dmfree(&p, queue );
    magma_dmfree(&v, queue );
    magma_dmfree(&s, queue );
    magma_dmfree(&t, queue );

    solver_par->info = info;
    return info;
}   /* magma_dbicgstab */
