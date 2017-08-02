/*
    -- MAGMA (version 2.2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date November 2016

       @author Hartwig Anzt

       @generated from sparse/src/ztfqmr.cpp, normal z -> s, Sun Nov 20 20:20:44 2016
*/

#include "magmasparse_internal.h"

#define RTOLERANCE     lapackf77_slamch( "E" )
#define ATOLERANCE     lapackf77_slamch( "E" )


/**
    Purpose
    -------

    Solves a system of linear equations
       A * X = B
    where A is a real matrix A.
    This is a GPU implementation of the transpose-free Quasi-Minimal Residual 
    method (TFQMR).

    Arguments
    ---------

    @param[in]
    A           magma_s_matrix
                input matrix A

    @param[in]
    b           magma_s_matrix
                RHS b

    @param[in,out]
    x           magma_s_matrix*
                solution approximation

    @param[in,out]
    solver_par  magma_s_solver_par*
                solver parameters

    @param[in]
    queue       magma_queue_t
                Queue to execute in.

    @ingroup magmasparse_sgesv
    ********************************************************************/

extern "C" magma_int_t
magma_stfqmr(
    magma_s_matrix A, magma_s_matrix b, magma_s_matrix *x,
    magma_s_solver_par *solver_par,
    magma_queue_t queue )
{
    magma_int_t info = MAGMA_NOTCONVERGED;

    // prepare solver feedback
    solver_par->solver = Magma_TFQMR;
    solver_par->numiter = 0;
    solver_par->spmv_count = 0;
    
    // local variables
    float c_zero = MAGMA_S_ZERO, c_one = MAGMA_S_ONE;
    // solver variables
    float nom0, r0,  res, nomb; //, normx, normd, normr_act;
    float rho = c_one, rho_l = c_one, eta = c_zero , c = c_zero , 
                        theta = c_zero , tau = c_zero, alpha = c_one, beta = c_zero,
                        sigma = c_zero;
    
    magma_int_t dofs = A.num_rows* b.num_cols;
    
    // magma_int_t stag = 0;

    // GPU workspace
    magma_s_matrix r={Magma_CSR}, r_tld={Magma_CSR}, pu_m={Magma_CSR},
                    d={Magma_CSR}, w={Magma_CSR}, v={Magma_CSR},
                    u_mp1={Magma_CSR}, u_m={Magma_CSR}, Au={Magma_CSR}, 
                    Ad={Magma_CSR}, Au_new={Magma_CSR};
    CHECK( magma_svinit( &r, Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_svinit( &u_mp1,Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_svinit( &r_tld,Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_svinit( &u_m, Magma_DEV, A.num_rows, b.num_cols, c_one, queue ));
    CHECK( magma_svinit( &pu_m, Magma_DEV, A.num_rows, b.num_cols, c_one, queue ));
    CHECK( magma_svinit( &v, Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_svinit( &d, Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_svinit( &w, Magma_DEV, A.num_rows, b.num_cols, c_one, queue ));
    CHECK( magma_svinit( &Ad, Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_svinit( &Au_new, Magma_DEV, A.num_rows, b.num_cols, c_zero, queue ));
    CHECK( magma_svinit( &Au, Magma_DEV, A.num_rows, b.num_cols, c_one, queue ));
    
    // solver setup
    CHECK(  magma_sresidualvec( A, b, *x, &r, &nom0, queue));
    solver_par->init_res = nom0;
    magma_scopy( dofs, r.dval, 1, r_tld.dval, 1, queue );   
    magma_scopy( dofs, r.dval, 1, w.dval, 1, queue );   
    magma_scopy( dofs, r.dval, 1, u_m.dval, 1, queue );   
    magma_scopy( dofs, u_m.dval, 1, pu_m.dval, 1, queue );   
    CHECK( magma_s_spmv( c_one, A, pu_m, c_zero, v, queue ));   // v = A u
    magma_scopy( dofs, v.dval, 1, Au.dval, 1, queue );  
    nomb = magma_snrm2( dofs, b.dval, 1, queue );
    if ( nomb == 0.0 ){
        nomb=1.0;
    }       
    if ( (r0 = nomb * solver_par->rtol) < ATOLERANCE ){
        r0 = ATOLERANCE;
    }
    solver_par->final_res = solver_par->init_res;
    solver_par->iter_res = solver_par->init_res;
    if ( solver_par->verbose > 0 ) {
        solver_par->res_vec[0] = (real_Double_t)nom0;
        solver_par->timing[0] = 0.0;
    }
    if ( nom0 < r0 ) {
        info = MAGMA_SUCCESS;
        goto cleanup;
    }

    tau = magma_ssqrt( magma_sdot( dofs, r.dval, 1, r_tld.dval, 1, queue ));
    rho = magma_sdot( dofs, r.dval, 1, r_tld.dval, 1, queue );
    rho_l = rho;
    
    //Chronometry
    real_Double_t tempo1, tempo2;
    tempo1 = magma_sync_wtime( queue );
    
    solver_par->numiter = 0;
    solver_par->spmv_count = 0;
    // start iteration
    do
    {
        solver_par->numiter++;
        if( solver_par->numiter%2 == 1 ){
            alpha = rho / magma_sdot( dofs, v.dval, 1, r_tld.dval, 1, queue );
            magma_scopy( dofs, u_m.dval, 1, u_mp1.dval, 1, queue );   
            magma_saxpy( dofs,  -alpha, v.dval, 1, u_mp1.dval, 1, queue );     // u_mp1 = u_m - alpha*v;
        }
        magma_saxpy( dofs,  -alpha, Au.dval, 1, w.dval, 1, queue );     // w = w - alpha*Au;
        sigma = theta * theta / alpha * eta;    
        magma_sscal( dofs, sigma, d.dval, 1, queue );    
        magma_saxpy( dofs, c_one, pu_m.dval, 1, d.dval, 1, queue );     // d = pu_m + sigma*d;
        magma_sscal( dofs, sigma, Ad.dval, 1, queue );         
        magma_saxpy( dofs, c_one, Au.dval, 1, Ad.dval, 1, queue );     // Ad = Au + sigma*Ad;

        
        theta = magma_ssqrt( magma_sdot(dofs, w.dval, 1, w.dval, 1, queue) ) / tau;
        c = c_one / magma_ssqrt( c_one + theta*theta );
        tau = tau * theta *c;
        eta = c * c * alpha;
        // normd = magma_snrm2( dofs, d.dval, 1, queue );
        // normx = magma_snrm2( dofs, x->dval, 1, queue );
        // 
        // 
        // if ( MAGMA_S_ABS(eta)*normd < 1e-15*normx ){
        //     stag = stag + 1;
        // } else {
        //     stag = 0;
        // }

        magma_saxpy( dofs, eta, d.dval, 1, x->dval, 1, queue );     // x = x + eta * d
        magma_saxpy( dofs, -eta, Ad.dval, 1, r.dval, 1, queue );     // r = r - eta * Ad
        res = magma_snrm2( dofs, r.dval, 1, queue );
        // normr_act = res;
        
        if ( solver_par->verbose > 0 ) {
            tempo2 = magma_sync_wtime( queue );
            if ( (solver_par->numiter)%solver_par->verbose == 0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) res;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }

        if ( res/nomb <= solver_par->rtol || res <= solver_par->atol ){
            info = MAGMA_SUCCESS;
            break;
        }
        
        // if (normr_act < normrmin){      //update minimal norm quantities
        //     normrmin = normr_act;
        //     //xmin = x;
        // }
        // if ( stag >= maxstagsteps ){      // 3 iterates are the same
        //     break;
        // }
    
    
        if( solver_par->numiter%2 == 0 ){
            rho = magma_sdot( dofs, w.dval, 1, r_tld.dval, 1, queue );
            beta = rho / rho_l;
            rho_l = rho;
            magma_scopy( dofs, w.dval, 1, u_mp1.dval, 1, queue );  
            magma_saxpy( dofs, beta, u_m.dval, 1, u_mp1.dval, 1, queue );     // u_mp1 = w + beta*u_m;
        }
              
        magma_scopy( dofs, u_mp1.dval, 1, pu_m.dval, 1, queue );  
        CHECK( magma_s_spmv( c_one, A, pu_m, c_zero, Au_new, queue )); // Au_new = A pu_m
        solver_par->spmv_count++;
        if( solver_par->numiter%2 == 0 ){
            magma_sscal( dofs, beta*beta, v.dval, 1, queue );                    
            magma_saxpy( dofs, beta, Au.dval, 1, v.dval, 1, queue );              
            magma_saxpy( dofs, c_one, Au_new.dval, 1, v.dval, 1, queue );      // v = Au_new + beta*(Au+beta*v);
        }
        magma_scopy( dofs, Au_new.dval, 1, Au.dval, 1, queue );  
        magma_scopy( dofs, u_mp1.dval, 1, u_m.dval, 1, queue );  
    }
    while ( solver_par->numiter+1 <= solver_par->maxiter );
    
    tempo2 = magma_sync_wtime( queue );
    solver_par->runtime = (real_Double_t) tempo2-tempo1;
    float residual;
    CHECK(  magma_sresidualvec( A, b, *x, &r, &residual, queue));
    solver_par->iter_res = res;
    solver_par->final_res = residual;

    if ( solver_par->numiter < solver_par->maxiter ) {
        info = MAGMA_SUCCESS;
    } else if ( solver_par->init_res > solver_par->final_res ) {
        if ( solver_par->verbose > 0 ) {
            if ( (solver_par->numiter)%solver_par->verbose == 0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) res;
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
            if ( (solver_par->numiter)%solver_par->verbose == 0 ) {
                solver_par->res_vec[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) res;
                solver_par->timing[(solver_par->numiter)/solver_par->verbose]
                        = (real_Double_t) tempo2-tempo1;
            }
        }
        info = MAGMA_DIVERGENCE;
    }
    
cleanup:
    magma_smfree(&r, queue );
    magma_smfree(&r_tld, queue );
    magma_smfree(&d, queue );
    magma_smfree(&w, queue );
    magma_smfree(&v, queue );
    magma_smfree(&pu_m, queue );
    magma_smfree(&u_m, queue );
    magma_smfree(&u_mp1, queue );
    magma_smfree(&d, queue );
    magma_smfree(&Au, queue );
    magma_smfree(&Au_new, queue );
    magma_smfree(&Ad, queue );
    
    solver_par->info = info;
    return info;
}   /* magma_stfqmr */
