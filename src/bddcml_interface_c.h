/* signatures of BDDCML interface functions */
/* meaning of the variables is described in src/bddcml_fortran_interface.f90 */
#ifndef bddcml_h
#define bddcml_h

/*
 * Macro for producing Fortran symbols based on 
 * useful for calling Fortran/C interoperability
 */
#include "f_symbol.h"

#define bddcml_init F_SYMBOL(bddcml_init, BDDCML_INIT)
void bddcml_init( int *nl, int *nsublev, int *lnsublev, int *nsub_loc_1, int *comm_init, int *verbose_level, int *numbase, 
                  int *just_direct_solve_int );

#define bddcml_upload_global_data F_SYMBOL(bddcml_upload_global_data, BDDCML_UPLOAD_GLOBAL_DATA)
void bddcml_upload_global_data( int *nelem, int *nnod, int *ndof, int *ndim, int *meshdim, 
                                int *inet, int *linet, int *nnet, int *lnnet, int *nndf, int *lnndf, 
                                double *xyz, int *lxyz1, int *lxyz2,
                                int *ifix, int *lifix, double *fixv, int *lfixv, double *rhs, int *lrhs, double *sol, int *lsol, int *idelm, 
                                int *neighbouring, int *load_division_int );

#define bddcml_upload_subdomain_data F_SYMBOL(bddcml_upload_subdomain_data, BDDCML_UPLOAD_SUBDOMAIN_DATA)
void bddcml_upload_subdomain_data( int *nelem, int *nnod, int *ndof, int *ndim, int *meshdim,
                                   int *isub, int *nelems, int *nnods, int *ndofs, 
                                   int *inet, int *linet, int *nnet, int *lnnet, int *nndf, int *lnndf, 
                                   int *isngn, int *lisngn, int *isvgvn, int *lisvgvn, int *isegn, int *lisegn, 
                                   double *xyz, int *lxyz1, int *lxyz2, 
                                   int *ifix, int *lifix, double *fixv, int *lfixv, 
                                   double *rhs, int *lrhs, int *is_rhs_complete, 
                                   double *sol, int *lsol, 
                                   int *matrixtype, int *i_sparse, int *j_sparse, double *a_sparse, int *la, int *is_assembled_int,
                                   double *user_constraints, int *luser_constraints1, int *luser_constraints2,
                                   double *element_data, int *lelement_data1, int *lelement_data2,
                                   double *dof_data, int *ldof_data );

#define bddcml_setup_preconditioner F_SYMBOL(bddcml_setup_preconditioner, BDDCML_SETUP_PRECONDITIONER)
void bddcml_setup_preconditioner( int *matrixtype, int *use_defaults_int,
                                  int *parallel_division_int, 
                                  int *use_arithmetic_constraints_int, 
                                  int *use_adaptive_constraints_int,
                                  int *use_user_constraints_int,
                                  int *weights_type );

#define bddcml_solve F_SYMBOL(bddcml_solve, BDDCML_SOLVE)
void bddcml_solve( int *comm_all, int *method, double *tol, int *maxit, int *ndecrmax, 
                   int *recycling_int, int *max_number_of_stored_vectors,
                   int *num_iter, int *converged_reason, double *condition_number);

#define bddcml_download_local_solution F_SYMBOL(bddcml_download_local_solution, BDDCML_DOWNLOAD_LOCAL_SOLUTION)
void bddcml_download_local_solution( int *isub, double *sols, int *lsols );

#define bddcml_download_global_solution F_SYMBOL(bddcml_download_global_solution, BDDCML_DOWNLOAD_GLOBAL_SOLUTION)
void bddcml_download_global_solution( double *sol, int *lsol );

#define bddcml_download_local_reactions F_SYMBOL(bddcml_download_local_reactions, BDDCML_DOWNLOAD_LOCAL_REACTIONS)
void bddcml_download_local_reactions( int *isub, double *reas, int *lreas );

#define bddcml_download_global_reactions F_SYMBOL(bddcml_download_global_reactions, BDDCML_DOWNLOAD_GLOBAL_REACTIONS)
void bddcml_download_global_reactions( double *rea, int *lrea );

#define bddcml_change_global_data F_SYMBOL(bddcml_change_global_data, BDDCML_CHANGE_GLOBAL_DATA)
void bddcml_change_global_data( int *ifix, int *lifix, double *fixv, int *lfixv, double *rhs, int *lrhs, double *sol, int *lsol);

#define bddcml_change_subdomain_data F_SYMBOL(bddcml_change_subdomain_data, BDDCML_CHANGE_SUBDOMAIN_DATA)
void bddcml_change_subdomain_data( int *isub, 
                                   int *ifix, int *lifix, double *fixv, int *lfixv, 
                                   double *rhs, int *lrhs, int *is_rhs_complete, 
                                   double *sol, int *lsol );

#define bddcml_setup_new_data F_SYMBOL(bddcml_setup_new_data, BDDCML_SETUP_NEW_DATA)
void bddcml_setup_new_data( );

#define bddcml_dotprod_subdomain F_SYMBOL(bddcml_dotprod_subdomain, BDDCML_DOTPROD_SUBDOMAIN)
void bddcml_dotprod_subdomain( int *isub, double *vec1, int *lvec1, double *vec2, int *lvec2, double *dotprod );

#define bddcml_finalize F_SYMBOL(bddcml_finalize, BDDCML_FINALIZE)
void bddcml_finalize( );

#endif
