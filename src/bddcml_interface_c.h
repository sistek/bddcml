/* signatures of BDDCML interface functions */
/* meaning of the variables is described in src/bddcml_fortran_interface.f90 */
#ifndef bddcml_h
#define bddcml_h

// macro for calling Fortran routines from C - changing the name as Fortran compiler does it
#ifndef F_SYMBOL

# if defined(UPPER) 
#  define F_SYMBOL(lower_case,upper_case) upper_case
# elif defined(Add_)
#  define F_SYMBOL(lower_case,upper_case) lower_case##_
# elif defined(Add__)
#  define F_SYMBOL(lower_case,upper_case) lower_case##__
# else
#  define F_SYMBOL(lower_case,upper_case) lower_case
# endif

# endif

#define bddcml_init F_SYMBOL(bddcml_init, BDDCML_INIT)
void bddcml_init( int *nl, int *nsublev, int *lnsublev, int *comm_init );

#define bddcml_upload_global_data F_SYMBOL(bddcml_upload_global_data, BDDCML_UPLOAD_GLOBAL_DATA)
void bddcml_upload_global_data( int *nelem, int *nnod, int *ndof,
                                int *numbase, int *inet, int *linet, int *nnet, int *lnnet, int *nndf, int *lnndf, 
                                double *xyz, int *lxyz1, int *lxyz2,
                                int *ifix, int *lifix, double *fixv, int *lfixv, double *rhs, int *lrhs, double *sol, int *lsol, int *idelm );

#define bddcml_upload_local_data F_SYMBOL(bddcml_upload_local_data, BDDCML_UPLOAD_LOCAL_DATA)
void bddcml_upload_local_data( int *nelem, int *nnod, int *ndof, int *ndim, 
                               int *isub, int *nelems, int *nnods, int *ndofs, 
                               int *numbase, int *inet, int *linet, int *nnet, int *lnnet, int *nndf, int *lnndf, 
                               int *isngn, int *lisngn, int *isvgvn, int *lisvgvn, int *isegn, int *lisegn, 
                               double *xyz, int *lxyz1, int *lxyz2, 
                               int *ifix, int *lifix, double *fixv, int *lfixv, 
                               double *rhs, int *lrhs, 
                               int *matrixtype, int *i_sparse, int *j_sparse, double *a_sparse, int *la, int *is_assembled_int );

#define bddcml_setup_preconditioner F_SYMBOL(bddcml_setup_preconditioner, BDDCML_SETUP_PRECONDITIONER)
void bddcml_setup_preconditioner( int *matrixtype, int *ndim, int *meshdim, int *neighbouring, 
                                  int *use_defaults_int, int *load_division_int,
                                  int *parallel_division_int, int *correct_division_int, int *parallel_neighbouring_int,
                                  int *parallel_globs_int, int *use_arithmetic_int, int *use_adaptive_int );

#define bddcml_solve F_SYMBOL(bddcml_solve, BDDCML_SOLVE)
void bddcml_solve( int *comm_all, int *method, double *tol, int *maxit, int *ndecrmax, 
                   int *num_iter, int *converged_reason, double *condition_number);

#define bddcml_download_local_solution F_SYMBOL(bddcml_download_local_solution, BDDCML_DOWNLOAD_LOCAL_SOLUTION)
void bddcml_download_local_solution( double *sols, int *lsols, double *norm_sol );

#define bddcml_download_global_solution F_SYMBOL(bddcml_download_global_solution, BDDCML_DOWNLOAD_GLOBAL_SOLUTION)
void bddcml_download_global_solution( double *sol, int *lsol );

#define bddcml_finalize F_SYMBOL(bddcml_finalize, BDDCML_FINALIZE)
void bddcml_finalize( );

#endif
