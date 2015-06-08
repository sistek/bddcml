/* signatures of BDDCML interface functions */
/* meaning of the variables is described in src/module_bddcml.f90 */
#ifndef bddcml_h
#define bddcml_h

void bddcml_init_c( int *nl, int *nsublev, int *lnsublev, int *nsub_loc_1, int *comm_init, int *verbose_level, int *numbase, 
                    int *just_direct_solve_int );

void bddcml_upload_global_data_c( int *nelem, int *nnod, int *ndof, int *ndim, int *meshdim, 
                                  int *inet, int *linet, int *nnet, int *lnnet, int *nndf, int *lnndf, 
                                  double *xyz, int *lxyz1, int *lxyz2,
                                  int *ifix, int *lifix, double *fixv, int *lfixv, double *rhs, int *lrhs, double *sol, int *lsol, int *idelm, 
                                  int *neighbouring, int *load_division_int );

void bddcml_upload_subdomain_data_c( int *nelem, int *nnod, int *ndof, int *ndim, int *meshdim,
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
                                     double *dof_data, int *ldof_data, int *find_components_int );

void bddcml_setup_preconditioner_c( int *matrixtype, int *use_defaults_int,
                                    int *parallel_division_int, 
                                    int *use_corner_constraints_int, 
                                    int *use_arithmetic_constraints_int, 
                                    int *use_adaptive_constraints_int,
                                    int *use_user_constraints_int,
                                    int *weights_type );

void bddcml_solve_c( int *comm_all, int *method, double *tol, int *maxit, int *ndecrmax, 
                     int *recycling_int, int *max_number_of_stored_vectors,
                     int *num_iter, int *converged_reason, double *condition_number);

void bddcml_download_local_solution_c( int *isub, double *sols, int *lsols );

void bddcml_download_global_solution_c( double *sol, int *lsol );

void bddcml_download_local_reactions_c( int *isub, double *reas, int *lreas );

void bddcml_download_global_reactions_c( double *rea, int *lrea );

void bddcml_change_global_data_c( int *ifix, int *lifix, double *fixv, int *lfixv, double *rhs, int *lrhs, double *sol, int *lsol);

void bddcml_change_subdomain_data_c( int *isub, 
                                     int *ifix, int *lifix, double *fixv, int *lfixv, 
                                     double *rhs, int *lrhs, int *is_rhs_complete, 
                                     double *sol, int *lsol );

void bddcml_setup_new_data_c( );

void bddcml_dotprod_subdomain_c( int *isub, double *vec1, int *lvec1, double *vec2, int *lvec2, double *dotprod );

void bddcml_finalize_c( );

#endif
