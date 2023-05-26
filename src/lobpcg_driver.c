/* BDDCML - Multilevel BDDC
 *
 * This program is a free software.
 * You can redistribute it and/or modify it under the terms of 
 * the GNU Lesser General Public License 
 * as published by the Free Software Foundation, 
 * either version 3 of the license, 
 * or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details
 * <http://www.gnu.org/copyleft/lesser.html>.
 *_______________________________________________________________*/

/* driver for running lobpcg code of Andrew Knyazev et al. */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/* lobpcg includes */
#include "lobpcg.h"
#include "matmultivec.h"
#include "multivector.h"
#include "multi_vector.h"
#include "pcg_multi.h"

/* function prototypes */
int lobpcg_dsygv(int *itype, char *jobz, char *uplo, int *n, double *a, int *lda, double *b, int *ldb,
                 double *w, double *work, int *lwork, int *info);
                    
int lobpcg_dpotrf(char *uplo, int *n, double *a, int *lda, int *info);

void lobpcg_mvecmult_f(int *n, double *x, int *lx, double *y, int *ly, int *idmatrix);
/* end prototypes */

void mvecmult_c (void * data, void * x_p, void * y_p, int idmatrix)
{
   serial_Multi_Vector *x = (serial_Multi_Vector *) x_p;
   serial_Multi_Vector *y = (serial_Multi_Vector *) y_p;
   double  *x_data; 
   double  *y_data;
   double * src;
   double * dest;
   int * x_active_ind;
   int * y_active_ind;
   int i;
   int size;
   int num_active_vectors;
   int lx;
   int ly;

   assert (x->size == y->size && x->num_active_vectors == y->num_active_vectors);
   assert (x->size>1);
   
   x_data = x->data;
   y_data = y->data;
   size = x->size;
   lx = size;
   ly = size;
   num_active_vectors = x->num_active_vectors;
   x_active_ind = x->active_indices;
   y_active_ind = y->active_indices;
   
   for(i=0; i<num_active_vectors; i++)
   {
      src = x_data + x_active_ind[i]*size;
      dest = y_data + y_active_ind[i]*size;

      /* call Fortran function */
      lobpcg_mvecmult_f(&size,src,&lx,dest,&ly,&idmatrix);
   }
   
}

void mvecmultA (void * data, void * x_p, void * y_p)
   /* auxiliary routine to multiply by matrix A */
{
   mvecmult_c (data, x_p, y_p, 1);
}

void mvecmultB (void * data, void * x_p, void * y_p)
   /* auxiliary routine to multiply by matrix B */
{
   mvecmult_c (data, x_p, y_p, 2);
}

void mvecmultM (void * data, void * x_p, void * y_p)
   /* auxiliary routine to multiply by matrix B */
{
   mvecmult_c (data, x_p, y_p, 5);
}

void mv_extract_values (serial_Multi_Vector * x_p, double * eigvec, int nvec)
{
   serial_Multi_Vector *x = (serial_Multi_Vector *) x_p;
   double  *x_data; 
   double * src;
   double * dest;
   int * x_active_ind;
   int i;
   int j;
   int size;
   int num_active_vectors;
   int idest;

   assert (x->size>1);
   
   x_data = x->data;
   size = x->size;
   num_active_vectors = x->num_active_vectors;
   x_active_ind = x->active_indices;

   assert(num_active_vectors == nvec);
   
   idest = 0;
   for(i=0; i<num_active_vectors; i++)
   {
      src = x_data + x_active_ind[i]*size;
      dest = eigvec + idest*size;

      /* copy particular eigenvector */
      for(j=0; j<size; j++)
      {
	 dest[j] = src[j];
      }
      idest = idest + 1;
   }
   
}

void mv_load_values (double * eigvec, int nvec, int size, serial_Multi_Vector * x_p)
{
   serial_Multi_Vector *x = (serial_Multi_Vector *) x_p;
   double  *x_data; 
   double * src;
   double * dest;
   int * x_active_ind;
   int i;
   int j;
   int idest;

   /* Perform checks */
   assert(x->size == size);
   assert(x->num_active_vectors == nvec);


   x_data = x->data;
   x_active_ind = x->active_indices;

   idest = 0;
   for(i=0; i<nvec; i++)
   {
      src  = eigvec + idest*size;
      dest = x_data + x_active_ind[i]*size;

      /* copy particular eigenvector */
      for(j=0; j<size; j++)
      {
	 dest[j] = src[j];
      }
      idest = idest + 1;
   }
   
}

/* Fortran callable subroutine arguments are passed by reference */
extern void lobpcg_driver(int *N, int *NVEC, double *TOL, int *MAXIT, int *VERBOSITY_LEVEL, int *USE_X_VALUES, int *LOBPCG_PRECONDITIONER, double *lambda, double *vec, int *ITERATIONS, int *IERR) 
{
   int n=*N; 
   int nvec=*NVEC; 
   double tol=*TOL;
   int maxit=*MAXIT; 
   int verbosity_level=*VERBOSITY_LEVEL; 
   int use_x_values=*USE_X_VALUES; 
   int lobpcg_preconditioner=*LOBPCG_PRECONDITIONER; 

 /* lobpcg data */
   serial_Multi_Vector * x;
   mv_MultiVectorPtr xx; 
   double * resid;
   int iterations;
   lobpcg_Tolerance lobpcg_tol;
   mv_InterfaceInterpreter ii;
   lobpcg_BLASLAPACKFunctions blap_fn;
   int ierr;

   void (*precfun)( void*, void*, void* );

 /* prepare data for LOBPCG */
  /* create multivector */
   x = serial_Multi_VectorCreate(n, nvec);
   serial_Multi_VectorInitialize(x);

   if (use_x_values == 1)
   {
      /* use predefined values in vec */
      mv_load_values (vec, nvec, n, x);
   }
   else if (use_x_values == 0)
   {
      /* fill it with random numbers */
      serial_Multi_VectorSetRandomValues(x, 1);
   }
   else
   {
      printf("\nLOBPCG_DRIVER: unknown value of use supplied values as start %d \n",use_x_values);
   }

  /* set correct pointer to preconditioner */
   if (lobpcg_preconditioner == 1)
   {
      /* use BDDC */
      precfun = mvecmultM;
   }
   else
   {
      /* use NULL */
      precfun = NULL;
   }


/* get memory for eigenvalues, eigenvalue history, residual norms, residual norms history */
   /* memory is already prepared for eigenvalues */

   /* request memory for resid. norms */
   resid = (double *)malloc(sizeof(double)*nvec);

   /* set tolerances */
   lobpcg_tol.absolute = tol;
   lobpcg_tol.relative = 1e-50;

/* setup interface interpreter and wrap around "x" another structure */

   SerialSetupInterpreter( &ii );
   xx = mv_MultiVectorWrap( &ii, x, 0);

/* set pointers to lapack functions */
   blap_fn.dpotrf = lobpcg_dpotrf;
   blap_fn.dsygv  = lobpcg_dsygv;

   /* execute lobpcg */
   ierr = lobpcg_solve_double( 
	  xx,
          NULL,
          mvecmultA,
          NULL,
          mvecmultB,
          NULL,
	  precfun, /* precfun is set above: for option with preconditioner, use mvecmultM here, without it, use NULL */
          NULL,
          blap_fn,
          lobpcg_tol,
          maxit,
          verbosity_level,
          &iterations,
          
          /* eigenvalues; "lambda_values" should point to array  containing <blocksize> doubles where <blocksi
          ze> is the width of multivector "blockVectorX" */
          lambda,
          
          /* eigenvalues history; a pointer to the entries of the  <blocksize>-by-(<maxIterations>+1) matrix s
          tored
          in  fortran-style. (i.e. column-wise) The matrix may be  a submatrix of a larger matrix, see next
          argument; If you don't need eigenvalues history, provide NULL in this entry */
          NULL,
          
          /* global height of the matrix (stored in fotran-style)  specified by previous argument */
          0,
          
          /* residual norms; argument should point to array of <blocksize> doubles */
          resid,
          
          /* residual norms history; a pointer to the entries of the  <blocksize>-by-(<maxIterations>+1) matri
          x
          stored in  fortran-style. (i.e. column-wise) The matrix may be  a submatrix of a larger matrix, see
          next
          argument If you don't need residual norms history, provide NULL in this entry */
          NULL,
          
          /* global height of the matrix (stored in fotran-style)  specified by previous argument */
          0
          );

/*   if (ierr)
     { 
        printf("LOBPCG exited with nonzero code: %d \n",ierr); 
     }
*/

/* print eigenvectors to file */
/*   serial_Multi_VectorPrint(x, "eigenvectors"); */
   mv_extract_values (x, vec, nvec);

   serial_Multi_VectorDestroy(x);
   mv_MultiVectorDestroy(xx);
   free(resid); 

   *ITERATIONS = iterations;
   *IERR       = ierr;
 
}

