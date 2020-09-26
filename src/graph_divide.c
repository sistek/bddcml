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

#include "metis.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

void graph_divide_c(int *numflag, int *graphtype, int *nvertex, int *xadj,
                    int *lxadj, int *adjncy, int *ladjncy, int *vwgt,
                    int *lvwgt, int *adjwgt, int *ladjwgt, int *nsub,
                    int *contiguous_clusters,
                    int *edgecut, int *part, int *lpart )
{

// Check version of METIS.
#if (METIS_VER_MAJOR < 5)
   printf("ERROR in GRAPH_DIVIDE_C: Unsupported version of METIS."
          " Upgrade to version >= 5.X.X \n");
   abort();
#endif

  // Check inputs.
  if (numflag == NULL || graphtype == NULL || nvertex == NULL || xadj == NULL ||
      lxadj == NULL || adjncy == NULL || ladjncy == NULL || vwgt == NULL ||
      lvwgt == NULL || adjwgt == NULL || ladjwgt == NULL || nsub == NULL ||
      edgecut == NULL || part == NULL || lpart == NULL) {
     printf("ERROR in GRAPH_DIVIDE_C: One or more required parameters is NULL. "
            "Aborting.\n");
     abort();
  }

  // Check that type of integers is compatible.
  if (sizeof(idx_t) != sizeof(int)) {
     printf("ERROR in GRAPH_DIVIDE_C: Wrong type of integers for METIS.\n");
     abort();
  }

  // Prepare the options array.
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);

  // Set contiguous subdomains.
  options[METIS_OPTION_CONTIG] = *contiguous_clusters; 

  // Use the given numbering (0 from C, 1 from Fortran).
  options[METIS_OPTION_NUMBERING] = *numflag;

  // Switch off debugging info.
  options[METIS_OPTION_DBGLVL] = 0;
  //options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO;

  // Allow load imbalance in division 20%
  options[METIS_OPTION_UFACTOR] = 200;

  // Initialize parts.
  for (int i = 0; i < *lpart; i++ ) {
     part[i] = *numflag;
  }

  // Number of constraints.
  int ncon = 1;

  // Error code.
  int ierr = METIS_OK;

  // Divide graph.
  if (*nsub == 0) {
     // Cannot divide any mesh into 0 number of subdomains.
     printf("ERROR in GRAPH_DIVIDE_C: Illegal number of subdomains %d, "
            " Aborting.\n", *nsub);
     abort();
  }
  else if (*nsub == 1) {
     // This is a simple case, just return at this point.
     *edgecut = 0;
  }
  else if (*nsub > 1 && *nsub <= 8) {
     // For small number of subdomains, call recursive bisection.
     ierr = METIS_PartGraphRecursive(nvertex,&ncon,xadj,adjncy,vwgt,NULL,adjwgt,
                                     nsub, NULL,NULL,options,edgecut,part);
  }
  else {
     // For larger number of subdomains, call k-way algorithm.
     ierr = METIS_PartGraphKway(nvertex,&ncon,xadj,adjncy,vwgt,NULL,adjwgt,
                                nsub, NULL,NULL,options,edgecut,part);
  }

  // Check error code.
  if (ierr != METIS_OK) {
     printf("ERROR in GRAPH_DIVIDE_C: Metis error %d . Aborting... \n ",
            ierr);
  }

  return;
}
