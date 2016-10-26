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

/*****************************************
* Wrapper of METIS_PartGraphRecursive and METIS_PartGraphKWay functions
******************************************/

void graph_divide_c( int *numflag, int *graphtype, int *nvertex, int *xadj, int *lxadj, int *adjncy, int *ladjncy, 
                     int *vwgt, int *lvwgt, int *adjwgt, int *ladjwgt, int *nsub, 
                     int *edgecut, int *part, int *lpart )
{
  int i;
  int *options;
  int wgtflag;

  /***********************************/
  /* Try and take care of bad inputs */
  /***********************************/
  if (numflag == NULL || graphtype == NULL || nvertex == NULL || 
      xadj == NULL || lxadj == NULL ||
      adjncy == NULL || ladjncy == NULL || vwgt == NULL || lvwgt == NULL || adjwgt == NULL || ladjwgt == NULL ||
      nsub == NULL || edgecut == NULL || part == NULL || lpart == NULL) {
     printf("ERROR in GRAPH_DIVIDE_C: One or more required parameters is NULL. Aborting.\n");
     abort();
  }


  /* prepare options */
#if (METIS_VER_MAJOR >= 5)
  if ( sizeof(idx_t) != sizeof(int) ) {
     printf("ERROR in GRAPH_DIVIDE_C: Wrong type of integers for METIS.\n");
     abort();
  }
  /*printf(" METIS >=5.0 recognized.\n");*/
  int ncon = 1;
  /*real_t ubvec[1];*/
  /*ubvec[0] = 1.001;*/
  /*int *options = NULL;*/
  options = malloc(METIS_NOPTIONS * sizeof(int));

  for (i = 0;i < METIS_NOPTIONS;i++) {
     options[i] = -1;
  }
  options[METIS_OPTION_OBJTYPE]   = METIS_OBJTYPE_CUT;
  options[METIS_OPTION_CTYPE]     = METIS_CTYPE_RM;
  options[METIS_OPTION_IPTYPE]    = METIS_IPTYPE_GROW;
  options[METIS_OPTION_RTYPE]     = METIS_RTYPE_GREEDY;
  options[METIS_OPTION_NCUTS]     = 1;
  options[METIS_OPTION_NSEPS]     = 1;
  options[METIS_OPTION_NUMBERING] = *numflag;
  options[METIS_OPTION_NITER]     = 10;
  options[METIS_OPTION_SEED]      = 12345;
  options[METIS_OPTION_MINCONN]   = 1;
  options[METIS_OPTION_CONTIG]    = 0;
  options[METIS_OPTION_COMPRESS]  = 0;
  options[METIS_OPTION_CCORDER]   = 0;
  options[METIS_OPTION_UFACTOR]   = 0;
  /*options[METIS_OPTION_DBGLVL]    = METIS_DBG_INFO;*/
  options[METIS_OPTION_DBGLVL]    = 0;
#else
  /*printf(" METIS < 5.0 recognized.\n");*/
  /* weights */
  if (*graphtype == 1) {
     wgtflag = 1;
  }
  else {
     wgtflag = 0;
  }
  options = malloc(8 * sizeof(int));
  for ( i = 0; i < 8; i++ ) {
     options[i] = 0;
  }
#endif

  /* Initialize parts */
  for ( i = 0; i < *lpart; i++ ) {
     part[i] = *numflag;
  }

  /* divide graph */
  if (*nsub == 0) {
     printf("ERROR in GRAPH_DIVIDE_C: Illegal number of subdomains %d,  Aborting.\n", *nsub);
     abort();
  }
  else if (*nsub == 1) {
     *edgecut = 0;
  }
  else if (*nsub > 1 && *nsub <= 8) {

#if (METIS_VER_MAJOR >= 5)
     options[METIS_OPTION_PTYPE]     = METIS_PTYPE_RB;
     options[METIS_OPTION_UFACTOR]   = 1;
     METIS_PartGraphRecursive(nvertex,&ncon,xadj,adjncy,vwgt,NULL,adjwgt,nsub,NULL,NULL,options,edgecut,part);
#else
     METIS_PartGraphRecursive(nvertex,xadj,adjncy,vwgt,adjwgt,&wgtflag,numflag,nsub,options,edgecut,part);
#endif
  }
  else {
#if (METIS_VER_MAJOR >= 5)
     options[METIS_OPTION_PTYPE]     = METIS_PTYPE_KWAY;
     options[METIS_OPTION_UFACTOR]   = 30;
     METIS_PartGraphKway(nvertex,&ncon,xadj,adjncy,vwgt,NULL,adjwgt,nsub,NULL,NULL,options,edgecut,part);
#else
     METIS_PartGraphKway(nvertex,xadj,adjncy,vwgt,adjwgt,&wgtflag,numflag,nsub,options,edgecut,part);
#endif
  }

  free(options);

  return;
}
