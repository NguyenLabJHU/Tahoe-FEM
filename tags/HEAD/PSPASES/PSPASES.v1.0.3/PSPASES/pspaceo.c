/*****************************************************************************/
/*                                                                           */
/*   (C) Copyright IBM Corporation, 1997                                     */
/*   (C) Copyright Regents of the University of Minnesota, 1997              */
/*                                                                           */
/*   pspaceo.c                                                               */
/*                                                                           */
/*   Written by Mahesh Joshi, U of MN.                                       */
/*                                                                           */
/*****************************************************************************/
/*                                                                           */
/* This code is meant to be used solely for educational, research, and       */
/* benchmarking purposes by non-profit institutions and US government        */
/* agencies only.  Use by any other organization requires prior written      */
/* permission from both IBM Corporation and the University of Minnesota.     */
/* The software may not be sold or redistributed.  One may make copies       */
/* of the software or modify it for their use provided that the copies,      */
/* modified or otherwise, are not sold or distributed, are used under the    */
/* same terms and conditions, and this notice and any part of the source     */
/* code that follows this notice are not separated.                          */
/*                                                                           */
/* As unestablished research software, this code is provided on an           */
/* ``as is'' basis without warranty of any kind, either expressed or         */
/* implied, including but not limited to implied warranties of               */
/* merchantability and fitness for a particular purpose.  IBM does not       */
/* warrant that the functions contained in this software will meet the       */
/* user's requirements or that the operation of its routines will be         */
/* uninterrupted or error-free.  Acceptance and use of this program          */
/* constitutes the user's understanding that he/she will have no recourse    */
/* to IBM for any actual or consequential damages, including, but not        */
/* limited to, lost profits or savings, arising out of the use or inability  */
/* to use these libraries.  Even if the user informs IBM of the possibility  */
/* of such damages, IBM expects the user to accept the risk of any such      */
/* harm, or the user shall not attempt to use these libraries for any        */
/* purpose.                                                                  */
/*                                                                           */
/* The downloading, compiling, or executing any part of this software        */
/* constitutes an implicit agreement to these terms.  These terms and        */
/* conditions are subject to change at any time without prior notice.        */
/*                                                                           */
/*****************************************************************************/
/* $Id: pspaceo.c,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $ */
/*****************************************************************************/

#include "pspaces.h"

#ifndef FCMATCH
void PSPACEO(int *rowdista,int *aptrs,int *ainds,int *order,
             int *sizes,int *ioptions,MPI_Comm *pcomm) {
  Fpspaceo(rowdista,aptrs,ainds,order,sizes,ioptions,pcomm);
}
#endif

void Fpspaceo(int *rowdista,int *aptrs,int *ainds,int *order,
             int *sizes,int *ioptions,MPI_Comm *pcomm)
{

int myid,pp,dd,mynnodes,N,i,j,k,l,m;
int ntemp,*temparr1,*temparr2,serialorder;
MPI_Comm comm;
int *dbgsizes,dbgpp,dbgdd;

  MPI_Comm_dup((*pcomm),&comm);

  MPI_Comm_rank(comm,&myid);
  MPI_Comm_size(comm,&pp);

  dd = floor(log((double)pp)/log(2.0));

  if(!myid && ((1<<dd != pp) || pp < 2) ) {
    printf("The number of processors must be > 1, and a power of 2.\n");
    MPI_Abort(comm,0);
  }

  ntemp = pp+1+O_NOPTS;
  if(!(temparr1 = (int *)malloc(ntemp*sizeof(int)))) {
    printf("memory allocation error\n");
    MPI_Abort(comm,0);
  }
  for(i=0;i<pp+1;i++) temparr1[i] = rowdista[i];
  for(j=0;j<O_NOPTS;j++) temparr1[j+i] = ioptions[j];

  if(!myid) {
    if(!(temparr2 = (int *)malloc(ntemp*sizeof(int)))) {
      printf("memory allocation error\n");
      MPI_Abort(comm,0);
    }
  }
  MPI_Reduce(temparr1,temparr2,ntemp,MPI_INT,MPI_BXOR,0,comm);
  if(!myid) {
    for(i=0;i<ntemp;i++) {
      if(temparr2[i]) {
        printf("PSPACEO: Global parameters not matching on all processors!\n");
        if(i<pp+1)
          printf("Check rowdista[%d]!\n",i);
        else
          printf("Check ioptions[%d] !\n",i-(pp+1));
        MPI_Abort(comm,0);
      }
    }
    free(temparr2);
  }
  MPI_Barrier(comm);
  free(temparr1);

  N = rowdista[pp];

  serialorder = ioptions[3];
  if(serialorder!=1) serialorder = 0;

  mynnodes = rowdista[myid+1]-rowdista[myid];

  if(serialorder) {
    dbgdd = dd;  /* change this to a fixed value to obtain same ordering */
                 /* irrespective of number of processors.                */
    dbgpp = 1<<dbgdd;

    if(!(dbgsizes = (int *)malloc(2*dbgpp*sizeof(int)))) {
      printf("memory allocation error (dgbsizes=%d) \n",2*dbgpp);
      MPI_Abort(comm,0);
    }

    porder(rowdista,aptrs,ainds,order,dbgsizes,&myid,&pp,&serialorder,&dbgpp,
           &comm);

    l=0;
    for(j=0;j<pp;j++) 
      sizes[j] = 0;
    for(i=dbgdd-dd;i>=0;i--)
      for(j=0;j<pp;j++)
        for(k=0;k<(1<<i);k++)
          sizes[j] += dbgsizes[l++];
    for(j=pp;j<2*pp-1;j++)
      sizes[j] = dbgsizes[l++];

    free(dbgsizes);

  } else {

    i = max(pp*50,5000);
    if( N < i ) {
      serialorder = 1;
      if(!myid) 
	printf("  [WARNING] Small matrix. Using serial ordering.\n");
    }

    porder(rowdista,aptrs,ainds,order,sizes,&myid,&pp,&serialorder,&pp,
           &comm);

  }

  MPI_Comm_free(&comm);

}
