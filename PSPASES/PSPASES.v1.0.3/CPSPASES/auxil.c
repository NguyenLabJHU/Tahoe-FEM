/*****************************************************************************/
/*                                                                           */
/*   (C) Copyright IBM Corporation, 1997                                     */
/*   (C) Copyright Regents of the University of Minnesota, 1997              */
/*                                                                           */
/*   auxil.c                                                                 */
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
/* $Id: auxil.c,v 1.1 2004-12-10 20:26:44 paklein Exp $ */
/*****************************************************************************/

#include <pspaces.h>

void initrnd (int *iseed)
{
void srand48();
long jseed;

jseed = *iseed + 69805;
srand48(jseed);
return;
}

void mydrand48n (double *rnum)
{
double drand48();

*rnum = drand48();
}

#ifndef FCMATCH
void CHECKB_AX(int *rowdista,int *aptrs,int *ainds,double *avals,
               int *rowdistb,int *pnrhs,double *b,int *pldb,double *x,
               int *pldx,double *perr,MPI_Comm *pcomm)
{
  Fcheckb_ax(rowdista,aptrs,ainds,avals,rowdistb,pnrhs,b,pldb,x,pldx,perr,
             pcomm);
}
#endif

void Fcheckb_ax(int *rowdista,int *aptrs,int *ainds,double *avals,
                int *rowdistb,int *pnrhs,double *b,int *pldb,double *x,
                int *pldx,double *perr,MPI_Comm *pcomm)
{
int myid,pp,N;
MPI_Comm comm;

MPI_Comm_dup((*pcomm),&comm);

MPI_Comm_rank(comm,&myid);
MPI_Comm_size(comm,&pp);

N = rowdista[pp];
if(N!=rowdistb[pp]) {
  if(!myid) 
    printf("Dimensions of A and b do not match\n");
  MPI_Barrier(comm);
  MPI_Abort(comm,0);
}

db_ax(&N,rowdista,rowdistb,pnrhs,aptrs,ainds,avals,b,pldb,x,pldx,
      &myid,&pp,perr,&comm);

MPI_Comm_free(&comm);
}
