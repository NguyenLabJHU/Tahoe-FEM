/*****************************************************************************/
/*                                                                           */
/*   (C) Copyright IBM Corporation, 1997                                     */
/*   (C) Copyright Regents of the University of Minnesota, 1997              */
/*                                                                           */
/*   dpspacef.c                                                              */
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
/* $Id: dpspacef.c,v 1.1 2004-12-10 20:26:44 paklein Exp $ */
/*****************************************************************************/

#include "pspaces.h"

#ifndef FCMATCH
void DPSPACEF(int *rowdista,int *aptrs,int *ainds,double *avals,
            int *ioptions,double *doptions,long *pspcomm,MPI_Comm *pcomm) {
  Fdpspacef(rowdista,aptrs,ainds,avals,ioptions,doptions,pspcomm,pcomm);
}
#endif

void Fdpspacef(int *rowdista,int *aptrs,int *ainds,double *avals,
             int *ioptions,double *doptions,long *pspcomm,MPI_Comm *pcomm)
{

int *order,*sizes,*Irowdist,*wrkint,*ranmasks,*whichsnode;
int myid,pp,dd,mynnodes,N,i,j,k,l,m,mybeginleaf,pasize,lgblk;
int *tptrs,*tinds,*lptrs,*tainds,*paptrs,*painds,*linds,*sup,*supinds;
int *tlinds,*tsup,*tsupinds,*psupinds,*psup,*plinds,*tsind,*iptrs,*lc,*ifopts;
int ntemp,*temparr1,*temparr2,nnz_ga,serialorder,nrecvsizs,ilbd;
int supsize,supindsize,tlplindsize,nnzleafa,root,rbits,cbits,myidv,myidh;
int leaflindsize,lindsize,info,lindsptr,lvalsptr,supptr,supindsptr,iu;
int nnz,ns,lspace,blk,wsolvesize,taivsize,ilcsize,lcsize,checksymm,sortinds;
double *tavals,*pavals,*lvals,dfopts[7];
double *lbal,lopc,vl,limbalfac,memory_pscomm_i,memory_pscomm_d;
int *cinfo;
MPI_Comm *pmcomm,comm;
PTRS *ap;

int *dbgsizes,dbgpp,dbgdd;
double totsan,sanity;

  if(!(pmcomm = (MPI_Comm *)malloc(sizeof(MPI_Comm)))) {
    printf("[dpspacef.c:0] memory allocation error\n");
    MPI_Abort((*pcomm),0);
  }

  MPI_Comm_dup((*pcomm),pmcomm);

  comm = *pmcomm;

  MPI_Comm_rank(comm,&myid);
  MPI_Comm_size(comm,&pp);

  dd = floor(log((double)pp)/log(2.0));

  if(!myid && (1<<dd != pp || pp<2)) {
    printf("The number of processors must be > 1, and a power of 2.\n");
    MPI_Abort(comm,0);
  }

  ntemp = pp+1+F_NOPTS;
  if(!(temparr1 = (int *)malloc(ntemp*sizeof(int)))) {
    printf("[dpspacef.c:1] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  for(i=0;i<pp+1;i++) temparr1[i] = rowdista[i];
  for(j=0;j<F_NOPTS;j++) temparr1[j+i] = ioptions[j];

  if(!myid) {
    if(!(temparr2 = (int *)malloc(ntemp*sizeof(int)))) {
      printf("[dpspacef.c:2] memory allocation error\n");
      MPI_Abort(comm,0);
    }
  }
  MPI_Reduce(temparr1,temparr2,ntemp,MPI_INT,MPI_BXOR,0,comm);
  if(!myid) {
    for(i=0;i<ntemp;i++) {
      if(temparr2[i]) {
        printf("DPSPACEF: Global parameters not matching on all processors!\n");
      if(i<pp+1)
        printf("DPSPACEF: Check rowdista[%d] !\n",i);
      else
        printf("DPSPACEF: Check ioptions[%d] !\n",i-(pp+1));
        MPI_Abort(comm,0);
      }
    }
    free(temparr2);
  }
  MPI_Barrier(comm);
  free(temparr1);

  N = rowdista[pp];

  lgblk = 0;
  while((blk = 1<<lgblk) < ioptions[0]) lgblk++;
  ioptions[0] = blk;
  checksymm = ioptions[1];
  if(checksymm!=1) checksymm = 0;
  sortinds = ioptions[2];
  if(sortinds!=1) sortinds = 0;
  serialorder = ioptions[3];
  if(serialorder!=1) serialorder = 0;

  mynnodes = rowdista[myid+1]-rowdista[myid];

  memory_pscomm_i = 0.0;
  memory_pscomm_d = 0.0;

  memory_pscomm_i += (double)(pp+1);
  if(!(Irowdist = (int *)malloc((pp+1)*sizeof(int)))) {
    printf("[dpspacef.c:3] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  for(i=0;i<pp+1;i++) Irowdist[i] = rowdista[i];

  memory_pscomm_i += (double)mynnodes;
  if(!(order = (int *)malloc(mynnodes*sizeof(int)))) {
    printf("[dpspacef.c:4] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  if(!(sizes = (int *)malloc(2*pp*sizeof(int)))) {
    printf("[dpspacef.c:5] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  if(serialorder) {
    dbgdd = dd;  /* change this to a fixed value to obtain same ordering */
                 /* irrespective of number of processors.                */
    dbgpp = 1<<dbgdd;

    if(!(dbgsizes = (int *)malloc(2*dbgpp*sizeof(int)))) {
      printf("[dpspacef.c:6] memory allocation error (dgbsizes=%d) \n",2*dbgpp);
      MPI_Abort(comm,0);
    }

    porder(Irowdist,aptrs,ainds,order,dbgsizes,&myid,&pp,&serialorder,&dbgpp,
           pmcomm);

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

    porder(Irowdist,aptrs,ainds,order,sizes,&myid,&pp,&serialorder,&pp,
           pmcomm);

  }

  i = (11*pp > 2*N) ? 11*pp : 2*N;
  if(!(wrkint = (int *)malloc(i*sizeof(int)))) {
    printf("[dpspacef.c:7] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  memory_pscomm_i += (double)(3*N);
  if(!(tptrs = (int *)malloc(3*N*sizeof(int)))) {
    printf("[dpspacef.c:8] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  memory_pscomm_i += (double)N;
  if(!(tinds = (int *)malloc(N*sizeof(int)))) {
    printf("[dpspacef.c:9] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  memory_pscomm_i += (double)(3*N);
  if(!(lptrs = (int *)malloc(3*N*sizeof(int)))) {
    printf("[dpspacef.c:10] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  i = 5*(2*pp-1);
  memory_pscomm_i += (double)i;
  if(!(ranmasks = (int *)malloc(i*sizeof(int)))) {
    printf("[dpspacef.c:11] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  taivsize = aptrs[2*(mynnodes-1)]+aptrs[2*(mynnodes-1)+1]-1;
  if(!(tainds = (int *)malloc(taivsize*sizeof(int)))) {
    printf("[dpspacef.c:12] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  if(!(tavals = (double *)malloc(taivsize*sizeof(double)))) {
    printf("[dpspacef.c:13] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  if(!(whichsnode = (int *)malloc(mynnodes*sizeof(int)))) {
    printf("[dpspacef.c:14] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  emovea(&N,&dd,&pp,&lgblk,&myid,&mynnodes,rowdista,order,sizes,
       &mybeginleaf,&pasize,aptrs,ainds,avals,tainds,tavals,
       lptrs,tinds,tptrs,wrkint,ranmasks,whichsnode,&checksymm,&sortinds,
       pmcomm);

  if(!(pavals = (double *)malloc(pasize*sizeof(double)))) {
    printf("[dpspacef.c:15] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  if(!(painds = (int *)malloc(pasize*sizeof(int)))) {
    printf("[dpspacef.c:16] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  if(!(paptrs = (int *)malloc(2*N*sizeof(int)))) {
    printf("[dpspacef.c:17] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  pmovea(&N,&dd,&pp,&lgblk,&myid,&mynnodes,order,paptrs,painds,
       pavals,aptrs,tainds,tavals,wrkint,ranmasks,whichsnode,
       &(lptrs[N]),&nrecvsizs,pmcomm);

  free(whichsnode);
  free(tainds);
  free(tavals);

  gentree(&N,paptrs,painds,pavals,&mybeginleaf,sizes,lptrs,tptrs,tinds,
        &dd,&myid,&pp,wrkint,&pasize,&(lptrs[N]),&nrecvsizs,pmcomm);

  if(!(pavals = (double *)realloc(pavals,pasize*sizeof(double)))) {
    printf("[dpspacef.c:18] %d: memory reallocation error\n",myid);
    MPI_Abort(comm,0);
  }
  if(!(painds = (int *)realloc(painds,pasize*sizeof(int)))) {
    printf("[dpspacef.c:19] %d: memory reallocation error\n",myid);
    MPI_Abort(comm,0);
  }

  if(!(cinfo = (int *)malloc(N*sizeof(int)))) {
    printf("[dpspacef.c:20] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  for(i=0;i<N;cinfo[i++]=0);

  supsize = 4*(sizes[myid]+dd); 

  supindsize = sizes[2*pp-2]; 
  k = 1;
  m = supindsize;
  for(i=1;i<dd+1;i++) {
    j = myid>>(dd-i);
    l = 2*(pp-k)-2; 
    if(j&1) l=l+1; 
    m = m+sizes[l];
    supindsize += m; 
    k = 2*pp-1-l;
  }

  tlplindsize = supindsize-m;
  supindsize <<= 1;
  supindsize += 1;

  nnzleafa = 0;
  for(i=mybeginleaf;i<mybeginleaf+sizes[myid];i++) {
    nnzleafa += paptrs[2*i+1];
    cinfo[i] = 1;
  }

  if(!(sup = (int *)malloc(supsize*sizeof(int)))) {
    printf("[dpspacef.c:21] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  if(!(supinds = (int *)malloc(supindsize*sizeof(int)))) {
    printf("[dpspacef.c:22] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  root = N-1;

  rbits = dd>>1;
  cbits = dd-(dd>>1);

  myidv = 0;
  i = myid>>1;
  for(j=0;j<rbits;j++) {
    myidv |= (i&1)<<j;
    i>>=2;
  }

  myidh = 0;
  i = myid;
  for(j=0;j<cbits;j++) {
    myidh |= (i&1)<<j;
    i>>=2;
  }

  for(j=0;j<N;j++) tptrs[3*j+2] = 0;

  leaflindsize = (int)((float)nnzleafa*FILLFACTOR);
  lindsize = tlplindsize+leaflindsize;

  if(!(linds = (int *)malloc(lindsize*sizeof(int)))) {
    printf("[dpspacef.c:23] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  if(!(lbal = (double *)malloc(4*pp*sizeof(double)))) {
    printf("[dpspacef.c:24] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  for(i=0;i<2*pp-1;i++) lbal[i] = 0.0;

  info = 1;
  while(info) {
    leaflindsize += 2*nnzleafa;
    lindsize = tlplindsize+leaflindsize;

    if(!(linds = (int *)realloc(linds,lindsize*sizeof(int)))) {
      printf("[dpspacef.c:25] memory allocation error\n");
      MPI_Abort(comm,0);
    }

    parsymb(&root,paptrs,painds,tptrs,tinds,lptrs,linds,&lindsptr,
          &lvalsptr,&nnz,lbal,sizes,sup,supinds,&supptr,&supindsptr,
          &iu,&leaflindsize,&dd,&N,&myid,&lgblk,&myidv,&myidh,
          wrkint,&(wrkint[N]),&info,cinfo,pmcomm);

  }

  MPI_Allreduce(&pasize,&nnz_ga,1,MPI_INT,MPI_SUM,comm);
  doptions[2] = nnz_ga;
  doptions[3] = nnz;

  ilbd = 2*pp;
  MPI_Reduce(lbal,&(lbal[ilbd]),2*pp-1,MPI_DOUBLE,MPI_SUM,0,comm);

  if(!myid) {
    lopc = 0.0;
    for(i=0;i<2*pp-1;i++) lopc += lbal[ilbd+i];
    doptions[4] = lopc;

    j = ilbd+pp;
    k = pp;
    m = ilbd;
    for(i=1;i<=dd;i++) {
      ntemp = 1<<i;
      for(l=m;l<m+k-1;l+=2) {
        vl = lbal[l] > lbal[l+1] ? lbal[l] : lbal[l+1];
        lbal[j] = lbal[j]/(double)(ntemp) + vl;
        j++;
      }
      k >>= 1;
      m = j-k;
    }

    limbalfac = lbal[m]/(lopc/pp);
    doptions[5] = limbalfac;
  }

  free(sizes);
  free(lbal);

  ns = (supptr-1)/4;
  lspace = lvalsptr-1;
  lindsize = lindsptr-1;
  supsize = supptr-1;
  supindsize = supindsptr-1+iu;

  memory_pscomm_i += (double)lindsize;
  if(!(linds = (int *)realloc(linds,lindsize*sizeof(int)))) {
    printf("[dpspacef.c:26] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  memory_pscomm_i += (double)supsize;
  if(!(sup = (int *)realloc(sup,supsize*sizeof(int)))) {
    printf("[dpspacef.c:27] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  if(!(supinds = (int *)realloc(supinds,supindsize*sizeof(int)))) {
    printf("[dpspacef.c:28] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  memory_pscomm_i += (double)(supindsptr-1);
  if(!(tsind = (int *)malloc((supindsptr-1)*sizeof(int)))) {
    printf("[dpspacef.c:29] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  for(i=0;i<supindsptr-1;i++) tsind[i] = supinds[i];

  memory_pscomm_d += (double)lspace;
  if(!(lvals = (double *)malloc(lspace*sizeof(double)))) {
    printf("[dpspacef.c:30] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  memory_pscomm_i += (double)(2*N);
  if(!(iptrs = (int *)malloc(2*N*sizeof(int)))) {
    printf("[dpspacef.c:31] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  lcsize = ((nnz_ga+N)*3)/dd;
  ilcsize = lcsize;

  memory_pscomm_i += (double)ilcsize;
  if(!(lc = (int *)malloc(ilcsize*sizeof(int)))) {
    printf("[dpspacef.c:32] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  memory_pscomm_i += (double)IFOPTS_SIZE;
  if(!(ifopts = (int *)malloc(IFOPTS_SIZE*sizeof(int)))) {
    printf("[dpspacef.c:33] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  for(i=0;i<7;i++) dfopts[i] = 0.0;
  for(i=0;i<5;i++) ifopts[i] = 0;

  info = 0;
  parfact1(&N,paptrs,painds,pavals,lptrs,linds,lvals,tptrs,tinds,sup,wrkint,
         &(wrkint[N]),&root,&dd,&lgblk,&blk,&myid,cinfo,supinds,&supindsptr,
         dfopts,ifopts,lc,iptrs,&lcsize,&wsolvesize,&info,pmcomm);

  if(info) {
    printf("Non Positive Definite Matrix found!\n");
    MPI_Abort(comm,0);
  }

  if(lcsize<ilcsize) {
    memory_pscomm_i -= (double)(ilcsize-lcsize);
    if(!lcsize) lcsize = 1;
    if(!(lc = (int *)realloc(lc,sizeof(int)*lcsize))) {
      printf("[dpspacef.c:34] memory allocation error\n");
      MPI_Abort(comm,0);
    }
  }

  doptions[0] = memory_pscomm_i*sizeof(int)+memory_pscomm_d*sizeof(double);
  doptions[1] = 0.0;

  ifopts[0] = N;
  ifopts[1] = dd;
  ifopts[2] = lgblk;
  ifopts[3] = wsolvesize;
  ifopts[4] = ns;
  ifopts[5] = supindsptr-1;
  ifopts[6] = myid;
  ifopts[7] = myidh;
  ifopts[8] = myidv;
  ifopts[9] = supsize;
  ifopts[10] = nnz;
  ifopts[11] = pasize;
  ifopts[12] = lspace;

  if(!(ap = (PTRS *)malloc(sizeof(PTRS)))) {
    printf("[dpspacef.c:35] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  ap->Ylink = NULL;
  ap->howManyNfacts = -1;

  ap->rowdist = Irowdist;
  ap->order = order;
  ap->ranmasks = ranmasks;

  /*** following are specifically for DPSPACEN ***/
  ap->wrkint = NULL;
  ap->cinfo = NULL;
  ap->paptrs = NULL;
  ap->painds = NULL;
  ap->supinds = NULL;
  /***********************************************/

  ap->lptrs = lptrs;
  ap->linds = linds;
  ap->tptrs = tptrs;
  ap->tinds = tinds;
  ap->sup = sup;
  ap->tsind = tsind;
  ap->lc = lc;
  ap->iptrs = iptrs;
  ap->ifopts = ifopts;
  ap->lvals = lvals;
  ap->pmcomm = pmcomm;
  
  *pspcomm = check_handle((long)ap,0);

/*
compmysan(&sanity,&N,lptrs,lvals,linds,cinfo);

MPI_Reduce(&sanity,&totsan,1,MPI_DOUBLE,MPI_SUM,0,comm);
if(!myid) {
  printf(" L Sanity = %20.15le\n",totsan/(double)N);
}
*/

  free(wrkint);
  free(cinfo);
  free(paptrs);
  free(painds);
  free(supinds);
  free(pavals);

}
