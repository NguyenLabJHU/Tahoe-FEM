/*****************************************************************************/
/*                                                                           */
/*   (C) Copyright IBM Corporation, 1997                                     */
/*   (C) Copyright Regents of the University of Minnesota, 1997              */
/*                                                                           */
/*   pspacey.c                                                              */
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
/* $Id: pspacey.c,v 1.2 2005-01-15 07:37:34 paklein Exp $ */
/*****************************************************************************/

#include "pspaces.h"

#ifndef FCMATCH
void PSPACEY(int *rowdista,int *aptrs,int *ainds,int *order,int *sizes,
             int *ioptions,double *doptions,long *pspcomm,MPI_Comm *pcomm) {
  Fpspacey(rowdista,aptrs,ainds,order,sizes,ioptions,doptions,pspcomm,pcomm);
}
#endif

void Fpspacey(int *rowdista,int *aptrs,int *ainds,int *order,int *sizes,
             int *ioptions,double *doptions,long *pspcomm,MPI_Comm *pcomm)
{

int *Iorder,*Irowdist,*wrkint,*ranmasks,*whichsnode,*wa1;
int myid,pp,dd,mynnodes,N,i,j,k,l,m,mybeginleaf,pasize,lgblk,IWN;
int *tptrs,*tinds,*lptrs,*tainds,*paptrs,*painds,*linds,*sup,*supinds;
int *tsind,*cinfo,*iptrs,*lc,*ifopts,maxnzpercol,nsend,nrecvsizs;
int ntemp,*temparr1,*temparr2,nnz_ga;
int wsize0,wsize1,ibuflen,dbuflen,iwspace,node,stakptr,nptr;
int supsize,supindsize,tlplindsize,nnzleafa,root,rbits,cbits,myidv,myidh;
int leaflindsize,lindsize,info,lindsptr,lvalsptr,supptr,supindsptr,iu;
int nnz,ns,lspace,blk,wsolvesize,taivsize,ilcsize,lcsize,checksymm,sortinds;
double memory_pscomm,parfact_memreq,moveav_memreq,dpspacen_memreq;
double pparfact_memreq,*lbal,lopc,vl,limbalfac;
int ilbd;
MPI_Comm *pmcomm,comm;
PTRS *ap;
int power2;

  if(!(pmcomm = (MPI_Comm *)malloc(sizeof(MPI_Comm)))) {
    printf("[pspacey.c:0] memory allocation error\n");
    MPI_Abort((*pcomm),0);
  }

  MPI_Comm_dup((*pcomm),pmcomm);

  comm = *pmcomm;

  MPI_Comm_rank(comm,&myid);
  MPI_Comm_size(comm,&pp);

#if 0
  dd = floor(log((double)pp)/log(2.0));
  if(!myid && (1<<dd != pp || pp<2)) {
    printf("The number of processors must be > 1, and a power of 2.\n");
    MPI_Abort(comm,0);
  }
#endif

	power2 = 2;
	dd = 0;
	while (pp != power2) {
		if (power2 > pp) {
			printf("The number of processors must be > 1, and a power of 2.\n");
    		MPI_Abort(comm,0);
		}
		dd++;
		power2 *= 2;
	}

  ntemp = pp+1+Y_NOPTS;
  if(!(temparr1 = (int *)malloc(ntemp*sizeof(int)))) {
    printf("[pspacey.c:1] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  for(i=0;i<pp+1;i++) temparr1[i] = rowdista[i];
  for(j=0;j<Y_NOPTS;j++) temparr1[j+i] = ioptions[j];

  if(!myid) {
    if(!(temparr2 = (int *)malloc(ntemp*sizeof(int)))) {
      printf("[pspacey.c:2] memory allocation error\n");
      MPI_Abort(comm,0);
    }
  }
  MPI_Reduce(temparr1,temparr2,ntemp,MPI_INT,MPI_BXOR,0,comm);
  if(!myid) {
    for(i=0;i<ntemp;i++) {
      if(temparr2[i]) {
        printf("PSPACEY: Global parameters not matching on all processors!\n");
      if(i<pp+1)
        printf("Check rowdista[%d] !\n",i);
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

  lgblk = 0;
  while((blk = 1<<lgblk) < ioptions[0]) lgblk++;
  ioptions[0] = blk;
  checksymm = ioptions[1];
  if(checksymm!=1) checksymm = 0;
  sortinds = ioptions[2];
  if(sortinds!=1) sortinds = 0;

  mynnodes = rowdista[myid+1]-rowdista[myid];

  memory_pscomm = 0.0;

  i = 11*pp+2*N+10;
  IWN = 11*pp;
  memory_pscomm += (double)i;
  if(!(wrkint = (int *)malloc(i*sizeof(int)))) {
    printf("[pspacey.c:3] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  memory_pscomm += (double)(3*N);
  if(!(tptrs = (int *)malloc(3*N*sizeof(int)))) {
    printf("[pspacey.c:4] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  memory_pscomm += (double)N;
  if(!(tinds = (int *)malloc(N*sizeof(int)))) {
    printf("[pspacey.c:5] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  memory_pscomm += (double)(3*N);
  if(!(lptrs = (int *)malloc(3*N*sizeof(int)))) {
    printf("[pspacey.c:6] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  i = 5*(2*pp-1);
  memory_pscomm += (double)i;
  if(!(ranmasks = (int *)malloc(i*sizeof(int)))) {
    printf("[pspacey.c:7] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  taivsize = aptrs[2*(mynnodes-1)]+aptrs[2*(mynnodes-1)+1]-1;
  if(!(tainds = (int *)malloc(taivsize*sizeof(int)))) {
    printf("[pspacey.c:8] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  if(!(whichsnode = (int *)malloc(mynnodes*sizeof(int)))) {
    printf("[pspacey.c:9] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  premovea(&N,&dd,&pp,&lgblk,&myid,&mynnodes,rowdista,order,sizes,
       &mybeginleaf,&pasize,aptrs,ainds,tainds,lptrs,tinds,tptrs,
       wrkint,ranmasks,whichsnode,&maxnzpercol,&checksymm,&sortinds,
       pmcomm);

  if(!(painds = (int *)malloc(pasize*sizeof(int)))) {
    printf("[pspacey.c:10] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  memory_pscomm += (double)(2*N);
  if(!(paptrs = (int *)malloc(2*N*sizeof(int)))) {
    printf("[pspacey.c:11] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  moveai(&N,&dd,&pp,&lgblk,&myid,&mynnodes,order,paptrs,&(lptrs[N]),
       painds,aptrs,tainds,wrkint,ranmasks,whichsnode,&nsend,&nrecvsizs,
       pmcomm);

  free(whichsnode);
  free(tainds);

  ygentree(&N,paptrs,painds,&mybeginleaf,sizes,lptrs,tptrs,tinds,
        &dd,&myid,&pp,&(wrkint[IWN]),&pasize,&(lptrs[N]),&nrecvsizs,pmcomm);

  memory_pscomm += (double)pasize;
  if(!(painds = (int *)realloc(painds,pasize*sizeof(int)))) {
    printf("%d: memory reallocation error\n",myid);
    MPI_Abort(comm,0);
  }

  memory_pscomm += (double)N;
  if(!(cinfo = (int *)malloc(N*sizeof(int)))) {
    printf("[pspacey.c:12] memory allocation error\n");
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
    printf("[pspacey.c:13] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  if(!(supinds = (int *)malloc(supindsize*sizeof(int)))) {
    printf("[pspacey.c:14] memory allocation error\n");
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
    printf("[pspacey.c:15] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  if(!(lbal = (double *)malloc(4*pp*sizeof(double)))) {
    printf("[pspacey.c:16] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  for(i=0;i<2*pp-1;i++) lbal[i] = 0.0;

  info = 1;
  while(info) {
    leaflindsize += 2*nnzleafa;
    lindsize = tlplindsize+leaflindsize;

    if(!(linds = (int *)realloc(linds,lindsize*sizeof(int)))) {
      printf("[pspacey.c:17] memory reallocation error\n");
      MPI_Abort(comm,0);
    }

    parsymb(&root,paptrs,painds,tptrs,tinds,lptrs,linds,&lindsptr,
          &lvalsptr,&nnz,lbal,sizes,sup,supinds,&supptr,&supindsptr,
          &iu,&leaflindsize,&dd,&N,&myid,&lgblk,&myidv,&myidh,
          &(wrkint[IWN]),&(wrkint[IWN+N]),&info,cinfo,pmcomm);

  }

  MPI_Allreduce(&pasize,&nnz_ga,1,MPI_INT,MPI_SUM,comm);

  ilbd = 2*pp;
  MPI_Reduce(lbal,&(lbal[ilbd]),2*pp-1,MPI_DOUBLE,MPI_SUM,0,comm);

  if(!myid) {

    doptions[2] = (double)nnz_ga;
    doptions[3] = (double)nnz;

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

  free(lbal);

  ns = (supptr-1)/4;
  lspace = lvalsptr-1;
  lindsize = lindsptr-1;
  supsize = supptr-1;
  supindsize = supindsptr-1+iu;

  memory_pscomm += (double)lindsize;
  if(!(linds = (int *)realloc(linds,lindsize*sizeof(int)))) {
    printf("[pspacey.c:18] memory reallocation error\n");
    MPI_Abort(comm,0);
  }
  memory_pscomm += (double)supsize;
  if(!(sup = (int *)realloc(sup,supsize*sizeof(int)))) {
    printf("[pspacey.c:19] memory reallocation error\n");
    MPI_Abort(comm,0);
  }
  memory_pscomm += (double)supindsize;
  if(!(supinds = (int *)realloc(supinds,supindsize*sizeof(int)))) {
    printf("[pspacey.c:20] memory reallocation error\n");
    MPI_Abort(comm,0);
  }
  memory_pscomm += (double)(supindsptr-1);
  if(!(tsind = (int *)malloc((supindsptr-1)*sizeof(int)))) {
    printf("[pspacey.c:21] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  for(i=0;i<supindsptr-1;i++) tsind[i] = supinds[i];

  memory_pscomm += (double)(2*N);
  if(!(iptrs = (int *)malloc(2*N*sizeof(int)))) {
    printf("[pspacey.c:22] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  lcsize = ((nnz_ga+N)*3)/dd;
  ilcsize = lcsize;

  memory_pscomm += (double)ilcsize;
  if(!(lc = (int *)malloc(ilcsize*sizeof(int)))) {
    printf("[pspacey.c:23] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  memory_pscomm += (double)IFOPTS_SIZE;
  if(!(ifopts = (int *)malloc(IFOPTS_SIZE*sizeof(int)))) {
    printf("[pspacey.c:24] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  if(!(wa1 = (int *)malloc(sizes[myid]*sizeof(int)))) {
    printf("[pspacey.c:25] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  eparfact1(&N,paptrs,painds,lptrs,linds,tptrs,tinds,sup,&(wrkint[IWN]),
         &(wrkint[IWN+N]),&root,&dd,&lgblk,&myid,supinds,&(wrkint[IWN+N+N]),
         &wsize0,&wsize1,&ibuflen,&dbuflen,&iwspace,&node,&stakptr,&nptr,
         lc,iptrs,&lcsize,&wsolvesize,wa1,pmcomm);

  free(wa1);

  if(lcsize<ilcsize) {
    memory_pscomm -= (double)(ilcsize-lcsize);
    if(!lcsize) lcsize = 1;
    if(!(lc = (int *)realloc(lc,sizeof(int)*lcsize))) {
      printf("[pspacey.c:26] memory reallocation error\n");
      MPI_Abort(comm,0);
    }
  }

  memory_pscomm += (double)mynnodes;
  if(!(Iorder = (int *)malloc(mynnodes*sizeof(int)))) {
    printf("[pspacey.c:27] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  memory_pscomm += (double)(pp+1);
  if(!(Irowdist = (int *)malloc((pp+1)*sizeof(int)))) {
    printf("[pspacey.c:28] memory allocation error\n");
    MPI_Abort(comm,0);
  }
  for(i=0;i<mynnodes;i++) Iorder[i] = order[i];
  for(i=0;i<pp+1;i++) Irowdist[i] = rowdista[i];

  parfact_memreq = (double)(pasize*sizeof(double));      /* pavals */

  moveav_memreq = (double)(nsend*sizeof(double));         /* sendvals */
  i = aptrs[2*(mynnodes-1)]+aptrs[2*(mynnodes-1)+1]-1;
  moveav_memreq += (double)(2*i*sizeof(int));             /* tainds */
  moveav_memreq += (double)(mynnodes*sizeof(int));        /* whichsnode */

  dpspacen_memreq = (double)(lspace*sizeof(double));      /* lvals */
  dpspacen_memreq += (double)(dbuflen*2*sizeof(double));  /* dbuf_s */
  dpspacen_memreq += (double)(ibuflen*2*sizeof(int));     /* ibuf_s */
  dpspacen_memreq += (double)(iwspace*sizeof(double));    /* wmem */
  dpspacen_memreq += (double)(N*sizeof(int));             /* locinds */

  parfact_memreq += (moveav_memreq > dpspacen_memreq) ? 
                     moveav_memreq : dpspacen_memreq;

  doptions[0] = memory_pscomm*(double)(sizeof(int));
  doptions[1] = parfact_memreq;

  doptions[12] = (double)(lspace*sizeof(double));
  doptions[13] = (double)(dbuflen*2*sizeof(double));
  doptions[14] = (double)(ibuflen*2*sizeof(int));
  doptions[15] = (double)(iwspace*sizeof(double));

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
  ifopts[11] = lspace;
  ifopts[12] = wsize0;
  ifopts[13] = wsize1;
  ifopts[14] = ibuflen;
  ifopts[15] = dbuflen;
  ifopts[16] = iwspace;
  ifopts[17] = node;
  ifopts[18] = stakptr;
  ifopts[19] = nptr;
  ifopts[20] = IWN;
  ifopts[21] = maxnzpercol;
  ifopts[22] = pasize;

  if(!(ap = (PTRS *)malloc(sizeof(PTRS)))) {
    printf("[pspacey.c:29] memory allocation error\n");
    MPI_Abort(comm,0);
  }

  ap->Ylink = NULL;
  ap->howManyNfacts = 0;

  ap->rowdist = Irowdist;
  ap->order = Iorder;
  ap->ranmasks = ranmasks;

  ap->wrkint = wrkint;
  ap->cinfo = cinfo;
  ap->painds = painds;
  ap->paptrs = paptrs;
  ap->supinds = supinds;

  ap->lptrs = lptrs;
  ap->linds = linds;
  ap->tptrs = tptrs;
  ap->tinds = tinds;
  ap->sup = sup;
  ap->tsind = tsind;
  ap->lc = lc;
  ap->iptrs = iptrs;
  ap->ifopts = ifopts;
  ap->lvals = NULL;
  ap->pmcomm = pmcomm;
  
  *pspcomm = check_handle((long)ap,0);
  if( *pspcomm < 0 ) {
    if(!myid) {
      printf("Unable to create new PSPASES communicator, all %d slots\n",
	     MAX_OPEN_PSPCOMMS);
      printf("are taken up. Please free some of those communicators using\n");
      printf("PSPACEC.\n");
    }
    MPI_Barrier(comm);
    MPI_Abort(comm,0);
  }

}
