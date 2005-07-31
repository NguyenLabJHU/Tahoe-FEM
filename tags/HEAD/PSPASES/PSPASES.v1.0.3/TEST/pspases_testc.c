
/***** This code reads in a matrix file in .bin format, and tests PSPASES */
/***** functionality for various parameters and solution paths.           */
/***** 					        written by: Mahesh Joshi  */
/*****   $Id: pspases_testc.c,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <cpspases.h>

#ifdef INT8
typedef short int4;
#else
typedef int int4;
#endif

#define DPSPACEN_REPS	3

void main(int argc, char *argv[]) 
{

FILE *fp;
int4 *aptrs4,*ainds4,*taptrs,*tainds;
int *aptrs,*ainds,*rowdist,*rowdistbx,*aisizes;
double *avals,*b,*x,*tavals,*tb;
const int itag=1;
int pp,myid,dd,blk,nrhs,nbuf[10],msglen,Npp,Npp1,asize,lsize,i,j,k,l,size,proc;
int br,mynb,ldb,ldx,mxasize,maxm,checksymm,sortinds,serialorder,mynnodes;
int4 N;
int nown,ioasize,m,options[16],runopt,beginOpt,endOpt;
long pspcomm,pspcommY,pspcommF,pspcommN[DPSPACEN_REPS],jseed;
MPI_Status mpistat;
double doptions[16],emax;
MPI_Comm comm;
int *order,*sizes,clean_option;
double time0,Ftime,Otime,Ytime,Ntime,Ttime,Ltime;
double *memPM,*memFM,maxPM,minPM,sumPM,maxFM,minFM,sumFM,maxTM,minTM;

MPI_Init(&argc,&argv);

MPI_Comm_dup(MPI_COMM_WORLD,&comm);
MPI_Comm_size(comm,&pp);
MPI_Comm_rank(comm,&myid);

if(argc<2) {
  if(!myid) {
   printf("Usage: %s <filename> [<blk> [<nrhs> [<br> [<symm> [<sort> [<oopt> [<ropt>] ]]]]]]\n",
             argv[0]);
   printf("   <filename> : contains Full Symmetric A in binary SPASES format\n");
   printf("\n Optional Parameters --\n");
   printf("   <blk>  : used in block cyclic distribution of factor. (default=64)\n");
   printf("   <nrhs> : number of right hand sides to solve. (default=1)\n");
   printf("   <br>   : blocking factor for right hand side. (default=nrhs)\n");
   printf("   <symm> : if 1, checks symmetry of input matrix. (default=0)\n");
   printf("   <sort> : if 1, sorts the indices in ainds. (default=0)\n");
   printf("   <oopt> : if 1, uses serial ordering. (default=0 for ParMetis)\n");
   printf("   <ropt> : if 0, DPSPACEF+DPSPACET.\n");
   printf("            if 1, PSPACEO+PSPACEY+DPSPACEN+DPSPACET.\n");
   printf("            if 2, PSPACEO+PSPACEY.\n");
   printf("            (default: 2 -> 1 -> 0)\n");
  }
  MPI_Barrier(comm);
  MPI_Abort(comm,0);
}

if(!myid) {
  blk = 0;
  nrhs = 0;
  br = 0;
  checksymm = 0;
  sortinds = 0;
  serialorder = 0;
  runopt = -1;
  if(argc > 2) {
    blk = atoi(argv[2]);
    if(argc > 3) {
      nrhs = atoi(argv[3]);
      if(argc > 4) {
        br = atoi(argv[4]);
        if(argc > 5) {
          checksymm = atoi(argv[5]); 
          if(argc > 6) {
            sortinds = atoi(argv[6]);
            if(argc > 7) {
              serialorder = atoi(argv[7]);
	      if(argc > 8) {
                runopt = atoi(argv[8]);
	      }
            }
          }
        } 
      } 
    }
  }
  if(!blk) blk = 64;
  if(!nrhs) nrhs = 1;
  if(!br) br = nrhs;
  if(runopt < 0 || runopt > 2) runopt = -1;

  fp = fopen(argv[1],"r");
  if(!fp) {
    printf("Unable to open %s\n",argv[1]);
    MPI_Abort(comm,0);
  }

  fread(&N,sizeof(int4),1,fp);
  nbuf[0] = N;
  nbuf[1] = nrhs;
  nbuf[2] = blk;
  nbuf[3] = br;
  nbuf[4] = checksymm;
  nbuf[5] = sortinds;
  nbuf[6] = serialorder;
  nbuf[7] = runopt;

  printf("\n\n* Testing %s\n",argv[1]);
  printf("Dimension = %d, #processors = %d\n",N,pp);
  printf("Parameters: blk  = %d, nrhs = %d, br = %d\n",blk,nrhs,br);
  printf("            symm = %d, sort = %d, serialorder = %d\n\n",
	   checksymm,sortinds,serialorder);
}
MPI_Barrier(comm);

MPI_Bcast(nbuf,10,MPI_INT,0,comm);

if(myid) {
  N = nbuf[0];
  nrhs = nbuf[1];
  blk = nbuf[2];
  br = nbuf[3];
  checksymm = nbuf[4];
  sortinds = nbuf[5];
  serialorder = nbuf[6];
  runopt = nbuf[7];
}

if(!(rowdist=(int *)malloc((pp+1)*sizeof(int)))) {
  printf("memory allocation failure\n");
  MPI_Abort(comm,0);
}

Npp = N/pp;
Npp1 = (N+pp-1)/pp;
if(N - Npp1*(pp-1) > 0 && (Npp1*pp - N < N - Npp*pp)) Npp = Npp1;

for(i=0;i<pp;i++) rowdist[i] = i*Npp;
rowdist[pp] = N;

m = rowdist[myid+1]-rowdist[myid];
maxm = Npp > N-(pp-1)*Npp ? Npp : N-(pp-1)*Npp;

if(!(aptrs4=(int4 *)malloc(2*m*sizeof(int4)))) {
  printf("memory allocation failure\n");
  MPI_Abort(comm,0);
}

if(!(aisizes = (int *)malloc(sizeof(int)*(pp+1)))) {
  printf("memory allocation failure\n");
  MPI_Abort(comm,0);
}

/* read aptrs */
if(!myid) {
  fread(aptrs4,sizeof(int4),2*m,fp);
  aisizes[0] = aptrs4[2*(m-1)]+aptrs4[2*m-1]-aptrs4[0];
  mxasize = aisizes[0];

  if(!(taptrs=(int4 *)malloc(sizeof(int4)*maxm*2))) {
    printf("memory allocation failure\n");
    MPI_Abort(comm,0);
  }
  for(proc=1;proc<pp;proc++) {
    m = rowdist[proc+1]-rowdist[proc];
    fread(taptrs,sizeof(int4),2*m,fp);
    aisizes[proc] = taptrs[2*(m-1)]+taptrs[2*m-1]-taptrs[0];
    mxasize = (mxasize > aisizes[proc]) ? mxasize : aisizes[proc];
    for(i=2;i<2*m;i+=2) taptrs[i] -= taptrs[0]-1;
    taptrs[0] = 1;
    MPI_Send(taptrs,sizeof(int4)*2*m,MPI_BYTE,proc,itag,comm);
  }
  free(taptrs);
}
else {
  MPI_Recv(aptrs4,sizeof(int4)*2*m,MPI_BYTE,0,itag,comm,&mpistat);
}

MPI_Bcast(aisizes,sizeof(int)*(pp+1),MPI_BYTE,0,comm);
asize = aisizes[myid];

if(!(ainds4=(int4 *)malloc(asize*sizeof(int4)))) {
  printf("memory allocation failure\n");
  MPI_Abort(comm,0);
}

/* read ainds */
if(!myid) {
  fread(ainds4,sizeof(int4),asize,fp);
  if(!(tainds=(int4 *)malloc(sizeof(int4)*mxasize))) {
    printf("memory allocation failure\n");
    MPI_Abort(comm,0);
  }
  for(proc=1;proc<pp;proc++) {
    fread(tainds,sizeof(int4),aisizes[proc],fp);
    MPI_Send(tainds,sizeof(int4)*aisizes[proc],MPI_BYTE,proc,itag,comm);
  }
  free(tainds);
}
else {
  MPI_Recv(ainds4,sizeof(int4)*asize,MPI_BYTE,0,itag,comm,&mpistat);
}

/* read avals */
if(!(avals=(double *)malloc(asize*sizeof(double)))) {
  printf("memory allocation failure\n");
  MPI_Abort(comm,0);
}
if(!myid) {
  fread(avals,sizeof(double),asize,fp);
  if(!(tavals=(double *)malloc(sizeof(double)*mxasize))) {
    printf("memory allocation failure\n");
    MPI_Abort(comm,0);
  }
  for(proc=1;proc<pp;proc++) {
    fread(tavals,sizeof(double),aisizes[proc],fp);
    MPI_Send(tavals,sizeof(double)*aisizes[proc],MPI_BYTE,proc,itag,comm);
  }
  free(tavals);
  fclose(fp);
}
else {
  MPI_Recv(avals,sizeof(double)*asize,MPI_BYTE,0,itag,comm,&mpistat);
}

#ifdef INT8
m = rowdist[myid+1]-rowdist[myid];
if(!(aptrs=(int *)malloc(sizeof(int)*2*m))) {
  printf("memory allocation failure\n");
  MPI_Abort(comm,0);
}

if(!(ainds=(int *)malloc(sizeof(int)*asize))) {
  printf("memory allocation failure\n");
  MPI_Abort(comm,0);
}

for(i=0;i<2*m;i++) {
  aptrs[i] = aptrs4[i];
}
for(i=0;i<asize;i++) {
  ainds[i] = ainds4[i];
}
free(ainds4);
free(aptrs4);
#else
aptrs = aptrs4;
ainds = ainds4;
#endif

options[0] = blk;        /* block size */
options[1] = checksymm;  /* check symmetry : 1: yes, 0: no (Default:0) */
options[2] = sortinds;   /* sort indices : 1: yes, 0: no (Default:0) */
options[3] = serialorder;/* use serial ordering : 1: yes, 0: no (Default:0) */

if(runopt == -1) {
  beginOpt = 2;
  endOpt = 0;
} else {
  beginOpt = runopt;
  endOpt = runopt;
}

for(runopt=beginOpt;runopt>=endOpt;runopt--) {

MPI_Barrier(comm);
if(!myid) {
  switch (runopt) {

  case 0 :
    printf("\n-------> Testing DPSPACEF+DPSPACET\n\n");
    break;
  case 1 :
    printf("\n-------> Testing PSPACEO+PSPACEY+DPSPACEN+DSPACET\n\n");
    break;
  case 2 :
    printf("\n-------> Testing PSPACEO+PSPACEY\n\n");
    break;
  }
}

if(runopt == 0) {

  MPI_Barrier(comm);
  if(!myid) printf("calling DPSPACEF (Ordering+sYmbolic+Numerical)..\n");

  MPI_Barrier(comm);
  time0 = MPI_Wtime();

  DPSPACEF(rowdist,aptrs,ainds,avals,options,doptions,&pspcommF,&comm);

  MPI_Barrier(comm);
  Ftime = MPI_Wtime()-time0;

} else {

  mynnodes = rowdist[myid+1]-rowdist[myid];
  if(!(order = (int *)malloc(sizeof(int)*mynnodes))) {
    printf("memory allocation failure\n");
    MPI_Abort(comm,0);
  }
  if(!(sizes = (int *)malloc(sizeof(int)*(2*pp)))) {
    printf("memory allocation failure\n");
    MPI_Abort(comm,0);
  }

  MPI_Barrier(comm);
  if(!myid) printf("calling  PSPACEO (Compute fill-reducing Ordering)..\n");

  MPI_Barrier(comm);
  time0 = MPI_Wtime();

  PSPACEO(rowdist,aptrs,ainds,order,sizes,options,&comm);

  MPI_Barrier(comm);
  Otime = MPI_Wtime()-time0;

  MPI_Barrier(comm);
  if(!myid) printf("calling  PSPACEY (sYmbolic Factorization)..\n");

  MPI_Barrier(comm);
  time0 = MPI_Wtime();

  PSPACEY(rowdist,aptrs,ainds,order,sizes,options,doptions,&pspcommY,&comm);

  MPI_Barrier(comm);
  Ytime = MPI_Wtime()-time0;

  free(order);
  free(sizes);

}

if(!myid) {
  if(!(memPM = (double *)malloc(pp*sizeof(double)))) {
    printf("memory allocation failure\n");
    MPI_Abort(comm,0);
  }
  if(!(memFM = (double *)malloc(pp*sizeof(double)))) {
    printf("memory allocation failure\n");
    MPI_Abort(comm,0);
  }
}

doptions[0] /= (double)(1024*1024);
doptions[1] /= (double)(1024*1024);
MPI_Gather(&(doptions[0]),1,MPI_DOUBLE,memPM,1,MPI_DOUBLE,0,comm);
MPI_Gather(&(doptions[1]),1,MPI_DOUBLE,memFM,1,MPI_DOUBLE,0,comm);

MPI_Barrier(comm);
if(!myid) {
  minPM = memPM[0];
  maxPM = memPM[0];
  sumPM = memPM[0];
  for(i=1;i<pp;i++) {
    if(memPM[i] < minPM) minPM = memPM[i];
    if(memPM[i] > maxPM) maxPM = memPM[i];
    sumPM += memPM[i];
  }

  printf("\n");

  printf("                                     Max       Min       Sum\n");
  printf("  Memory Consumed by PSPCOMM : %9.3lf %9.3lf %9.3lf MB\n",
     maxPM,minPM,sumPM);

  if(runopt != 0) {
    minFM = memFM[0];
    maxFM = memFM[0];
    sumFM = memFM[0];
    minTM = memFM[0]+memPM[0];
    maxTM = memFM[0]+memPM[0];

    for(i=1;i<pp;i++) {
      if(memFM[i] < minFM) minFM = memFM[i];
      if(memFM[i] > maxFM) maxFM = memFM[i];
      sumFM += memFM[i];
      if(memFM[i]+memPM[i] < minTM) minTM = memFM[i]+memPM[i];
      if(memFM[i]+memPM[i] > maxTM) maxTM = memFM[i]+memPM[i];
    }

    printf("  More Factor Memory Needed  : %9.3lf %9.3lf %9.3lf MB\n",
	maxFM,minFM,sumFM);

    printf("  Total Factor Memory Needed : %9.3lf %9.3lf %9.3lf MB\n",
	maxTM,minTM,sumFM+sumPM);
  }
  printf("\n");

  free(memPM);
  free(memFM);

  printf("  NNZ_Lower_Triangle_of_A = %10.4le  NNZ_L          = %10.4le\n",
          doptions[2],doptions[3]);
  printf("  Tree_Opcount_Imbalance  = %10.3lf  Factor_Opcount = %10.4le\n",
          doptions[5],doptions[4]);
  printf("\n");

  if(runopt == 0) {
    Ltime = Ftime;
    printf("  DPSPACEF Time (Ftime)      = %10.3lf sec\n",Ftime);
    printf("  Factor_Opcount / Ftime     = %10.3lf MFLOPS\n",
		doptions[4]/Ftime*1.e-6);
  } else {
    Ltime = Otime+Ytime;
    printf("  Order Time                 = %10.3lf sec\n",Otime);
    printf("  sYmbolic Time              = %10.3lf sec\n",Ytime);
  }

  printf("\n");
}

if(runopt == 1) {
  MPI_Barrier(comm);
  if(!myid) printf("calling DPSPACEN (Numerical Factorization)..\n");

  MPI_Barrier(comm);
  time0 = MPI_Wtime();

  for(i=0;i<DPSPACEN_REPS;i++) {
    DPSPACEN(rowdist,aptrs,ainds,avals,&pspcommY,&(pspcommN[i]),&comm);
  }

  MPI_Barrier(comm);
  Ntime = (MPI_Wtime()-time0)/(double)DPSPACEN_REPS;

  if(!myid) {
    Ltime += Ntime;
    printf("\n");
    printf("  DPSPACEN Time              = %10.3lf sec\n",Ntime);
    printf("  Numerical Factor Perf      = %10.3lf MFLOPS\n",
		doptions[4]/Ntime*1.e-6);
    printf("\n");
  }

}

if(runopt == 1) {

  MPI_Barrier(comm);
  if(!myid) printf("calling PSPACEC (option=1) on pspcommY..\n");
  clean_option = 1;
  PSPACEC(&pspcommY,&clean_option);

  for(i=1; i<DPSPACEN_REPS; i++) {
    MPI_Barrier(comm);
    if(!myid) printf("calling PSPACEC (option=0) on pspcommN[%d]..\n",i);
    clean_option = 0;
    PSPACEC(&(pspcommN[i]),&clean_option);
  }

  if(!myid) printf("\n");

} else if( runopt == 2) {

  MPI_Barrier(comm);
  if(!myid) printf("calling PSPACEC (option=0) on pspcommY..\n");
  clean_option = 0;
  PSPACEC(&pspcommY,&clean_option);

}

if(runopt != 2) {

  if(!(rowdistbx = (int *)malloc(sizeof(int)*(pp+1)))) {
    printf("memory allocation failure\n");
    MPI_Abort(comm,0);
  }

  j = 4;
  for(i=0;i<pp;i++) rowdistbx[i] = i*(Npp-Npp/j);
  rowdistbx[pp] = N;
  m = N-(pp-1)*(Npp-Npp/j);
  maxm = m > (Npp-Npp/j) ? m : (Npp-Npp/j);

  mynb = rowdistbx[myid+1]-rowdistbx[myid];

  ldb = mynb;
  if(!(b = (double *)malloc(sizeof(double)*ldb*nrhs))) {
    printf("memory allocation failure\n");
    MPI_Abort(comm,0);
  }

  if(!myid) {
    m = rowdistbx[1];
    jseed = 13 + 69805;
    srand48(jseed);
    for(j=0;j<nrhs;j++) {
      k = j*ldb;
      for(i=k;i<k+m;b[i++] = 2.0*(drand48()-0.5));
    }
    if(!(tb = (double *)malloc(maxm*sizeof(double)))) {
      printf("memory allocation failure\n");
      MPI_Abort(comm,0);
    }
    for(proc=1;proc<pp;proc++) {
      m = rowdistbx[proc+1]-rowdistbx[proc];
      for(j=0;j<nrhs;j++) {
        for(i=0;i<m;tb[i++]=2.0*(drand48()-0.5));
          MPI_Send(tb,m*sizeof(double),MPI_BYTE,proc,itag,comm);
      }
    }
    free(tb);
  } else {
    for(j=0;j<nrhs;j++) {
      k = j*ldb;
      MPI_Recv(&(b[k]),mynb*sizeof(double),MPI_BYTE,0,itag,comm,&mpistat);
    }
  }

  ldx = mynb;
  if(!(x = (double *)malloc(sizeof(double)*ldx*nrhs))) {
    printf("memory allocation failure\n");
    MPI_Abort(comm,0);
  }

  MPI_Barrier(comm);
  if(!myid) printf("calling DPSPACET (Triangular Systems Solution)..\n");

  if( runopt == 0 ) {
    pspcomm = pspcommF;
  } else {
    pspcomm = pspcommN[0];
  }

  options[0] = br; /* blocking factor on nrhs. */

  MPI_Barrier(comm);
  time0 = MPI_Wtime();

  DPSPACET(rowdistbx,&nrhs,b,&ldb,x,&ldx,options,&pspcomm,&comm);

  MPI_Barrier(comm);
  Ttime = MPI_Wtime()-time0;

  if(!myid) {
    printf("\n");
    printf("  DPSPACET Time              = %10.3lf sec\n",Ttime);
    printf("  Triangular Solution Perf   = %10.3lf MFLOPS\n",
  	    doptions[3]*(double)nrhs*4.e-6/Ttime);
    printf("\n");
    Ltime += Ttime;
    printf("  Total Solver Time          = %10.3lf sec\n",Ltime);
    printf("\n");
  }

  if( runopt == 0 ) {

    MPI_Barrier(comm);
    if(!myid) printf("calling PSPACEC (option=0) on pspcommF..\n");
    clean_option = 0;
    PSPACEC(&pspcommF,&clean_option);

  } else {

    MPI_Barrier(comm);
    if(!myid) printf("calling PSPACEC (option=0) on pspcommN[0]..\n");
    clean_option = 0;
    PSPACEC(&(pspcommN[0]),&clean_option);

    MPI_Barrier(comm);
    if(!myid) printf("calling PSPACEC (option=0) on pspcommY..\n");
    clean_option = 0;
    PSPACEC(&pspcommY,&clean_option);
  }

  MPI_Barrier(comm);
  if(!myid) printf("all stages of PSPASES are done! checking B-AX ..\n");

  CHECKB_AX(rowdist,aptrs,ainds,avals,rowdistbx,&nrhs,b,&ldb,x,&ldx,
            &emax,&comm);

  MPI_Barrier(comm);
  if(!myid) printf("max |B - AX| = %20.14le\n",emax);

  free(b);
  free(x);
  free(rowdistbx);
}

}

free(aptrs);
free(ainds);
free(avals);
free(rowdist);

MPI_Finalize();
}
