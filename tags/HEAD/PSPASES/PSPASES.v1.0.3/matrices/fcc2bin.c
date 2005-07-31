
/**** This code converts a matrix file  in .fcc format to */
/**** .bin format file. */
/****   $Id: fcc2bin.c,v 1.1.1.1 2004-10-07 16:05:25 paklein Exp $ */

#include <stdio.h>
#include <stdlib.h>

#ifdef INT8
typedef short int4;
#else
typedef int int4;
#endif

#define DEFAULTNZPERCOL 100

void main(int argc,char *argv[]) {

FILE *fpin,*fpout;
int4 N;
int4 *aptrs,*ainds,*pi;
double *avals,*pd;
int asize,pai,i,rownnz,j,pap;
int NZPERCOL;

if(argc < 3) {
  printf("Usage: %s <fcc_input_file> <bin_output_file> <NZPERCOL> \n",argv[0]);
  printf("       <NZPERCOL> : avg number of nonzeros per column or row\n"); 
  printf("                    (default = %d).\n",DEFAULTNZPERCOL);
  exit(0);
}

if(!(fpin = fopen(argv[1],"r"))) {
 printf("Unable to open file %s\n",argv[1]);
 exit(0);
}

if(argc > 3) {
  NZPERCOL = atoi(argv[3]);
} else {
  printf("Using default value (=%d) for NZPERCOL.\n",DEFAULTNZPERCOL);
  printf("    You may specify different value on command line.\n");
  printf("    Usage: %s <fcc_input_file> <bin_output_file> <NZPERCOL> \n",argv[0]);
  printf("           <NZPERCOL> : avg number of nonzeros per column or row\n"); 
  NZPERCOL = DEFAULTNZPERCOL;
}

printf("\nUsing integer of size %d bytes for encoding in binary.\n\n",sizeof(int4));

fscanf(fpin,"%d",&N);
asize = N*NZPERCOL;

if(!(aptrs = (int4 *)malloc(2*N*sizeof(int4)))) {
  printf("memory allocation failure\n");
  exit(0);
}

if(!(ainds = (int4 *)malloc(asize*sizeof(int4)))) {
  printf("memory allocation failure\n");
  exit(0);
}

if(!(avals = (double *)malloc(asize*sizeof(double)))) {
  printf("memory allocation failure\n");
  exit(0);
}

pai = 0;
pap = 0;
for(i=0;i<N;i++) {
  fscanf(fpin,"%d",&rownnz);
  if(rownnz+pai >= asize) {
    printf("Insufficient memory allocated for ainds/vals (at row %d/%d)\n",i,N);
    printf("Try to increase NZPERCOL value.\n");
    exit(0);
  }
  for(j=pai;j<pai+rownnz;j++) {
    fscanf(fpin,"%lf",avals+j);
    fscanf(fpin,"%d",ainds+j);
  }
  aptrs[pap] = pai+1;
  aptrs[pap+1] = rownnz;
  pap += 2;
  pai += rownnz;
}
asize = pai;

fclose(fpin);

if(!(fpout = fopen(argv[2],"w"))) {
 printf("Unable to open file %s\n",argv[2]);
 exit(0);
}

fwrite(&N,sizeof(int4),1,fpout);
fwrite(aptrs,sizeof(int4),2*N,fpout);
fwrite(ainds,sizeof(int4),asize,fpout);
fwrite(avals,sizeof(double),asize,fpout);

}
