/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * wavefrontK.c 
 *
 * This file contains code for the initial directed diffusion at the coarsest
 * graph
 *
 * Started 5/19/97, Kirk, George
 *
 * $Id: diffutil.c,v 1.1.1.1 2004-10-07 16:05:25 paklein Exp $
 *
 */

#include <parmetis.h>

/*************************************************************************
*  This function computes the edge-cut of a serial graph.
**************************************************************************/
int ComputeSerialEdgeCut(GraphType *graph)
{
  int i, j;
  int cut = 0;

  for (i=0; i<graph->nvtxs; i++) {
    for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
      if (graph->where[i] != graph->where[graph->adjncy[j]])
        cut += graph->adjwgt[j];
  }
  graph->mincut = cut/2;

  return graph->mincut;
}

/*************************************************************************
*  This function computes the TotalV of a serial graph.
**************************************************************************/
int ComputeSerialTotalV(GraphType *graph, idxtype *home)
{
  int i;
  int totalv = 0;

  for (i=0; i<graph->nvtxs; i++) 
    if (graph->where[i] != home[i])
      totalv += graph->vwgt[i];

  return totalv;
}

/*************************************************************************
*  This function computes the TotalV of a distributed graph.
**************************************************************************/
int ComputeParallelTotalV(CtrlType *ctrl, GraphType *graph, idxtype *home)
{
  int i;
  int totalv = 0;

  for (i=0; i<graph->nvtxs; i++)
    if (graph->where[i] != home[i])
      totalv += graph->vwgt[i];

  return GlobalSESum(ctrl, totalv);
}


/*************************************************************************
*  This function computes the balance of a distributed graph.
**************************************************************************/
float ComputeParallelBalance(CtrlType *ctrl, GraphType *graph)
{
  return 1.0*ctrl->nparts*graph->gpwgts[idxamax(ctrl->nparts, graph->gpwgts)]/
	(1.0*idxsum(ctrl->nparts, graph->gpwgts));
}


/*************************************************************************
*  This function compares the values of two referenced integers.
**************************************************************************/
int GreaterThan(const void *first, const void *second)
{
  if ((*(int*)(first)) > (*(int*)(second)))
    return 1;
  
  if ((*(int*)(first)) < (*(int*)(second)))
    return -1;

  return 0;
}



/*************************************************************************
* This function implements the CG solver used during the directed diffusion
**************************************************************************/
void ConjGrad(int n, idxtype *rowptr, idxtype *colind, float *values, float *b, float *x, float tol, float *workspace)
{
  int i, j, k;
  float *p, *r, *q, *z, *M; 
  float alpha, beta, rho, rho_1, error, bnrm2;

/*
  printf("%d %d %d %d %d | %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", rowptr[0], rowptr[1],
     rowptr[2], rowptr[3], rowptr[4], colind[0], colind[1], colind[2], colind[3],
     colind[4], colind[5], colind[6], colind[7], colind[8], colind[9], colind[10],
     colind[11], colind[12], colind[13], colind[14], colind[15]);
  printf("%2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f\n", 
     values[0], values[1], values[2], values[3], values[4], values[5], values[6], 
     values[7], values[8], values[9], values[10], values[11], values[12], values[13], 
     values[14], values[15]);
*/

  /* Initial Setup */   
  p = workspace;
  r = workspace + n;
  q = workspace + 2*n;
  z = workspace + 3*n;
  M = workspace + 4*n;

  for (i=0; i<n; i++) {
    x[i] = 0.0;
    if (values[rowptr[i]] != 0.0)
      M[i] = 1.0/values[rowptr[i]];
    else
      M[i] = 0.0;
  }

  /* r = b - Ax */
  mvMult(n, rowptr, colind, values, x, r);
  for (i=0; i<n; i++)
    r[i] = b[i]-r[i];

  bnrm2 = snorm2(n, b);
  if (bnrm2 > 0.0) {
    error = snorm2(n, r) / bnrm2;

    if (error > tol) { 
      /* Begin Iterations */
      for (k=0; k<n; k++) {
        for (i=0; i<n; i++)
          z[i] = r[i]*M[i];

        rho = sdot(n, r, z);

        if (k == 0)
          scopy(n, z, p);
        else {
          beta = rho/rho_1;
          for (i=0; i<n; i++)
            p[i] = z[i] + beta*p[i]; 
        }

        mvMult(n, rowptr, colind, values, p, q); /* q = A*p */

        alpha = rho/sdot(n, p, q);
        saxpy(n, alpha, p, x);    /* x = x + alpha*p */
        saxpy(n, -alpha, q, r);   /* r = r - alpha*q */
        error = snorm2(n, r) / bnrm2;
        if (error < tol)
          break;

        rho_1 = rho;
      }
    }
  }

  /*
  printf("%f %f %f %f => %f %f %f %f\n", b[0], b[1], b[2], b[3], x[0], x[1], x[2], x[3]);
  */
}


/*************************************************************************
* This function implements a sparse matrix-vector multiplication
**************************************************************************/
void mvMult(int nvtxs, idxtype *rowptr, idxtype *colind, float *values, float *x, float *y)
{
  int i, j;

  for (i=0; i<nvtxs; i++) {
    y[i] = 0.0;
    for (j=rowptr[i]; j<rowptr[i+1]; j++) 
      y[i] += values[j]*x[colind[j]];
  }
}
