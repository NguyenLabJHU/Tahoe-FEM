/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * wavefront.c 
 *
 * This file contains code for the initial directed diffusion at the coarsest
 * graph
 *
 * Started 5/19/97, Kirk, George
 *
 * $Id: wave.c,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $
 *
 */

#include <parmetis.h>

/*************************************************************************
* This function performs a k-way directed diffusion
**************************************************************************/
int KWay_WavefrontDiffuser(CtrlType *ctrl, GraphType *graph, int nparts, float flowFactor)
{
  int ii, i, j, k, l, nvtxs, from, to, done, totalv;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *vsize, *where, *home, *pwgts, *perm;
  idxtype *transfer, *tmpvec, *trvec, *wspace, *rowptr, *colind;
  int swaps, first, second, third, fromflow, toflow;
  float balance, mean, *load, *values;
  PQueueType q;
  MatrixType lmat;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vsize = graph->vsize;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;

  srand(ctrl->mype);
  srand48(ctrl->mype);
  PQueueInit(&q, nparts);

  home = idxmalloc(nvtxs, "KWay_WavefrontDiffuser: home");
  pwgts = idxsmalloc(nparts, 0, "KWay_WavefrontDiffuser: pwgts");
  trvec = idxmalloc(nparts, "KWay_WavefrontDiffuser: trvec");
  load = fmalloc(nparts, "KWay_WavefrontDiffuser: load");
  perm = idxmalloc(nvtxs, "KWay_WavefrontDiffuser: perm");

  rowptr = lmat.rowptr = idxmalloc(nparts+1, "KWay_WavefrontDiffuser: lmat.rowptr");
  colind = lmat.colind = idxmalloc(graph->nedges, "KWay_WavefrontDiffuser: lmat.colind");
  transfer = lmat.transfer = idxmalloc(amax(graph->nedges, graph->nedges*sizeof(float)/sizeof(idxtype)), "KWay_WavefrontDiffuser: lmat.transfer");
  values = lmat.values = (float *)lmat.transfer;
  wspace = idxmalloc(amax(2*nparts+1+nvtxs,6*nparts*sizeof(float)/sizeof(idxtype)), "KWay_WavefrontDiffuser: wspace");
  tmpvec = idxmalloc(nparts, "KWay_WavefrontDiffuser: transfer");

  for (i=0; i<nvtxs; i++) {
    pwgts[where[i]] += vwgt[i];
    home[i] = where[i];
    perm[i] = i;
  }
  mean = idxsum(nparts, pwgts)/nparts;

  for (done=0, l=amin(nparts/2, NGD_PASSES); l>0; l--) {
    /* Set-up and solve the diffusion equation */
    PQueueReset(&q);
    for (j=0; j<nparts; j++) 
      load[j] = (float)(pwgts[j]) - mean;

    SolveDiffusionEquation(graph, &lmat, nparts, load, trvec, wspace);

    for (i=0; i<nparts; i++) 
      PQueueInsert(&q, i, trvec[i]);

    FastRandomPermute(nvtxs, perm, 0);

    for (ii=0; ii<nvtxs; ii++) {
      i = perm[ii];
      from = where[i];

      if (pwgts[from] == vwgt[i] || (from != q.perm[0] && from != q.perm[1] && from != q.perm[2] && from == home[i]))
        continue;

      /* Scatter the sparse transfer row into the dense tmpvec row */
      for (j=rowptr[from]+1; j<rowptr[from+1]; j++)
        tmpvec[colind[j]] = transfer[j];

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        to = where[adjncy[j]];
        if (from != to) {
          if (tmpvec[to] > (flowFactor * vwgt[i])) {
            tmpvec[to] -= vwgt[i];
            where[i] = to;
            INC_DEC(pwgts[to], pwgts[from], vwgt[i]);
            INC_DEC(trvec[to], trvec[from], vwgt[i]);
            PQueueUpdate(&q, from, trvec[from]+vwgt[i], trvec[from]);
            PQueueUpdate(&q, to, trvec[to]-vwgt[i], trvec[to]);
            break;
          }
        }
      }

      /* Gather the dense tmpvec row into the sparse transfer row */
      for (j=rowptr[from]+1; j<rowptr[from+1]; j++)
        transfer[j] = tmpvec[colind[j]];

    }

    balance = (1.0*pwgts[idxamax(nparts, pwgts)])/(1.0*mean);

    if (balance < UNBALANCE_FRACTION + .05)
      done = 1;

    if (GlobalSESum(ctrl, done) > 0)
      break;
  }

  graph->mincut = ComputeSerialEdgeCut(graph);
  totalv = ComputeSerialTotalV(graph, home);

  /* printf("%d %d %5.3f\n", graph->mincut, totalv, balance); */

  GKfree(&home, &pwgts, &trvec, &load, &perm, &rowptr, &colind, &transfer, &wspace, &tmpvec, LTERM);

  return (int)(balance * (float)(nparts*mean)) + totalv;
}


/*************************************************************************
* This function performs a k-way directed diffusion
**************************************************************************/
int SortedKWay_WD(CtrlType *ctrl, GraphType *graph, int nparts)
{
  int ii, i, j, k, l, nvtxs, from, to, edge, ewgt, done, max1, max2, max3, totalv;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *vsize;
  idxtype *where, *home, *pwgts, *ed, *perm;
  idxtype *transfer, *tmpvec, *trvec, *wspace;
  float sum, mean, *load, *workspace;
  idxtype *rowptr, *colind;
  float balance, *values;
  MatrixType lmat;
  KeyValueType *cand;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vsize = graph->vsize;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;

  srand(ctrl->mype);
  srand48(ctrl->mype);

  home = idxmalloc(nvtxs, "SortedDiffuser: home");
  perm = idxmalloc(nvtxs, "SortedDiffuser: perm");
  ed = idxsmalloc(nvtxs, 0, "SortedDiffuser: ed");
  pwgts = idxsmalloc(nparts, 0, "SortedDiffuser: pwgts");
  trvec = idxmalloc(nparts, "SortedDiffuser: trvec");
  load = fmalloc(nparts, "SortedDiffuser: load");
  cand = (KeyValueType *)GKmalloc(nvtxs*sizeof(KeyValueType), "SortedKWay_WD: cand");

  rowptr = lmat.rowptr = idxmalloc(nparts+1, "KWay_WavefrontDiffuser: lmat.rowptr");
  colind = lmat.colind = idxmalloc(graph->nedges, "KWay_WavefrontDiffuser: lmat.colind");
  transfer = lmat.transfer = idxmalloc(amax(graph->nedges, graph->nedges*sizeof(float)/sizeof(idxtype)), "KWay_WavefrontDiffuser: lmat.transfer");
  values = lmat.values = (float *)lmat.transfer;
  wspace = idxmalloc(amax(2*nparts+1+nvtxs,6*nparts*sizeof(float)/sizeof(idxtype)), "KWay_WavefrontDiffuser: wspace");
  tmpvec = idxmalloc(nparts, "SortedDiffuser: tmpvec");


  for (i=0; i<nvtxs; i++) {
    pwgts[where[i]] += vwgt[i];
    home[i] = where[i];

    /* Compute the ed */
    for (j=xadj[i]; j<xadj[i+1]; j++) 
      ed[i] += (where[i] != where[adjncy[j]] ? adjwgt[j] : 0);
  }
  mean = idxsum(nparts, pwgts)/nparts;


  /*----------------------------------------------------------------- 
   * Perform a number of diffusion iterations                     
   *----------------------------------------------------------------- */
  for (done=0, l=amin(nparts/2, NGD_PASSES); l>0; l--) {
    FastRandomPermute(nvtxs, perm, 1);
    for (j=nvtxs, k=0, ii=0; ii<nvtxs; ii++) {
      i = perm[ii];
      if (ed[i] != 0) {
        cand[k].key = -ed[i];
        cand[k++].val = i;
      }
      else {
        cand[--j].key = 0;
        cand[j].val = i;
      }
    }
    ikeysort(k, cand);
  
    /* Set-up and solve the diffusion equations */
    for (i=0; i<nparts; i++) 
      load[i] = (float)(pwgts[i]) - mean;

    SolveDiffusionEquation(graph, &lmat, nparts, load, trvec, wspace);

    /* Determine the three largest partitions */
    for (max1 = max2 = max3 = 0, i=1; i<nparts; i++) {
      if (trvec[max1] < trvec[i]) {
        max3 = max2;
        max2 = max1;
        max1 = i;
      }
      else if (trvec[max2] < trvec[i]) {
        max3 = max2;
        max2 = i;
      } 
      else if (trvec[max3] < trvec[i]) {
        max3 = i;
      }
    }

    /* Get into the loop and start moving vertices */
    for (ii=0; ii<nvtxs; ii++) {
      i = cand[ii].val;
      from = where[i];

      if (pwgts[from] == vwgt[i] || (from != max1 && from != max2 && from != max3 && from == home[i]))
        continue;

      /* Scatter the sparse transfer row into the dense tmpvec row */
      for (j=rowptr[from]+1; j<rowptr[from+1]; j++)
        tmpvec[colind[j]] = transfer[j];

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        to = where[adjncy[j]];
        if (from != to) {
          if (tmpvec[to] > .25*vwgt[i]) {
            tmpvec[to] -= vwgt[i];
            INC_DEC(pwgts[to], pwgts[from], vwgt[i]);
            where[i] = to;

            /* Update all external degrees */
            ed[i] = 0;
            for (k=xadj[i]; k<xadj[i+1]; k++) {
              edge = adjncy[k];
              ed[i] += (to != where[edge] ? adjwgt[k] : 0);

              if (where[edge] == from) 
                ed[edge] += adjwgt[k];
              if (where[edge] == to) 
                ed[edge] -= adjwgt[k];
            }
            break;
          }
        }
      }

      /* gather the dense tmpvec row into the sparse transfer row */
      for (j=rowptr[from]+1; j<rowptr[from+1]; j++)
        transfer[j] = tmpvec[colind[j]];

    }

    balance = (1.0*pwgts[idxamax(nparts, pwgts)])/(1.0*mean);

    if (balance < UNBALANCE_FRACTION + .05)
      done = 1;

    if (GlobalSESum(ctrl, done) > 0)
      break;
  }

  graph->mincut = ComputeSerialEdgeCut(graph);
  totalv = ComputeSerialTotalV(graph, home);

  /* printf("%d %d %5.3f\n", graph->mincut, totalv, balance); */

  GKfree(&home, &perm, &ed, &pwgts, &trvec, &load, &cand, &rowptr, &colind, &transfer, &wspace, &tmpvec, LTERM);

  return (int)(balance * (float)(nparts*mean)) + totalv;
}




/*************************************************************************
* This function sets up the Laplacian matrix
**************************************************************************/
void SolveDiffusionEquation(GraphType *graph, MatrixType *lmat, int nparts, float *load, 
       idxtype *trvec, idxtype *wspace)
{
  int i, ii, j, jj, k, l, nvtxs;
  idxtype *xadj, *adjncy, *where;
  idxtype *pcounts, *perm, *marker;
  idxtype *rowptr, *colind, *transfer;
  float *values, *fwspace, *solution;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;

  rowptr = lmat->rowptr;
  colind = lmat->colind;
  values = lmat->values;
  transfer = lmat->transfer;

  pcounts = idxset(nparts, 0, wspace);
  marker = idxset(nparts, -1, wspace+nparts+1);
  perm = wspace+2*nparts+1;

  for (i=0; i<nvtxs; i++)
    pcounts[where[i]]++;
  MAKECSR(i, nparts, pcounts);

  for (i=0; i<nvtxs; i++) 
    perm[pcounts[where[i]]++] = i;

  for (i=nparts; i>0; i--) 
    pcounts[i] = pcounts[i-1];
  pcounts[0] = 0;


  /* Construct the matrix */
  rowptr[0] = k = 0;
  for (ii=0; ii<nparts; ii++) {
    colind[k++] = ii;
    marker[ii] = ii;

    for (jj=pcounts[ii]; jj<pcounts[ii+1]; jj++) {
      i = perm[jj];
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        l = where[adjncy[j]];
        if (marker[l] != ii) {
          colind[k] = l;
          values[k++] = -1;
          marker[l] = ii;
        }
      }
    }
    values[rowptr[ii]] = k-rowptr[ii]-1;
    rowptr[ii+1] = k;
  }

  fwspace = (float*)wspace;
  solution = fwspace;
  ConjGrad(nparts, rowptr, colind, values, load, solution, .001, fwspace+nparts);

  /* Create the trvec array */
  idxset(nparts, 0, trvec);
  for (j=0; j<nparts; j++) {
    for (k=rowptr[j]+1; k<rowptr[j+1]; k++) {
      if (solution[j] > solution[colind[k]]) {
        transfer[k] = (int) (solution[j] - solution[colind[k]]);
        trvec[j] += transfer[k];
        trvec[colind[k]] -= transfer[k];
      }
    }
  }

} 


