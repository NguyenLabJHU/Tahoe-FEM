/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * initpart.c
 *
 * This file contains code that performs log(p) parallel multilevel
 * recursive bissection
 *
 * Started 3/4/96
 * George
 *
 * $Id: initpart.c,v 1.1.1.1 2004-10-07 16:05:25 paklein Exp $
 */

#include <parmetis.h>


#define DEBUG_IPART_



/*************************************************************************
* This function is the entry point of the initial partitioning algorithm.
* This algorithm assembles the graph to all the processors and preceed
* serially.
**************************************************************************/
void InitPartition(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace, int context)
{
  int i, j, lpecut[2], gpecut[2];
  int numflag=0, wgtflag = 3, moptions[10];
  idxtype *vtxdist, *part;
  GraphType *agraph;
  int *sendcounts, *displs;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));

  agraph = AssembleGraph(ctrl, graph, wspace, context);
  part = idxmalloc(agraph->nvtxs, "InitPartition: part");

/*
  if (ctrl->mype == 0) 
    WriteMetisGraph(agraph->nvtxs, agraph->xadj, agraph->adjncy, agraph->vwgt, agraph->adjwgt);
*/

  MALLOC_CHECK(NULL);

  moptions[0] = 1;
  moptions[1] = 3;
  moptions[2] = 1;
  moptions[3] = 1;
  moptions[4] = 0; /*(ctrl->mype == 0 ? 15 : 0);*/
  moptions[7] = ctrl->mype;
  METIS_PartGraphKway2(&agraph->nvtxs, agraph->xadj, agraph->adjncy, agraph->vwgt, agraph->adjwgt, 
        &wgtflag, &numflag, &ctrl->nparts, moptions, &lpecut[0], part);

  MALLOC_CHECK(NULL);

  /* Determine which PE got the minimum cut */
  MPI_Comm_rank(ctrl->gcomm, lpecut+1);
  MPI_Allreduce(lpecut, gpecut, 1, MPI_2INT, MPI_MINLOC, ctrl->gcomm);

  /* myprintf(ctrl, "Mincut: %d, GMincut: %d\n", lpecut[0], gpecut[0]); */

  if (lpecut[1] == gpecut[1] && gpecut[1] != 0) 
    MPI_Send((void *)part, agraph->nvtxs, IDX_DATATYPE, 0, 1, ctrl->gcomm);
  if (lpecut[1] == 0 && gpecut[1] != 0)
    MPI_Recv((void *)part, agraph->nvtxs, IDX_DATATYPE, gpecut[1], 1, ctrl->gcomm, &ctrl->status);
  
  if (context == 1) {
    /* The minimum PE performs the Scatter */
    vtxdist = graph->vtxdist;
    graph->where = idxmalloc(graph->nvtxs+graph->nrecv, "InitPartition: where");

    sendcounts = imalloc(ctrl->npes, "InitPartitionNew: sendcounts");
    displs = imalloc(ctrl->npes, "InitPartitionNew: displs");

    for (i=0; i<ctrl->npes; i++) {
      sendcounts[i] = vtxdist[i+1]-vtxdist[i];
      displs[i] = vtxdist[i];
    }

    MPI_Scatterv((void *)part, sendcounts, displs, IDX_DATATYPE, 
                 (void *)graph->where, graph->nvtxs, IDX_DATATYPE, 0, ctrl->comm);

    GKfree(&sendcounts, &displs, LTERM);
  }

  free(part);
  FreeGraph(agraph);

  MALLOC_CHECK(NULL);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));

}




/*************************************************************************
* This function is the entry point of the initial partition algorithm
* that does recursive bissection.
* This algorithm assembles the graph to all the processors and preceeds
* by parallelizing the recursive bisection step.
**************************************************************************/
void InitPartition_RB(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace, int context)
{
  int i, j, mype, npes, gnvtxs, mypartnum;
  idxtype *part, *gwhere0, *gwhere1;
  GraphType *agraph;
  int lnparts, fpart, fpe, lnpes; 
  int dpwgt, twoparts=2, numflag = 0, wgtflag = 3, moptions[10], edgecut;
  float tpwgts[2];

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));

  agraph = AssembleGraph(ctrl, graph, wspace, context);
  part = idxmalloc(agraph->nvtxs, "InitPartition: part");

  MPI_Comm_rank(ctrl->gcomm, &mype);
  MPI_Comm_size(ctrl->gcomm, &npes);

  gnvtxs = agraph->nvtxs;
  dpwgt = idxsum(gnvtxs, agraph->vwgt)/ctrl->nparts;

  gwhere0 = idxsmalloc(gnvtxs, 0, "InitPartition: gwhere0");
  gwhere1 = idxmalloc(gnvtxs, "InitPartition: gwhere1");


  /* Go into the recursive bisection */
  moptions[0] = 0;
  lnparts = ctrl->nparts;
  fpart = fpe = 0;
  lnpes = npes;
  while (lnpes > 1 && lnparts > 1) {
    tpwgts[0] = (1.0*(lnparts>>1))/(1.0*lnparts);
    tpwgts[1] = 1.0 - tpwgts[0];

    METIS_WPartGraphRecursive(&agraph->nvtxs, agraph->xadj, agraph->adjncy, agraph->vwgt, 
          agraph->adjwgt, &wgtflag, &numflag, &twoparts, tpwgts, moptions, &edgecut, part);

    if (mype < fpe+lnpes/2) {  /* I'm picking the left branch */
      KeepPart(ctrl, agraph, wspace, part, 0);
      lnpes = lnpes/2;
      lnparts = lnparts/2;
    }
    else {
      KeepPart(ctrl, agraph, wspace, part, 1);
      fpart = fpart + lnparts/2;
      fpe = fpe + lnpes/2;
      lnpes = lnpes - lnpes/2;
      lnparts = lnparts - lnparts/2;
    }
  }

  if (lnparts == 1) { /* Take care the case in which npes is greater than or equal to nparts */
    if (mype == fpe) { /* Only the first process will assign labels (for the reduction to work) */
      for (i=0; i<agraph->nvtxs; i++) 
        gwhere0[agraph->label[i]] = fpart;
    }
  }
  else { /* Take care the case in which npes is smaller than nparts */
    METIS_PartGraphKway(&agraph->nvtxs, agraph->xadj, agraph->adjncy, agraph->vwgt, 
          agraph->adjwgt, &wgtflag, &numflag, &lnparts, moptions, &edgecut, part);

    for (i=0; i<agraph->nvtxs; i++) 
      gwhere0[agraph->label[i]] = fpart + part[i];
  }


  MPI_Allreduce((void *)gwhere0, (void *)gwhere1, gnvtxs, IDX_DATATYPE, MPI_SUM, ctrl->gcomm);


  if (context == 1) {  /* The active processors set the where information */
    graph->where = idxmalloc(graph->nvtxs+graph->nrecv, "InitPartition: where");
    idxcopy(graph->nvtxs, gwhere1+graph->vtxdist[ctrl->mype], graph->where);
  }

  FreeGraph(agraph);
  GKfree(&gwhere0, &gwhere1, &part, LTERM);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));

}



/*************************************************************************
* This function assembles the graph into a single processor
**************************************************************************/
GraphType *AssembleGraph(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace, int context)
{
  int i, j, k, l, gnvtxs, nvtxs, gnedges, nedges, firstvtx, gsize;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *vtxdist, *edgedist, *imap;
  idxtype *axadj, *aadjncy, *aadjwgt, *avwgt, *alabel;
  idxtype *mygraph, *ggraph;
  int *recvcounts, *displs, mysize;
  GraphType *agraph;

  if (context == 1) {
    gnvtxs = graph->gnvtxs;
    nvtxs = graph->nvtxs;
    nedges = graph->xadj[nvtxs];
    xadj = graph->xadj;
    vwgt = graph->vwgt;
    adjncy = graph->adjncy;
    adjwgt = graph->adjwgt;
    vtxdist = graph->vtxdist;
    imap = graph->imap;

    /* Determine the # of idxtype to receive from each processor */
    recvcounts = imalloc(ctrl->npes, "AssembleGraph: recvcounts");
    mysize = 2*nvtxs + 2*nedges;
    MPI_Allgather((void *)(&mysize), 1, MPI_INT, (void *)recvcounts, 1, MPI_INT, ctrl->comm);
  
    displs = imalloc(ctrl->npes+1, "AssembleGraph: displs");
    displs[0] = 0;
    for (i=1; i<ctrl->npes+1; i++) 
      displs[i] = displs[i-1] + recvcounts[i-1];

    /* Construct the one-array storage format of the assembled graph */
    mygraph = (mysize <= wspace->maxcore ? wspace->core : idxmalloc(mysize, "AssembleGraph: mygraph"));
    for (k=i=0; i<nvtxs; i++) {
      mygraph[k++] = xadj[i+1]-xadj[i];
      mygraph[k++] = vwgt[i];
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        mygraph[k++] = imap[adjncy[j]];
        mygraph[k++] = adjwgt[j];
      }
    }
    ASSERT(ctrl, mysize == k);

    /* Assemble the entire graph */
    gsize = displs[ctrl->npes]+2;
    ggraph = (gsize <= wspace->maxcore-mysize ? wspace->core+mysize : idxmalloc(gsize, "AssembleGraph: ggraph"));
    ggraph[0] = gnvtxs;
    ggraph[1] = (displs[ctrl->npes] - 2*gnvtxs)/2;
    MPI_Gatherv((void *)mygraph, mysize, IDX_DATATYPE, (void *)(ggraph+2), 
                   recvcounts, displs, IDX_DATATYPE, 0, ctrl->comm);

    GKfree(&recvcounts, &displs, LTERM);
    if (mysize > wspace->maxcore)
      free(mygraph);
  }

  MPI_Bcast((void *)&gsize, 1, MPI_INT, 0, ctrl->gcomm);

  /* printf("Gsize: %d, %d\n", gsize, context); */

  if (context == 0) 
    ggraph = (gsize <= wspace->maxcore ? wspace->core : idxmalloc(gsize, "AssembleGraph: ggraph"));

  MPI_Bcast((void *)ggraph, gsize, IDX_DATATYPE, 0, ctrl->gcomm);

  agraph = CreateGraph();
  agraph->maxvwgt = graph->maxvwgt;
  agraph->nvtxs = gnvtxs = ggraph[0];
  agraph->nedges = gnedges = ggraph[1];

  /* Allocate memory for the assembled graph */
  axadj = agraph->xadj = idxmalloc(gnvtxs+1, "AssembleGraph: axadj");
  avwgt = agraph->vwgt = idxmalloc(gnvtxs, "AssembleGraph: avwgt");
  aadjncy = agraph->adjncy = idxmalloc(gnedges, "AssembleGraph: adjncy");
  aadjwgt = agraph->adjwgt = idxmalloc(gnedges, "AssembleGraph: adjwgt");
  alabel = agraph->label = idxmalloc(gnvtxs, "AssembleGraph: alabel");

  for (k=2, j=i=0; i<gnvtxs; i++) {
    axadj[i] = ggraph[k++];
    avwgt[i] = ggraph[k++];
    for (l=0; l<axadj[i]; l++) {
      aadjncy[j] = ggraph[k++];
      aadjwgt[j] = ggraph[k++];
      j++;
    }
  }

  /* Now fix up the received graph */
  MAKECSR(i, gnvtxs, axadj);

  for (i=0; i<gnvtxs; i++)
    alabel[i] = i;

  if (context == 1 && gsize > wspace->maxcore-mysize)
    free(ggraph);
  else if (context == 0 && gsize > wspace->maxcore)
    free(ggraph);

/*
  if (mysize > wspace->nlarge || displs[ctrl->npes] > 2*wspace->nlarge)
    myprintf(ctrl, "Had to allocate memory %d %d %d\n", mysize, displs[ctrl->npes], wspace->nlarge);
*/

  return agraph;
}




/*************************************************************************
* This function keeps one parts
**************************************************************************/
void KeepPart(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace, idxtype *part, int mypart)
{
  int i, j, k, l, nvtxs, mynvtxs, mynedges;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *label;
  idxtype *rename;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  label = graph->label;

  rename = idxmalloc(nvtxs, "KeepPart: rename");
  
  for (mynvtxs=0, i=0; i<nvtxs; i++) {
    if (part[i] == mypart)
      rename[i] = mynvtxs++;
  }

  for (mynvtxs=0, mynedges=0, j=xadj[0], i=0; i<nvtxs; i++) {
    if (part[i] == mypart) {
      for (; j<xadj[i+1]; j++) {
        k = adjncy[j];
        if (part[k] == mypart) {
          adjncy[mynedges] = rename[k];
          adjwgt[mynedges++] = adjwgt[j];
        }
      }
      j = xadj[i+1];  /* Save xadj[i+1] for later use */

      vwgt[mynvtxs] = vwgt[i];
      label[mynvtxs] = label[i];
      xadj[++mynvtxs] = mynedges;

    }
    else {
      j = xadj[i+1];  /* Save xadj[i+1] for later use */
    }
  }

  graph->nvtxs = mynvtxs;
  graph->nedges = mynedges;

  free(rename);
}



