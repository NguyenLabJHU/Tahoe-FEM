/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * initdiff.c 
 *
 * This file contains code for the initial directed diffusion at the coarsest
 * graph
 *
 * Started 5/19/97, Kirk, George
 *
 * $Id: initdiff.c,v 1.1.1.1 2004-10-07 16:05:25 paklein Exp $
 *
 */

#include <parmetis.h>


/*************************************************************************
* This function is the entry point of the directed diffusion algorithm
* for the coarsest graph.  This function assembles the graph to all the 
* processors and preceed serially.
**************************************************************************/
void InitDiffusion(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i, j, lpecut[2], gpecut[2];
  idxtype *vtxdist;
  GraphType *agraph;
  int *sendcounts, *displs;
  int rating;

  ASSERT(ctrl, ctrl->npes == ctrl->nparts);

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));

  agraph = AssembleAdaptiveGraph(ctrl, graph, wspace);

  /* myprintf(ctrl, "Assembled graph: %d %d\n", agraph->nvtxs, agraph->nedges); */
  MPI_Barrier(ctrl->comm);

  vtxdist = graph->vtxdist;

  agraph->where = idxmalloc(agraph->nvtxs, "InitDiffusion: agraph->where");
  for (i=0; i<ctrl->npes; i++) {
    for (j=vtxdist[i]; j<vtxdist[i+1]; j++)
      agraph->where[j] = i;
  }

  switch (ctrl->mype) {
    case 0:
      rating = SortedKWay_WD(ctrl, agraph, ctrl->nparts);
      break;
    case 1:
      /* rating = SortedKWay_WD(ctrl, agraph, ctrl->nparts); */
      rating = KWay_WavefrontDiffuser(ctrl, agraph, ctrl->nparts, .65);
      break;
    case 2:
      rating = KWay_WavefrontDiffuser(ctrl, agraph, ctrl->nparts, .75);
      break;
    case 3:
      rating = KWay_WavefrontDiffuser(ctrl, agraph, ctrl->nparts, .50);
      break;
    default:
      rating = KWay_WavefrontDiffuser(ctrl, agraph, ctrl->nparts, .25);
      break;
  }

  /* Determine which PE got the minimum cut */
  lpecut[0] = rating;
  MPI_Comm_rank(ctrl->comm, lpecut+1);
  MPI_Allreduce(lpecut, gpecut, 1, MPI_2INT, MPI_MINLOC, ctrl->gcomm);

  /* myprintf(ctrl, "Mincut: %d, (%d %d)\n", agraph->mincut, rating, gpecut[0]); */

  graph->where = idxmalloc(graph->nvtxs+graph->nrecv, "InitPartition: where");
  sendcounts = imalloc(ctrl->npes, "InitPartitionNew: sendcounts");
  displs = imalloc(ctrl->npes, "InitPartitionNew: displs");

  for (i=0; i<ctrl->npes; i++) {
    sendcounts[i] = vtxdist[i+1]-vtxdist[i];
    displs[i] = vtxdist[i];
  }

  MPI_Scatterv((void *)agraph->where, sendcounts, displs, IDX_DATATYPE, 
               (void *)graph->where, graph->nvtxs, IDX_DATATYPE, gpecut[1], ctrl->comm);

  GKfree(&sendcounts, &displs, LTERM);

  FreeGraph(agraph);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));

}



/*************************************************************************
* This function assembles the graph into a single processor
**************************************************************************/
GraphType *AssembleAdaptiveGraph(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i, j, k, l, gnvtxs, nvtxs, gnedges, nedges, firstvtx, gsize;
  idxtype *xadj, *vwgt, *vsize, *adjncy, *adjwgt, *vtxdist, *edgedist, *imap;
  idxtype *axadj, *aadjncy, *aadjwgt, *avwgt, *avsize, *alabel;
  idxtype *mygraph, *ggraph;
  int *recvcounts, *displs, mysize;
  GraphType *agraph;

  gnvtxs = graph->gnvtxs;
  nvtxs = graph->nvtxs;
  nedges = graph->xadj[nvtxs];
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vsize = graph->vsize;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  vtxdist = graph->vtxdist;
  imap = graph->imap;

  /* Determine the # of idxtype to receive from each processor */
  recvcounts = imalloc(ctrl->npes, "AssembleGraph: recvcounts");
  mysize = 3*nvtxs + 2*nedges;
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
    mygraph[k++] = vsize[i];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      mygraph[k++] = imap[adjncy[j]];
      mygraph[k++] = adjwgt[j];
    }
  }
  ASSERT(ctrl, mysize == k);

  /* Assemble the entire graph */
  gsize = displs[ctrl->npes];
  ggraph = (gsize <= wspace->maxcore-mysize ? wspace->core+mysize : idxmalloc(gsize, "AssembleGraph: ggraph"));
  MPI_Allgatherv((void *)mygraph, mysize, IDX_DATATYPE, (void *)ggraph, recvcounts, displs, IDX_DATATYPE, ctrl->comm);

  /* MPI_Bcast((void *)ggraph, gsize, IDX_DATATYPE, 0, ctrl->comm); */

  GKfree(&recvcounts, &displs, LTERM);
  if (mysize > wspace->maxcore)
    free(mygraph);

  agraph = CreateGraph();
  agraph->maxvwgt = graph->maxvwgt;
  agraph->nvtxs = gnvtxs;
  agraph->nedges = gnedges = (gsize-3*gnvtxs)/2;

  /* Allocate memory for the assembled graph */
  axadj = agraph->xadj = idxmalloc(gnvtxs+1, "AssembleGraph: axadj");
  avwgt = agraph->vwgt = idxmalloc(gnvtxs, "AssembleGraph: avwgt");
  avsize = agraph->vsize = idxmalloc(gnvtxs, "AssembleGraph: avsize");
  aadjncy = agraph->adjncy = idxmalloc(gnedges, "AssembleGraph: adjncy");
  aadjwgt = agraph->adjwgt = idxmalloc(gnedges, "AssembleGraph: adjwgt");
  alabel = agraph->label = idxmalloc(gnvtxs, "AssembleGraph: alabel");

  for (k=j=i=0; i<gnvtxs; i++) {
    axadj[i] = ggraph[k++];
    avwgt[i] = ggraph[k++];
    avsize[i] = ggraph[k++];
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

  if (gsize > wspace->maxcore-mysize)
    free(ggraph);

/*
  if (mysize > wspace->nlarge || displs[ctrl->npes] > 2*wspace->nlarge)
    myprintf(ctrl, "Had to allocate memory %d %d %d\n", mysize, displs[ctrl->npes], wspace->nlarge);
*/

  return agraph;
}

