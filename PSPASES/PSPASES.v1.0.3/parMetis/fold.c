/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * fold.c
 *
 * This file contains routines that deal with folding and unfolding a graph
 *
 * Started 11/19/96
 * George
 *
 * $Id: fold.c,v 1.1.1.1 2004-10-07 16:05:25 paklein Exp $
 *
 */

#include <parmetis.h>


/*************************************************************************
* This function checks if there is enough memory to do folding
**************************************************************************/
int EnoughMemory(CtrlType *ctrl, int mynedges, int nlarge)
{
  int othernedges, val;

  if (ctrl->mype%2 == 0) {
    MPI_Recv((void *)&othernedges, 1, MPI_INT, ctrl->mype+1, 1, ctrl->comm, &ctrl->status);
    val = (mynedges+othernedges <= nlarge ? 0 : 1);
  }
  else {
    MPI_Send((void *)&mynedges, 1, MPI_INT, ctrl->mype-1, 1, ctrl->comm);
    val = 0;
  }

  return (GlobalSESum(ctrl, val) == 0);
}


/*************************************************************************
* This function folds the graph to the even processors. 
**************************************************************************/
GraphType *FoldGraph(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i, j, k, nvtxs, even, npes=ctrl->npes, fnpes=npes/2, fnvtxs, fnedges, fmype;
  idxtype *vtxdist, *fldxadj, *fldvwgt, *fldvtxdist;
  GraphType *fldgraph;

  even = ctrl->mype%2 == 0;
  fmype = ctrl->mype/2;

  nvtxs = graph->nvtxs;
  vtxdist = graph->vtxdist;

  fldgraph = NULL;

  /* Construct the fldvtxdist */
  if (even) {
    fldgraph = CreateGraph();
    fldgraph->maxvwgt = graph->maxvwgt;

    fldvtxdist = fldgraph->vtxdist = idxmalloc(fnpes+1, "FoldGraph: fldvtxdist");
    for (i=0; i<=fnpes; i++)
      fldvtxdist[i] = vtxdist[2*i];

    fnvtxs = fldvtxdist[fmype+1]-fldvtxdist[fmype];

    fldxadj = fldgraph->xadj = idxmalloc(fnvtxs+1, "FoldGraph: fldxadj");
    fldvwgt = fldgraph->vwgt = idxmalloc(fnvtxs, "FoldGraph: fldxvwgt");

    for (i=0; i<nvtxs; i++)
      fldxadj[i] = graph->xadj[i+1]-graph->xadj[i];
    idxcopy(graph->nvtxs, graph->vwgt, fldvwgt);

    MPI_Recv((void *)(fldxadj+graph->nvtxs), fnvtxs-nvtxs+1, IDX_DATATYPE, ctrl->mype+1, 1, 
             ctrl->comm, &ctrl->status);
    MPI_Recv((void *)(fldvwgt+graph->nvtxs), fnvtxs-nvtxs, IDX_DATATYPE, ctrl->mype+1, 1, 
             ctrl->comm, &ctrl->status);

    for (i=0; i<fnvtxs-nvtxs; i++) 
      fldxadj[nvtxs+i] = fldxadj[nvtxs+i+1]-fldxadj[nvtxs+i];
    MAKECSR(i, fnvtxs, fldxadj);

    fnedges = fldxadj[fnvtxs];
    fldgraph->adjncy = idxmalloc(fnedges, "FoldGraph: fldadjncy");
    fldgraph->adjwgt = idxmalloc(fnedges, "FoldGraph: fldadjwgt");

    idxcopy(graph->nedges, graph->adjncy, fldgraph->adjncy);
    idxcopy(graph->nedges, graph->adjwgt, fldgraph->adjwgt);

    MPI_Recv((void *)(fldgraph->adjncy+graph->nedges), fnedges-graph->nedges, IDX_DATATYPE, 
             ctrl->mype+1, 1, ctrl->comm, &ctrl->status);
    MPI_Recv((void *)(fldgraph->adjwgt+graph->nedges), fnedges-graph->nedges, IDX_DATATYPE, 
             ctrl->mype+1, 1, ctrl->comm, &ctrl->status);

    fldgraph->level = graph->level + 1;
    fldgraph->nvtxs = fnvtxs;
    fldgraph->nedges = fnedges;
    fldgraph->gnvtxs = graph->gnvtxs;

    if (fnedges > wspace->nlarge/2)
      printf("FNEDGES is higher that original! %d %d\n", fnedges, wspace->nlarge/2);
  }
  else {
    MPI_Send((void *)graph->xadj, graph->nvtxs+1, IDX_DATATYPE, ctrl->mype-1, 1, ctrl->comm);
    MPI_Send((void *)graph->vwgt, graph->nvtxs, IDX_DATATYPE, ctrl->mype-1, 1, ctrl->comm);
    MPI_Send((void *)graph->adjncy, graph->nedges, IDX_DATATYPE, ctrl->mype-1, 1, ctrl->comm);
    MPI_Send((void *)graph->adjwgt, graph->nedges, IDX_DATATYPE, ctrl->mype-1, 1, ctrl->comm);
  }

  return fldgraph;
}


/*************************************************************************
* This function folds the graph to the even processors. 
**************************************************************************/
void UnFoldGraph(CtrlType *ctrl, GraphType *graph, GraphType *fldgraph, WorkSpaceType *wspace)
{
  int even;

  even = ctrl->mype%2 == 0;

  graph->where = idxmalloc(graph->nvtxs, "UnFoldGraph: graph->where");

  /* Fillin graph->where properly */
  if (even) {
    idxcopy(graph->nvtxs, fldgraph->where, graph->where);
    MPI_Send((void *)(fldgraph->where+graph->nvtxs), fldgraph->nvtxs-graph->nvtxs, IDX_DATATYPE, 
             ctrl->mype+1, 1, ctrl->comm);

    FreeGraph(fldgraph); 
  }
  else {
    MPI_Recv((void *)graph->where, graph->nvtxs, IDX_DATATYPE, ctrl->mype-1, 1, ctrl->comm, &ctrl->status);
  }

}


