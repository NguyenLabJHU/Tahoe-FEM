/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * memory.c
 *
 * This file contains routines that deal with memory allocation
 *
 * Started 2/24/96
 * George
 *
 * $Id: memory.c,v 1.1.1.1 2004-10-07 16:05:25 paklein Exp $
 *
 */

#include <parmetis.h>


/*************************************************************************
* This function allocate various pools of memory
**************************************************************************/
void PreAllocateMemory(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i;

  wspace->nlarge = 2*graph->nedges;

  wspace->maxcore = 8*graph->nedges;
  wspace->core = idxmalloc(wspace->maxcore, "PreAllocateMemory: wspace->core");

  wspace->pairs = (KeyValueType *)wspace->core;
  wspace->indices = (idxtype *)(wspace->pairs + wspace->nlarge);
  wspace->degrees = (EdgeType *)(wspace->indices + wspace->nlarge);


  wspace->pv1 = idxmalloc(ctrl->nparts+ctrl->npes+1, "PreAllocateMemory: wspace->pv?");
  wspace->pv2 = idxmalloc(ctrl->nparts+ctrl->npes+1, "PreAllocateMemory: wspace->pv?");
  wspace->pv3 = idxmalloc(ctrl->nparts+ctrl->npes+1, "PreAllocateMemory: wspace->pv?");
  wspace->pv4 = idxmalloc(ctrl->nparts+ctrl->npes+1, "PreAllocateMemory: wspace->pv?");

  wspace->pepairs1 = (KeyValueType *)GKmalloc(sizeof(KeyValueType)*(ctrl->nparts+ctrl->npes+1), "PreAllocateMemory: wspace->pepairs?");
  wspace->pepairs2 = (KeyValueType *)GKmalloc(sizeof(KeyValueType)*(ctrl->nparts+ctrl->npes+1), "PreAllocateMemory: wspace->pepairs?");

}


/*************************************************************************
* This function de-allocate various pools of memory
**************************************************************************/
void FreeWSpace(WorkSpaceType *wspace)
{

  GKfree(&wspace->core, 
         &wspace->pv1, 
         &wspace->pv2, 
         &wspace->pv3,
         &wspace->pv4, 
         &wspace->pepairs1, 
         &wspace->pepairs2, 
         LTERM);
}


/*************************************************************************
* This function de-allocates memory allocated for the control structures
**************************************************************************/
void FreeCtrl(CtrlType *ctrl)
{
  MPI_Comm_free(&(ctrl->gcomm));
}


/*************************************************************************
* This function creates a CoarseGraphType data structure and initializes
* the various fields
**************************************************************************/
GraphType *CreateGraph(void)
{
  GraphType *graph;

  graph = (GraphType *)GKmalloc(sizeof(GraphType), "CreateCoarseGraph: graph");

  InitGraph(graph);

  return graph;
}


/*************************************************************************
* This function creates a CoarseGraphType data structure and initializes
* the various fields
**************************************************************************/
void InitGraph(GraphType *graph) 
{
  graph->gnvtxs = graph->nvtxs = graph->nedges = graph->nsep = -1;
  graph->nnbrs = graph->nrecv = graph->nsend = graph->nlocal = -1;
  graph->xadj = graph->vwgt = graph->vsize = graph->adjncy = graph->adjwgt = NULL;
  graph->vtxdist = NULL;
  graph->match = graph->cmap = NULL;
  graph->label = NULL;

  graph->peind = NULL;
  graph->sendptr = graph->sendind = graph->recvptr = graph->recvind = NULL;
  graph->imap = NULL;
  graph->pexadj = graph->peadjncy = graph->peadjloc = NULL;
  graph->lperm = NULL;

  graph->slens = graph->rlens = NULL;
  graph->rcand = NULL;

  graph->where = graph->lpwgts = graph->gpwgts = NULL;
  graph->rinfo = NULL;

  graph->nrinfo = NULL;
  graph->sepind = NULL;

  graph->coarser = graph->finer = NULL;

}

/*************************************************************************
* This function deallocates any memory stored in a graph
**************************************************************************/
void FreeGraph(GraphType *graph) 
{

  GKfree(&graph->xadj, 
         &graph->vwgt,
         &graph->vsize,
         &graph->adjncy,
         &graph->adjwgt,
         &graph->vtxdist, 
         &graph->match, 
         &graph->cmap, 
         &graph->lperm, 
         &graph->label, 
         &graph->where, 
         &graph->rinfo, 
         &graph->nrinfo, 
         &graph->sepind,
         &graph->lpwgts, 
         &graph->gpwgts, 
         &graph->peind, 
         &graph->sendptr, 
         &graph->sendind, 
         &graph->recvptr, 
         &graph->recvind, 
         &graph->imap,
         &graph->rlens,
         &graph->slens,
         &graph->rcand,
         &graph->pexadj,
         &graph->peadjncy,
         &graph->peadjloc,
         LTERM);

  free(graph);
}


/*************************************************************************
* This function deallocates any memory stored in a graph
**************************************************************************/
void FreeGraphContent(GraphType *graph) 
{

  GKfree(&graph->xadj, 
         &graph->vwgt,
         &graph->vsize,
         &graph->adjncy,
         &graph->adjwgt,
         &graph->vtxdist, 
         &graph->match, 
         &graph->cmap, 
         &graph->lperm, 
         &graph->where, 
         &graph->label,
         &graph->rinfo, 
         &graph->nrinfo, 
         &graph->lpwgts, 
         &graph->gpwgts, 
         &graph->peind, 
         &graph->sendptr, 
         &graph->sendind, 
         &graph->recvptr, 
         &graph->recvind, 
         &graph->imap,
         &graph->rlens,
         &graph->slens,
         &graph->rcand,
         &graph->pexadj,
         &graph->peadjncy,
         &graph->peadjloc,
         LTERM);

}


/*************************************************************************
* This function deallocates any memory stored in a graph
**************************************************************************/
void FreeInitialGraphAndRemap(GraphType *graph, int wgtflag) 
{
  int i, j, k, nvtxs, nedges;
  idxtype *adjncy, *imap;

  nedges = graph->nedges;
  adjncy = graph->adjncy;
  imap = graph->imap;

  if (imap != NULL) {
    for (i=0; i<nedges; i++)
      adjncy[i] = imap[adjncy[i]];  /* Apply local to global transformation */
  }

  /* Free Metis's things */
  GKfree(&graph->match, 
         &graph->cmap, 
         &graph->lperm, 
         &graph->where, 
         &graph->label, 
         &graph->rinfo, 
         &graph->nrinfo, 
         &graph->lpwgts, 
         &graph->gpwgts, 
         &graph->sepind,
         &graph->peind, 
         &graph->sendptr, 
         &graph->sendind, 
         &graph->recvptr, 
         &graph->recvind, 
         &graph->imap,
         &graph->rlens,
         &graph->slens,
         &graph->rcand,
         &graph->pexadj,
         &graph->peadjncy,
         &graph->peadjloc,
         LTERM);

  if ((wgtflag&2) == 0) 
    GKfree(&graph->vwgt, &graph->vsize, LTERM);
  if ((wgtflag&1) == 0) 
    GKfree(&graph->adjwgt, LTERM);

  free(graph);
}
