/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * drivers.c
 *
 * This file contains the driving routines for the various parallel
 * multilevel partitioning and repartitioning algorithms
 *
 * Started 11/19/96
 * George
 *
 * $Id: drivers.c,v 1.1.1.1 2004-10-07 16:05:25 paklein Exp $
 *
 */

#include <parmetis.h>



/*************************************************************************
* This function is the driver to the cold-start global partitioning
* algorithm.
**************************************************************************/
void Global_Partition(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{

  SetUp(ctrl, graph, wspace);

  IFSET(ctrl->dbglvl, DBG_PROGRESS, rprintf(ctrl, "[%6d %8d %5d %5d][%d][%d,%d]\n", 
        graph->gnvtxs, GlobalSESum(ctrl, graph->nedges), GlobalSEMin(ctrl, graph->nvtxs),
        GlobalSEMax(ctrl, graph->nvtxs), ctrl->CoarsenTo, graph->maxvwgt,
        GlobalSEMax(ctrl, graph->vwgt[idxamax(graph->nvtxs, graph->vwgt)])));  

  if (graph->gnvtxs < 1.3*ctrl->CoarsenTo || (graph->finer != NULL && graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)) {
    /* Done with coarsening. Find a partition */
    switch (ctrl->ipart) {
      case IPART_SER:
        InitPartition(ctrl, graph, wspace, 1);
        break;
      case IPART_RB:
        InitPartition_RB(ctrl, graph, wspace, 1);
        break;
      default:
        errexit("Unknown IPART type!\n");
    }

    if (graph->finer == NULL) { /* Do that only of no-coarsening took place */
      ComputePartitionParams(ctrl, graph, wspace);
      KWayRefine(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
    }
  }
  else { /* Coarsen it and the partition it */
    GlobalMatch_HEM(ctrl, graph, wspace);

    Global_Partition(ctrl, graph->coarser, wspace);

    ProjectPartition(ctrl, graph, wspace);
    ComputePartitionParams(ctrl, graph, wspace);
    KWayRefine(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
  }

}



/*************************************************************************
* This function is the driver to the cold-start global partitioning
* algorithm.
**************************************************************************/
void FldGlobal_Partition(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace, int level)
{
  GraphType *fldgraph;
  MPI_Comm newcomm, oldcomm;
  int gnedges;


  IFSET(ctrl->dbglvl, DBG_PROGRESS, rprintf(ctrl, "[%6d %8d %5d %5d][%d][%d,%d]\n", 
        graph->gnvtxs, GlobalSESum(ctrl, graph->nedges), GlobalSEMin(ctrl, graph->nvtxs),
        GlobalSEMax(ctrl, graph->nvtxs), ctrl->CoarsenTo, graph->maxvwgt,
        GlobalSEMax(ctrl, graph->vwgt[idxamax(graph->nvtxs, graph->vwgt)])));  

  if (graph->gnvtxs < 1.3*ctrl->CoarsenTo || (graph->finer != NULL && graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)) {
    /* Done with coarsening. Find a partition */
    SetUp(ctrl, graph, wspace);

    switch (ctrl->ipart) {
      case IPART_SER:
        InitPartition(ctrl, graph, wspace, 1);
        break;
      case IPART_RB:
        InitPartition_RB(ctrl, graph, wspace, 1);
        break;
      default:
        errexit("Unknown IPART type!\n");
    }

    if (graph->finer == NULL) { /* Do that only of no-coarsening took place */
      ComputePartitionParams(ctrl, graph, wspace);
      KWayRefine(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
    }
  }
  else { /* Coarsen it and then partition it */
    gnedges = GlobalSESum(ctrl, graph->nedges)/(ctrl->npes*ctrl->npes);
    if (gnedges-ctrl->foldf < (int)(0.3*ctrl->foldf) && EnoughMemory(ctrl, graph->nedges, wspace->nlarge>>1)) {
      IFSET(ctrl->dbglvl, DBG_PROGRESS, rprintf(ctrl, "Folding: %d\n", gnedges));
      fldgraph = FoldGraph(ctrl, graph, wspace);

      MPI_Comm_split(ctrl->comm, ctrl->mype%2 == 0, 0, &newcomm);

      if (ctrl->mype%2 == 0) { /* Even-numbered processors proceed only */
        oldcomm = ctrl->comm;
        ctrl->comm = newcomm;
        ctrl->npes = ctrl->npes>>1;
        ctrl->mype = ctrl->mype>>1;

        if (ctrl->npes%2 == 0 && ctrl->npes > 2)  /* Even # of processors and more than 2 */
          FldGlobal_Partition(ctrl, fldgraph, wspace, level+1);
        else
          Global_Partition(ctrl, fldgraph, wspace);

        ctrl->comm = oldcomm;
        ctrl->npes = ctrl->npes<<1;
        ctrl->mype = ctrl->mype<<1;
      }
      else { /* The rest just go and participate in IPART */
        switch (ctrl->ipart) {
          case IPART_SER:
            InitPartition(ctrl, graph, wspace, 0);
            break;
          case IPART_RB:
            InitPartition_RB(ctrl, graph, wspace, 0);
            break;
          default:
            errexit("Unknown IPART type!\n");
        }
      }

      MPI_Comm_free(&newcomm);

      UnFoldGraph(ctrl, graph, fldgraph, wspace);
    }
    else { /* No folding, business as usual */
      SetUp(ctrl, graph, wspace);
    
      GlobalMatch_HEM(ctrl, graph, wspace);

      FldGlobal_Partition(ctrl, graph->coarser, wspace, level);

      ProjectPartition(ctrl, graph, wspace);
      ComputePartitionParams(ctrl, graph, wspace);
      KWayRefine(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
    }
  }
}




/*************************************************************************
* This function is the driver for the partition refinement mode of ParMETIS
**************************************************************************/
void Refine_Partition(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i;

  SetUp(ctrl, graph, wspace);

  IFSET(ctrl->dbglvl, DBG_PROGRESS, rprintf(ctrl, "[%6d %8d %5d %5d][%d][%d,%d]\n", 
        graph->gnvtxs, GlobalSESum(ctrl, graph->nedges), GlobalSEMin(ctrl, graph->nvtxs),
        GlobalSEMax(ctrl, graph->nvtxs), ctrl->CoarsenTo, graph->maxvwgt,
        GlobalSEMax(ctrl, graph->vwgt[idxamax(graph->nvtxs, graph->vwgt)])));  

  if (graph->gnvtxs < 1.3*ctrl->CoarsenTo || (graph->finer != NULL && graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)) {
    /* Set the initial partition */
    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));
    graph->where = idxmalloc(graph->nvtxs+graph->nrecv, "Refine_Partition: graph->where");
    for (i=0; i<graph->nvtxs; i++)
      graph->where[i] = ctrl->mype;
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));

    if (graph->finer == NULL) { /* Do that only of no-coarsening took place */
      ComputePartitionParams(ctrl, graph, wspace);
      KWayRefine(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
    }
  }
  else { /* Coarsen it and the partition it */
    LocalMatch_HEM(ctrl, graph, wspace);

    Refine_Partition(ctrl, graph->coarser, wspace);

    ProjectPartition(ctrl, graph, wspace);
    ComputePartitionParams(ctrl, graph, wspace);
    KWayRefine(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
  }

}


/*************************************************************************
* This function is the driver for the adaptive refinement mode of ParMETIS
**************************************************************************/
void AdaptiveUndirected_Partition(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i;

  SetUp(ctrl, graph, wspace);

  IFSET(ctrl->dbglvl, DBG_PROGRESS, rprintf(ctrl, "[%6d %8d %5d %5d][%d][%d,%d]\n", 
        graph->gnvtxs, GlobalSESum(ctrl, graph->nedges), GlobalSEMin(ctrl, graph->nvtxs),
        GlobalSEMax(ctrl, graph->nvtxs), ctrl->CoarsenTo, graph->maxvwgt,
        GlobalSEMax(ctrl, graph->vwgt[idxamax(graph->nvtxs, graph->vwgt)])));  

  if (graph->gnvtxs < 1.3*ctrl->CoarsenTo || (graph->finer != NULL && graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)) {
    /* Set the initial partition */
    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));
    graph->where = idxmalloc(graph->nvtxs+graph->nrecv, "Adaptive_Partition: graph->where");
    for (i=0; i<graph->nvtxs; i++)
      graph->where[i] = ctrl->mype;
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));

    ComputePartitionParams(ctrl, graph, wspace);
    KWayAdaptiveRefineClean(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);

  }
  else { /* Coarsen it and the partition it */
    LocalMatch_HEM(ctrl, graph, wspace);

    AdaptiveUndirected_Partition(ctrl, graph->coarser, wspace);

    ProjectPartition(ctrl, graph, wspace);
    ComputePartitionParams(ctrl, graph, wspace);
    if (1.0*ctrl->nparts*graph->gpwgts[idxamax(ctrl->nparts, graph->gpwgts)]/(1.0*idxsum(ctrl->nparts, graph->gpwgts)) - UNBALANCE_FRACTION > 0.004)
      KWayAdaptiveRefineClean(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
    else
      KWayRefineClean(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
  }
}



/*************************************************************************
* This function is the driver for the adaptive refinement mode of ParMETIS
**************************************************************************/
void AdaptiveDirected_Partition(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i;

  SetUp(ctrl, graph, wspace);

  IFSET(ctrl->dbglvl, DBG_PROGRESS, rprintf(ctrl, "[%6d %8d %5d %5d][%d][%d,%d]\n", 
        graph->gnvtxs, GlobalSESum(ctrl, graph->nedges), GlobalSEMin(ctrl, graph->nvtxs),
        GlobalSEMax(ctrl, graph->nvtxs), ctrl->CoarsenTo, graph->maxvwgt,
        GlobalSEMax(ctrl, graph->vwgt[idxamax(graph->nvtxs, graph->vwgt)])));  

  if (graph->gnvtxs < 1.3*ctrl->CoarsenTo || (graph->finer != NULL && graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)) {
    /* Compute the directed diffusion */
    InitDiffusion(ctrl, graph, wspace);

    if (graph->finer == NULL) /* Do that only of no-coarsening took place */
      ComputePartitionParams(ctrl, graph, wspace);
  }
  else { /* Coarsen it and the partition it */
    LocalMatch_HEM(ctrl, graph, wspace);

    AdaptiveDirected_Partition(ctrl, graph->coarser, wspace);

    ProjectPartition(ctrl, graph, wspace);
    ComputePartitionParams(ctrl, graph, wspace);
    if (1.0*ctrl->nparts*graph->gpwgts[idxamax(ctrl->nparts, graph->gpwgts)]/(1.0*idxsum(ctrl->nparts, graph->gpwgts)) - UNBALANCE_FRACTION > 0.004)
      KWayAdaptiveRefineClean(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
    else
      KWayRefineClean(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
  }
}





/*************************************************************************
* This function is the driver for the adaptive refinement mode of ParMETIS
**************************************************************************/
void AdaptiveReMap_Partition(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i;

  SetUp(ctrl, graph, wspace);

  IFSET(ctrl->dbglvl, DBG_PROGRESS, rprintf(ctrl, "[%6d %8d %5d %5d][%d][%d,%d]\n", 
        graph->gnvtxs, GlobalSESum(ctrl, graph->nedges), GlobalSEMin(ctrl, graph->nvtxs),
        GlobalSEMax(ctrl, graph->nvtxs), ctrl->CoarsenTo, graph->maxvwgt,
        GlobalSEMax(ctrl, graph->vwgt[idxamax(graph->nvtxs, graph->vwgt)])));  

  if (graph->gnvtxs < 1.3*ctrl->CoarsenTo || (graph->finer != NULL && graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)) {
    /* Done with coarsening. Find a partition */
    switch (ctrl->ipart) {
      case IPART_SER:
        InitPartition(ctrl, graph, wspace, 1);
        break;
      case IPART_RB:
        InitPartition_RB(ctrl, graph, wspace, 1);
        break;
      default:
        errexit("Unknown IPART type!\n");
    }

    ComputePartitionParams(ctrl, graph, wspace);
    KWayRefine(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
    ReMapGraph(ctrl, graph, 1, wspace);
    GKfree(&graph->rinfo, &graph->lpwgts, &graph->gpwgts, LTERM);

    if (graph->finer == NULL) { /* Do that only of no-coarsening took place */
      ComputePartitionParams(ctrl, graph, wspace);
      if (1.0*ctrl->nparts*graph->gpwgts[idxamax(ctrl->nparts, graph->gpwgts)]/(1.0*idxsum(ctrl->nparts, graph->gpwgts)) - UNBALANCE_FRACTION > 0.004)
        KWayAdaptiveRefineClean2(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
    }
  }
  else { /* Coarsen it and the partition it */
    LocalMatch_HEM(ctrl, graph, wspace);

    AdaptiveReMap_Partition(ctrl, graph->coarser, wspace);

    ProjectPartition(ctrl, graph, wspace);
    ComputePartitionParams(ctrl, graph, wspace);
    if (1.0*ctrl->nparts*graph->gpwgts[idxamax(ctrl->nparts, graph->gpwgts)]/(1.0*idxsum(ctrl->nparts, graph->gpwgts)) - UNBALANCE_FRACTION > 0.004)
      KWayAdaptiveRefineClean2(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
    else
      KWayRefineClean(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
  }
}



/*************************************************************************
* This function is the driver for the partition refinement mode of ParMETIS
**************************************************************************/
void Order_Partition(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i;

  SetUp(ctrl, graph, wspace);

  IFSET(ctrl->dbglvl, DBG_PROGRESS, rprintf(ctrl, "[%6d %8d %5d %5d][%d][%d,%d]\n", 
        graph->gnvtxs, GlobalSESum(ctrl, graph->nedges), GlobalSEMin(ctrl, graph->nvtxs),
        GlobalSEMax(ctrl, graph->nvtxs), ctrl->CoarsenTo, graph->maxvwgt,
        GlobalSEMax(ctrl, graph->vwgt[idxamax(graph->nvtxs, graph->vwgt)])));  

  if (graph->gnvtxs < 1.3*ctrl->CoarsenTo || (graph->finer != NULL && graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)) {
    /* Compute the initial npart-way multisection */
    InitMultisection(ctrl, graph, wspace);

    if (graph->finer == NULL) { /* Do that only of no-coarsening took place */
      ComputeNodePartitionParams(ctrl, graph, wspace);
      KWayNodeRefine(ctrl, graph, wspace, 2*NGR_PASSES, ORDER_UNBALANCE_FRACTION);
    }
  }
  else { /* Coarsen it and the partition it */
    LocalMatch_HEM(ctrl, graph, wspace); 

    Order_Partition(ctrl, graph->coarser, wspace);

    ProjectPartition(ctrl, graph, wspace);
    ComputeNodePartitionParams(ctrl, graph, wspace);
    KWayNodeRefine(ctrl, graph, wspace, 2*NGR_PASSES, ORDER_UNBALANCE_FRACTION);
  }
}


/*************************************************************************
* This function returns the smaller number of vertices stored by any
* processor.
**************************************************************************/
int SmallerSubGraph(CtrlType *ctrl, GraphType *graph)
{
  int i, min;

  min = graph->vtxdist[1]-graph->vtxdist[0];
  for (i=1; i<ctrl->npes; i++) {
    if (min > graph->vtxdist[i+1]-graph->vtxdist[i])
      min = graph->vtxdist[i+1]-graph->vtxdist[i];
  }

  return min;
}
