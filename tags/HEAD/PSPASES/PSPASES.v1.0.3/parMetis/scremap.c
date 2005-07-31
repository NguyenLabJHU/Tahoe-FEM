/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * scremap.c
 *
 * This is the entry point of parallel scratch-remap repartitioning routines
 *
 * Started 10/19/96
 * George
 *
 * $Id: scremap.c,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $
 *
 */

#include <parmetis.h>




/***********************************************************************************
* This function is the entry point of the parallel repartitioning algortihm that
* partitions from scratch. 
************************************************************************************/
void ParMETIS_RepartRemap(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, 
       idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options,
       int *edgecut, idxtype *part, MPI_Comm *comm)
{
  int i, j, k, min, max, tvwgt;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph;
  int nmoved, maxin, maxout;

  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);

  /* Take care the npes = 1 case */
  if (npes == 1) { 
    idxset(vtxdist[1], 0, part);
    *edgecut = 0;
    return;
  }

  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  SetUpCtrl(&ctrl, npes, options, *comm);
  ctrl.CoarsenTo = amax(1000, amin(vtxdist[npes]+1, 25*npes));
  ctrl.ipart = IPART_RB;

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt, *wgtflag);
  graph->vsize = idxsmalloc(graph->nvtxs, 1, "Par_KMetis: vsize");

  PreAllocateMemory(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  Global_Partition(&ctrl, graph, &wspace);
  ReMapGraph(&ctrl, graph, 0, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));

  idxcopy(graph->nvtxs, graph->where, part);
  *edgecut = graph->mincut;

  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_INFO, ComputeMoveStatistics(&ctrl, graph, &nmoved, &maxin, &maxout));
  IFSET(ctrl.dbglvl, DBG_INFO, tvwgt = GlobalSESum(&ctrl, idxsum(graph->nvtxs, graph->vwgt)));
  IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "Final Cut: %6d \tBalance: %6.3f \nNMoved: %d %d %d %d [%d %d %d %d]\n", 
          graph->mincut, 1.0*npes*graph->gpwgts[idxamax(npes, graph->gpwgts)]/(1.0*tvwgt),
          nmoved, maxin, maxout, maxin+maxout, npes, graph->gpwgts[idxamax(npes, graph->gpwgts)], tvwgt, graph->gnvtxs));

  GKfree(&graph->vsize, LTERM);
  FreeInitialGraphAndRemap(graph, *wgtflag);
  FreeWSpace(&wspace);
  FreeCtrl(&ctrl);

  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

}




/***********************************************************************************
* This function is the entry point of the parallel repartitioning algortihm that
* partitions from scratch. It is for interanl use only!
************************************************************************************/
void ParMETIS_RepartMLRemap(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, 
       idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options,
       int *edgecut, idxtype *part, MPI_Comm *comm)
{
  int i, j, k, min, max, tvwgt;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph;
  int nmoved, maxin, maxout;

  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);

  /* Take care the npes = 1 case */
  if (npes == 1) { 
    idxset(vtxdist[1], 0, part);
    *edgecut = 0;
    return;
  }

  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  SetUpCtrl(&ctrl, npes, options, *comm);
  ctrl.ipart = IPART_RB;
  ctrl.CoarsenTo = amax(1000, amin(vtxdist[npes]+1, 25*npes));

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt, *wgtflag);
  graph->vsize = idxsmalloc(graph->nvtxs, 1, "Par_KMetis: vsize");

  PreAllocateMemory(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  AdaptiveReMap_Partition(&ctrl, graph, &wspace);
  ReMapGraph(&ctrl, graph, 0, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));

  idxcopy(graph->nvtxs, graph->where, part);
  *edgecut = graph->mincut;

  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_INFO, ComputeMoveStatistics(&ctrl, graph, &nmoved, &maxin, &maxout));
  IFSET(ctrl.dbglvl, DBG_INFO, tvwgt = GlobalSESum(&ctrl, idxsum(graph->nvtxs, graph->vwgt)));
  IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "Final Cut: %6d \tBalance: %6.3f \nNMoved: %d %d %d %d [%d %d %d %d]\n", 
          graph->mincut, 1.0*npes*graph->gpwgts[idxamax(npes, graph->gpwgts)]/(1.0*tvwgt),
          nmoved, maxin, maxout, maxin+maxout, npes, graph->gpwgts[idxamax(npes, graph->gpwgts)], tvwgt, graph->gnvtxs));

  GKfree(&graph->vsize, LTERM);
  FreeInitialGraphAndRemap(graph, *wgtflag);
  FreeWSpace(&wspace);
  FreeCtrl(&ctrl);

  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

}

