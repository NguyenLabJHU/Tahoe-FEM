/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * diffuse.c
 *
 * This is the entry point of parallel difussive repartitioning routines
 *
 * Started 10/19/96
 * George
 *
 * $Id: diffuse.c,v 1.1.1.1 2004-10-07 16:05:25 paklein Exp $
 *
 */

#include <parmetis.h>



/***********************************************************************************
* This function is the entry point of the parallel multilevel local diffusion
* algorithm. It uses parallel undirected diffusion followed by adaptive k-way 
* refinement. This function utilizes local coarsening.
************************************************************************************/
void ParMETIS_RepartLDiffusion(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, 
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

  if (npes == 1) { /* Take care the npes = 1 case */
    idxset(vtxdist[1], 0, part);
    *edgecut = 0;
    return;
  }

  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  SetUpCtrl(&ctrl, npes, options, *comm);
  ctrl.CoarsenTo = amin(vtxdist[npes]+1, 70*npes);

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt, *wgtflag);
  graph->vsize = idxsmalloc(graph->nvtxs, 1, "Par_KMetis: vsize");

  PreAllocateMemory(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  AdaptiveUndirected_Partition(&ctrl, graph, &wspace);
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
* This function is the entry point of the parallel multilevel global diffusion
* algorithm. It uses parallel undirected diffusion followed by adaptive k-way 
* refinement. This function utilizes local coarsening.
************************************************************************************/
void ParMETIS_RepartGDiffusion(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, 
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

  if (npes == 1) { /* Take care the npes = 1 case */
    idxset(vtxdist[1], 0, part);
    *edgecut = 0;
    return;
  }

  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  SetUpCtrl(&ctrl, npes, options, *comm);
  ctrl.CoarsenTo = amin(vtxdist[npes]+1, 50*npes);

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt, *wgtflag);
  graph->vsize = idxsmalloc(graph->nvtxs, 1, "Par_KMetis: vsize");

  PreAllocateMemory(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  AdaptiveDirected_Partition(&ctrl, graph, &wspace);
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
* This function is the entry point of the parallel multilevel undirected diffusion
* algorithm. It uses parallel undirected diffusion followed by adaptive k-way 
* refinement. This function utilizes local coarsening.
************************************************************************************/
void PARUAMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
               idxtype *part, int *options, MPI_Comm comm)
{
  int wgtflag, numflag, edgecut, newoptions[5];

  newoptions[0] = 1;
  newoptions[OPTION_IPART] = options[2];
  newoptions[OPTION_FOLDF] = options[1];
  newoptions[OPTION_DBGLVL] = options[4];

  numflag = options[3];
  wgtflag = (vwgt == NULL ? 0 : 2) + (adjwgt == NULL ? 0 : 1);

  ParMETIS_RepartLDiffusion(vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, 
     newoptions, &edgecut, part, &comm);

  options[0] = edgecut;

}



/***********************************************************************************
* This function is the entry point of the parallel multilevel directed diffusion
* algorithm. It uses parallel undirected diffusion followed by adaptive k-way 
* refinement. This function utilizes local coarsening.
************************************************************************************/
void PARDAMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
               idxtype *part, int *options, MPI_Comm comm)
{
  int wgtflag, numflag, edgecut, newoptions[5];

  newoptions[0] = 1;
  newoptions[OPTION_IPART] = options[2];
  newoptions[OPTION_FOLDF] = options[1];
  newoptions[OPTION_DBGLVL] = options[4];

  numflag = options[3];
  wgtflag = (vwgt == NULL ? 0 : 2) + (adjwgt == NULL ? 0 : 1);

  ParMETIS_RepartGDiffusion(vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, 
     newoptions, &edgecut, part, &comm);

  options[0] = edgecut;

}

