/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * gmetis.c
 *
 * This is the entry point of parallel geometry based partitioning 
 * routines
 *
 * Started 10/19/96
 * George
 *
 * $Id: gmetis.c,v 1.1.1.1 2004-10-07 16:05:25 paklein Exp $
 *
 */

#include <parmetis.h>




/***********************************************************************************
* This function is the entry point of the parallel kmetis algorithm that uses
* coordinates to compute an initial graph distribution.
************************************************************************************/
void ParMETIS_PartGeomKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
       idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, int *nparts, 
       int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  int i, j, k, min, max, tvwgt;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph, *mgraph;

  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);

  if (*nparts == 1) { /* Take care the npes = 1 case */
    idxset(vtxdist[1], 0, part);
    *edgecut = 0;
    return;
  }
  if (npes == 1 && *nparts > 1) { /* Call METIS to do the partitioning */
    int moptions[10];

    moptions[0] = 0;
    METIS_PartGraphKway(&vtxdist[1], xadj, adjncy, vwgt, adjwgt, wgtflag, numflag,
                        nparts, moptions, edgecut, part);
    
    return;
  }

  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  SetUpCtrl(&ctrl, npes, options, *comm);
  ctrl.CoarsenTo = amin(vtxdist[npes]+1, 25*amax(npes, *nparts));

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt, *wgtflag);

  PreAllocateMemory(&ctrl, graph, &wspace);

  /*=================================================================
   * Compute the initial npes-way partitioning geometric partitioning 
   =================================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  Coordinate_Partition(&ctrl, graph, *ndims, xyz, 1, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));

  /*=================================================================
   * Move the graph according to the partitioning
   =================================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.MoveTmr));

  mgraph = MoveGraph(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.MoveTmr));

  if (ctrl.dbglvl&DBG_INFO) {
    tvwgt = GlobalSESum(&ctrl, idxsum(graph->nvtxs, graph->vwgt));
    ComputePartitionParams(&ctrl, graph, &wspace);
    rprintf(&ctrl, "XYZ Cut: %6d \tBalance: %6.3f [%d %d %d][%d %d %d %d]\n", 
          graph->mincut, 1.0*npes*graph->gpwgts[idxamax(npes, graph->gpwgts)]/(1.0*tvwgt), 
          graph->gpwgts[idxamax(npes, graph->gpwgts)], tvwgt, graph->gnvtxs,
          GlobalSEMax(&ctrl, graph->nrecv), GlobalSESum(&ctrl, graph->nrecv), 
          GlobalSEMax(&ctrl, graph->nsend), GlobalSESum(&ctrl, graph->nsend));
  }

  /*=======================================================
   * Now compute the partition of the moved graph
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  FreeWSpace(&wspace);

  ctrl.nparts = *nparts;
  PreAllocateMemory(&ctrl, mgraph, &wspace);

  mgraph->maxvwgt = graph->maxvwgt;

  if (npes%2 == 0 && npes > 2)
    FldGlobal_Partition(&ctrl, mgraph, &wspace, 0);
  else
    Global_Partition(&ctrl, mgraph, &wspace);

  ctrl.nparts = npes;

  /* Invert the ordering back to the original graph */
  ProjectInfoBack(&ctrl, graph, part, mgraph->where, &wspace);

  *edgecut = mgraph->mincut;

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));

  IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "Final %d-way Cut: %6d \tBalance: %6.3f [%d %d %d][%d %d %d %d]\n", 
          *nparts, mgraph->mincut, 1.0*(*nparts)*mgraph->gpwgts[idxamax(*nparts, mgraph->gpwgts)]/(1.0*tvwgt), 
          mgraph->gpwgts[idxamax(*nparts, mgraph->gpwgts)], tvwgt, mgraph->gnvtxs, 
          GlobalSEMax(&ctrl, mgraph->nrecv), GlobalSESum(&ctrl, mgraph->nrecv), 
          GlobalSEMax(&ctrl, mgraph->nsend), GlobalSESum(&ctrl, mgraph->nsend)));

  FreeGraph(mgraph);
  FreeInitialGraphAndRemap(graph, *wgtflag);
  FreeWSpace(&wspace);
  FreeCtrl(&ctrl);

  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

}




/***********************************************************************************
* This function is the entry point of the parallel rmetis algorithm that uses
* coordinates to compute an initial graph distribution.
************************************************************************************/
void ParMETIS_PartGeomRefine(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, 
       idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, 
       float *xyz, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  int i, j, k, min, max, tvwgt;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph, *mgraph;

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
  ctrl.CoarsenTo = amin(vtxdist[npes]+1, 25*npes);

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt, *wgtflag);

  PreAllocateMemory(&ctrl, graph, &wspace);

  /*=======================================================
   * Compute the initial geometric partitioning 
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  Coordinate_Partition(&ctrl, graph, *ndims, xyz, 1, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));

  /*=======================================================
   * Move the graph according to the partitioning
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.MoveTmr));

  mgraph = MoveGraph(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.MoveTmr));

  if (ctrl.dbglvl&DBG_INFO) {
    tvwgt = GlobalSESum(&ctrl, idxsum(graph->nvtxs, graph->vwgt));
    ComputePartitionParams(&ctrl, graph, &wspace);
    rprintf(&ctrl, "XYZ Cut: %6d \tBalance: %6.3f [%d %d %d][%d %d %d %d]\n", 
          graph->mincut, 1.0*npes*graph->gpwgts[idxamax(npes, graph->gpwgts)]/(1.0*tvwgt), 
          graph->gpwgts[idxamax(npes, graph->gpwgts)], tvwgt, graph->gnvtxs,
          GlobalSEMax(&ctrl, graph->nrecv), GlobalSESum(&ctrl, graph->nrecv), 
          GlobalSEMax(&ctrl, graph->nsend), GlobalSESum(&ctrl, graph->nsend));
  }

  /*=======================================================
   * Now compute the partition of the moved graph
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  FreeWSpace(&wspace);
  PreAllocateMemory(&ctrl, mgraph, &wspace);

  mgraph->maxvwgt = graph->maxvwgt;

  Refine_Partition(&ctrl, mgraph, &wspace);

  /* Invert the ordering back to the original graph */
  ProjectInfoBack(&ctrl, graph, part, mgraph->where, &wspace);

  *edgecut = mgraph->mincut;

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));

  IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "Final Cut: %6d \tBalance: %6.3f [%d %d %d][%d %d %d %d]\n", 
          mgraph->mincut, 1.0*npes*mgraph->gpwgts[idxamax(npes, mgraph->gpwgts)]/(1.0*tvwgt), 
          mgraph->gpwgts[idxamax(npes, mgraph->gpwgts)], tvwgt, mgraph->gnvtxs, 
          GlobalSEMax(&ctrl, mgraph->nrecv), GlobalSESum(&ctrl, mgraph->nrecv), 
          GlobalSEMax(&ctrl, mgraph->nsend), GlobalSESum(&ctrl, mgraph->nsend)));

  FreeGraph(mgraph);
  FreeInitialGraphAndRemap(graph, *wgtflag);
  FreeWSpace(&wspace);
  FreeCtrl(&ctrl);

  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

}




/***********************************************************************************
* This function is the entry point of the parallel ordering algorithm.
* This function assumes that the graph is already nice partitioned among the 
* processors and then proceeds to perform recursive bisection.
************************************************************************************/
void ParMETIS_PartGeom(idxtype *vtxdist, int *ndims, float *xyz, idxtype *part, MPI_Comm *comm)
{
  int i, j, k, min, max, tvwgt, npes, mype, nvtxs, firstvtx, options[10];
  idxtype *xadj, *adjncy;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph;

  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);

  if (npes == 1) { /* Take care the npes = 1 case */
    idxset(vtxdist[1], 0, part);
    return;
  }
  
  /* Setup a fake graph to allow the rest of the code to work unchanged */
  options[0] = 0;

  nvtxs = vtxdist[mype+1]-vtxdist[mype];
  firstvtx = vtxdist[mype];
  xadj = idxmalloc(nvtxs+1, "ParMETIS_PartGeom: xadj");
  adjncy = idxmalloc(nvtxs, "ParMETIS_PartGeom: adjncy");
  for (i=0; i<nvtxs; i++) {
    xadj[i] = i;
    adjncy[i] = firstvtx + (i+1)%nvtxs;
  }
  xadj[nvtxs] = nvtxs;


  /* Proceed with the rest of the code */
  SetUpCtrl(&ctrl, npes, options, *comm);
  ctrl.CoarsenTo = amin(vtxdist[npes]+1, 25*npes);

  graph = SetUpGraph(&ctrl, vtxdist, xadj, NULL, adjncy, NULL, 0);

  PreAllocateMemory(&ctrl, graph, &wspace);

  /*=======================================================
   * Compute the initial geometric partitioning 
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  Coordinate_Partition(&ctrl, graph, *ndims, xyz, 0, &wspace);

  idxcopy(graph->nvtxs, graph->where, part);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));

  FreeInitialGraphAndRemap(graph, 0);
  FreeWSpace(&wspace);
  FreeCtrl(&ctrl);

  GKfree(&xadj, &adjncy, LTERM);
}


/***********************************************************************************
* This function is the entry point of the parallel kmetis algorithm that uses
* coordinates to compute an initial graph distribution.
************************************************************************************/
void PARGKMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
                int ndims, float *xyz, idxtype *part, int *options, MPI_Comm comm)
{
  int npes, wgtflag, numflag, edgecut, newoptions[5];

  MPI_Comm_size(comm, &npes);

  newoptions[0] = 1;
  newoptions[OPTION_IPART] = options[2];
  newoptions[OPTION_FOLDF] = options[1];
  newoptions[OPTION_DBGLVL] = options[4];

  numflag = options[3];
  wgtflag = (vwgt == NULL ? 0 : 2) + (adjwgt == NULL ? 0 : 1);

  ParMETIS_PartGeomKway(vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, 
     &ndims, xyz, &npes, newoptions, &edgecut, part, &comm);

  options[0] = edgecut;

}


/***********************************************************************************
* This function is the entry point of the parallel rmetis algorithm that uses
* coordinates to compute an initial graph distribution.
************************************************************************************/
void PARGRMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
                int ndims, float *xyz, idxtype *part, int *options, MPI_Comm comm)
{
  int wgtflag, numflag, edgecut, newoptions[5];

  newoptions[0] = 1;
  newoptions[OPTION_IPART] = options[2];
  newoptions[OPTION_FOLDF] = options[1];
  newoptions[OPTION_DBGLVL] = options[4];

  numflag = options[3];
  wgtflag = (vwgt == NULL ? 0 : 2) + (adjwgt == NULL ? 0 : 1);

  ParMETIS_PartGeomRefine(vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, 
     &ndims, xyz, newoptions, &edgecut, part, &comm);

  options[0] = edgecut;

}

/***********************************************************************************
* This function is the entry point of the parallel ordering algorithm.
* This function assumes that the graph is already nice partitioned among the 
* processors and then proceeds to perform recursive bisection.
************************************************************************************/
void PARGMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int ndims, float *xyz, 
               idxtype *part, int *options, MPI_Comm comm)
{

  ParMETIS_PartGeom(vtxdist, &ndims, xyz, part, &comm);

  options[0] = -1;

}
