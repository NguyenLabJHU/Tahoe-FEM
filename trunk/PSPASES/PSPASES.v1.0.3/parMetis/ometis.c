/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kmetis.c
 *
 * This is the entry point of parallel ordering
 *
 * Started 10/19/96
 * George
 *
 * $Id: ometis.c,v 1.1.1.1 2004-10-07 16:05:25 paklein Exp $
 *
 */

#include <parmetis.h>




/***********************************************************************************
* This function is the entry point of the parallel ordering algorithm.
* This function assumes that the graph is already nice partitioned among the 
* processors and then proceeds to perform recursive bisection.
************************************************************************************/
void ParMETIS_NodeND(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag,
        int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm)
{
  int i, j, k, min, max, nparts;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph, *mgraph;
  idxtype *morder;

  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);

  if (!ispow2(npes)) {
    if (mype == 0)
      printf("Error: The number of processors must be a power of 2!\n");
    return;
  }

  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, order, npes, mype, 1);

  SetUpCtrl(&ctrl, npes, options, *comm);
  ctrl.ipart = IPART_RB;
  ctrl.CoarsenTo = amin(vtxdist[npes]+1, 25*npes);

  graph = SetUpGraph(&ctrl, vtxdist, xadj, NULL, adjncy, NULL, 0);

  PreAllocateMemory(&ctrl, graph, &wspace);

  /*=======================================================
   * Compute the initial k-way partitioning 
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  if (npes > 2)
    FldGlobal_Partition(&ctrl, graph, &wspace, 0);
  else
    Global_Partition(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));

  /*=======================================================
   * Move the graph according to the partitioning
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.MoveTmr));

  MALLOC_CHECK(NULL);
  mgraph = MoveGraph(&ctrl, graph, &wspace);
  MALLOC_CHECK(NULL);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.MoveTmr));

  /*=======================================================
   * Now compute an ordering of the moved graph
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  FreeWSpace(&wspace);
  PreAllocateMemory(&ctrl, mgraph, &wspace);

  ctrl.ipart = options[OPTION_IPART];
  ctrl.CoarsenTo = amin(vtxdist[npes]+1, amax(20*npes, 1000));
  mgraph->maxvwgt = graph->maxvwgt;

  morder = idxmalloc(mgraph->nvtxs, "PAROMETIS: morder");
  MultilevelOrder(&ctrl, mgraph, morder, sizes, &wspace);

  MALLOC_CHECK(NULL);

  /* Invert the ordering back to the original graph */
  ProjectInfoBack(&ctrl, graph, order, morder, &wspace);

  MALLOC_CHECK(NULL);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));

  free(morder);
  FreeGraph(mgraph);
  FreeInitialGraphAndRemap(graph, 0);
  FreeWSpace(&wspace);
  FreeCtrl(&ctrl);

  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, order, npes, mype, 0);

  MALLOC_CHECK(NULL);
}


/***********************************************************************************
* This function is the entry point of the parallel ordering algorithm.
* This function assumes that the graph is already nice partitioned among the 
* processors and then proceeds to perform recursive bisection.
************************************************************************************/
void PAROMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
                idxtype *order, idxtype *sizes, int *options, MPI_Comm comm)
{
  int numflag, newoptions[5];

  newoptions[0] = 1;
  newoptions[OPTION_IPART] = options[2];
  newoptions[OPTION_FOLDF] = options[1];
  newoptions[OPTION_DBGLVL] = options[4];

  numflag = options[3];

  ParMETIS_NodeND(vtxdist, xadj, adjncy, &numflag, newoptions, order, sizes, &comm);

  options[0] = -1;

}

