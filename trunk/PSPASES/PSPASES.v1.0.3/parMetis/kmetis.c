/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kmetis.c
 *
 * This is the entry point of PARMETIS_PartGraphKway
 *
 * Started 10/19/96
 * George
 *
 * $Id: kmetis.c,v 1.1.1.1 2004-10-07 16:05:25 paklein Exp $
 *
 */

#include <parmetis.h>


/***********************************************************************************
* This function is the entry point of the parallel k-way multilevel partitionioner. 
* This function assumes nothing about the graph distribution.
* It is the general case.
************************************************************************************/
void ParMETIS_PartKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
       idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, 
       idxtype *part, MPI_Comm *comm)
{
  int i, j, k, min, max, tvwgt;
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph;

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

  SetUpCtrl(&ctrl, *nparts, options, *comm);
  ctrl.CoarsenTo = amin(vtxdist[npes]+1, 25*amax(npes, *nparts));

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt, *wgtflag);

  PreAllocateMemory(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  if (npes%2 == 0 && npes > 2)
    FldGlobal_Partition(&ctrl, graph, &wspace, 0);
  else
    Global_Partition(&ctrl, graph, &wspace);

  ReMapGraph(&ctrl, graph, 0, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));

  idxcopy(graph->nvtxs, graph->where, part);
  *edgecut = graph->mincut;

  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_INFO, tvwgt = GlobalSESum(&ctrl, idxsum(graph->nvtxs, graph->vwgt)));
  IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "Final %d-way Cut: %6d \tBalance: %6.3f [%d %d %d]\n", 
          *nparts, graph->mincut, 1.0*(*nparts)*graph->gpwgts[idxamax(*nparts, graph->gpwgts)]/(1.0*tvwgt), 
          graph->gpwgts[idxamax(*nparts, graph->gpwgts)], tvwgt, graph->gnvtxs));

  FreeInitialGraphAndRemap(graph, *wgtflag);
  FreeWSpace(&wspace);
  FreeCtrl(&ctrl);

  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

}





/***********************************************************************************
* This function is the entry point of the parallel k-way multilevel partitionioner. 
* This function assumes nothing about the graph distribution.
* It is the general case.
************************************************************************************/
void PARKMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
               idxtype *part, int *options, MPI_Comm comm)
{
  int wgtflag, numflag, edgecut, newoptions[5];
  int npes;

  MPI_Comm_size(comm, &npes);

  newoptions[0] = 1;
  newoptions[OPTION_IPART] = options[2];
  newoptions[OPTION_FOLDF] = options[1];
  newoptions[OPTION_DBGLVL] = options[4];

  numflag = options[3];
  wgtflag = (vwgt == NULL ? 0 : 2) + (adjwgt == NULL ? 0 : 1);

  ParMETIS_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, &npes,
           newoptions, &edgecut, part, &comm);

  options[0] = edgecut;

}


