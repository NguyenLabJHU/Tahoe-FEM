/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * proto.h
 *
 * This file contains header files
 *
 * Started 10/19/95
 * George
 *
 * $Id: proto.h,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $
 *
 */

/* kmetis.c */
void ParMETIS_PartKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void PARKMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);


/* rmetis.c */
void ParMETIS_RefineKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void PARRMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);


/* diffuse.c */
void ParMETIS_RepartLDiffusion(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_RepartGDiffusion(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void PARUAMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);
void PARDAMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);


/* scremap.c */
void ParMETIS_RepartRemap(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_RepartMLRemap(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);


/* ometis.c */
void ParMETIS_NodeND(idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, MPI_Comm *);
void PAROMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, MPI_Comm);


/* gmetis.c */
void ParMETIS_PartGraphGeomKway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_PartGraphGeomRefine(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int *, int *, int *, float *, int *, int *, idxtype *, MPI_Comm *);
void ParMETIS_PartGraphGeom(idxtype *, int *, float *, idxtype *, MPI_Comm *);
void PARGKMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, float *, idxtype *, int *, MPI_Comm);
void PARGRMETIS(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int, float *, idxtype *, int *, MPI_Comm);
void PARGMETIS(idxtype *, idxtype *, idxtype *, int, float *, idxtype *, int *, MPI_Comm);


/* pspases.c */
void ParMETIS_SerialNodeND(idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, MPI_Comm *);
GraphType *AssembleEntireGraph(CtrlType *, idxtype *, idxtype *, idxtype *);


/* debug.c */
void PrintVector(CtrlType *, int, int, idxtype *, char *);
void PrintVector2(CtrlType *, int, int, idxtype *, char *);
void PrintPairs(CtrlType *, int, KeyValueType *, char *);
void PrintGraph(CtrlType *, GraphType *);
void PrintGraph2(CtrlType *, GraphType *);
void PrintSetUpInfo(CtrlType *ctrl, GraphType *graph);
void PrintTransferedGraphs(CtrlType *, int, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void WriteMetisGraph(int, idxtype *, idxtype *, idxtype *, idxtype *);


/* comm.c */
void CommInterfaceData(CtrlType *, GraphType *, idxtype *, idxtype *, idxtype *);
void CommChangedInterfaceData(CtrlType *, GraphType *, int, idxtype *, idxtype *, KeyValueType *, KeyValueType *, idxtype *);
int GlobalSEMax(CtrlType *, int);
double GlobalSEMaxDouble(CtrlType *, double);
int GlobalSEMin(CtrlType *, int);
int GlobalSESum(CtrlType *, int);


/* io.c */
void ReadGraph(GraphType *, char *, MPI_Comm);
void ReadGraph3(GraphType *, MPI_Comm);
void ReadPartitionedGraph(GraphType *, char *, MPI_Comm);
float *ReadCoordinates(GraphType *, char *, int, MPI_Comm);
void WritePVector(char *, idxtype *, idxtype *, MPI_Comm);


/* util.c */
void errexit(char *,...);
void myprintf(CtrlType *, char *f_str,...);
void rprintf(CtrlType *, char *f_str,...);
#ifndef DMALLOC
int *imalloc(int, char *);
idxtype *idxmalloc(int, char *);
float *fmalloc(int, char *);
int *ismalloc(int, int, char *);
idxtype *idxsmalloc(int, idxtype, char *);
void *GKmalloc(int, char *);
#endif
/*void GKfree(void **,...); */
int *iset(int n, int val, int *x);
idxtype * idxset(int n, idxtype val, idxtype *x);
int idxamax(int n, idxtype *x);
int idxamin(int n, idxtype *x);
int idxasum(int n, idxtype *x);
float snorm2(int, float *);
float sdot(int n, float *, float *);
void saxpy(int, float, float *, float *);
void ikeyvalsort_org(int, KeyValueType *);
int IncKeyValueCmp(const void *, const void *);
void dkeyvalsort(int, KeyValueType *);
int DecKeyValueCmp(const void *, const void *);
int BSearch(int, idxtype *, int);
void RandomPermute(int, idxtype *, int);
void FastRandomPermute(int, idxtype *, int);
double drand48();
void srand48(long);
int ispow2(int);
int log2(int);
void BucketSortKeysDec(int, int, idxtype *, idxtype *);


/* qsort_special.c */
void iidxsort(int, idxtype *);
void iintsort(int, int *);
void ikeysort(int, KeyValueType *);
void ikeyvalsort(int, KeyValueType *);


/* memory.c */
void PreAllocateMemory(CtrlType *, GraphType *, WorkSpaceType *);
void FreeWSpace(WorkSpaceType *);
void FreeCtrl(CtrlType *);
GraphType *CreateGraph(void);
void InitGraph(GraphType *);
void FreeGraph(GraphType *);
void FreeInitialGraphAndRemap(GraphType *, int);


/* grsetup.c */
GraphType *SetUpGraph(CtrlType *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, int);
void SetUpCtrl(CtrlType *ctrl, int, int *, MPI_Comm);
void ChangeNumbering(idxtype *, idxtype *, idxtype *, idxtype *, int, int, int);
void GraphRandomPermute(GraphType *);
void ComputeMoveStatistics(CtrlType *, GraphType *, int *, int *, int *);


/* timer.c */
void InitTimers(CtrlType *);
void PrintTimingInfo(CtrlType *);
void PrintTimer(CtrlType *, timer, char *);


/* setup.c */
void SetUp(CtrlType *, GraphType *, WorkSpaceType *);
int Home_PE(int, int, idxtype *, int);


/* coarsen.c */
void GlobalMatch_HEM(CtrlType *, GraphType *, WorkSpaceType *);
void Global_CreateCoarseGraph(CtrlType *, GraphType *, WorkSpaceType *, int);
void LocalMatch_HEM(CtrlType *, GraphType *, WorkSpaceType *);
void Local_CreateCoarseGraph(CtrlType *, GraphType *, WorkSpaceType *, int);


/* edge_refine.c */
void ProjectPartition(CtrlType *, GraphType *, WorkSpaceType *);
void ComputePartitionParams(CtrlType *, GraphType *, WorkSpaceType *);
void ComputeCut(CtrlType *, GraphType *, WorkSpaceType *);
void KWayRefine(CtrlType *, GraphType *, WorkSpaceType *, int, float);
void KWayRefineClean(CtrlType *, GraphType *, WorkSpaceType *, int, float);
void KWayAdaptiveRefineClean(CtrlType *, GraphType *, WorkSpaceType *, int, float);
void KWayAdaptiveRefineClean2(CtrlType *, GraphType *, WorkSpaceType *, int, float);


/* node_refine.c */
void ComputeNodePartitionParams0(CtrlType *, GraphType *, WorkSpaceType *);
void ComputeNodePartitionParams(CtrlType *, GraphType *, WorkSpaceType *);
void KWayNodeRefine0(CtrlType *, GraphType *, WorkSpaceType *, int, float);
void KWayNodeRefine(CtrlType *, GraphType *, WorkSpaceType *, int, float);
void KWayNodeRefine2(CtrlType *, GraphType *, WorkSpaceType *, int, float);
void PrintNodeBalanceInfo(CtrlType *, int, idxtype *, idxtype *, idxtype *, int);


/* initpart.c */
void InitPartition(CtrlType *, GraphType *, WorkSpaceType *, int);
void InitPartition_RB(CtrlType *, GraphType *, WorkSpaceType *, int);
GraphType *AssembleGraph(CtrlType *, GraphType *, WorkSpaceType *, int);
void KeepPart(CtrlType *, GraphType *, WorkSpaceType *, idxtype *, int);


/* initmsection.c */
void InitMultisection(CtrlType *, GraphType *, WorkSpaceType *);
GraphType *AssembleMultisectedGraph(CtrlType *, GraphType *, WorkSpaceType *);


/* drivers.c */
void Global_Partition(CtrlType *, GraphType *, WorkSpaceType *);
void FldGlobal_Partition(CtrlType *, GraphType *, WorkSpaceType *, int);
void Refine_Partition(CtrlType *, GraphType *, WorkSpaceType *);
void AdaptiveUndirected_Partition(CtrlType *, GraphType *, WorkSpaceType *);
void AdaptiveDirected_Partition(CtrlType *, GraphType *, WorkSpaceType *);
void AdaptiveReMap_Partition(CtrlType *, GraphType *, WorkSpaceType *);
void Order_Partition(CtrlType *, GraphType *, WorkSpaceType *);
int SmallerSubGraph(CtrlType *, GraphType *);


/* fold.c */
int EnoughMemory(CtrlType *, int, int);
GraphType *FoldGraph(CtrlType *, GraphType *, WorkSpaceType *);
void UnFoldGraph(CtrlType *, GraphType *, GraphType *, WorkSpaceType *);



/* move.c */
GraphType *MoveGraph(CtrlType *, GraphType *, WorkSpaceType *);
void ProjectInfoBack(CtrlType *, GraphType *, idxtype *, idxtype *, WorkSpaceType *);
void FindVtxPerm(CtrlType *, GraphType *, idxtype *, WorkSpaceType *);


/* initdiff.c */
void InitDiffusion(CtrlType *, GraphType *, WorkSpaceType *);
GraphType *AssembleAdaptiveGraph(CtrlType *, GraphType *, WorkSpaceType *);


/* test.c */
void TestParMetis(char *, MPI_Comm);
void TestAdaptiveMETIS(idxtype *, idxtype *, idxtype *, idxtype *, int *, int, MPI_Comm);
void AdaptGraph(GraphType *, int, MPI_Comm);
void AdaptGraph2(GraphType *, int, MPI_Comm);
int ComputeRealCut(idxtype *, idxtype *, char *, MPI_Comm);
int ComputeRealCut2(idxtype *, idxtype *, idxtype *, idxtype *, char *, MPI_Comm);
void TestMoveGraph(GraphType *, GraphType *, idxtype *, MPI_Comm);


/* order.c */
void MultilevelOrder(CtrlType *, GraphType *, idxtype *, idxtype *, WorkSpaceType *);
void LabelSeparators(CtrlType *, GraphType *, idxtype *, idxtype *, idxtype *, idxtype *, WorkSpaceType *);
void CompactGraph(CtrlType *, GraphType *, idxtype *, WorkSpaceType *);
void LocalOrder(CtrlType *, GraphType *, idxtype *, int, WorkSpaceType *);
void LocalNDOrder(CtrlType *, GraphType *, idxtype *, int, WorkSpaceType *);


/* pqueue.c */
void PQueueInit(PQueueType *, int);
void PQueueReset(PQueueType *);
void PQueueFree(PQueueType *);
int PQueueInsert(PQueueType *, int, int);
int PQueueUpdate(PQueueType *, int, int, int);
int PQueueDelete(PQueueType *, int);
int PQueueGetMax(PQueueType *);
int PQueueSeeMax(PQueueType *);
int PQueueCheck(PQueueType *);


/* xyzpart.c */
void Coordinate_Partition(CtrlType *, GraphType *, int, float *, int, WorkSpaceType *);
void PartSort(CtrlType *, GraphType *, KeyValueType *, WorkSpaceType *);


/* remap.c */
void ReMapGraph(CtrlType *, GraphType *, int, WorkSpaceType *);
void ComputeTotalVReMap(CtrlType *, idxtype *, idxtype *, WorkSpaceType *);
void ComputeTotalVReMap1(CtrlType *, idxtype *, idxtype *, WorkSpaceType *);


/* wave.c */
int KWay_WavefrontDiffuser(CtrlType *, GraphType *, int, float);
int SortedKWay_WD(CtrlType *, GraphType *, int);
void SolveDiffusionEquation(GraphType *, MatrixType *, int, float *, idxtype *, idxtype *);
int SetUpLaplace(GraphType *, int, int, idxtype *, idxtype *, float *, idxtype *);
void ConjGrad(int, idxtype *, idxtype *, float *, float *, float *, float, float *);
void mvMult(int, idxtype *, idxtype *, float *, float *, float *);


/* diffutil.c */
int ComputeSerialEdgeCut(GraphType *);
int ComputeSerialTotalV(GraphType *, idxtype *);
int ComputeParallelTotalV(CtrlType *, GraphType *, idxtype *);
float ComputeParallelBalance(CtrlType *, GraphType *);
int GreaterThan(const void *, const void *);

/* tio.c */
void ReadTestGraph(GraphType *, char *, MPI_Comm);
float *ReadTestCoordinates(GraphType *, char *, int, MPI_Comm);
void ReadMetisGraph(char *, int *, idxtype **, idxtype **);
