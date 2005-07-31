/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * rename.h
 *
 * This file contains renaming #defines to remove any conflicts that the
 * library may have with the users functions.
 *
 * Started 7/17/97
 * George
 *
 * $Id: rename.h,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $
 */

/* coarsen.c */
#define GlobalMatch_HEM GlobalMatch_HEM__
#define Global_CreateCoarseGraph Global_CreateCoarseGraph__
#define LocalMatch_HEM LocalMatch_HEM__
#define Local_CreateCoarseGraph Local_CreateCoarseGraph__

/* comm.c */
#define CommInterfaceData CommInterfaceData__
#define CommChangedInterfaceData CommChangedInterfaceData__
#define GlobalSEMax GlobalSEMax__
#define GlobalSEMaxDouble GlobalSEMaxDouble__
#define GlobalSEMin GlobalSEMin__
#define GlobalSESum GlobalSESum__

/* debug.c */
#define PrintVector PrintVector__
#define PrintVector2 PrintVector2__
#define PrintPairs PrintPairs__
#define PrintGraph PrintGraph__
#define PrintGraph2 PrintGraph2__
#define PrintSetUpInfo PrintSetUpInfo__
#define PrintTransferedGraphs PrintTransferedGraphs__

/* diffuse.c */

/* diffutil.c */
#define ComputeSerialEdgeCut ComputeSerialEdgeCut__
#define ComputeSerialTotalV ComputeSerialTotalV__
#define ComputeParallelTotalV ComputeParallelTotalV__
#define ComputeParallelBalance ComputeParallelBalance__
#define GreaterThan GreaterThan__
#define ConjGrad ConjGrad__
#define mvMult mvMult__

/* drivers.c */
#define Global_Partition Global_Partition__
#define FldGlobal_Partition FldGlobal_Partition__
#define Refine_Partition Refine_Partition__
#define AdaptiveUndirected_Partition AdaptiveUndirected_Partition__
#define AdaptiveDirected_Partition AdaptiveDirected_Partition__
#define AdaptiveReMap_Partition AdaptiveReMap_Partition__
#define Order_Partition Order_Partition__
#define SmallerSubGraph SmallerSubGraph__

/* edge_refine.c */
#define ProjectPartition ProjectPartition__
#define ComputePartitionParams ComputePartitionParams__
#define ComputeCut ComputeCut__
#define KWayRefine KWayRefine__
#define KWayRefineClean KWayRefineClean__
#define KWayAdaptiveRefineClean KWayAdaptiveRefineClean__
#define KWayAdaptiveRefineClean2 KWayAdaptiveRefineClean2__

/* fold.c */
#define EnoughMemory EnoughMemory__
#define FoldGraph FoldGraph__
#define UnFoldGraph UnFoldGraph__

/* frename.c */

/* gmetis.c */

/* grsetup.c */
#define SetUpGraph SetUpGraph__
#define SetUpCtrl SetUpCtrl__
#define ChangeNumbering ChangeNumbering__
#define GraphRandomPermute GraphRandomPermute__
#define ComputeMoveStatistics ComputeMoveStatistics__

/* iidxsort.c */
#define iidxsort iidxsort__

/* iintsort.c */
#define iintsort iintsort__

/* ikeysort.c */
#define ikeysort ikeysort__

/* ikeyvalsort.c */
#define ikeyvalsort ikeyvalsort__

/* initdiff.c */
#define InitDiffusion InitDiffusion__
#define AssembleAdaptiveGraph AssembleAdaptiveGraph__


/* initmsection.c */
#define InitMultisection InitMultisection__
#define AssembleMultisectedGraph AssembleMultisectedGraph__

/* initpart.c */
#define InitPartition InitPartition__
#define InitPartition_RB InitPartition_RB__
#define AssembleGraph AssembleGraph__
#define KeepPart KeepPart__

/* kmetis.c */

/* memory.c */
#define PreAllocateMemory PreAllocateMemory__
#define FreeWSpace FreeWSpace__
#define FreeCtrl FreeCtrl__
#define CreateGraph CreateGraph__
#define InitGraph InitGraph__
#define FreeGraph FreeGraph__
#define FreeGraphContent FreeGraphContent__
#define FreeInitialGraphAndRemap FreeInitialGraphAndRemap__

/* move.c */
#define MoveGraph MoveGraph__
#define ProjectInfoBack ProjectInfoBack__
#define FindVtxPerm FindVtxPerm__


/* node_refine.c */
#define ComputeNodePartitionParams ComputeNodePartitionParams__
#define KWayNodeRefine KWayNodeRefine__
#define PrintNodeBalanceInfo PrintNodeBalanceInfo__

/* ometis.c */

/* order.c */
#define MultilevelOrder MultilevelOrder__
#define LabelSeparators LabelSeparators__
#define CompactGraph CompactGraph__
#define LocalNDOrder LocalNDOrder__

/* pqueue.c */
#define PQueueInit PQueueInit__
#define PQueueReset PQueueReset__
#define PQueueFree PQueueFree__
#define PQueueInsert PQueueInsert__
#define PQueueUpdate PQueueUpdate__
#define PQueueDelete PQueueDelete__
#define PQueueGetMax PQueueGetMax__
#define PQueueSeeMax PQueueSeeMax__
#define PQueueCheck PQueueCheck__

/* pspases.c */
#define AssembleEntireGraph AssembleEntireGraph__

/* remap.c */
#define ReMapGraph ReMapGraph__
#define ComputeTotalVReMap ComputeTotalVReMap__
#define ComputeTotalVReMap1 ComputeTotalVReMap1__

/* rmetis.c */

/* scremap.c */


/* setup.c */
#define SetUp SetUp__

/* timer.c */
#define InitTimers InitTimers__
#define PrintTimingInfo PrintTimingInfo__
#define PrintTimer PrintTimer__

/* util.c */
#define errexit errexit__
#define myprintf myprintf__
#define rprintf rprintf__
#ifndef DMALLOC
#define imalloc imalloc__
#define idxmalloc idxmalloc__
#define fmalloc fmalloc__
#define ismalloc ismalloc__
#define idxsmalloc idxsmalloc__
#define GKmalloc GKmalloc__
#endif
#define GKfree GKfree__
#define iset iset__
#define idxset idxset__
#define idxamax idxamax__
#define idxamin idxamin__
#define idxasum idxasum__
#define charsum charsum__
#define isum isum__
#define snorm2 snorm2__
#define sdot sdot__
#define saxpy saxpy__
#define ikeyvalsort_org ikeyvalsort_org__
#define IncKeyValueCmp IncKeyValueCmp__
#define dkeyvalsort dkeyvalsort__
#define DecKeyValueCmp DecKeyValueCmp__
#define BSearch BSearch__
#define RandomPermute RandomPermute__
#define FastRandomPermute FastRandomPermute__
#define ispow2 ispow2__
#define log2 log2_

/* wave.c */
#define KWay_WavefrontDiffuser KWay_WavefrontDiffuser__
#define SortedKWay_WD SortedKWay_WD__
#define SolveDiffusionEquation SolveDiffusionEquation__

/* xyzpart.c */
#define Coordinate_Partition Coordinate_Partition__
#define PartSort PartSort__
