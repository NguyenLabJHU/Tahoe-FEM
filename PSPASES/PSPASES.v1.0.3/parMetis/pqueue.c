/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * pqueue.c
 *
 * This file contains functions for manipulating the priority queues.
 *
 * Started 9/2/94
 * George
 *
 * $Id: pqueue.c,v 1.1.1.1 2004-10-07 16:05:26 paklein Exp $
 *
 */

#include "parmetis.h"


/*************************************************************************
* This function initializes the klpartdef data structures
**************************************************************************/
void PQueueInit(PQueueType *queue, int maxnnodes)
{
  queue->nnodes = 0;
  queue->maxnnodes = maxnnodes;

  queue->iperm = idxsmalloc(maxnnodes, -1, "PQueueInit: iperm");
  queue->perm = idxmalloc(maxnnodes, "PQueueInit: perm");
  queue->values = idxmalloc(maxnnodes, "PQueueInit: values");

}


/*************************************************************************
* This function resets the buckets
**************************************************************************/
void PQueueReset(PQueueType *queue)
{
  queue->nnodes = 0;
  idxset(queue->maxnnodes, -1, queue->iperm);
}


/*************************************************************************
* This function frees the buckets
**************************************************************************/
void PQueueFree(PQueueType *queue)
{
  queue->maxnnodes = 0;

  GKfree(&queue->iperm, &queue->perm, &queue->values, LTERM);
}


/*************************************************************************
* This function adds a node of certain gain into a partition
**************************************************************************/
int PQueueInsert(PQueueType *queue, int node, int gain)
{
  int i, j;
  idxtype *iperm, *perm, *values;

  iperm = queue->iperm;
  perm = queue->perm;
  values = queue->values;

  ASSERTSP(iperm[node] == -1, ("\t%d %d\n", node, iperm[node]));

  i = queue->nnodes++;
  while (i > 0) {
    j = (i-1)/2;
    if (values[j] < gain) {
      values[i] = values[j];
      perm[i] = perm[j];
      iperm[perm[i]] = i;
      i = j;
    }
    else 
      break;
  }
  ASSERTS(i >= 0);
  values[i] = gain;
  perm[i] = node;
  iperm[node] = i;

  ASSERTS(PQueueCheck(queue));

  return 0;
}


/*************************************************************************
* This function deletes a node from a partition and reinserts it with
* an updated gain
**************************************************************************/
int PQueueUpdate(PQueueType *queue, int node, int oldgain, int newgain)
{
  int snode, tmp, i, j;
  idxtype *iperm, *perm, *values;

  if (oldgain == newgain) 
    return 0;

  iperm = queue->iperm;
  perm = queue->perm;
  values = queue->values;

  i = iperm[node];
  ASSERTSP(perm[i] == node, ("\t%d %d %d\n", node, iperm[node], perm[iperm[node]]));

  if (oldgain < newgain) { /* Filter-up */
    while (i > 0) {
      j = (i-1)/2;
      if (values[j] < newgain) {
        values[i] = values[j];
        perm[i] = perm[j];
        iperm[perm[i]] = i;
        i = j;
      }
      else 
        break;
    }
  }
  else { /* Filter down */
    while ((j=2*i+1) < queue->nnodes) {
      if (values[j] > newgain) {
        if (j+1 < queue->nnodes && values[j+1] > values[j])
          j = j+1;
        values[i] = values[j];
        perm[i] = perm[j];
        iperm[perm[i]] = i;
        i = j;
      }
      else if (j+1 < queue->nnodes && values[j+1] > newgain) {
        j = j+1;
        values[i] = values[j];
        perm[i] = perm[j];
        iperm[perm[i]] = i;
        i = j;
      }
      else
        break;
    }
  }

  values[i] = newgain;
  perm[i] = node;
  iperm[node] = i;

  ASSERTS(PQueueCheck(queue));

  return 0;
}



/*************************************************************************
* This function deletes a node from a priority queue 
**************************************************************************/
int PQueueDelete(PQueueType *queue, int node)
{
  int snode, tmp, i, j, newgain, oldgain;
  idxtype *iperm, *perm, *values;

  if (queue->nnodes == 0)
    return -1;

  iperm = queue->iperm;
  perm = queue->perm;
  values = queue->values;

  i = iperm[node];
  iperm[node] = -1;
  ASSERTS(perm[i] == node);

  queue->nnodes--; 
  node = perm[queue->nnodes];
  newgain = values[queue->nnodes];
  oldgain = values[i];

  if (queue->nnodes > 0) {
    if (oldgain < newgain) { /* Filter-up */
      while (i > 0) {
        j = (i-1)/2;
        if (values[j] < newgain) {
          values[i] = values[j];
          perm[i] = perm[j];
          iperm[perm[i]] = i;
          i = j;
        }
        else 
          break;
      }
    }
    else { /* Filter down */
      while ((j=2*i+1) < queue->nnodes) {
        if (values[j] > newgain) {
          if (j+1 < queue->nnodes && values[j+1] > values[j])
            j = j+1;
          values[i] = values[j];
          perm[i] = perm[j];
          iperm[perm[i]] = i;
          i = j;
        }
        else if (j+1 < queue->nnodes && values[j+1] > newgain) {
          j = j+1;
          values[i] = values[j];
          perm[i] = perm[j];
          iperm[perm[i]] = i;
          i = j;
        }
        else
          break;
      }
    }

    values[i] = newgain;
    perm[i] = node;
    iperm[node] = i;
  }

  ASSERTS(PQueueCheck(queue));

  return 0;
}


/*************************************************************************
* This function returns the vertex with the largest gain from a partition
* and removes the node from the bucket list
**************************************************************************/
int PQueueGetMax(PQueueType *queue)
{
  int vtx, i, j, gain, node;
  idxtype *iperm, *perm, *values;

  if (queue->nnodes == 0)
    return -1;

  iperm = queue->iperm;
  perm = queue->perm;
  values = queue->values;

  queue->nnodes--;
  vtx = perm[0];
  iperm[vtx] = -1;

  if ((i = queue->nnodes) > 0) {
    gain = values[i];
    node = perm[i];
    i = 0;
    while ((j=2*i+1) < queue->nnodes) {
      if (values[j] > gain) {
        if (j+1 < queue->nnodes && values[j+1] > values[j])
          j = j+1;
        values[i] = values[j];
        perm[i] = perm[j];
        iperm[perm[i]] = i;
        i = j;
      }
      else if (j+1 < queue->nnodes && values[j+1] > gain) {
        j = j+1;
        values[i] = values[j];
        perm[i] = perm[j];
        iperm[perm[i]] = i;
        i = j;
      }
      else
        break;
    }

    values[i] = gain;
    perm[i] = node;
    iperm[node] = i;
  }

  ASSERTS(PQueueCheck(queue));

  return vtx;
}
      

/*************************************************************************
* This function returns the vertex with the largest gain from a partition
**************************************************************************/
int PQueueSeeMax(PQueueType *queue)
{

  if (queue->nnodes == 0)
    return -1;
  else
    return queue->perm[0];

}
      


/*************************************************************************
* This functions checks the consistency of the heap
**************************************************************************/
int PQueueCheck(PQueueType *queue)
{
  int i, j, nnodes;

  nnodes = queue->nnodes;
  if (nnodes == 0)
    return 1;

  ASSERTSP(queue->iperm[queue->perm[0]] == 0, ("%d %d\n", queue->perm[0], queue->iperm[queue->perm[0]]));

  for (i=1; i<nnodes; i++) {
    /* printf("%d %d %d\n", queue->perm[i], queue->iperm[queue->perm[i]], i); */
    ASSERTSP(queue->iperm[queue->perm[i]] == i, ("%d %d %d %d\n", nnodes, i, queue->perm[i], queue->iperm[queue->perm[i]])); 
    ASSERTS(queue->values[i] <= queue->values[(i-1)/2]);
  }
  for (i=1; i<nnodes; i++)
    ASSERTS(queue->values[i] <= queue->values[0]);

  return 1;
}
