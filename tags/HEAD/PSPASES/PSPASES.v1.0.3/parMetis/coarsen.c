/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * coarsen.c
 *
 * This file contains code that finds a matching and performs the coarsening
 *
 * Started 2/22/96
 * George
 *
 * $Id: coarsen.c,v 1.1.1.1 2004-10-07 16:05:25 paklein Exp $
 *
 */

#include <parmetis.h>


/*************************************************************************
* This function finds a matching
**************************************************************************/
void GlobalMatch_HEM(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i, ii, j, k, l, nnbrs, nvtxs, cnvtxs, firstvtx, lastvtx, maxi, nkept, 
      otherlastvtx, nrequests, nchanged, pass, nmatched, wside;
  idxtype *xadj, *ladjncy, *adjwgt, *vwgt, *vtxdist;
  idxtype *match, *rmatch, *smatch;
  idxtype *peind, *sendptr, *recvptr;
  idxtype *perm, *iperm, *nperm, *changed;
  int *nreqs_pe;
  KeyValueType *match_requests, *match_granted, *pe_requests;

  graph->match_type = MATCH_GLOBAL;

  IFSET(ctrl->dbglvl, DBG_TIME, MPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  ladjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  vwgt = graph->vwgt;

  vtxdist = graph->vtxdist;
  firstvtx = vtxdist[ctrl->mype];
  lastvtx = vtxdist[ctrl->mype+1];

  match = graph->match = idxsmalloc(nvtxs+graph->nrecv, UNMATCHED, "HEM_Match: match");

  nnbrs = graph->nnbrs;
  peind = graph->peind;
  sendptr = graph->sendptr;
  recvptr = graph->recvptr;

  /* Use wspace->indices as the tmp space for matching info of the boundary
   * vertices that are sent and received */
  rmatch = match + nvtxs;
  smatch = wspace->indices;
  changed = smatch+graph->nsend;

  /* Use wspace->indices as the tmp space for match requests of the boundary
   * vertices that are sent and received */
  match_requests = wspace->pairs;
  match_granted = match_requests + graph->nsend;

  nreqs_pe = ismalloc(nnbrs, 0, "Match_HEM: nreqs_pe");

  nkept = graph->gnvtxs/ctrl->npes - nvtxs;

  perm = (idxtype *)wspace->degrees;
  iperm = perm + nvtxs;
  FastRandomPermute(nvtxs, perm, 1);
  for (i=0; i<nvtxs; i++)
    iperm[perm[i]] = i;

  nperm = iperm + nvtxs;
  for (i=0; i<nnbrs; i++)
    nperm[i] = i;

  /*************************************************************
   * Go now and find a matching by doing multiple iterations
   *************************************************************/
  /* First nullify the heavy vertices */
  for (nchanged=i=0; i<nvtxs; i++) {
    if (vwgt[i] > graph->maxvwgt) {
      match[i] = TOO_HEAVY;
      nchanged++;
    }
  }
  if (GlobalSESum(ctrl, nchanged) > 0) {
    IFSET(ctrl->dbglvl, DBG_PROGRESS, rprintf(ctrl, "We found %d heavy vertices!\n", GlobalSESum(ctrl, nchanged)));
    CommInterfaceData(ctrl, graph, match, smatch, rmatch);
  }


  for (nmatched=pass=0; pass<4; pass++) {
#ifdef DEBUG_MATCH
    PrintVector2(ctrl, nvtxs, firstvtx, match, "Match0"); 
#endif

    wside = (graph->level+pass)%2;
    nchanged = nrequests = 0;
    for (ii=nmatched; ii<nvtxs; ii++) {
      i = perm[ii];
      if (match[i] == UNMATCHED) {  /* Unmatched */
        maxi = -1;

        /* Find a heavy-edge matching */
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          k = ladjncy[j];
          if (match[k] == UNMATCHED && 
              (maxi == -1 || adjwgt[maxi] < adjwgt[j] || (adjwgt[maxi] == adjwgt[j] && k < nvtxs)))
            maxi = j;
        }

        if (maxi != -1) {
          k = ladjncy[maxi];
          if (k < nvtxs) { /* Take care the local vertices first */
            /* Here we give preference the local matching by granding it right away */
            if (i <= k) {
              match[i] = firstvtx+k + KEEP_BIT;
              match[k] = firstvtx+i;
            }
            else {
              match[i] = firstvtx+k;
              match[k] = firstvtx+i + KEEP_BIT;
            }
            changed[nchanged++] = i;
            changed[nchanged++] = k;
          }
          else { /* Take care any remote boundary vertices */
            match[k] = MAYBE_MATCHED;
            /* Alternate among which vertices will issue the requests */
            if ((wside ==0 && firstvtx+i < graph->imap[k]) || (wside == 1 && firstvtx+i > graph->imap[k])) { 
              match[i] = MAYBE_MATCHED;
              match_requests[nrequests].key = graph->imap[k];
              match_requests[nrequests].val = firstvtx+i;
              nrequests++;
            }
          }
        }
      }
    }


#ifdef DEBUG_MATCH
    PrintVector2(ctrl, nvtxs, firstvtx, match, "Match1");
    myprintf(ctrl, "[c: %2d] Nlocal: %d, Nrequests: %d\n", c, nlocal, nrequests);
#endif


    /***********************************************************
    * Exchange the match_requests, requests for me are stored in
    * match_granted 
    ************************************************************/
    /* Issue the receives first. Note that from each PE can receive a maximum
       of the interface node that it needs to send it in the case of a mat-vec */
    for (i=0; i<nnbrs; i++) {
      MPI_Irecv((void *)(match_granted+recvptr[i]), 2*(recvptr[i+1]-recvptr[i]), IDX_DATATYPE,
                peind[i], 1, ctrl->comm, ctrl->rreq+i);
    }

    /* Issue the sends next. This needs some work */
    ikeysort(nrequests, match_requests);
    for (j=i=0; i<nnbrs; i++) {
      otherlastvtx = vtxdist[peind[i]+1];
      for (k=j; k<nrequests && match_requests[k].key < otherlastvtx; k++);
      MPI_Isend((void *)(match_requests+j), 2*(k-j), IDX_DATATYPE, peind[i], 1, ctrl->comm, ctrl->sreq+i);
      j = k;
    }

    /* OK, now get into the loop waiting for the operations to finish */
    MPI_Waitall(nnbrs, ctrl->rreq, ctrl->statuses);
    for (i=0; i<nnbrs; i++) {
      MPI_Get_count(ctrl->statuses+i, IDX_DATATYPE, nreqs_pe+i);
      nreqs_pe[i] = nreqs_pe[i]/2;  /* Adjust for pairs of IDX_DATATYPE */
    }
    MPI_Waitall(nnbrs, ctrl->sreq, ctrl->statuses);


    /***********************************************************
    * Now, go and service the requests that you received in 
    * match_granted 
    ************************************************************/
    RandomPermute(nnbrs, nperm, 0);
    for (ii=0; ii<nnbrs; ii++) {
      i = nperm[ii];
      pe_requests = match_granted+recvptr[i];
      for (j=0; j<nreqs_pe[i]; j++) {
        k = pe_requests[j].key;
        ASSERTP(ctrl, k >= firstvtx && k < lastvtx, (ctrl, "%d %d %d %d %d\n", firstvtx, lastvtx, k, j, peind[i]));
        /* myprintf(ctrl, "Requesting a match %d %d\n", pe_requests[j].key, pe_requests[j].val); */
        if (match[k-firstvtx] == UNMATCHED) { /* Bingo, lets grand this request */
          changed[nchanged++] = k-firstvtx;
          if (nkept >= 0) { /* Flip a coin for who gets it */
            match[k-firstvtx] = pe_requests[j].val + KEEP_BIT;
            nkept--;
          }
          else {
            match[k-firstvtx] = pe_requests[j].val;
            pe_requests[j].key += KEEP_BIT;
            nkept++;
          }
          /* myprintf(ctrl, "Request from pe:%d (%d %d) granded!\n", peind[i], pe_requests[j].val, pe_requests[j].key); */ 
        }
        else { /* We are not granding the request */
          /* myprintf(ctrl, "Request from pe:%d (%d %d) not granded!\n", peind[i], pe_requests[j].val, pe_requests[j].key); */ 
          pe_requests[j].key = UNMATCHED;
        }
      }
    }


    /***********************************************************
    * Exchange the match_granded information. It is stored in
    * match_requests 
    ************************************************************/
    /* Issue the receives first. Note that from each PE can receive a maximum
       of the interface node that it needs to send during the case of a mat-vec */
    for (i=0; i<nnbrs; i++) {
      MPI_Irecv((void *)(match_requests+sendptr[i]), 2*(sendptr[i+1]-sendptr[i]), IDX_DATATYPE,
                peind[i], 1, ctrl->comm, ctrl->rreq+i);
    }

    /* Issue the sends next. */
    for (i=0; i<nnbrs; i++) {
      MPI_Isend((void *)(match_granted+recvptr[i]), 2*nreqs_pe[i], IDX_DATATYPE, 
                peind[i], 1, ctrl->comm, ctrl->sreq+i);
    }

    /* OK, now get into the loop waiting for the operations to finish */
    MPI_Waitall(nnbrs, ctrl->rreq, ctrl->statuses);
    for (i=0; i<nnbrs; i++) {
      MPI_Get_count(ctrl->statuses+i, IDX_DATATYPE, nreqs_pe+i);
      nreqs_pe[i] = nreqs_pe[i]/2;  /* Adjust for pairs of IDX_DATATYPE */
    }
    MPI_Waitall(nnbrs, ctrl->sreq, ctrl->statuses);


    /***********************************************************
    * Now, go and through the match_requests and update local
    * match information for the matchings that were granded.
    ************************************************************/
    for (i=0; i<nnbrs; i++) {
      pe_requests = match_requests+sendptr[i];
      for (j=0; j<nreqs_pe[i]; j++) {
        match[pe_requests[j].val-firstvtx] = pe_requests[j].key;
        if (pe_requests[j].key != UNMATCHED)
          changed[nchanged++] = pe_requests[j].val-firstvtx;
      }
    }

    for (i=0; i<nchanged; i++) {
      ii = iperm[changed[i]];
      perm[ii] = perm[nmatched];
      iperm[perm[nmatched]] = ii;
      nmatched++;
    }

    CommChangedInterfaceData(ctrl, graph, nchanged, changed, match, match_requests, match_granted, wspace->pv4);
  }

  /* Traverse the vertices and those that were unmatched, match them with themselves */
  cnvtxs = 0;
  for (i=0; i<nvtxs; i++) {
    if (match[i] == UNMATCHED || match[i] == TOO_HEAVY) {
      match[i] = (firstvtx+i) + KEEP_BIT;
      cnvtxs++;
    }
    else if (match[i] >= KEEP_BIT) {  /* A matched vertex which I get to keep */
      cnvtxs++;
    }
  }

  if (ctrl->dbglvl&DBG_MATCHINFO) {
    PrintVector2(ctrl, nvtxs, firstvtx, match, "Match");
    myprintf(ctrl, "Cnvtxs: %d\n", cnvtxs);
    rprintf(ctrl, "Done with matching...\n");
  }

  IFSET(ctrl->dbglvl, DBG_TIME, MPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  GKfree(&nreqs_pe, LTERM);

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ContractTmr));
  Global_CreateCoarseGraph(ctrl, graph, wspace, cnvtxs);
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ContractTmr));

}



/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void Global_CreateCoarseGraph(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace, int cnvtxs)
{
  int i, j, k, l, ii, jj, kk, ll, nnbrs, nvtxs, nedges;
  int firstvtx, lastvtx, cfirstvtx, clastvtx, otherlastvtx;
  int npes=ctrl->npes, mype=ctrl->mype;
  int cnedges, nedgesm, nsend, nrecv, nkeepsize, nrecvsize, nsendsize, v, u;
  idxtype *xadj, *vwgt, *ladjncy, *adjwgt, *vtxdist;
  idxtype *match, *cmap, *rcmap, *scmap;
  idxtype *cxadj, *cvwgt, *cadjncy, *cadjwgt, *cvtxdist;
  idxtype *rsizes, *ssizes, *rlens, *slens, *rgraph, *sgraph, *perm;
  idxtype *peind, *recvptr, *recvind;
  GraphType *cgraph;
  KeyValueType *scand, *rcand;
  int mask=(1<<13)-1, htable[8192], htableidx[8192];

  nvtxs = graph->nvtxs;

  vtxdist = graph->vtxdist;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  ladjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  match = graph->match;

  firstvtx = vtxdist[mype];
  lastvtx = vtxdist[mype+1];

  cmap = graph->cmap = idxmalloc(nvtxs+graph->nrecv, "CreateCoarseGraph: cmap");

  nnbrs = graph->nnbrs;
  peind = graph->peind;
  recvind = graph->recvind;
  recvptr = graph->recvptr;

  /* Use wspace->indices as the tmp space for map of the boundary
   * vertices that are sent and received */
  scmap = wspace->indices;
  rcmap = cmap + nvtxs;


  /* Initialize the coarser graph */
  cgraph = CreateGraph();
  cgraph->nvtxs = cnvtxs;
  cgraph->maxvwgt = graph->maxvwgt;
  cgraph->level = graph->level+1;

  cgraph->finer = graph;
  graph->coarser = cgraph;



  /*************************************************************
  * Obtain the vtxdist of the coarser graph 
  **************************************************************/
  cvtxdist = cgraph->vtxdist = idxmalloc(npes+1, "CreateCoarseGraph: cvtxdist");
  cvtxdist[npes] = cnvtxs;  /* Use last position in the cvtxdist as a temp buffer */

  MPI_Allgather((void *)(cvtxdist+npes), 1, IDX_DATATYPE, (void *)cvtxdist, 1, IDX_DATATYPE, ctrl->comm);

  MAKECSR(i, npes, cvtxdist);

  cgraph->gnvtxs = cvtxdist[npes];

#ifdef DEBUG_CONTRACT
  PrintVector(ctrl, npes+1, 0, cvtxdist, "cvtxdist");
#endif


  /*************************************************************
  * Construct the cmap vector 
  **************************************************************/
  cfirstvtx = cvtxdist[mype];
  clastvtx = cvtxdist[mype+1];

  /* Create the cmap of what you know so far locally */
  cnvtxs = 0;
  for (i=0; i<nvtxs; i++) {
    if (match[i] >= KEEP_BIT) {
      k = match[i] - KEEP_BIT;
      if (k>=firstvtx && k<firstvtx+i)
        continue;  /* Both (i,k) are local and i has been matched via the (k,i) side */

      cmap[i] = cfirstvtx + cnvtxs++;
      if (k != firstvtx+i && (k>=firstvtx && k<lastvtx)) { /* I'm matched locally */
        cmap[k-firstvtx] = cmap[i];
        match[k-firstvtx] += KEEP_BIT;  /* Add the KEEP_BIT to simplify coding */
      }
    }
  }
  ASSERT(ctrl, cnvtxs == clastvtx-cfirstvtx);

  CommInterfaceData(ctrl, graph, cmap, scmap, rcmap);

  /* Update the cmap of the locally stored vertices that will go away. 
   * The remote processor assigned cmap for them */
  for (i=0; i<nvtxs; i++) {
    if (match[i] < KEEP_BIT) { /* Only vertices that go away satisfy this*/
      cmap[i] = rcmap[BSearch(graph->nrecv, recvind, match[i])];
    }
  }

  CommInterfaceData(ctrl, graph, cmap, scmap, rcmap);


#ifdef DEBUG_CONTRACT
  PrintVector(ctrl, nvtxs, firstvtx, cmap, "Cmap");
#endif


  /*************************************************************
  * Determine how many adjcency lists you need to send/receive.
  **************************************************************/
  /* Use wspace->pairs as the tmp space for the boundary vertices that are sent and received */
  scand = wspace->pairs;
  rcand = graph->rcand = (KeyValueType *)GKmalloc(recvptr[nnbrs]*sizeof(KeyValueType), "CreateCoarseGraph: rcand");

  nkeepsize = nsend = nrecv = 0;
  for (i=0; i<nvtxs; i++) {
    if (match[i] < KEEP_BIT) { /* This is going away */
      scand[nsend].key = match[i];
      scand[nsend].val = i;
      nsend++;
    }
    else {
      nkeepsize += (xadj[i+1]-xadj[i]);

      k = match[i]-KEEP_BIT;
      if (k<firstvtx || k>=lastvtx) { /* This is comming from afar */
        rcand[nrecv].key = k;
        rcand[nrecv].val = cmap[i] - cfirstvtx;  /* Set it for use during the partition projection */
        ASSERT(ctrl, rcand[nrecv].val>=0 && rcand[nrecv].val<cnvtxs);
        nrecv++;
      }
    }
  }


#ifdef DEBUG_CONTRACT
  PrintPairs(ctrl, nsend, scand, "scand");
  PrintPairs(ctrl, nrecv, rcand, "rcand");
#endif

  /***************************************************************
  * Determine how many lists and their sizes  you will send and 
  * received for each of the neighboring PEs
  ****************************************************************/
  rsizes = wspace->pv1;
  ssizes = wspace->pv2;
  idxset(nnbrs, 0, ssizes);
  idxset(nnbrs, 0, rsizes);
  rlens = graph->rlens = idxmalloc(nnbrs+1, "CreateCoarseGraph: graph->rlens");
  slens = graph->slens = idxmalloc(nnbrs+1, "CreateCoarseGraph: graph->slens");

  /* Take care the sending data first */
  ikeyvalsort(nsend, scand);
  slens[0] = 0;
  for (k=i=0; i<nnbrs; i++) {
    otherlastvtx = vtxdist[peind[i]+1];
    for (; k<nsend && scand[k].key < otherlastvtx; k++)
      ssizes[i] += (xadj[scand[k].val+1]-xadj[scand[k].val]);
    slens[i+1] = k;
  }

  /* Take care the receiving data next. You cannot yet determine the rsizes[] */
  ikeyvalsort(nrecv, rcand);
  rlens[0] = 0;
  for (k=i=0; i<nnbrs; i++) {
    otherlastvtx = vtxdist[peind[i]+1];
    for (; k<nrecv && rcand[k].key < otherlastvtx; k++);
    rlens[i+1] = k;
  }

#ifdef DEBUG_CONTRACT
  PrintVector(ctrl, nnbrs+1, 0, slens, "slens");
  PrintVector(ctrl, nnbrs+1, 0, rlens, "rlens");
#endif

  /***************************************************************
  * Exchange size information
  ****************************************************************/
  /* Issue the receives first. */
  for (i=0; i<nnbrs; i++) {
    if (rlens[i+1]-rlens[i] > 0)  /* Issue a receive only if you are getting something */
      MPI_Irecv((void *)(rsizes+i), 1, IDX_DATATYPE, peind[i], 1, ctrl->comm, ctrl->rreq+i);
  }

  /* Take care the sending data next */
  for (i=0; i<nnbrs; i++) {
    if (slens[i+1]-slens[i] > 0)  /* Issue a send only if you are sending something */
      MPI_Isend((void *)(ssizes+i), 1, IDX_DATATYPE, peind[i], 1, ctrl->comm, ctrl->sreq+i);
  }

  /* OK, now get into the loop waiting for the operations to finish */
  for (i=0; i<nnbrs; i++) {
    if (rlens[i+1]-rlens[i] > 0)  
      MPI_Wait(ctrl->rreq+i, &ctrl->status);
  }
  for (i=0; i<nnbrs; i++) {
    if (slens[i+1]-slens[i] > 0)  
      MPI_Wait(ctrl->sreq+i, &ctrl->status);
  }


#ifdef DEBUG_CONTRACT
  PrintVector(ctrl, nnbrs, 0, rsizes, "rsizes");
  PrintVector(ctrl, nnbrs, 0, ssizes, "ssizes");
#endif

  /*************************************************************
  * Allocate memory for received/sent graphs and start sending 
  * and receiving data.
  * rgraph and sgraph is a different data structure than CSR
  * to fasilitate single message exchange.
  **************************************************************/
  nrecvsize = idxsum(nnbrs, rsizes);
  nsendsize = idxsum(nnbrs, ssizes);
  if (3*(nrecv+nsend) + 2*(nrecvsize+nsendsize) <= wspace->nlarge) {  
    rgraph = (idxtype *)wspace->degrees;
    sgraph = rgraph + 3*nrecv+2*nrecvsize;
  }
  else {
    rgraph = idxmalloc(3*nrecv+2*nrecvsize, "CreateCoarseGraph: rgraph");
    sgraph = idxmalloc(3*nsend+2*nsendsize, "CreateCoarseGraph: sgraph");
  }


  /* Deal with the received portion first */
  for (l=i=0; i<nnbrs; i++) {
    if (rlens[i+1]-rlens[i] > 0) {  /* Issue a receive only if you are getting something */
      /* myprintf(ctrl, "Receiveing from pe:%d, %d lists of size %d\n", peind[i], rlens[i+1]-rlens[i], rsizes[i]); */
      MPI_Irecv((void *)(rgraph+l), 3*(rlens[i+1]-rlens[i])+2*rsizes[i], IDX_DATATYPE, peind[i], 1, ctrl->comm, ctrl->rreq+i);
      l += 3*(rlens[i+1]-rlens[i])+2*rsizes[i];
    }
  }


  /* Deal with the sent portion now */
  for (ll=l=i=0; i<nnbrs; i++) {
    if (slens[i+1]-slens[i] > 0) {  /* Issue a send only if you are sending something */
      for (k=slens[i]; k<slens[i+1]; k++) {
        ii = scand[k].val;
        sgraph[ll++] = firstvtx+ii;
        sgraph[ll++] = xadj[ii+1]-xadj[ii];
        sgraph[ll++] = vwgt[ii];
        for (jj=xadj[ii]; jj<xadj[ii+1]; jj++) {
          sgraph[ll++] = cmap[ladjncy[jj]];
          sgraph[ll++] = adjwgt[jj];
        }
      }

      ASSERT(ctrl, ll-l == 3*(slens[i+1]-slens[i])+2*ssizes[i]);

      /* myprintf(ctrl, "Sending to pe:%d, %d lists of size %d\n", peind[i], slens[i+1]-slens[i], ssizes[i]); */
      MPI_Isend((void *)(sgraph+l), ll-l, IDX_DATATYPE, peind[i], 1, ctrl->comm, ctrl->sreq+i);
      l = ll;
    }
  }

  /* OK, now get into the loop waiting for the operations to finish */
  for (i=0; i<nnbrs; i++) {
    if (rlens[i+1]-rlens[i] > 0)  
      MPI_Wait(ctrl->rreq+i, &ctrl->status);
  }
  for (i=0; i<nnbrs; i++) {
    if (slens[i+1]-slens[i] > 0)  
      MPI_Wait(ctrl->sreq+i, &ctrl->status);
  }


#ifdef DEBUG_CONTRACT
  rprintf(ctrl, "Graphs were sent!\n");
  PrintTransferedGraphs(ctrl, nnbrs, peind, slens, rlens, sgraph, rgraph);
#endif

  /*************************************************************
  * Setup the mapping from indices returned by BSearch to 
  * those that are actually stored 
  **************************************************************/
  perm = idxsmalloc(recvptr[nnbrs], -1, "CreateCoarseGraph: perm");
  for (j=i=0; i<nrecv; i++) {
    /* myprintf(ctrl, "For received vertex %d, set perm[%d]=%d\n", rgraph[j], BSearch(graph->nrecv, recvind, rgraph[j]), j+1); */
    perm[BSearch(graph->nrecv, recvind, rgraph[j])] = j+1;
    j += 3+2*rgraph[j+1];
  }


  /*************************************************************
  * Finally, create the coarser graph
  **************************************************************/
  /* Allocate memory for the coarser graph, and fire up coarsening */
  cxadj = cgraph->xadj = idxmalloc(cnvtxs+1, "CreateCoarserGraph: cxadj");
  cvwgt = cgraph->vwgt = idxmalloc(cnvtxs, "CreateCoarserGraph: cvwgt");
  cadjncy = idxmalloc(2*(nkeepsize+nrecvsize), "CreateCoarserGraph: cadjncy");
  cadjwgt = cadjncy + nkeepsize+nrecvsize;

  iset(8192, -1, htable);

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    if (match[i] >= KEEP_BIT) {
      v = firstvtx+i; 
      u = match[i]-KEEP_BIT;

      if (u>=firstvtx && u<lastvtx && v > u) 
        continue;  /* I have already collapsed it as (u,v) */

      /* Collapse the v vertex first, which you know that is local */
      cvwgt[cnvtxs] = vwgt[i];
      nedges = 0;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        k = cmap[ladjncy[j]];
        if (k != cfirstvtx+cnvtxs) {  /* If this is not an internal edge */
          l = k&mask;
          if (htable[l] == -1) { /* Seeing this for first time */
            htable[l] = k;
            htableidx[l] = cnedges+nedges;
            cadjncy[cnedges+nedges] = k; 
            cadjwgt[cnedges+nedges++] = adjwgt[j];
          }
          else if (htable[l] == k) {
            cadjwgt[htableidx[l]] += adjwgt[j];
          }
          else { /* Now you have to go and do a search. Expensive case */
            for (l=0; l<nedges; l++) {
              if (cadjncy[cnedges+l] == k)
                break;
            }
            if (l < nedges) {
              cadjwgt[cnedges+l] += adjwgt[j];
            }
            else {
              cadjncy[cnedges+nedges] = k; 
              cadjwgt[cnedges+nedges++] = adjwgt[j];
            }
          }
        }
      }

      /* Collapse the u vertex next */
      if (v != u) { 
        if (u>=firstvtx && u<lastvtx) { /* Local vertex */
          u -= firstvtx;
          cvwgt[cnvtxs] += vwgt[u];

          for (j=xadj[u]; j<xadj[u+1]; j++) {
            k = cmap[ladjncy[j]];
            if (k != cfirstvtx+cnvtxs) {  /* If this is not an internal edge */
              l = k&mask;
              if (htable[l] == -1) { /* Seeing this for first time */
                htable[l] = k;
                htableidx[l] = cnedges+nedges;
                cadjncy[cnedges+nedges] = k; 
                cadjwgt[cnedges+nedges++] = adjwgt[j];
              }
              else if (htable[l] == k) {
                cadjwgt[htableidx[l]] += adjwgt[j];
              }
              else { /* Now you have to go and do a search. Expensive case */
                for (l=0; l<nedges; l++) {
                  if (cadjncy[cnedges+l] == k)
                    break;
                }
                if (l < nedges) {
                  cadjwgt[cnedges+l] += adjwgt[j];
                }
                else {
                  cadjncy[cnedges+nedges] = k; 
                  cadjwgt[cnedges+nedges++] = adjwgt[j];
                }
              }
            }
          }
        }
        else { /* Remote vertex */
          u = perm[BSearch(graph->nrecv, recvind, u)];
          cvwgt[cnvtxs] += rgraph[u+1];  /* Remember that the +1 stores the vertex weight */

          for (j=0; j<rgraph[u]; j++) {
            k = rgraph[u+2+2*j];
            if (k != cfirstvtx+cnvtxs) {  /* If this is not an internal edge */
              l = k&mask;
              if (htable[l] == -1) { /* Seeing this for first time */
                htable[l] = k;
                htableidx[l] = cnedges+nedges;
                cadjncy[cnedges+nedges] = k; 
                cadjwgt[cnedges+nedges++] = rgraph[u+2+2*j+1];
              }
              else if (htable[l] == k) {
                cadjwgt[htableidx[l]] += rgraph[u+2+2*j+1];
              }
              else { /* Now you have to go and do a search. Expensive case */
                for (l=0; l<nedges; l++) {
                  if (cadjncy[cnedges+l] == k)
                    break;
                }
                if (l < nedges) {
                  cadjwgt[cnedges+l] += rgraph[u+2+2*j+1];
                }
                else {
                  cadjncy[cnedges+nedges] = k; 
                  cadjwgt[cnedges+nedges++] = rgraph[u+2+2*j+1];
                }
              }
            }
          }
        }
      }

      cnedges += nedges;
      for (j=cxadj[cnvtxs]; j<cnedges; j++)
        htable[cadjncy[j]&mask] = -1;  /* reset the htable */
      cxadj[++cnvtxs] = cnedges;
    }
  }

  cgraph->nedges = cnedges;

  cgraph->adjncy = idxmalloc(cnedges, "CreateCoarserGraph: cadjncy");
  cgraph->adjwgt = idxmalloc(cnedges, "CreateCoarserGraph: cadjwgt");
  idxcopy(cnedges, cadjncy, cgraph->adjncy);
  idxcopy(cnedges, cadjwgt, cgraph->adjwgt);
  free(cadjncy);

  free(perm);

  if (rgraph != (idxtype *)wspace->degrees) 
    GKfree(&rgraph, &sgraph, LTERM);

}



/*************************************************************************
* This function finds a HEM matching between local vertices only
**************************************************************************/
void LocalMatch_HEM(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i, ii, j, k, l, nvtxs, cnvtxs, firstvtx, maxi; 
  idxtype *xadj, *ladjncy, *adjwgt, *vtxdist;
  idxtype *perm, *match;

  graph->match_type = MATCH_LOCAL;

  IFSET(ctrl->dbglvl, DBG_TIME, MPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  ladjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  vtxdist = graph->vtxdist;
  firstvtx = vtxdist[ctrl->mype];

  match = graph->match = idxmalloc(nvtxs+graph->nrecv, "HEM_Match: match");
  idxset(nvtxs, UNMATCHED, match);
  idxset(graph->nrecv, 0, match+nvtxs);  /* Easy way to handle remote vertices */


  /*************************************************************
   * Go now and find a local matching 
   *************************************************************/
  perm = wspace->indices;
  FastRandomPermute(nvtxs, perm, 1);
  cnvtxs = 0;
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];
    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxi = -1;

      /* Find a heavy-edge matching, if the weight of the vertex is OK */
      if (graph->vwgt[i] < graph->maxvwgt) {
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          if (match[ladjncy[j]] == UNMATCHED && graph->vwgt[ladjncy[j]] <= graph->maxvwgt &&
              (maxi == -1 || adjwgt[maxi] < adjwgt[j]))
            maxi = j;
        }
      }

      if (maxi != -1) {
        k = ladjncy[maxi];
        if (i <= k) {
          match[i] = firstvtx+k + KEEP_BIT;
          match[k] = firstvtx+i;
        }
        else {
          match[i] = firstvtx+k;
          match[k] = firstvtx+i + KEEP_BIT;
        }
      }
      else {
        match[i] = (firstvtx+i) + KEEP_BIT;
      }
      cnvtxs++;
    }
  }

  CommInterfaceData(ctrl, graph, match, wspace->indices, match+nvtxs);


#ifdef DEBUG_MATCH
  PrintVector2(ctrl, nvtxs, firstvtx, match, "Match1");
#endif


  if (ctrl->dbglvl&DBG_MATCHINFO) {
    PrintVector2(ctrl, nvtxs, firstvtx, match, "Match");
    myprintf(ctrl, "Cnvtxs: %d\n", cnvtxs);
    rprintf(ctrl, "Done with matching...\n");
  }

  IFSET(ctrl->dbglvl, DBG_TIME, MPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));


  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ContractTmr));
  Local_CreateCoarseGraph(ctrl, graph, wspace, cnvtxs);
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ContractTmr));

}





/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void Local_CreateCoarseGraph(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace, int cnvtxs)
{
  int i, j, k, l, ii, jj, kk, ll, nnbrs, nvtxs, nedges, firstvtx, cfirstvtx;
  int npes=ctrl->npes, mype=ctrl->mype;
  int cnedges, v, u;
  idxtype *xadj, *vwgt, *vsize, *ladjncy, *adjwgt, *vtxdist, *where;
  idxtype *match, *cmap, *rcmap, *scmap;
  idxtype *cxadj, *cvwgt, *cvsize, *cadjncy, *cadjwgt, *cvtxdist, *cwhere;
  GraphType *cgraph;
  int mask=(1<<13)-1, htable[8192], htableidx[8192];

  nvtxs = graph->nvtxs;

  vtxdist = graph->vtxdist;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vsize = graph->vsize;
  ladjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;

  match = graph->match;

  firstvtx = vtxdist[mype];

  cmap = graph->cmap = idxmalloc(nvtxs+graph->nrecv, "CreateCoarseGraph: cmap");

  /* Initialize the coarser graph */
  cgraph = CreateGraph();
  cgraph->nvtxs = cnvtxs;
  cgraph->maxvwgt = graph->maxvwgt;
  cgraph->level = graph->level+1;

  cgraph->finer = graph;
  graph->coarser = cgraph;


  /*************************************************************
  * Obtain the vtxdist of the coarser graph 
  **************************************************************/
  cvtxdist = cgraph->vtxdist = idxmalloc(npes+1, "CreateCoarseGraph: cvtxdist");
  cvtxdist[npes] = cnvtxs;  /* Use last position in the cvtxdist as a temp buffer */

  MPI_Allgather((void *)(cvtxdist+npes), 1, IDX_DATATYPE, (void *)cvtxdist, 1, IDX_DATATYPE, ctrl->comm);

  MAKECSR(i, npes, cvtxdist);

  cgraph->gnvtxs = cvtxdist[npes];

#ifdef DEBUG_CONTRACT
  PrintVector(ctrl, npes+1, 0, cvtxdist, "cvtxdist");
#endif


  /*************************************************************
  * Construct the cmap vector 
  **************************************************************/
  cfirstvtx = cvtxdist[mype];

  /* Create the cmap of what you know so far locally */
  cnvtxs = 0;
  for (i=0; i<nvtxs; i++) {
    if (match[i] >= KEEP_BIT) {
      k = match[i] - KEEP_BIT;
      if (k<firstvtx+i)
        continue;  /* i has been matched via the (k,i) side */

      cmap[i] = cfirstvtx + cnvtxs++;
      if (k != firstvtx+i) {
        cmap[k-firstvtx] = cmap[i];
        match[k-firstvtx] += KEEP_BIT;  /* Add the KEEP_BIT to simplify coding */
      }
    }
  }

  CommInterfaceData(ctrl, graph, cmap, wspace->indices, cmap+nvtxs);


#ifdef DEBUG_CONTRACT
  PrintVector(ctrl, nvtxs, firstvtx, cmap, "Cmap");
#endif



  /*************************************************************
  * Finally, create the coarser graph
  **************************************************************/
  /* Allocate memory for the coarser graph, and fire up coarsening */
  cxadj = cgraph->xadj = idxmalloc(cnvtxs+1, "CreateCoarserGraph: cxadj");
  cvwgt = cgraph->vwgt = idxmalloc(cnvtxs, "CreateCoarserGraph: cvwgt");
  if (vsize != NULL)
    cvsize = cgraph->vsize = idxmalloc(cnvtxs, "CreateCoarserGraph: cvsize");
  if (where != NULL)
    cwhere = cgraph->where = idxmalloc(cnvtxs, "CreateCoarserGraph: cwhere");
  cadjncy = idxmalloc(2*graph->nedges, "CreateCoarserGraph: cadjncy");
  cadjwgt = cadjncy+graph->nedges;

  iset(8192, -1, htable);

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = firstvtx+i; 
    u = match[i]-KEEP_BIT;

    if (v > u) 
      continue;  /* I have already collapsed it as (u,v) */

    /* Collapse the v vertex first, which you know that is local */
    cvwgt[cnvtxs] = vwgt[i];
    if (vsize != NULL)
      cvsize[cnvtxs] = vsize[i];
    if (where != NULL)
      cwhere[cnvtxs] = where[i];
    nedges = 0;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = cmap[ladjncy[j]];
      if (k != cfirstvtx+cnvtxs) {  /* If this is not an internal edge */
        l = k&mask;
        if (htable[l] == -1) { /* Seeing this for first time */
          htable[l] = k;
          htableidx[l] = cnedges+nedges;
          cadjncy[cnedges+nedges] = k; 
          cadjwgt[cnedges+nedges++] = adjwgt[j];
        }
        else if (htable[l] == k) {
          cadjwgt[htableidx[l]] += adjwgt[j];
        }
        else { /* Now you have to go and do a search. Expensive case */
          for (l=0; l<nedges; l++) {
            if (cadjncy[cnedges+l] == k)
              break;
          }
          if (l < nedges) {
            cadjwgt[cnedges+l] += adjwgt[j];
          }
          else {
            cadjncy[cnedges+nedges] = k; 
            cadjwgt[cnedges+nedges++] = adjwgt[j];
          }
        }
      }
    }

    /* Collapse the u vertex next */
    if (v != u) { 
      u -= firstvtx;
      cvwgt[cnvtxs] += vwgt[u];
      if (vsize != NULL)
        cvsize[cnvtxs] += vsize[u];
      if (where != NULL && cwhere[cnvtxs] != where[u])
        myprintf(ctrl, "Something went wrong with the where local matching! %d %d\n", cwhere[cnvtxs], where[u]);

      for (j=xadj[u]; j<xadj[u+1]; j++) {
        k = cmap[ladjncy[j]];
        if (k != cfirstvtx+cnvtxs) {  /* If this is not an internal edge */
          l = k&mask;
          if (htable[l] == -1) { /* Seeing this for first time */
            htable[l] = k;
            htableidx[l] = cnedges+nedges;
            cadjncy[cnedges+nedges] = k; 
            cadjwgt[cnedges+nedges++] = adjwgt[j];
          }
          else if (htable[l] == k) {
            cadjwgt[htableidx[l]] += adjwgt[j];
          }
          else { /* Now you have to go and do a search. Expensive case */
            for (l=0; l<nedges; l++) {
              if (cadjncy[cnedges+l] == k)
                break;
            }
            if (l < nedges) {
              cadjwgt[cnedges+l] += adjwgt[j];
            }
            else {
              cadjncy[cnedges+nedges] = k; 
              cadjwgt[cnedges+nedges++] = adjwgt[j];
            }
          }
        }
      }
    }

    cnedges += nedges;
    for (j=cxadj[cnvtxs]; j<cnedges; j++)
      htable[cadjncy[j]&mask] = -1;  /* reset the htable */
    cxadj[++cnvtxs] = cnedges;
  }

  cgraph->nedges = cnedges;

  cgraph->adjncy = idxmalloc(cnedges, "CreateCoarserGraph: cadjncy");
  cgraph->adjwgt = idxmalloc(cnedges, "CreateCoarserGraph: cadjwgt");
  idxcopy(cnedges, cadjncy, cgraph->adjncy);
  idxcopy(cnedges, cadjwgt, cgraph->adjwgt);
  GKfree(&cadjncy, &graph->where, LTERM); /* Note that graph->where works fine even if it is NULL */

}


