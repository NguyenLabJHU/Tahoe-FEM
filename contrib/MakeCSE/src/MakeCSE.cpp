/* $Id: MakeCSE.cpp,v 1.3 2002-10-08 20:51:50 paklein Exp $ */
// for debugging, this prints extra information
//#define _PRINT_DEBUG_
#include "MakeCSE.h"
#include "MakeCSE_IOManager.h"
#include "MakeCSE_FEManager.h"

const int kPrint = 40000;

using namespace Tahoe;

/* constructor */
MakeCSE::MakeCSE (ostream& log, GlobalEdgeFinderT& Edger) :
  out (log),
  fPrintUpdate (false),
  theEdger (&Edger)
{
  ArrayT<StringT> cnames (4);
  cnames[0] = "initial setup";
  cnames[1] = "collect facets";
  cnames[2] = "connectivity renumbering";
  cnames[3] = "map node sets";
  cClock.Set (cnames);
}

MakeCSE::~MakeCSE (void)
{
}

void MakeCSE::Initialize (MakeCSE_IOManager& theInput, MakeCSE_FEManager& FEM, int comments)
{
  if (comments == 1) fPrintUpdate = true;
  fNumStartElements = theEdger->TotalElements();

  cClock.Initial ();
  SetFE (FEM);
  SetInput (theInput);
  cClock.Sum (0);
}

void MakeCSE::Create (void)
{
  // Reconstruct Inverse connects using CSE's
  cClock.Initial ();
  theEdger->ResetInvConnects(theNodes->NumNodes());

  // index over the number of side sets used to insert at
  cClock.Sum(0);
  RenumberFaceNodes ();
  cClock.Sum (2);

  // map node sets
  theEdger->ResetInvConnects(theNodes->NumNodes());
  theNodes->MapNodeSets (fSurface1Facets, *theEdger);
  cClock.Sum (3);

  // collect data for use with contact algorithms
  if (fContact.Length() > 0) CollectSurfaceData ();

  // report mass less nodes
  CollectMassLessNodes ();

  cClock.Print (cout);
  cClock.Print (out);
  PrintControlEnd (out);
  PrintControlEnd (cout);

  if (fNumStartElements > 50000)
    cout << "\n Writing Output File." << endl << endl;
}

/************ Private ***************************/
/************ Private ***************************/
/************ Private ***************************/

void MakeCSE::SetFE (MakeCSE_FEManager& FEM)
{
  // set up link to regular elements
  int num = FEM.NumRegularGroups ();
  theElements.Allocate (num);
  for (int i=0; i < num; i++)
    theElements[i] = FEM.ElementGroup (i);

  // set link to CS elements
  int numcse = FEM.NumCSEGroups ();
  theCSEs.Allocate (numcse);
  for (int ic=0; ic < numcse; ic++)
    theCSEs[ic] = FEM.ElementGroup (num + ic);
    
  // set up link to nodes
  theNodes = FEM.NodeManager();
  fNumStartNodes = theNodes->NumNodes();
}

void MakeCSE::SetInput (MakeCSE_IOManager& theInput)
{
  // grab input data, collect individual facets later
  out << "\n M a k e   C S E   D a t a :\n" << '\n';
  iArrayT facetdata, zonedata, boundarydata;
  theInput.InputData (facetdata, MakeCSE_IOManager::kFacet);
  theInput.InputData (zonedata, MakeCSE_IOManager::kZone);
  theInput.InputData (boundarydata, MakeCSE_IOManager::kBoundary);

  int nummeth = 0;
  if (boundarydata.Length() > 0) nummeth++;
  if (zonedata.Length() > 0) nummeth++;
  if (facetdata.Length() > 0) nummeth++;
  if (nummeth > 1)
    cout << "\n\n *** Warning, you are using more than one method. ***\n\n";

  // collect data
  cClock.Initial();
  if (facetdata.Length() > 0) CollectFacets (theInput, facetdata);
  if (zonedata.Length() > 0) CollectZones (theInput, zonedata);
  if (boundarydata.Length() > 0) CollectBoundaries (boundarydata);
  cClock.Sum (1);

  // remove single nodes from potential split node list
  CollectSingleNodes (theInput);

  // set up contact data ID values
  InitializeContact (theInput);

  out.flush ();
}

void MakeCSE::InitializeContact (MakeCSE_IOManager& theInput)
{
  theInput.InputData (fContact, MakeCSE_IOManager::kContactData);

  // collect existing side set ID's
  fSSetID = 1;
  for (int e=0; e < theElements.Length(); e++)
    for (int s=0; s < theElements[e]->NumSideSets(); s++)
      {
	int id = theElements[e]->SideSetID(s);
	while (fSSetID <= id) fSSetID++;
      }

  // collect existing node set ID's
  fNSetID = 1;
  for (int n=0; n < theNodes->NumNodeSets(); n++)
    {
      int id = theNodes->NodeSetID(n);
      while (fNSetID <= id) fNSetID++;
    }

  out << " CSE Block ID's to prep for contact surfaces . . = "
      << fContact.Length() << '\n';
  fContact.WriteWrapped (out, 5);
}

void MakeCSE::CollectFacets (MakeCSE_IOManager& theInput, const iArrayT& facetdata)
{
  cout << "\n Collecting Facet Data . . ." << endl;
  out << " Side Sets Specified by user. .  . . . . . . . . = "
      << facetdata.Length()/3 << '\n';
  if (fPrintUpdate)
    {
      out << " SideSet   Facet   CSEID\n";
      facetdata.WriteWrapped (out, 3);
    }

  iArray2DT set;
  int cs = theEdger->TotalElements();
  for (int i=0; i < facetdata.Length(); i += 3)
    {
      int group = theEdger->ElementGroup (facetdata [i + 1]);
      int csegroup = theEdger->ElementGroup (facetdata [i + 2]);
      int cselemgroup = csegroup - theElements.Length();
      InitializeSet (group, cselemgroup);

      // read data, read as local, change to my global, 
      // which might differ depending on the order of element block storage
      StringT sid;
      sid.Append (facetdata[i]);
      const iArray2DT temp = theInput.SideSet (sid);
      bool local = theInput.IsSideSetLocal (sid);
      if (local)
	theInput.SideSetLocalToGlobal (sid, temp, set);
      else
	set = temp;

      // check set validity
      if (!theElements[group]->CheckSideSet (set))
	{
	  cout << "MakeCSE::CollectFacets, side set " << facetdata[i] 
	       << " fails CheckSideSet for element group id: " 
	       << facetdata [i + 1] << endl;
	  throw eBadInputValue;
	}

      // convert to global numbering 
      int *elem = set.Pointer();
      for (int j=0; j < set.MajorDim(); j++, elem += 2)
	*elem = theEdger->GlobalElement (*elem, group);

      // initialize set data
      AddElements (set.MajorDim(), cselemgroup);
      
      // insert cse
      int *pelem = set.Pointer();
      int *pface = set.Pointer(1);
      for (int e=0; e < set.MajorDim(); e++, pelem += 2, pface += 2, cs++)
	{
	  InitializeFacet (*pelem, *pface, group, cs, cselemgroup);
	  if ((set.MajorDim() > kPrint && (e+1)%kPrint == 0) ||
	      e+1 == set.MajorDim())
	    cout << "   " << e+1 << " done initializing of " << set.MajorDim() << " facets " << endl;
	}
    }
}

void MakeCSE::CollectSingleNodes (MakeCSE_IOManager& theInput)
{
  iArrayT setids;
  theInput.InputData (setids, MakeCSE_IOManager::kSingleNodes);
  iAutoArrayT nodes;
  iArrayT nodeset;
  for (int i=0; i < setids.Length(); i++)
    {
      StringT sids;
      sids.Append (setids[i]);
      nodeset = theInput.NodeSet (sids);
      nodes.AppendUnique (nodeset);
    }
  out << " Nodes specified by user not to be split . . . . = "
      << nodes.Length() << '\n';
  if (fPrintUpdate)
    {
      for (int p=0; p < nodes.Length(); p++)
	{
	  out << setw (8) << nodes[p] + 1;
	  if ((p+1)%6 == 0) out << '\n';
	}
      out << '\n';
    }
  
  RemoveSingleNodes (nodes);
}

void MakeCSE::CollectZones (MakeCSE_IOManager& theInput, const iArrayT& zonedata)
{
  cout << "\n Collecting Zone Data . . ." << endl;
  out << " Element Group Zones . . . . . . . . . . . . . . = "
      << zonedata.Length()/2 << '\n';
  if (fPrintUpdate)
    {
      out << "  BlockID OutputID\n";
      zonedata.WriteWrapped (out, 2);
    }

  iArrayT elemids (zonedata.Length()/2), cseids (zonedata.Length()/2);
  for (int y=0, z=0; y < elemids.Length(); y++, z+= 2)
    {
      elemids[y] = zonedata [z];
      cseids[y] = zonedata [z+1];
    }

  // collect zone facets plus those facets shared boundary facets
  // but only shared boundary facets from blocks with higher group numbers
  iAutoArrayT boundarynodes;
  iArray2DT set;
  int cs = theEdger->TotalElements();
  for (int i=0; i < elemids.Length(); i++)
    {
      int group = theEdger->ElementGroup (elemids[i]);
      int cselemgroup = theEdger->ElementGroup (cseids[i]) - theElements.Length();
      InitializeSet (group, cselemgroup);
      theEdger->ZoneFacets (elemids[i], elemids, set, boundarynodes);

      // insert cse
      int *pelem = set.Pointer();
      int *pface = set.Pointer(1);
      AddElements (set.MajorDim(), cselemgroup);
      for (int e=0; e < set.MajorDim(); e++, pelem += 2, pface += 2, cs++)
	{
	  InitializeFacet (*pelem, *pface, group, cs, cselemgroup);
	  if ((set.MajorDim() > kPrint && (e+1)%kPrint == 0) ||
	      e+1 == set.MajorDim())
	    cout << "   " << e+1 << " done initializing of " << set.MajorDim() << " facets " << endl;
	}

      if (elemids.Length() < 10) 
	cout << "  Done with Zone " << elemids[i] << endl;
    }

  // based on user input, allow some boundary nodes to be split
  int zoneedgetype = kSingleZE;
  theInput.InputData (zoneedgetype, MakeCSE_IOManager::kZoneEdge);
  out << " Zone Edge Node Splitting Method . . . . . . . . = "
      << zoneedgetype << "\n";

  switch (zoneedgetype)
    {
    case kSingleZE: break; // don't modify boundary node list
    case kDoubleZE: boundarynodes.Free(); break; // remove all nodes from list
    case kMixSingZE:
      {
	boundarynodes.Free();
	iArrayT nodesets;
	theInput.InputData (nodesets, MakeCSE_IOManager::kZoneEdgeNSets);
	for (int i=0; i < nodesets.Length(); i++)
	  {
	    iArrayT temp;
	    StringT sids;
	    sids.Append (nodesets[i]);
	    temp = theInput.NodeSet (sids);
	    boundarynodes.AppendUnique (temp);
	  }
	out << " Zone Edge Node Sets not to be Split . . . . . . = "
	    << nodesets << "\n";
	break;
      }
    case kMixDoubZE:
      {
	iArrayT nodesets;
	theInput.InputData (nodesets, MakeCSE_IOManager::kZoneEdgeNSets);
	for (int i=0; i < nodesets.Length(); i++)
	  {
	    iArrayT temp;
	    StringT sids;
	    sids.Append (nodesets[i]);
	    temp = theInput.NodeSet (sids);
	    out << " Zone Edge Node Sets to be Split . . . . . . . . . = "
		<< nodesets << "\n";

	    int *pt = temp.Pointer();
	    for (int j=0; j < temp.Length(); j++, pt++)
	      {
		int dex = boundarynodes.PositionOf (*pt);
		if (dex > -1) boundarynodes.DeleteAt (dex);
	      }
	  }
	break;
      }
    }

  // remove boundary nodes from split list
  out << " Nodes on edge of zone, not to be split. . . . . = "
      << boundarynodes.Length() << "\n";
  if (fPrintUpdate)
    {
      for (int p=0; p < boundarynodes.Length(); p++)
	{
	  out << setw (8) << boundarynodes[p] + 1;
	  if ((p+1)%6 == 0) out << '\n';
	}
      out << '\n';
    }

  RemoveSingleNodes (boundarynodes);  

  cout << "  Done with Zone Data" << endl;
}

void MakeCSE::CollectBoundaries (const iArrayT& boundarydata)
{
  cout << "\n Collecting Boundary Data . . ." << endl;
  out << " Element Group Boundaries. . . . . . . . . . . . = "
      << boundarydata.Length()/4 << '\n';
  if (fPrintUpdate)
    {
      out << "  BlockID1 BlkID2Min BlkID2Max  OutputID\n";
      boundarydata.WriteWrapped (out, 4);
    }

  iArray2DT set;
  int cs = theEdger->TotalElements();
  for (int i=0; i < boundarydata.Length(); i += 4)
    {
      int groupid = boundarydata [i];
      int start = boundarydata [i + 1];
      int stop  = boundarydata [i + 2];
      // only compare to higher groups
      if (start < groupid) start = groupid + 1;
      if (stop >= start)
	{
	  int group = theEdger->ElementGroup (groupid);
	  int cselemgroup = theEdger->ElementGroup (boundarydata [i+3]) - theElements.Length();
	  InitializeSet (group, cselemgroup);

	  iArrayT bordergroupids (stop - start + 1);
	  for (int c=0; c < bordergroupids.Length(); c++)
	    bordergroupids[c] = c + start;
	  
	  theEdger->BoundaryFacets (groupid, bordergroupids, set);

	  // insert cse
	  int *pelem = set.Pointer();
	  int *pface = set.Pointer(1);
	  AddElements (set.MajorDim(), cselemgroup);
	  for (int e=0; e < set.MajorDim(); e++, pelem += 2, pface += 2, cs++)
	    {
	      InitializeFacet (*pelem, *pface, group, cs, cselemgroup);
	      if ((set.MajorDim() > kPrint && (e+1)%kPrint == 0) ||
	      e+1 == set.MajorDim())
		cout << "   " << e+1 << " done initializing of " << set.MajorDim() << " facets " << endl;
	    }
	}

      if (boundarydata.Length()/4 < 10) 
	cout << "  Done with Boundary " 
	     << boundarydata[i] << " " << boundarydata [i+1] << endl;
    }
  cout << "  Done with Boundary Data" << endl;
}

void MakeCSE::InitializeSet (int group, int cselemgroup)
{
  // has the group been initialized ?
  GeometryT::CodeT cgeocode = theCSEs[cselemgroup]->GeometryCode();
  if (cgeocode == GeometryT::kNone)
    {
      int numfacenodes = theElements[group]->NumFaceNodes (0);
      GeometryT::CodeT code = theElements[group]->GeometryCode();
      theCSEs[cselemgroup]->Initialize (code, numfacenodes);
    }
}

void MakeCSE::AddElements (int num, int cselemgroup)
{
  int csegroup = cselemgroup + theElements.Length();
  int numfaces = theCSEs[cselemgroup]->NumElemFaces();
  // allocate space for these CS elements
  theCSEs[cselemgroup]->AddElements (num);
  theEdger->AddElements (num, numfaces, csegroup);
}

void MakeCSE::InitializeFacet (int elem, int face, int group, int cs, int cselemgroup)
{
  // find face nodes
  int local = theEdger->LocalElement (elem, group);
  iArrayT facenodes;
  theElements[group]->FaceNodes (local, face, facenodes);

  // collect list of facet nodes to examine for connectivity renumbering
  fPotentialSplitNodes.AppendUnique (facenodes);

  // save surface 1 facets for later
  fSurface1Facets.Append (elem);
  fSurface1Facets.Append (face);

  // find current neighbor
  int neighbor, neighborfacet;
  theEdger->NeighborFacet (elem, face, neighbor, neighborfacet);
      
  // save surface 2 facets for later
  fSurface2Facets.Append (neighbor);
  fSurface2Facets.Append (neighborfacet);

  // set cs element nodes
  int csegroup = cselemgroup + theElements.Length();
  int localcs = theEdger->LocalElement (cs, csegroup);
  theCSEs[cselemgroup]->SetNodes (localcs, facenodes);
  
  // cross link neighbor data
  iArrayT csfaces;
  theCSEs[cselemgroup]->CSElemFaces (csfaces);
  theEdger->SetNeighbor (cs, csfaces[0], elem, face);
  theEdger->SetNeighbor (cs, csfaces[1], neighbor, neighborfacet);
}

void MakeCSE::RemoveSingleNodes (const ArrayT<int>& nodes)
{
  int *ps = nodes.Pointer();
  for (int i=0; i < nodes.Length(); i++, ps++)
    {
      int dex = fPotentialSplitNodes.PositionOf (*ps);
      if (dex > -1) fPotentialSplitNodes.DeleteAt (dex);
    }
}

void MakeCSE::RenumberFaceNodes (void)
{
  cout << "\n Potential Split Nodes " << fPotentialSplitNodes.Length() << "\n";
  iArrayT checkelems (theEdger->TotalElements() - fNumStartElements);
  int *node = fPotentialSplitNodes.Pointer();
  int num = fPotentialSplitNodes.Length(), freq;
  if (num > kPrint*10) freq = kPrint/2;
  else if (num > kPrint) freq = kPrint/10;
  else freq = kPrint/2/10;

  for (int n=0; n < num; n++, node++)
    {
      // gather elements using this node
      iAutoArrayT hit_elems;
      theEdger->HitElements (*node, hit_elems);
      checkelems = 0;

#ifdef _PRINT_DEBUG_
      cout << "\nnode: " << *node << endl;
      cout << "   hit:\n" << hit_elems;
#endif

      while (hit_elems.Length() > 0)
	{
	  iAutoArrayT elems, faces;

	  // create list of element facets to renumber
	  MakeList (*node, elems, faces, hit_elems);

	  // subtract list facets from hit_elems
	  ReduceList (hit_elems, elems, checkelems);

	  // renumber list facets, except the last time
	  if (hit_elems.Length() > 0) ReNumber (*node, elems, faces);
	}

      // tell status to user
      if ((n%freq == 0 && n > 0) || n == num - 1)
	{
	  cClock.Sum (2);
	  cClock.Print (cout);
	  cout << setw(5) << n + 1 << " Done of " << num << endl;
	}
    }
}

void MakeCSE::MakeList (int node, iAutoArrayT& elems, iAutoArrayT& faces, iAutoArrayT& hit_elems)
{
  // just want to look at regular elements
  int e;
  for (e=0; e < hit_elems.Length(); e++)
    if (hit_elems[e] < fNumStartElements) break;

#ifdef _PRINT_DEBUG_
  cout << "     e:" << e << endl;;
#endif
	
  int local, group;
  iArrayT faceswithnode;
  if (e == hit_elems.Length())
    {
      // only have CSEs left, must be massless
      fNoMassNodes.AppendUnique (node);
      
      // append all CS elements to the list
      elems.Append (hit_elems);
      for (int f=0; f < hit_elems.Length(); f++)
	{
	  theEdger->LocalElement (hit_elems[f], local, group);
	  int csegroup = group - theElements.Length();
	  theCSEs[csegroup]->FacesWithNode (local, node, faceswithnode);
	  faces.Append (faceswithnode[0]);
	}
    }
  else
    {
      // find facets on one side of the split node
      theEdger->LocalElement (hit_elems[e], local, group);
      theElements[group]->FacesWithNode (local, node, faceswithnode);
      elems.Append (hit_elems[e]);
      faces.Append (faceswithnode[0]);
      FindNeighbors (hit_elems[e], faceswithnode[0], node, elems, faces);
    }
}

void MakeCSE::FindNeighbors (int elm, int face, int node, iAutoArrayT& elemneighbors, iAutoArrayT& faceneighbors)
{
  // collect list of neighboring element faces using this node 
  iAutoArrayT nelems, nfaces;
  CheckNeighbor (elm, face, node, nelems, nfaces);

  // have we seen this neighboring element before
  for (int n=0; n < nelems.Length(); n++)
    if (!elemneighbors.HasValue (nelems[n]) && nelems[n] > -1)
      {
	// append new face 
	elemneighbors.Append (nelems[n]);
	faceneighbors.Append (nfaces[n]);
	
	// check this neighbor for its neighbors if not a CSE
	if (nelems[n] < fNumStartElements)
	  FindNeighbors (nelems[n], nfaces[n], node, elemneighbors, faceneighbors);
      }
}

/* collects faces of neighboring elements that have this node */
/* return 1, if the face with the node is an external face */
void MakeCSE::CheckNeighbor (int eglobal, int face, int node, iAutoArrayT& nelems, iAutoArrayT& nfaces)
{
  iArrayT facenodes, faces;
  nelems.Free();
  nfaces.Free();

  int elocal = MakeCSE_FEManager::kNotSet, g = MakeCSE_FEManager::kNotSet;
  int e3global = MakeCSE_FEManager::kNotSet;
  int f3 = MakeCSE_FEManager::kNotSet;

  // collect faces using this node from the given element
  theEdger->LocalElement (eglobal, elocal, g);
  theElements[g]->FacesWithNode (elocal, node, faces);

  // examine face neighbors
  int g1 = theElements[g]->GroupNumber();
  int *facei = faces.Pointer();
  for (int i=0; i < faces.Length(); i++, facei++)
    {
      theEdger->NeighborFacet (eglobal, *facei, e3global, f3);

      // only append neighbors that are not CSE
      if (e3global > -1)
	{
	  nelems.Append (e3global);
	  nfaces.Append (f3);
	}
    }
}

void MakeCSE::ReduceList (iAutoArrayT& hit_elems, const iAutoArrayT& elems, const iArrayT& checkelems) const
{
  for (int u=0; u < elems.Length(); u++)
    {
      int dex = hit_elems.PositionOf (elems[u]);
      if (dex > -1)
	{
	  // remove regular elements
	  if (elems[u] < fNumStartElements) hit_elems.DeleteAt (dex);

	  // only remove CSE at first sight
	  //if all regular neighbors are in list
	  else if (checkelems[elems[u] - fNumStartElements]++ == 0)
	    {
	      int local, group;
	      theEdger->LocalElement (elems[u], local, group);
	      int csegroup = group - theElements.Length();
	      int numf = theCSEs[csegroup]->NumElemFaces();
	      int count = 0, neighbor, neighborfacet;
	      for (int c=0; c < numf; c++)
		{
		  theEdger->NeighborFacet (elems[u], c, neighbor, neighborfacet);
		  if ((neighbor < fNumStartElements && neighbor > -1) &&
		      elems.HasValue (neighbor)) count++;
		}
	      if (count == 2) hit_elems.DeleteAt (dex);
	    }
	  // always remove from list the 2nd time you see a CSE
	  else
	    hit_elems.DeleteAt (dex);
	}
    }
}

int MakeCSE::ReNumber (int node, const ArrayT<int>& elems, const ArrayT<int>& faces)
{
  int newnode = theNodes->AddDuplicateCoord (node);
  
  iArrayT facenodes;
  int g, elocal;
  for (int i=0; i < elems.Length(); i++)
    {
#ifdef _PRINT_DEBUG_
      if (i==0) cout << "    newnode: " << newnode << endl;
      cout << "      " << elems[i] << " " << faces[i] << endl;
#endif

      theEdger->LocalElement (elems[i], elocal, g);

      if (g < theElements.Length())
	theElements[g]->ResetOneFaceNode (elocal, faces[i], node, newnode);
      else 
	theCSEs[g - theElements.Length()]->ResetOneFaceNode (elocal, faces[i], node, newnode);

    }
  return newnode;
}

void MakeCSE::CollectMassLessNodes (void)
{
  if (fNoMassNodes.Length() > 0)
    {
      RemoveRepeats (fNoMassNodes);
      out  << "\n Number of Massless Nodes Found. . . . . . . . . = "
	   << fNoMassNodes.Length() << '\n';
      cout << "\n Number of Massless Nodes Found. . . . . . . . . = "
	   << fNoMassNodes.Length() << '\n';
      theNodes->AddNodeSet (fNSetID++, fNoMassNodes, NodeManagerPrimitive::kSurface2);
    }
}

void MakeCSE::CollectSurfaceData (void)
{
  int *pelem = fSurface1Facets.Pointer();
  int *pface = fSurface1Facets.Pointer(1);
  int *pelem2 = fSurface2Facets.Pointer();
  int *pface2 = fSurface2Facets.Pointer(1);
  iAutoArrayT nodes1, nodes2;
  ArrayT<iAutoArrayT> faces1 (theElements.Length());
  ArrayT<iAutoArrayT> faces2 (theElements.Length());
  iArrayT facenodes, csfaces;
  for (int i=0; i < fSurface1Facets.Length()/2; i++)
    {
      int g1 = theEdger->WhichGroup (*pelem);
      int local = theEdger->LocalElement (*pelem, g1);
      int cs, csfacet;
      theEdger->NeighborFacet (*pelem, *pface, cs, csfacet);
      int cg = theEdger->WhichGroup (cs);
      int cgid = theCSEs[cg - theElements.Length()]->GroupNumber();

      if (fContact.HasValue (cgid))
	{
	  // surface 1 facets
	  faces1[g1].Append (local);
	  faces1[g1].Append (*pface);
	  
	  // surface 1 nodes
	  theElements[g1]->FaceNodes (local, *pface, facenodes);
	  nodes1.Append (facenodes);

	  if (*pelem2 > -1)
	    {
	      // surface 2 facets
	      int g2 = theEdger->WhichGroup (*pelem2);
	      int local2 = theEdger->LocalElement (*pelem2, g2);
	      faces2[g2].Append (local2);
	      faces2[g2].Append (*pface2);
	    }

	  // surface 2 nodes
	  theCSEs[cg - theElements.Length()]->CSElemFaces (csfaces);
	  int other = csfacet;
	  if (csfacet == csfaces[0]) other = csfaces[1];
	  int cslocal = theEdger->LocalElement (cs, cg);
	  theCSEs[cg - theElements.Length()]->FaceNodes (cslocal, other, facenodes);
	  nodes2.Append (facenodes);
	}
      pelem += 2;
      pface += 2;
      pelem2 += 2;
      pface2 += 2;
    }
  
  out  << "\n Contact Data, Surface 1 Facets. . . \n";
  cout  << "\n Contact Data, Surface 1 Facets. . . \n";
  for (int j=0; j < theElements.Length(); j++)
    if (faces1[j].Length() > 0)
      {
	iArray2DT temp;
	temp.Set (faces1[j].Length()/2, 2, faces1[j].Pointer());
	theElements[j]->AddSideSet (fSSetID++, temp);
      }

  out  << "\n Contact Data, Surface 2 Facets. . .  \n";
  cout  << "\n Contact Data, Surface 2 Facets. . .  \n";
  for (int j2=0; j2 < theElements.Length(); j2++)
    if (faces2[j2].Length() > 0)
      {
	iArray2DT temp;
	temp.Set (faces2[j2].Length()/2, 2, faces2[j2].Pointer());
	theElements[j2]->AddSideSet (fSSetID++, temp);
      }

  out  << "\n Contact Data, Surface 1 Nodes . . .  \n";
  cout  << "\n Contact Data, Surface 1 Nodes . . .  \n";
  if (nodes1.Length() > 0)
    {
      RemoveRepeats (nodes1);
      theNodes->AddNodeSet (fNSetID++, nodes1, NodeManagerPrimitive::kSurface1);
    }

  out  << "\n Contact Data, Surface 2 Nodes . . . .\n";
  cout  << "\n Contact Data, Surface 2 Nodes . . . .\n";
  if (nodes2.Length() > 0)
    {
      RemoveRepeats (nodes2);
      theNodes->AddNodeSet (fNSetID++, nodes2, NodeManagerPrimitive::kSurface2);
    }
}

// maybe someday this will be added to iArrayT ?
void MakeCSE::RemoveRepeats (ArrayT<int>& n) const
{
      iArrayT nodes;
      nodes.Swap (n);
      nodes.SortAscending();

      // determine number of nodes
      int count = 1;
      for (int m=1; m < nodes.Length(); m++)
	if (nodes[m] != nodes[m-1]) count++;

      // collect nodes, only once
      n.Allocate (count);
      int *pnew = n.Pointer();
      int *pold = nodes.Pointer();
      *pnew++ = *pold++;
      for (int ni=1; ni < nodes.Length(); ni++, *pold++)
	if (*pold != nodes[ni-1]) 
	  *pnew++ = *pold;
}

void MakeCSE::PrintControlEnd (ostream& o)
{
  o << "\n Number of CSE's Created . . . . . . . . . . . . = " 
    << theEdger->TotalElements() - fNumStartElements << '\n';
  o << " Number of Nodes Added . . . . . . . . . . . . . = "
    << theNodes->NumNodes() - fNumStartNodes << '\n';
  o << " Total Number of Nodes . . . . . . . . . . . . . = "
    << theNodes->NumNodes() << '\n';
  o << "\n Done Creating CSE's " << endl;
}
