/* $Id: GlobalEdgeFinderT.cpp,v 1.9 2003-11-21 22:47:39 paklein Exp $ */
#include "GlobalEdgeFinderT.h"
#include "ExceptionT.h"
#include "MakeCSE_FEManager.h"

const int fill_ = GlobalEdgeFinderT::kNoNeighbor;
const int kPrint = 20000;

using namespace Tahoe;

/* constructor */
GlobalEdgeFinderT::GlobalEdgeFinderT(ostream& out) :
  fCurrent (false),
  log (out),
  fNumNodes (0),
  fNumRegular (0),
  fNumElements (0)
{
}

GlobalEdgeFinderT::~GlobalEdgeFinderT (void) {}

void GlobalEdgeFinderT::Initialize (MakeCSE_FEManager& FEM, int num_nodes)
{
  fNumNodes = num_nodes;

  // set up link to regular elements and create global element map 
  fNumRegular = FEM.NumRegularGroups ();
  int numcse = FEM.NumCSEGroups ();
  fNumElements = 0;
  theElements.Allocate (fNumRegular + numcse);
  fElementMap.Allocate (fNumRegular + numcse);
  fElementID.Allocate (fNumRegular + numcse);
  for (int i=0; i < fNumRegular + numcse; i++)
    {
      theElements[i] = FEM.ElementGroup (i);
      fElementID[i] = theElements[i]->GroupNumber();
      int num_elems = theElements[i]->NumElements();
      fElementMap[i].Allocate (num_elems);
      int *p = fElementMap[i].Pointer();
      for (int i2=0; i2 < num_elems; i2++)
	*p++ = i2 + fNumElements;
      fNumElements += num_elems;
    }

  // create reverse element map
  fRevElementMap.Allocate (fNumElements, 2);
  int *rev = fRevElementMap.Pointer();
  for (int p=0; p < fNumRegular + numcse; p++)
    {
      int num_elems = theElements[p]->NumElements();
      for (int g=0; g < num_elems; g++)
	{
	  *rev++ = p;
	  *rev++ = g;
	}
    }

  // allocate neighbor space
  fNeighbors.Allocate (fNumElements);
  fNeighborFacets.Allocate (fNumElements);
  for (int j=0, e = 0; j < fNumRegular + numcse; j++)
    {
      int numfacets = theElements[j]->NumElemFaces();
      int numelems = theElements[j]->NumElements();
      for (int k=0; k < numelems; k++, e++)
	{
	  fNeighbors[e].Allocate (numfacets);
	  fNeighborFacets[e].Allocate (numfacets);
	  fNeighbors[e] = kNoNeighbor;
	  fNeighborFacets[e] = kNoNeighbor;
	}
    }
}

// when you find out how many CSE are in a set for which element group
void GlobalEdgeFinderT::AddElements (int numelems, int numfaces, int group)
{
  // resize element map
  int length = fElementMap[group].Length();
  int numold = fNumElements;
  fElementMap[group].Resize (length + numelems, fill_);

  // resize neighbor data
  fNumElements += numelems;
  int numfacets = theElements[group]->NumElemFaces();
  int neighborlength = fNeighbors.Length();
  fNeighbors.Resize (fNumElements);
  fNeighborFacets.Resize (fNumElements);

  // initialize data
  int *p = fElementMap[group].Pointer(length);
  for (int k=numold; k < fNumElements; k++)
    {
      *p++ = k;      
      fNeighbors[k].Allocate (numfacets);
      fNeighborFacets[k].Allocate (numfacets);
      fNeighbors[k] = kNoNeighbor;
      fNeighborFacets[k] = kNoNeighbor;      
    }

  // resize reverse element map
  int revlength = fRevElementMap.MajorDim();
  fRevElementMap.Resize (revlength + numelems, fill_);
  int *rev = fRevElementMap(numold);
  for (int c=0; c < numelems; c++)
    {
      *rev++ = group;
      *rev++ = c + length;
    }
}

int GlobalEdgeFinderT::ElementGroup (const StringT& groupid) const
{
  int group = -1;
  fElementID.HasValue (groupid, group);
  if (group < 0) 
    {
      cout << "\nGlobalEdgeFinderT::ElementGroup unable to find " 
	   << groupid << endl;
      throw ExceptionT::kGeneralFail;
    }
  return group;
}

int GlobalEdgeFinderT::WhichGroup (int elem) const
{
  int group = CSEConstants::kNotSet;
  if (elem < fNumElements && elem > -1)
    group = fRevElementMap (elem, 0);
  if (group == CSEConstants::kNotSet && elem > -1)
    cout << "GlobalEdgeFinderT::WhichGroup, unable to find " << elem << endl; 
  return group;
}

void GlobalEdgeFinderT::NeighborFacet (int elem, int face, int& neighbor, int& neighborfacet)
{
  //cout << "Neighbor " << elem << " " << face << endl;
  //cout << fNeighbors[elem] << endl;
  if (fNeighbors[elem][face] == kNoNeighbor)
    {
      SetInverseConnects();

      // data about this element
      int elemlocal, elemgroup;
      iArrayT facenodes1;
      LocalElement (elem, elemlocal, elemgroup);
      theElements[elemgroup]->AbbrFaceNodes (elemlocal, face, facenodes1);
      facenodes1.SortAscending();

      // assemble list of elements neighboring this element facet
      iAutoArrayT hit_elems;
      HitElements (facenodes1, hit_elems);
      
      // search for matching facet
      int matchelem = CSEConstants::kNotSet, matchface = CSEConstants::kNotSet;
      int *elem2 = hit_elems.Pointer();
      iArrayT facenodes2;
      int local2, group2;
      for (int i=0; i < hit_elems.Length() && matchelem < 0; i++, elem2++)
	if (*elem2 != elem) // ignore self
	  {
	    LocalElement (*elem2, local2, group2);
	    for (int f=0; f < theElements[group2]->NumElemFaces() && matchelem < 0; f++)
	      {
		theElements[group2]->AbbrFaceNodes (local2, f, facenodes2);
		if (facenodes1.Length() == facenodes2.Length())
		  { 
		    bool match = true;
		    facenodes2.SortAscending();
		    for (int n=0; n < facenodes1.Length() && match; n++)
		      if (facenodes1[n] != facenodes2[n]) match = false;
		    if (match)
		      {
			matchelem = *elem2;
			matchface = f;
		      }
		  }
	      }
	  }
      //cout << "hit\n" << hit_elems << endl;
      //cout << "facenodes\n" << facenodes1 << endl;
      //cout << matchelem << " " << matchface << endl;

      // store findings
      SetNeighbor (elem, face, matchelem, matchface);
    }

  // return data
  neighbor = fNeighbors[elem][face];
  neighborfacet = fNeighborFacets[elem][face];
  
  if (neighbor == kNoNeighbor)
    {
      cout << "GEF::NeighborFacet no neighbor set " << elem << " " << face << endl;
      throw ExceptionT::kGeneralFail;
    }
}

void GlobalEdgeFinderT::SetNeighbor (int elem, int face, int neighbor, int neighborfacet)
{
  if (neighbor < 0) 
    {
      neighbor = kExteriorFacet;
      neighborfacet = kExteriorFacet;
    }
  else
    {
      fNeighbors [neighbor][neighborfacet] = elem;
      fNeighborFacets [neighbor][neighborfacet] = face;
    }

  if (elem > -1)
    {
      fNeighbors [elem][face] = neighbor;
      fNeighborFacets [elem][face] = neighborfacet;
    }
  //cout << "Set " << elem << " " << face << " " << neighbor << " " << neighborfacet << endl;
  //cout << fNeighbors[elem] << endl;
  //if (neighbor > -1) cout << fNeighbors[neighbor] << endl;
}

void GlobalEdgeFinderT::LocalElement (int global, int& local, int& group) const
{
  local = CSEConstants::kNotSet;
  group = CSEConstants::kNotSet;
  if (global < fNumElements && global > -1)
    {
      group = fRevElementMap (global, 0);
      local = fRevElementMap (global, 1);
    }
  if (local == CSEConstants::kNotSet)
    {
    cout << "GlobalEdgeFinderT::LocalElement, unable to convert " << global << " " << fNumElements  << endl;
    fRevElementMap.WriteNumbered (cout);
    }
}

int GlobalEdgeFinderT::LocalElement (int global, int group) const
{
  int g, local;
  LocalElement (global, local, g);
  if (g != group)
    cout << "GlobalEdgeFinderT::LocalElement, unable to match group "
	 << g << " " << group << endl;
  return local;
}

int GlobalEdgeFinderT::GlobalElement (int local, int group) const
{
  int global = CSEConstants::kNotSet;
  if (local < fElementMap[group].Length()) return fElementMap[group][local];
  if (global == CSEConstants::kNotSet)
    cout << "GlobalEdgeFinderT::GlobalElement, unable to convert " << local << " " << group << endl;
  return global;
}

void GlobalEdgeFinderT::ZoneFacets (const StringT& groupid, const sArrayT& zonegroupids, iArray2DT& sideset, iAutoArrayT& boundarynodes)
{
  int group = ElementGroup (groupid);
  iArrayT zonegroups (zonegroupids.Length());
  for (int z=0; z < zonegroupids.Length(); z++)
    zonegroups[z] = ElementGroup (zonegroupids[z]);

  int *elem = fElementMap[group].Pointer();
  iArrayT facenodes;
  iAutoArrayT data;
  for (int e=0; e < fElementMap[group].Length(); e++, *elem++)
    {
      //if ((fElementMap[group].Length() > kPrint && (e+1)%kPrint == 0) ||
      //  e+1 == fElementMap[group].Length())
      //cout << "   " << e+1 << " done of " << fElementMap[group].Length() << " elements " << endl;
      for (int f=0; f < fNeighbors[*elem].Length(); f++)
	{
	  // examine neighbor
	  int neighbor, neighborfacet;
	  NeighborFacet (*elem, f, neighbor, neighborfacet);

	  // only look at higher numbered elements to remove cross reference
	  int neighborgroup = WhichGroup (neighbor);
	  if (neighbor > *elem && zonegroups.HasValue(neighborgroup))
	    {
	      data.Append (*elem);
	      data.Append (f);
	    }
	  
	  // zone boundary facet, exclude CSE facets from previous sets
	  else if (neighborgroup < fNumRegular &&
		   !zonegroups.HasValue (neighborgroup))
	    {
	      int local = LocalElement (*elem, group);
	      theElements[group]->FaceNodes (local, f, facenodes);
	      boundarynodes.AppendUnique (facenodes);
	    }
	}
    }
  sideset.Allocate (data.Length()/2, 2);
  sideset.CopyPart (0, data, 0, data.Length());
}

void GlobalEdgeFinderT::BoundaryFacets (const StringT& groupid, const sArrayT& bordergroupids, iArray2DT& sideset)
{
  int group = ElementGroup (groupid);
  iArrayT bordergroups (bordergroupids.Length());
  for (int z=0; z < bordergroupids.Length(); z++)
    bordergroups[z] = ElementGroup (bordergroupids[z]);

  int *elem = fElementMap[group].Pointer();
  iArrayT facenodes;
  iAutoArrayT data;
  for (int e=0; e < fElementMap[group].Length(); e++, *elem++)
    {
      //if ((fElementMap[group].Length() > kPrint && (e+1)%kPrint == 0) ||
      //  e+1 == fElementMap[group].Length())
      //cout << "   " << e+1 << " done of " << fElementMap[group].Length() << " elements " << endl;
      for (int f=0; f < fNeighbors[*elem].Length(); f++)
        {
          // examine neighbor
          int neighbor, neighborfacet;
          NeighborFacet (*elem, f, neighbor, neighborfacet);
          
          // only look at groups numbered higher to remove cross reference
          int neighborgroup = WhichGroup (neighbor);
          if (bordergroups.HasValue(neighborgroup))
            {
              data.Append (*elem);
              data.Append (f);
            }
        }
    }
  sideset.Allocate (data.Length()/2, 2);
  sideset.CopyPart (0, data, 0, data.Length());
}

bool GlobalEdgeFinderT::HasNode (int node, const ArrayT<int>& facets)
{
  iAutoArrayT hit_elems;
  HitElements (node, hit_elems);
  int edex;
  int *hit = hit_elems.Pointer();
  for (int i=0; i < hit_elems.Length(); i++, hit++)
    {
      const int *pelem = facets.Pointer();
      const int *pface = facets.Pointer(1);
      for (int j=0; j < facets.Length()/2; j++, pelem += 2, pface += 2)
	if (*pelem == *hit)
	  return true;
    }
  return false;
}

bool GlobalEdgeFinderT::AnotherNode (int node, const ArrayT<int>& facets, const iArrayT& set)
{
  for (int i=0; i < set.Length(); i++)
    if (set[i] != node)
      if (HasNode (set[i], facets)) return true;
  return false;
}

void GlobalEdgeFinderT::HitElements (int node, AutoArrayT<int>& hit_elems)
{
  for (int l=0, *elems = fInvConnects[node].Pointer(); 
       l < fInvConnects[node].Length(); l++, *elems++)
    hit_elems.AppendUnique (*elems);
}

/**********************************************************************
 * Private
 **********************************************************************/

// must do with all facenodes, not just vertex nodes, so that the public
// HitElements function works properly
void GlobalEdgeFinderT::SetInverseConnects (void)
{
  if (fCurrent) return;
  fCurrent = true;

  fInvConnects.Allocate (fNumNodes);
  for (int i=0; i < fNumElements; i++)
    {
      int local, group;
      LocalElement (i, local, group);
      iArrayT nodes;
      theElements[group]->ElementNodes (local, nodes);
      for (int j=0; j < nodes.Length(); j++)
	    {
	      if (nodes[j] > -1) 
		fInvConnects[nodes[j]].AppendUnique (i);
	    }
    }
}


void GlobalEdgeFinderT::HitElements (iArrayT& facenodes1, AutoArrayT<int>& hit_elems)
{
  for (int j=0, *row = facenodes1.Pointer(); j < facenodes1.Length(); j++, row++)
    for (int l=0, *elems = fInvConnects[*row].Pointer(); l < fInvConnects[*row].Length(); l++, *elems++)
      hit_elems.AppendUnique (*elems);
}

