// file: NodeManagerPrimitive.cpp

// created: SAW 10/07/99

#include "NodeManagerPrimitive.h"

#include "FEManager.h"
#include "MakeCSEIOManager.h"
#include "GlobalEdgeFinderT.h"
#include "dArray2DT.h"
#include "dArrayT.h"

NodeManagerPrimitive::NodeManagerPrimitive (ostream& fMainOut, int comments, FEManager& FEM) :
  out (fMainOut),
  fPrintUpdate (comments),
  theBoss (&FEM)
{
}

void NodeManagerPrimitive::Initialize (IOManager& theInput)
{
  EchoCoordinates (theInput);
  EchoNodeSets (theInput);

  fNew2Old.Allocate (fCoordinates.MajorDim());
  for (int i=0; i < fNew2Old.Length(); i++)
    fNew2Old[i] = i;
}

// FUTURE: update coordinate map in GlobalData
int NodeManagerPrimitive::AddCoord (const dArrayT& p)
{
  int newnode = fCoordinates.MajorDim();
  double fill = 0;
  fCoordinates.Resize (newnode + 1, fill, true);
  fCoordinates.WriteRow(newnode, p);
  
  fNew2Old.Resize (newnode + 1, FEManager::kNotSet, true, false);
  fNew2Old[newnode] = newnode;

  return newnode;
}

int NodeManagerPrimitive::AddDuplicateCoord (const int oldnode)
{
  int newnode = fCoordinates.MajorDim();
  double fill = 0;
  fCoordinates.Resize (newnode + 1, fill, true);
  fCoordinates.CopyRowFromRow (newnode, oldnode);

  fNew2Old.Resize (newnode + 1, FEManager::kNotSet, true, false);
  fNew2Old[newnode] = oldnode;

  if (fSplitNodes.AppendUnique (oldnode))
    {
      int length = fSplitNodes.Length();
      fOld2New.Resize (fSplitNodes.Length(), true, false);
      fOld2New[length-1].AppendUnique (newnode);
    }
  else 
    {
      int dex = fSplitNodes.PositionOf (oldnode);
      fOld2New[dex].AppendUnique (newnode);
    }

  return newnode;
}

void NodeManagerPrimitive::AddNodeSet (int setID, const ArrayT<int>& nodes, int transfermethod)
{
  int dex;
  fNodeSetID.HasValue (setID, dex);
  if (dex > -1)
    {
      int num = nodes.Length();
      int length = fNodeSetData[dex].Length();
      fNodeSetData[dex].Resize (length + num, FEManager::kNotSet, true, false);
      fNodeSetData[dex].CopyPart (length, nodes, 0, num);
      RemoveRepeats (fNodeSetData[dex]);
      out << "            Added to Node Set. . . . . . . . . . = " 
	  << setID << '\n' << endl;
      cout << "            Added to Node Set. . . . . . . . . . = " 
	   << setID << '\n' << endl;
    }
  else
    {
      int length = fNodeSetID.Length();
      fNodeSetData.Resize (length + 1, true, false);
      fNodeSetID.Resize (length + 1, FEManager::kNotSet, true, false);
      fNodeSetData[length].Allocate (nodes.Length());
      fNodeSetData[length].CopyPart (0, nodes, 0, nodes.Length());
      fNodeSetID[length] = setID;

      int ml = fTransMethods.Length();
      fTransMethods.Resize (ml + 1, transfermethod, true, false);

      out << "            Added Node Set . . . . . . . . . . . = " 
	  << setID << '\n' << endl;
      cout << "            Added Node Set . . . . . . . . . . . = " 
	   << setID << '\n' << endl;
    }  
}

int NodeManagerPrimitive::OriginalNode (const int node) const
{
  int previous = fNew2Old[node];
  if (previous == node)
    return node;
  else 
    return OriginalNode (previous);
}

void NodeManagerPrimitive::MapNodeSets (const ArrayT<int>& surface1facets, GlobalEdgeFinderT &E)
{
  for (int set=0; set < fTransMethods.Length(); set++)
    {
      out << "\n   Mapping Node Set " << fNodeSetID[set];
      switch (fTransMethods[set])
	{
	case kSurface1:
	  {
	    out << " to Surface 1.\n";
	    SurfaceNodeSet (fNodeSetData[set], true, surface1facets, E);
	    break;
	  }
	case kSurface2:
	  {
	    out << " to Surface 2.\n";
	    SurfaceNodeSet (fNodeSetData[set], false, surface1facets, E);
	    break;
	  }
	  
	case kMap: // use only one node from one surface 
	  {
	    out << " to one side.\n";
	    MapNodeSet (fNodeSetData[set], surface1facets, E);
	    break;
	  }
	  
	case kSplit: // use all nodes on all surfaces
	  {
	    out << " by adding all split nodes.\n";
	    Split (fNodeSetData[set]);
	    break;
	  }
	}
      fNodeSetData[set].SortAscending();
    }	
}

// renumbers nodes by geometric location
void NodeManagerPrimitive::Renumber (int option, iArrayT& map)
{
  cout << "\n Renumbering Nodes " << endl;
  int num = fCoordinates.MajorDim();
  int dof = fCoordinates.MinorDim();
  iAutoArrayT newlist;
  iArrayT offset (fNumInitCoordinates);
  offset = 0;

  // determine if all nodes are to be renumbered
  int start = 1;
  if (option == FEManager::kRenumberAdded) 
    {
      start = fNumInitCoordinates;
      for (int j=0; j < fNumInitCoordinates; j++)
	newlist.Append (j);
    }
  else
    newlist.Append (0);

  // determine new numbering order
  double *next = fCoordinates.Pointer(dof);
  for (int i=start; i < num; i++, next += dof)
    {
      int insert = -1;
      int oldnode = OriginalNode (i); // see if duplicate coord

      if (oldnode >= i)
	{
	  int *compare = newlist.Pointer();
	  for (int j=0; j < newlist.Length(); j++, compare++)
	    for (int d=0; d < dof; d++)
	      if (insert < 0 && next[d] < fCoordinates (*compare, d))
		insert = j;
	}
      else 
	{
	  offset [oldnode]++;
	  insert = newlist.PositionOf (oldnode) + offset [oldnode];
	}

      if (insert > -1)
	newlist.InsertAt (i, insert);
      else
	newlist.Append (i);
    }

  // reorder coordinate list
  dArray2DT old (num, dof);
  old.Swap (fCoordinates);
  dArrayT p (dof);
  for (int m=0; m < num; m++)
    fCoordinates.WriteRow(m, old(newlist[m]));

  map.Allocate (num);
  for (int h=0; h < num; h++)
    map[newlist[h]] = h;

  // map node sets
  cout << "\n Renumbering Node Sets " << endl;
  for (int n=0; n < fNodeSetData.Length(); n++)
    {
      int *nex = fNodeSetData[n].Pointer();
      for (int j=0; j < fNodeSetData[n].Length(); j++, nex++)
	*nex = map [*nex];
      fNodeSetData[n].SortAscending();
    }
}

void NodeManagerPrimitive::RegisterOutput (IOManager& theIO)
{
  iArrayT nodemap (0);
  theIO.SetCoordinates (fCoordinates, &nodemap);

  iArrayT blocktonodesets;
  theIO.InputData (blocktonodesets, MakeCSEIOManager::kBlockToNode);

  iArrayT nodes;
  int setID = fNodeSetID.Max();
  setID++;
  for (int i=0; i < blocktonodesets.Length(); i++)
    {
      out  << "\n Creating Node Set from Element Group ID . . . . = "
	   << blocktonodesets[i] << '\n';
      cout  << "\n Creating Node Set from Element Group ID . . . . = "
	    << blocktonodesets[i] << '\n';
      theBoss->NodesUsed (blocktonodesets[i], nodes);
      AddNodeSet (setID++, nodes, kSplit);
    }

  for (int n=0; n < fNodeSetData.Length(); n++)
    theIO.AddNodeSet (fNodeSetData[n], fNodeSetID[n]);
}

/********** private *****************/

void NodeManagerPrimitive::EchoCoordinates (IOManager& theInput)
{
  iArrayT map;
  theInput.ReadCoordinates (fCoordinates, map);
  fNumInitCoordinates = fCoordinates.MajorDim();

  out << " Number of nodal points. . . . . . . . . . . . . = " 
      << fNumInitCoordinates << '\n';
  out << " Number of nodal degrees of freedom. . . . . . . = " 
      << fCoordinates.MinorDim() << endl << endl;
}

void NodeManagerPrimitive::EchoNodeSets (IOManager& theInput)
{
  /* read in nodes set that are to be transferred */
  iArrayT ids;
  theInput.InputData (ids, MakeCSEIOManager::kMapNodes);

  fNodeSetData.Allocate (ids.Length()/2);
  fTransMethods.Allocate (ids.Length()/2);
  fNodeSetID.Allocate (ids.Length()/2);

  out << "\n N o d e   S e t   D a t a :\n\n";
  out << " Number of Node Sets . . . . . . . . . . . . . . = " 
      << fNodeSetData.Length() << endl;

  /* read node sets */
  for (int i=0, j=0; i < fNodeSetData.Length(); i++, j += 2)
    {
      fTransMethods[i] = ids[j+1];
      theInput.ReadNodeSet (ids[j], fNodeSetData[i]);
      out << "  Node Set . . . . . . . . . . . . . . . . . . . = " 
	  << ids[j] << '\n';
      out << "   Number of Nodes in Set. . . . . . . . . . . . = "
	  << fNodeSetData[i].Length() << '\n'; 

      /* save node set id for output manager */
      fNodeSetID[i] = ids[j];
    }

  out << " Node Set Transfer Methods . . . . . . . . . . . = "
      << fTransMethods.Length() << '\n';
  ids.PrintWithFormat (out, 8, 0, 2);

  /* checks */
  if (fTransMethods.Length() > 0 && 
      (fTransMethods.Min() < kSurface1 || fTransMethods.Max() > kSplit))
    {
      cout << "Invalid Node Transfer Method" << endl;
      throw eBadInputValue;
    }
}

void NodeManagerPrimitive::SurfaceNodeSet (iArrayT& set, bool wantsurf1, const ArrayT<int>& surface1facets, GlobalEdgeFinderT& theEdger)
{
  int *node = set.Pointer();
  int num = set.Length();
  for (int n=0; n < num; n++, node++)
    {
      // assume fOld2New has length of 1 per node
      int dex = fSplitNodes.PositionOf (*node);
      if (dex > -1 && fOld2New[dex].Length() > 0)
	{
	  // see if a surface 1 facet uses that node
	  bool surf1hasnode = theEdger.HasNode (*node, surface1facets);

	  // modify data if needed
	  if ((surf1hasnode && !wantsurf1) ||
	      (!surf1hasnode && wantsurf1))
	    {
	      if (fPrintUpdate)
		{
		  out << "     Node " << *node << " was replaced by "
		      << fOld2New[dex][0] << ".\n";
		}
	      *node = fOld2New[dex][0];
	    }
	  else if (fPrintUpdate)
	    out << "     Node " << *node << " remains.\n";
	}
    }
}

void NodeManagerPrimitive::MapNodeSet (iArrayT& set, const ArrayT<int>& surface1facets, GlobalEdgeFinderT& theEdger)
{
  int *node = set.Pointer();
  int num = set.Length();
  iArrayT nodes;
  for (int n=0; n < num; n++, node++)
    {
      // assume fOld2New has length of 1 per node
      int dex = fSplitNodes.PositionOf (*node);
      if (dex > -1 && fOld2New[dex].Length() > 0)
	{
	  // see if the node is from surface 1
	  bool replace = false;
	  bool surf1has = theEdger.HasNode (*node, surface1facets);

	  if (surf1has)
	    // see if a surface 1 element uses another node from the set
	    if (!theEdger.AnotherNode (*node, surface1facets, set)) 
	      replace = true;
	  
	  // else the node is from surface 2, 
	  //and another node is probably there also, so do nothing

	  if (replace)
	    {
	      if (fPrintUpdate)
		{
		  out << "     Node " << *node << " was replaced by "
		      << fOld2New[dex][0] << ".\n";
		}
	      *node = fOld2New[dex][0];
	    }
	}
    }
}

void NodeManagerPrimitive::Split (iArrayT& set)
{
  int *node = set.Pointer();
  int num = set.Length();
  iAutoArrayT add;
  for (int n=0; n < num; n++, node++)
    {
      int dex = fSplitNodes.PositionOf (*node);
      if (dex > -1)
	add.AppendUnique (fOld2New[dex]);
    }
  
  out << "     Number of Nodes Added to Set. . . . . . . . = " 
      << add.Length() << "\n";
  if (fPrintUpdate)
    {
      for (int i=0; i < add.Length(); i++)
	{
	  out << setw (10) << add[i] + 1;
	  if ((i+1)%6 == 0 && i > 0) out << "\n";
	}
      out << "\n";
    }

  if (add.Length() > 0)
    {
      set.Resize (num + add.Length(), -1, true, false);
      set.CopyPart (num, add, 0, add.Length());
      set.SortAscending ();
    }
}

// maybe someday this will be added to iArrayT ?
void NodeManagerPrimitive::RemoveRepeats (ArrayT<int>& n) const
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

