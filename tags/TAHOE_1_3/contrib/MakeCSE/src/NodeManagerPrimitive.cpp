// $Id: NodeManagerPrimitive.cpp,v 1.9 2002-11-05 13:26:26 sawimme Exp $
// created: SAW 10/07/99
#include "NodeManagerPrimitive.h"

#include "ExceptionT.h"
#include "MakeCSE_IOManager.h"
#include "MakeCSE_FEManager.h"
#include "dArrayT.h"
#include "CSEConstants.h"
#include "OutputBaseT.h"

using namespace Tahoe;

NodeManagerPrimitive::NodeManagerPrimitive (ostream& fMainOut, bool comments, MakeCSE_FEManager& FEM) :
  out (fMainOut),
  fPrintUpdate (comments),
  theBoss (&FEM)
{
}

void NodeManagerPrimitive::Initialize (ModelManagerT& model, MakeCSE_IOManager& theInput)
{
  EchoCoordinates (model);
  EchoNodeSets (model, theInput);

  fNew2Old.Allocate (fCoordinates.MajorDim());
  for (int i=0; i < fNew2Old.Length(); i++)
    fNew2Old[i] = i;
}

// FUTURE: update coordinate map in GlobalData
int NodeManagerPrimitive::AddCoord (const dArrayT& p)
{
  int newnode = fCoordinates.MajorDim();
  double fill = 0;
  fCoordinates.Resize (newnode + 1, fill);
  fCoordinates.SetRow(newnode, p);
  
  fNew2Old.Resize (newnode + 1, CSEConstants::kNotSet);
  fNew2Old[newnode] = newnode;

  return newnode;
}

int NodeManagerPrimitive::AddDuplicateCoord (const int oldnode)
{
  int newnode = fCoordinates.MajorDim();
  double fill = 0;
  fCoordinates.Resize (newnode + 1);
  fCoordinates.CopyRowFromRow (newnode, oldnode);

  fNew2Old.Resize (newnode + 1, CSEConstants::kNotSet);
  fNew2Old[newnode] = oldnode;

  if (fSplitNodes.AppendUnique (oldnode))
    {
      int length = fSplitNodes.Length();
      fOld2New.Resize (fSplitNodes.Length());
      fOld2New[length-1].AppendUnique (newnode);
    }
  else 
    {
      int dex = fSplitNodes.PositionOf (oldnode);
      fOld2New[dex].AppendUnique (newnode);
    }

  return newnode;
}

void NodeManagerPrimitive::AddNodeSet (const StringT& setID, const ArrayT<int>& nodes, CSEConstants::NodeMapMethodT transfermethod)
{
  int dex;
  fNodeSetID.HasValue (setID, dex);
  if (dex > -1)
    {
      int num = nodes.Length();
      int length = fNodeSetData[dex].Length();
      fNodeSetData[dex].Resize (length + num, CSEConstants::kNotSet);
      fNodeSetData[dex].CopyPart (length, nodes, 0, num);
      RemoveRepeats (fNodeSetData[dex]);
      out << "            Added to Node Set. . . . . . . . . . = " 
	  << setID << '\n' << endl;
    }
  else
    {
      int length = fNodeSetID.Length();
      fNodeSetData.Resize (length + 1);
      fNodeSetID.Resize (length + 1, "");
      fNodeSetData[length].Allocate (nodes.Length());
      fNodeSetData[length].CopyPart (0, nodes, 0, nodes.Length());
      fNodeSetID[length] = setID;

      int ml = fTransMethods.Length();
      fTransMethods.Resize (ml + 1, transfermethod);

      out << "            Added Node Set . . . . . . . . . . . = " 
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
	case CSEConstants::kSurface1:
	  {
	    out << " to Surface 1.\n";
	    SurfaceNodeSet (fNodeSetData[set], true, surface1facets, E);
	    break;
	  }
	case CSEConstants::kSurface2:
	  {
	    out << " to Surface 2.\n";
	    SurfaceNodeSet (fNodeSetData[set], false, surface1facets, E);
	    break;
	  }
	  
	case CSEConstants::kMap: // use only one node from one surface 
	  {
	    out << " to one side.\n";
	    MapNodeSet (fNodeSetData[set], surface1facets, E);
	    break;
	  }
	  
	case CSEConstants::kSplit: // use all nodes on all surfaces
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
void NodeManagerPrimitive::Renumber (CSEConstants::RenumberMethodT option, iArrayT& map)
{
  cout << "\n Renumbering Nodes " << endl;
  int num = fCoordinates.MajorDim();
  int dof = fCoordinates.MinorDim();
  iAutoArrayT newlist;
  iArrayT offset (fNumInitCoordinates);
  offset = 0;

  // determine if all nodes are to be renumbered
  int start = 1;
  if (option == CSEConstants::kRenumberAdded) 
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
    fCoordinates.SetRow(m, old(newlist[m]));

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

void NodeManagerPrimitive::RegisterOutput (OutputBaseT& output, MakeCSE_IOManager& input)
{
  // FUTURE: carry over node tags
  fOutputNodeMap.Dimension (fCoordinates.MajorDim());
  fOutputNodeMap.SetValueToPosition ();
  fOutputNodeMap ++;

  output.SetCoordinates (fCoordinates, &fOutputNodeMap);

  sArrayT blocktonodesets;
  input.BlockToNode (blocktonodesets);

  int nsetid = 1;
  
  for (int g=0; g < fNodeSetID.Length(); g++)
    if  (nsetid < atoi (fNodeSetID[g])) 
      nsetid = atoi (fNodeSetID[g]) + 1;

  iArrayT nodes;
  for (int i=0; i < blocktonodesets.Length(); i++)
    {
      StringT name;
      name.Append (nsetid + 1);
      out  << "\n Creating Node Set from Element Group ID . . . . = "
	   << blocktonodesets[i] << '\n';
      theBoss->NodesUsed (blocktonodesets[i], nodes);
      if (nodes.Length() > 0)
	AddNodeSet (name, nodes, CSEConstants::kSplit);
      else
	out << "\n     No nodes found...\n";
    }

  for (int n=0; n < fNodeSetData.Length(); n++)
    output.AddNodeSet (fNodeSetData[n], fNodeSetID[n]);
}

/********** private *****************/

void NodeManagerPrimitive::EchoCoordinates (ModelManagerT& theInput)
{
  iArrayT map;
  fCoordinates = theInput.Coordinates ();
  theInput.AllNodeIDs (map);
  fNumInitCoordinates = fCoordinates.MajorDim();

  out << " Number of nodal points. . . . . . . . . . . . . = " 
      << fNumInitCoordinates << '\n';
  out << " Number of nodal degrees of freedom. . . . . . . = " 
      << fCoordinates.MinorDim() << endl << endl;
}

void NodeManagerPrimitive::EchoNodeSets (ModelManagerT& model, MakeCSE_IOManager& theInput)
{
  /* read in nodes set that are to be transferred */
  theInput.NodeSetsMapped (fNodeSetID, fTransMethods);
  fNodeSetData.Allocate (fNodeSetID.Length());

  out << "\n N o d e   S e t   D a t a :\n\n";
  out << " Number of Node Sets . . . . . . . . . . . . . . = " 
      << fNodeSetData.Length() << endl;

  /* read node sets */
  for (int i=0; i < fNodeSetData.Length(); i++)
    {
      fNodeSetData[i] = model.NodeSet (fNodeSetID[i]);
      out << "  Node Set . . . . . . . . . . . . . . . . . . . = " 
	  << fNodeSetID[i] << '\n';
      out << "   Number of Nodes in Set. . . . . . . . . . . . = "
	  << fNodeSetData[i].Length() << '\n'; 
    }

  out << " Node Set Transfer Methods . . . . . . . . . . . = "
      << fTransMethods.Length() << '\n';
  out << "  Method  Set ID\n";
  for (int g=0; g < fNodeSetID.Length(); g++)
    out << setw (kIntWidth) << fTransMethods[g] << setw (kIntWidth) 
	<< fNodeSetID[g] << '\n';

  /* checks */
  bool okay = true;
  if (fTransMethods.Length() < 0) okay = false;
  for (int t=0; t < fTransMethods.Length() && okay; t++)
    if (fTransMethods[t] < CSEConstants::kSurface1 || fTransMethods[t] > CSEConstants::kSplit)
      okay = false;
  if (!okay)
    {
      cout << "Invalid Node Transfer Method" << endl;
      throw ExceptionT::kBadInputValue;
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
      set.Resize (num + add.Length(), -1);
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

