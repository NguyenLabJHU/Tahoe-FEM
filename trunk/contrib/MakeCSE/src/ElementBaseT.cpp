// file: ElementBaseT.cpp

// created: SAW 10/06/99

#include "ElementBaseT.h"
#include "FEManager.h"
#include "OutputSetT.h"
#include "GeometryBaseT.h"
#include "TriT.h"
#include "QuadT.h"
#include "HexahedronT.h"
#include "TetrahedronT.h"
#include "PentahedronT.h"

using namespace Tahoe;

ElementBaseT::ElementBaseT (ostream& fMainOut, int ID) :
  out (fMainOut),
  fNumElemNodes (0),
  fGeometryCode (GeometryT::kNone),
  fGroupID (ID),
  fOutputID (-1)
{
}

ElementBaseT::~ElementBaseT (void)
{
}

// use this initialize for preexisting data
void ElementBaseT::Initialize (MakeCSEIOManager& theInput)
{
  EchoConnectivity (theInput);
  EchoSideSets (theInput);
}

void ElementBaseT::Initialize (GeometryT::CodeT geocode, int numnodes)
{
  fGeometryCode = geocode;
  fNumElemNodes = numnodes;
  SetFace ();

  out << "\n  Initializing Element Group ID. . . . . . . . . = " 
      << fGroupID << '\n';
  PrintControlData ();
  out << "\n";
}

void ElementBaseT::AddElements (int numelems)
{
  int num = fNodeNums.MajorDim();
  // reallocate space;
  if (fNodeNums.Length() == 0)
    {
      fNodeNums.Allocate (numelems, fNumElemNodes);
      fNodeNums = FEManager::kNotSet;
    }
  else
    fNodeNums.Resize (num + numelems, FEManager::kNotSet);
}

void ElementBaseT::SetNodes (int e1local, const iArrayT& nodes)
{
  if (nodes.Length() != fNumElemNodes)
    {
      cout << "ElementBaseT::Cannot set nodes, wrong length" << endl;
      cout << nodes.Length() << " " << fNumElemNodes << endl;
      throw eSizeMismatch;
    }

  if (!IsElementValid (e1local)) throw eSizeMismatch;
  fNodeNums.SetRow (e1local, nodes);
}

void ElementBaseT::FacesWithNode (int e1local, int node, iArrayT& faces) const
{
  if (!IsElementValid (e1local)) throw eSizeMismatch;
  iAutoArrayT f;
  int *conn = fNodeNums(e1local);
  for (int i=0; i < fNumElemNodes; i++)
    if (*conn++ == node)
	f.AppendUnique (fRevFacetNodes[i]);

  iArrayT temp (f.Length(), f.Pointer());
  faces = temp;
}

bool ElementBaseT::FaceHasNode (int e1local, int f1, int node) const
{
  if (!IsElementValid (e1local) || !IsFaceValid (f1)) 
    { 
      cout << "FaceHasNode" << endl;
      throw eSizeMismatch;
    }
  int *pfN = fFacetNodes[f1].Pointer();
  for (int n=0; n < fFacetNodes[f1].Length(); n++, pfN++)
    if (fNodeNums (e1local, *pfN) == node)
      return true;
  return false;
}

void ElementBaseT::ResetOneFaceNode (int e1local, int f1, int oldnode, int newnode)
{
  if (!IsElementValid (e1local) || !IsFaceValid (f1)) 
    { 
      cout << "ResetOneFaceNode" << endl;
      throw eSizeMismatch;
    }
  int reset = 0;
  int *pfN = fFacetNodes[f1].Pointer();
  for (int n=0; n < fFacetNodes[f1].Length(); n++, pfN++)
    if (*pfN > -1)
      if (fNodeNums (e1local, *pfN) == oldnode)
	{
	  fNodeNums (e1local, *pfN) = newnode;
	  reset++;
	}
  
  if (reset == 0)
    {
      cout << "\n\nResetOneFaceNode failed to find node. " << endl;
      cout << fGroupID << " " << e1local << " " << f1 
	   << " " << oldnode << " " << newnode << endl;
      fNodeNums.PrintRow (e1local, cout);
      throw eGeneralFail;
    }
  
  else if (reset > 1)
    {
      cout << "\n\nResetOneFaceNode reset more than one node. " << endl;
      cout << fGroupID << " " << e1local << " " << f1 
	   << " " << oldnode << " " << newnode << endl;
      fNodeNums.PrintRow (e1local, cout);
      throw eGeneralFail;
    }
}

void ElementBaseT::ElementNodes (int e1local, iArrayT& nodes) const
{
  if (!IsElementValid (e1local)) 
    { 
      cout << "ElementNodes" << endl;
      throw eSizeMismatch;
    }
  nodes.Set (fNodeNums.MinorDim(), fNodeNums(e1local));
}

void ElementBaseT::FaceNodes (int e1local, int f1, iArrayT& nodes) const
{
  if (!IsElementValid (e1local) || !IsFaceValid (f1)) 
    { 
      cout << "FaceNodes" << endl;
      throw eSizeMismatch;
    }
  int *pfN = fFacetNodes[f1].Pointer();
  nodes.Allocate (fFacetNodes[f1].Length());
  for (int i=0; i <fFacetNodes[f1].Length(); i++)
    if (*pfN > -1) // account for ragged array (penta)
      nodes [i] = fNodeNums (e1local, *pfN++);
    else
      nodes [i] = FEManager::kNotSet;
}

void ElementBaseT::AbbrFaceNodes (int e1local, int f1, iArrayT& nodes) const
{
  if (!IsElementValid (e1local) || !IsFaceValid (f1)) 
    { 
      cout << "AbbrFaceNodes" << endl;
      throw eSizeMismatch;
    }
  iArrayT& vertexfacenodes = fVertexFaceNodes[f1];

  nodes.Allocate (vertexfacenodes.Length());
  for (int i=0; i < nodes.Length(); i++)
    if (vertexfacenodes [i] > -1) // account for ragged array (penta)
      nodes [i] = fNodeNums (e1local, vertexfacenodes[i]);
    else
      nodes [i] = FEManager::kNotSet;
}

bool ElementBaseT::CheckSideSet (const iArray2DT& sides) const
{
  int elem = NumElements();
  int face = NumElemFaces();
  int *s = sides.Pointer();
  for (int i=0; i < sides.MajorDim(); i++)
    {
      if (*s > elem || *s < 0) 
	{
	  cout << "ElementBaseT::CheckSideSet: element out of range"
	       << *s << " " << elem << " " << fGroupID << endl;
	  return false;
	}
      s++;
      if (*s > face || *s < 0)
	{
	  cout << "ElementBaseT::CheckSideSet: facet out of range"
	       << *s << " " << face << " " << fGroupID << endl;
	  return false;
	}
      s++;
    }
  return true;
}

void ElementBaseT::AddSideSet (int setID, const iArray2DT& sides)
{
  int dex;
  fSideSetID.HasValue (setID, dex);
  if (dex > -1)
    {
      int length = fSideSetData[dex].MajorDim();
      int num = sides.MajorDim();
      fSideSetData[dex].Resize (length + num, FEManager::kNotSet);
      fSideSetData[dex].BlockRowCopyAt (sides, length);
    }
  else
    {
      int length = fSideSetID.Length();
      fSideSetData.Resize (length + 1);
      fSideSetID.Resize (length + 1, FEManager::kNotSet);
      fSideSetData[length] = sides;
      fSideSetID[length] = setID;
      out << "\n  Element Group ID . . . . . . . . . . . . . . . = " 
	  << fGroupID << '\n';
      out << "    Added Side Set . . . . . . . . . . . . . . . = " 
	  << setID << endl;
      cout << "\n  Element Group ID . . . . . . . . . . . . . . . = " 
	   << fGroupID << '\n';
      cout << "    Added Side Set . . . . . . . . . . . . . . . = " 
	   << setID << endl;
    }
}

void ElementBaseT::Renumber (const iArrayT& map)
{
  // renumber connectivity
  int *n = fNodeNums.Pointer();
  for (int i=0; i < fNodeNums.Length(); i++, n++)
    *n = map [*n];
}


void ElementBaseT::Connectivities (AutoArrayT<const iArray2DT*>& conn, iAutoArrayT& geocodes, iAutoArrayT& change)
{
  conn.Append (&fNodeNums);
  geocodes.Append (fGeometryCode);
  change.Append (0);
}

void ElementBaseT::NodesUsed (iArrayT& nodes_used) const
{
	/* compressed number range */
	int min   = fNodeNums.Min();
	int range = fNodeNums.Max() - min + 1; 

	/* local map */
	iArrayT node_map(range);

	/* determine used nodes */
	node_map = 0;
	for (int i = 0; i < fNodeNums.Length(); i++)
		node_map[fNodeNums[i] - min] = 1;

	/* collect list */
	nodes_used.Allocate(node_map.Count(1));
	int dex = 0; 
	int*  p = node_map.Pointer();
	for (int j = 0; j < node_map.Length(); j++)
		if (*p++ == 1) nodes_used[dex++] = j + min;
}

void ElementBaseT::RegisterOutput (MakeCSEIOManager& theIO)
{
  ArrayT<StringT> n_labels;
  ArrayT<StringT> e_labels;
  //GenerateOutputLabels (codes, n_labels, e_labels);

  bool changing = false;
  //OutputSetT output_set (fGroupID, fGeometryCode, fNodeNums,
  //			 n_labels, e_labels, changing);

  //fOutputID = theIO.AddElementSet (output_set);
  
  //for (int s=0; s < fSideSetData.Length(); s++)
    //theIO.AddSideSet (fSideSetData[s], fSideSetID[s], fOutputID);
}

void ElementBaseT::WriteOutput (MakeCSEIOManager& theIO, IOBaseT::OutputModeT mode) const
{
  // send variable data
  /*iArrayT codes;
  SetOutputCodes (mode, codes);
  int num_out = codes.Sum();

  dArray2DT group_n_values;
  dArray2DT group_e_values (0,0);

  if (num_out > 0)
    {
      ;
    }

    theIO.WriteOutput (fOutputID, group_n_values, group_e_values);*/
}

bool ElementBaseT::IsElementValid (int e1local) const 
{
  if (e1local >= fNodeNums.MajorDim() || e1local < 0)
    {
      cout << "\n\nElementBaseT::IsElementValid, not valid\n";
      cout << "   element = " << e1local 
	   << "\n   num of elements in group = " << fNodeNums.MajorDim();
      cout << "\n   Element Group ID = " << fGroupID << endl;
      return false;
    }
  return true;
}

bool ElementBaseT::IsFaceValid (int face) const
{
  if (face >= fFacetNodes.Length() || face < 0)
    {
      cout << "\n\nElementBaseT::IsFaceValid, not valid\n";
      cout << "   face = " << face+1 
	   << "\n   number of allowed faces = " << fFacetNodes.Length();
      cout << "\n   Element Group ID = " << fGroupID << endl;
      return false;
    }
  return true;
}

// *********** PROTECTED *************

void ElementBaseT::EchoConnectivity (MakeCSEIOManager& theInput)
{
  ReadConnectivity (theInput, fGeometryCode, fNodeNums);
  InitializeConnectivity ();
}

void ElementBaseT::ReadConnectivity (MakeCSEIOManager& theInput, GeometryT::CodeT& geocode, iArray2DT& conn) const
{
  iArrayT map;
  StringT id;
  id.Append (fGroupID);
  geocode = theInput.ElementGroupGeometry (id);
  conn = theInput.ElementGroup (id);
  theInput.ElementIDs (id, map);
}

void ElementBaseT::InitializeConnectivity (void)
{
  fNumElemNodes = fNodeNums.MinorDim();
  SetFace ();

  out << "\n  Element Group ID . . . . . . . . . . . . . . . = " 
      << fGroupID << '\n';
  PrintControlData ();
}

void ElementBaseT::EchoSideSets (MakeCSEIOManager& theInput)
{
  ReadSideSetData (theInput, fSideSetData);
  CheckAllSideSets ();
}

void ElementBaseT::ReadSideSetData (MakeCSEIOManager& theInput, ArrayT<iArray2DT>& Data)
{
  /* read in side sets that are to be transferred */
  iArrayT sides;
  theInput.InputData (sides, MakeCSEIOManager::kCopySide);
  iAutoArrayT ids;

  /* list side sets in this element group */
  for (int s=0; s < sides.Length(); s += 2)
    if (sides[s+1] == fGroupID)
      ids.Append (sides[s]);

  Data.Allocate (ids.Length());
  out << "   Number of Side Sets . . . . . . . . . . . . . = " 
      << Data.Length() << endl;

  /* store side set id for output manager */
  fSideSetID.Allocate (ids.Length());

  /* read side sets */
  for (int i=0; i < Data.Length(); i++)
    {
      StringT sid;
      sid.Append (ids[i]);
      const iArray2DT temp = theInput.SideSet (sid);
      bool local = theInput.IsSideSetLocal (sid);
      if (local)
	theInput.SideSetLocalToGlobal (sid, temp, Data[i]);
      else
	Data[i] = temp;
      out << "    Side Set . . . . . . . . . . . . . . . . . . = " 
	  << ids[i] << '\n';
      out << "     Number of Facets in Set . . . . . . . . . . = "
	  << Data[i].MajorDim() << '\n'; 
      fSideSetID[i] = ids[i];
    }
}

void ElementBaseT::CheckAllSideSets (void)
{
  for (int i=0; i < fSideSetData.Length(); i++)
    {
      if (!CheckSideSet (fSideSetData[i]))
	{
	  cout << "ElementBaseT::CheckAllSideSets, side set " << i
	       << "\nfails CheckSideSet for element group id: " 
	       << fGroupID << " " << NumElements() << " " << NumElemFaces() << endl;
	  fSideSetData[i].WriteNumbered(cout);
	  throw eBadInputValue;
	}
    }    
}

void ElementBaseT::SetFace (void)
{
  GeometryBaseT *geo;
  switch (fGeometryCode)
    {
    case GeometryT::kTriangle:
      geo = new TriT (fNumElemNodes);
      break;

    case GeometryT::kQuadrilateral:
      geo = new QuadT (fNumElemNodes);
      break;

    case GeometryT::kHexahedron:
      geo = new HexahedronT (fNumElemNodes);
      break;

    case GeometryT::kTetrahedron:
      geo = new TetrahedronT (fNumElemNodes);
      break;
      
    case GeometryT::kPentahedron:
      geo = new PentahedronT (fNumElemNodes);
      break;

    default:
      {
	cout << "\n\n ElementBaseT::SetFace does not like GeoCode " 
	     << fGeometryCode << " " << fNumElemNodes << endl;
	throw eBadInputValue;
      }
    }
  
  // determine number of nodes on each facet and facet geometry codes
  iArrayT num_nodes;
  ArrayT<GeometryT::CodeT> fFacetCodes;
  geo->FacetGeometry (fFacetCodes, num_nodes);
       
  // create fFaceNodes, map between element nodes and face nodes
  fFacetNodes.Allocate (num_nodes.Length());
  fRevFacetNodes.Allocate (fNumElemNodes);
  for (int f=0; f < num_nodes.Length(); f++)
    {
      geo->NodesOnFacet (f, fFacetNodes[f]);
      int *fN = fFacetNodes[f].Pointer();
      for (int j=0; j < fFacetNodes[f].Length(); j++)
	fRevFacetNodes[*fN++].Append (f);
    }
       
  // vertexfacenodes, for use with EdgeFinderT and higher oder elements
  fVertexFaceNodes.Allocate (num_nodes.Length());
  iArrayT num (num_nodes.Length());
  if (fGeometryCode != GeometryT::kPentahedron)
    {
      iArray2DT temp;
      geo->NeighborNodeMap (temp);
      num = temp.MinorDim();
    }
  else
    {
      num = 4;
      num[0] = 3;
      num[1] = 3;
    }
  for (int k=0; k < num_nodes.Length(); k++)
    fVertexFaceNodes[k].Set (num[k], fFacetNodes[k].Pointer());
    
  delete geo;
}

void ElementBaseT::PrintControlData (void) const
{
  out << "   Number of Elements. . . . . . . . . . . . . . = "
      << fNodeNums.MajorDim() << '\n'; 
  out << "   Number of Element Nodes . . . . . . . . . . . = "
      << fNumElemNodes << '\n'; 
  out << "   Geometry Code . . . . . . . . . . . . . . . . = "
      << fGeometryCode << "\n"; 
}

