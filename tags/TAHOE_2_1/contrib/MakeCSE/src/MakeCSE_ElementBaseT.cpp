// $Id: MakeCSE_ElementBaseT.cpp,v 1.9 2003-11-21 22:47:39 paklein Exp $
// created: SAW 10/06/99
#include "MakeCSE_ElementBaseT.h"

#include "ExceptionT.h"
#include "OutputSetT.h"
#include "GeometryBaseT.h"
#include "TriT.h"
#include "QuadT.h"
#include "HexahedronT.h"
#include "TetrahedronT.h"
#include "PentahedronT.h"
#include "CSEConstants.h"
#include "ModelManagerT.h"
#include "MakeCSE_IOManager.h"
#include "OutputBaseT.h"

using namespace Tahoe;

MakeCSE_ElementBaseT::MakeCSE_ElementBaseT (ostream& fMainOut, const StringT& ID) :
  out (fMainOut),
  fNumElemNodes (0),
  fGeometryCode (GeometryT::kNone),
  fGroupID (ID)
{
}

MakeCSE_ElementBaseT::~MakeCSE_ElementBaseT (void)
{
}

// use this initialize for preexisting data
void MakeCSE_ElementBaseT::Initialize (ModelManagerT& model, MakeCSE_IOManager& theInput)
{
  EchoConnectivity (model);
  EchoSideSets (model, theInput);
}

void MakeCSE_ElementBaseT::Initialize (GeometryT::CodeT geocode, int numnodes)
{
  fGeometryCode = geocode;
  fNumElemNodes = numnodes;
  SetFace ();

  out << "\n  Initializing Element Group ID. . . . . . . . . = " 
      << fGroupID << '\n';
  PrintControlData ();
  out << "\n";
}

void MakeCSE_ElementBaseT::AddElements (int numelems)
{
  int num = fNodeNums.MajorDim();
  // reallocate space;
  if (fNodeNums.Length() == 0)
    {
      fNodeNums.Allocate (numelems, fNumElemNodes);
      fNodeNums = CSEConstants::kNotSet;
    }
  else
    fNodeNums.Resize (num + numelems, CSEConstants::kNotSet);
}

void MakeCSE_ElementBaseT::SetNodes (int e1local, const iArrayT& nodes)
{
  if (nodes.Length() != fNumElemNodes)
    {
      cout << "MakeCSE_ElementBaseT::Cannot set nodes, wrong length" << endl;
      cout << nodes.Length() << " " << fNumElemNodes << endl;
      throw ExceptionT::kSizeMismatch;
    }

  if (!IsElementValid (e1local)) throw ExceptionT::kSizeMismatch;
  fNodeNums.SetRow (e1local, nodes);
}

void MakeCSE_ElementBaseT::FacesWithNode (int e1local, int node, iArrayT& faces) const
{
  if (!IsElementValid (e1local)) throw ExceptionT::kSizeMismatch;
  iAutoArrayT f;
  const int *conn = fNodeNums(e1local);
  for (int i=0; i < fNumElemNodes; i++)
    if (*conn++ == node)
	f.AppendUnique (fRevFacetNodes[i]);

  iArrayT temp (f.Length(), f.Pointer());
  faces = temp;
}

bool MakeCSE_ElementBaseT::FaceHasNode (int e1local, int f1, int node) const
{
  if (!IsElementValid (e1local) || !IsFaceValid (f1)) 
    { 
      cout << "FaceHasNode" << endl;
      throw ExceptionT::kSizeMismatch;
    }
  const int *pfN = fFacetNodes[f1].Pointer();
  for (int n=0; n < fFacetNodes[f1].Length(); n++, pfN++)
    if (fNodeNums (e1local, *pfN) == node)
      return true;
  return false;
}

void MakeCSE_ElementBaseT::ResetOneFaceNode (int e1local, int f1, int oldnode, int newnode)
{
  if (!IsElementValid (e1local) || !IsFaceValid (f1)) 
    { 
      cout << "ResetOneFaceNode" << endl;
      throw ExceptionT::kSizeMismatch;
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
      throw ExceptionT::kGeneralFail;
    }
  
  else if (reset > 1)
    {
      cout << "\n\nResetOneFaceNode reset more than one node. " << endl;
      cout << fGroupID << " " << e1local << " " << f1 
	   << " " << oldnode << " " << newnode << endl;
      fNodeNums.PrintRow (e1local, cout);
      throw ExceptionT::kGeneralFail;
    }
}

void MakeCSE_ElementBaseT::ElementNodes (int e1local, iArrayT& nodes) const
{
  if (!IsElementValid (e1local)) 
    { 
      cout << "ElementNodes" << endl;
      throw ExceptionT::kSizeMismatch;
    }
  nodes.Alias(fNodeNums.MinorDim(), fNodeNums(e1local));
}

void MakeCSE_ElementBaseT::FaceNodes (int e1local, int f1, iArrayT& nodes) const
{
  if (!IsElementValid (e1local) || !IsFaceValid (f1)) 
    { 
      cout << "FaceNodes" << endl;
      throw ExceptionT::kSizeMismatch;
    }
  const int *pfN = fFacetNodes[f1].Pointer();
  nodes.Allocate (fFacetNodes[f1].Length());
  for (int i=0; i <fFacetNodes[f1].Length(); i++)
    if (*pfN > -1) // account for ragged array (penta)
      nodes [i] = fNodeNums (e1local, *pfN++);
    else
      nodes [i] = CSEConstants::kNotSet;
}

void MakeCSE_ElementBaseT::AbbrFaceNodes (int e1local, int f1, iArrayT& nodes) const
{
  if (!IsElementValid (e1local) || !IsFaceValid (f1)) 
    { 
      cout << "AbbrFaceNodes" << endl;
      throw ExceptionT::kSizeMismatch;
    }
  const iArrayT& vertexfacenodes = fVertexFaceNodes[f1];

  nodes.Allocate (vertexfacenodes.Length());
  for (int i=0; i < nodes.Length(); i++)
    if (vertexfacenodes [i] > -1) // account for ragged array (penta)
      nodes [i] = fNodeNums (e1local, vertexfacenodes[i]);
    else
      nodes [i] = CSEConstants::kNotSet;
}

bool MakeCSE_ElementBaseT::CheckSideSet (const iArray2DT& sides) const
{
  int elem = NumElements();
  int face = NumElemFaces();
  const int *s = sides.Pointer();
  for (int i=0; i < sides.MajorDim(); i++)
    {
      if (*s > elem || *s < 0) 
	{
	  cout << "MakeCSE_ElementBaseT::CheckSideSet: element out of range  "
	       << *s << " " << elem << " " << fGroupID << endl;
	  return false;
	}
      s++;
      if (*s > face || *s < 0)
	{
	  cout << "MakeCSE_ElementBaseT::CheckSideSet: facet out of range  "
	       << *s << " " << face << " " << fGroupID << endl;
	  return false;
	}
      s++;
    }
  return true;
}

void MakeCSE_ElementBaseT::AddSideSet (const StringT& setID, const iArray2DT& sides)
{
  int dex;
  fSideSetID.HasValue (setID, dex);
  if (dex > -1)
    {
      int length = fSideSetData[dex].MajorDim();
      int num = sides.MajorDim();
      fSideSetData[dex].Resize (length + num, CSEConstants::kNotSet);
      fSideSetData[dex].BlockRowCopyAt (sides, length);
    }
  else
    {
      int length = fSideSetID.Length();
      fSideSetData.Resize (length + 1);
      fSideSetID.Resize (length + 1, "");
      fSideSetData[length] = sides;
      fSideSetID[length] = setID;
      out << "\n  Element Group ID . . . . . . . . . . . . . . . = " 
	  << fGroupID << '\n';
      out << "    Added Side Set . . . . . . . . . . . . . . . = " 
	  << setID << endl;
    }
}

void MakeCSE_ElementBaseT::Renumber (const iArrayT& map)
{
  // renumber connectivity
  int *n = fNodeNums.Pointer();
  for (int i=0; i < fNodeNums.Length(); i++, n++)
    *n = map [*n];
}

void MakeCSE_ElementBaseT::NodesUsed (iArrayT& nodes_used) const
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

void MakeCSE_ElementBaseT::RegisterOutput (OutputBaseT& output)
{
  ArrayT<StringT> n_labels;
  ArrayT<StringT> e_labels;
  bool changing = false;
  ArrayT<StringT> block_ID (1);
  ArrayT<const iArray2DT*> conns (1);

  block_ID[0] = fGroupID;
  conns[0] = &fNodeNums;
  OutputSetT output_set (fGeometryCode, block_ID, conns, n_labels, e_labels, changing);

  int outID = output.AddElementSet (output_set);
  
  for (int s=0; s < fSideSetData.Length(); s++)
    output.AddSideSet (fSideSetData[s], fSideSetID[s], fGroupID);
}

bool MakeCSE_ElementBaseT::IsElementValid (int e1local) const 
{
  if (e1local >= fNodeNums.MajorDim() || e1local < 0)
    {
      cout << "\n\nMakeCSE_ElementBaseT::IsElementValid, not valid\n";
      cout << "   element = " << e1local 
	   << "\n   num of elements in group = " << fNodeNums.MajorDim();
      cout << "\n   Element Group ID = " << fGroupID << endl;
      return false;
    }
  return true;
}

bool MakeCSE_ElementBaseT::IsFaceValid (int face) const
{
  if (face >= fFacetNodes.Length() || face < 0)
    {
      cout << "\n\nMakeCSE_ElementBaseT::IsFaceValid, not valid\n";
      cout << "   face = " << face+1 
	   << "\n   number of allowed faces = " << fFacetNodes.Length();
      cout << "\n   Element Group ID = " << fGroupID << endl;
      return false;
    }
  return true;
}

// *********** PROTECTED *************

void MakeCSE_ElementBaseT::EchoConnectivity (ModelManagerT& theInput)
{
  ReadConnectivity (theInput, fGeometryCode, fNodeNums);
  InitializeConnectivity ();
}

void MakeCSE_ElementBaseT::ReadConnectivity (ModelManagerT& theInput, GeometryT::CodeT& geocode, iArray2DT& conn) const
{
  iArrayT map;
  StringT id;
  id.Append (fGroupID);
  geocode = theInput.ElementGroupGeometry (id);
  conn = theInput.ElementGroup (id);
  map.Dimension (conn.MajorDim());
  theInput.ElementIDs (id, map);
}

void MakeCSE_ElementBaseT::InitializeConnectivity (void)
{
  fNumElemNodes = fNodeNums.MinorDim();
  SetFace ();

  out << "\n  Element Group ID . . . . . . . . . . . . . . . = " 
      << fGroupID << '\n';
  PrintControlData ();
}

void MakeCSE_ElementBaseT::EchoSideSets (ModelManagerT& model, MakeCSE_IOManager& theInput)
{
  ReadSideSetData (model, theInput, fSideSetData);
  CheckAllSideSets ();
}

void MakeCSE_ElementBaseT::ReadSideSetData (ModelManagerT& model, MakeCSE_IOManager& theInput, ArrayT<iArray2DT>& Data)
{
  /* read in side sets that are to be transferred */
  sArrayT sides;
  theInput.SideSetsMapped (sides);
  AutoArrayT<StringT> ids;

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
      const iArray2DT temp = model.SideSet (ids[i]);
      bool local = model.IsSideSetLocal (ids[i]);
      const StringT elgroupid = model.SideSetGroupID (ids[i]);
      //if (local)
      //model.SideSetLocalToGlobal (elgroupid, temp, Data[i]);
      //else
	Data[i] = temp;
      out << "    Side Set . . . . . . . . . . . . . . . . . . = " 
	  << ids[i] << '\n';
      out << "     Number of Facets in Set . . . . . . . . . . = "
	  << Data[i].MajorDim() << '\n'; 
      fSideSetID[i] = ids[i];
     }
}

void MakeCSE_ElementBaseT::CheckAllSideSets (void)
{
  for (int i=0; i < fSideSetData.Length(); i++)
    {
      if (!CheckSideSet (fSideSetData[i]))
	{
	  cout << "MakeCSE_ElementBaseT::CheckAllSideSets, side set " << fSideSetID[i]
	       << "\nfails CheckSideSet for element group id: " 
	       << fGroupID << " " << NumElements() << " " << NumElemFaces() << endl;
	  fSideSetData[i].WriteNumbered(cout);
	  throw ExceptionT::kBadInputValue;
	}
    }    
}

void MakeCSE_ElementBaseT::SetFace (void)
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
	cout << "\n\n MakeCSE_ElementBaseT::SetFace does not like GeoCode " 
	     << fGeometryCode << " " << fNumElemNodes << endl;
	throw ExceptionT::kBadInputValue;
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

void MakeCSE_ElementBaseT::PrintControlData (void) const
{
  out << "   Number of Elements. . . . . . . . . . . . . . = "
      << fNodeNums.MajorDim() << '\n'; 
  out << "   Number of Element Nodes . . . . . . . . . . . = "
      << fNumElemNodes << '\n'; 
  out << "   Geometry Code . . . . . . . . . . . . . . . . = "
      << fGeometryCode << "\n"; 
}

