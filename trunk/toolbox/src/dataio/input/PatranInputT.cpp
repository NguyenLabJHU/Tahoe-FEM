/* $Id: PatranInputT.cpp,v 1.2 2001-08-07 23:11:54 paklein Exp $ */
/* created: sawimme July 2001 */

#include "PatranInputT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"

PatranInputT::PatranInputT (ostream& out) :
  InputBaseT (out),
  fPatran (out)
{
}

void PatranInputT::Open (const StringT& file)
{
  fPatran.OpenRead (file);
}

void PatranInputT::Close (void)
{
}

void PatranInputT::ElementGroupNames (ArrayT<StringT>& groupnames)
{
  int count = 0, numelems, numelemnodes;
  int numcomps = fPatran.NumNamedComponents ();
  ArrayT<StringT> names (numcomps);
  if (!fPatran.NamedComponents (names)) throw eGeneralFail;
  for (int i=0; i < numcomps; i++)
    {
      if (!fPatran.ReadElementBlockDims (names[i], numelems, numelemnodes)) 
	throw eGeneralFail;
      if (numelems > 0)
	  groupnames[count++] = names[i];
    }
}

void PatranInputT::NodeSetNames (ArrayT<StringT>& nodenames)
{
  int count = 0;
  int numcomps = fPatran.NumNamedComponents ();
  ArrayT<StringT> names (numcomps);
  iArrayT nodes;
  if (!fPatran.NamedComponents (names)) throw eGeneralFail;
  for (int i=0; i < numcomps; i++)
    {
      if (!fPatran.ReadNodeSet (names[i], nodes)) throw eGeneralFail;
      if (nodes.Length() > 0)
	nodenames[count++] = names[i];
    }
}

int  PatranInputT::NumElementGroups (void)
{
  int count = 0, numelems, numelemnodes;
  int numcomps = fPatran.NumNamedComponents ();
  ArrayT<StringT> names (numcomps);
  if (!fPatran.NamedComponents (names)) throw eGeneralFail;
  for (int i=0; i < numcomps; i++)
    {
      if (!fPatran.ReadElementBlockDims (names[i], numelems, numelemnodes)) 
	return false;
      if (numelems > 0)
	count++;
    }
  return count;
}

int  PatranInputT::NumNodeSets (void)
{
  int count = 0;
  int numcomps = fPatran.NumNamedComponents ();
  ArrayT<StringT> names (numcomps);
  iArrayT nodes;
  if (!fPatran.NamedComponents (names)) throw eGeneralFail;
  for (int i=0; i < numcomps; i++)
    {
      if (!fPatran.ReadNodeSet (names[i], nodes)) throw eGeneralFail;
      if (nodes.Length() > 0)
	count++;
    }
  return count;
}

void PatranInputT::ReadNodeMap (iArrayT& nodemap)
{
  if (!fPatran.ReadGlobalNodeMap (nodemap)) throw eGeneralFail;
  nodemap += -1;

  /* later adjust to be consecutive numbering?? */
}

void PatranInputT::ReadCoordinates (dArray2DT& coords)
{
  if (!fPatran.ReadCoordinates (coords, coords.MinorDim()))
    throw eGeneralFail;
}

void PatranInputT::ReadCoordinates (dArray2DT& coords, iArrayT& nodemap)
{
  ReadCoordinates (coords);
  ReadNodeMap (nodemap);
}

int PatranInputT::NumElements (StringT& name)
{
  int num, numnodes;
  if (!fPatran.ReadElementBlockDims (name, num, numnodes))
    throw eGeneralFail;
  return num;
}

int PatranInputT::NumElementNodes (StringT& name)
{
  int num, numnodes;
  if (!fPatran.ReadElementBlockDims (name, num, numnodes))
    throw eGeneralFail;
  return numnodes;  
}

void PatranInputT::ReadAllElementMap (iArrayT& elemmap)
{
  cout << "\n\n PatranInputT::Not programmed to read all element map\n\n";
  elemmap = -1;
}

void PatranInputT::ReadGlobalElementMap (StringT& name, iArrayT& elemmap)
{
  int namedtype;
  if (!fPatran.ReadElementSet (name, namedtype, elemmap))
    throw eGeneralFail;

  /* someday, convert from discontinuous to continuous numbering */

  elemmap += -1;
}

void PatranInputT::ReadConnectivity (StringT& name, iArray2DT& connects)
{
  int namedtype;
  if (!fPatran.ReadConnectivity (name, namedtype, connects))
    throw eGeneralFail;

  /* someday, convert from discontinuous to continuous numbering */

  connects += -1;
}

void PatranInputT::ReadGeometryCode (StringT& name, GeometryT::CodeT& code)
{
  iArrayT elems;
  int namedtype;
  if (!fPatran.ReadElementSet (name, namedtype, elems))
    throw eGeneralFail;

  SetCode (namedtype, code);
}


int PatranInputT::NumNodesInSet (StringT& name)
{
  int num;
  if (!fPatran.NumNodesInSet (name, num)) throw eGeneralFail;
  return num;
}

void PatranInputT::ReadNodeSet (StringT& name, iArrayT& nodes)
{
  if (!fPatran.ReadNodeSet (name, nodes)) throw eGeneralFail;

  /* someday, convert from discontinuous to continuous numbering */

  nodes += -1;
}

int PatranInputT::NumSidesInSet (StringT& anme)
{
#pragma unused(anme)
  cout << "\n\n PatranInputT::Not programmed to read side sets\n\n";
  return 0;
}

int PatranInputT::SideSetGroupIndex (StringT& name)
{
#pragma unused(name)
  cout << "\n\n PatranInputT::Not programmed to read side sets\n\n";
  return 0;
}

void PatranInputT::ReadSideSetLocal (StringT& name, iArray2DT& sides)
{
#pragma unused(name)
#pragma unused(sides)
  cout << "\n\n PatranInputT::Not programmed to read side sets\n\n";
  throw eGeneralFail;
}

void PatranInputT::ReadSideSetGlobal (StringT& name, iArray2DT& sides)
{
#pragma unused(name)
#pragma unused(sides)
  cout << "\n\n PatranInputT::Not programmed to read side sets\n\n";
  throw eGeneralFail;
}


/**************** PRIVATE *******************/

void PatranInputT::SetCode (int namedtype, GeometryT::CodeT& code) const
{
  switch (namedtype)
    {
    case PatranT::kNCLine:
    case PatranT::kNCLine2:
    case PatranT::kNCLine3:
      {
	code = GeometryT::kLine;
	break;
      }
    case PatranT::kNCTriangle: 
    case PatranT::kNCTriangle2:
    case PatranT::kNCTriangle3:
      {
	code = GeometryT::kTriangle;
	break;
      }
    case PatranT::kNCQuad: 
    case PatranT::kNCQuad2:
    case PatranT::kNCQuad3:
       {
	 code = GeometryT::kQuadrilateral;
	 break;
       }
    case PatranT::kNCTet:
    case PatranT::kNCTet2:
    case PatranT::kNCTet3:
      {
	code = GeometryT::kTetrahedron;
	break;
      }
    case PatranT::kNCWedge: 
    case PatranT::kNCWedge2:
    case PatranT::kNCWedge3:
      {
	code = GeometryT::kPentahedron;
	break;
      }
    case PatranT::kNCHex:
    case PatranT::kNCHex2:
    case PatranT::kNCHex3:
      {
	code = GeometryT::kHexahedron;
	break;
      }
    default:
      code = GeometryT::kNone;
    }
}
